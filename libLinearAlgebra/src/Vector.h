// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Utils/Array.h"
#include "Utils/Parallel.h"
#include "Simd/SimdBlock.h"

#include "Scalar.h"

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++                             Generic vectors.                             +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

template<typename T, alg_loc_t TLoc>
using alg_vector_base = alg_array<T>;

/// @defgroup alg_vector Vectors, linear algebra and some other operations
/// @{

///
/// Generic column vector.
///
/// @f$ \vec x = [ x_0, x_1, ..., x_{N-1} ]^T @f$
/// @see https://en.wikipedia.org/wiki/Vector_space
///
template<typename T, alg_loc_t TLoc = alg_loc_host>
class alg_vector : public alg_vector_base<T, TLoc>
{
public:

	/// Initializes a zero vector:
	/// @f$ x_i := 0 @f$
	/// @{
	ALG_INLINE explicit alg_vector()
		: alg_vector_base<T, TLoc>() {}
	ALG_INLINE explicit alg_vector(const size_t n)
		: alg_vector_base<T, TLoc>(n) {}
	/// @}

	/// Initializes a vector with all elements set to some value:
	/// @f$ x_i := t @f$
	ALG_INLINE explicit alg_vector(const size_t n, const T& t)
		: alg_vector_base<T, TLoc>(n, t) {}

	/// Initializes a vector with all elements are copied from some array:
	/// @f$ x_i := d_i @f$
	ALG_INLINE explicit alg_vector(const alg_array<T>& data)
		: alg_vector_base<T, TLoc>(data) {}

	/// Initializes a vector with all elements are copied from some array slice:
	/// @f$ x_i := d_{o+i*s} @f$
	ALG_INLINE explicit alg_vector(const alg_array<T>& data, const size_t offset, const size_t size, const size_t stride = 1)
		: alg_vector_base<T, TLoc>(data, offset, size, stride) {}

	/// Initializes x vector with initial data:
	/// @f$ x_i := d_i @f$
	ALG_INLINE alg_vector(const std::initializer_list<T>& data)
		: alg_vector_base<T, TLoc>(data) {}

public:
	
	/// @f$ r := x_i @f$ 
	/// @{
	ALG_INLINE T& operator() (const size_t i) { return (*this)[i]; }
	ALG_INLINE const T& operator() (const size_t i) const
	{
		return const_cast<const T&>(const_cast<alg_vector<T>&>(*this)(i));
	}
	/// @}

};	// class alg_vector<T>

template<typename T>
using alg_vector_host_t = alg_vector<T, alg_loc_host>;
template<typename T>
using alg_vector_device_t = alg_vector<T, alg_loc_device>;

///
/// Component type.
///
template<typename>
struct alg_component;
template<typename T, alg_loc_t TLoc>
struct alg_component<alg_vector<T, TLoc>> { using type = T; };

template<typename V>
using alg_component_t = typename alg_component<V>::type;

///
/// Related scalar type.
///
template<>
struct alg_related_scalar<alg_vector<alg_float32_t>> { using type = alg_float32_t; };
template<>
struct alg_related_scalar<alg_vector<alg_float64_t>> { using type = alg_float64_t; };
template<typename T, alg_loc_t TLoc>
struct alg_related_scalar<alg_vector<T, TLoc>> { using type = typename alg_related_scalar<T>::type; };

// ------------------------------------------------------------------------------------ //
// -----                          Vector kernel evaluators.                       ----- //
// ------------------------------------------------------------------------------------ //

///
/// Vector kernel evaluators.
///
template<typename TEnable, alg_loc_t TLoc, typename T, typename... Tk>
struct alg_eval_kernel_impl;
template<alg_loc_t TLoc, typename T, typename... Tk>
using alg_eval_kernel_impl_t = alg_eval_kernel_impl<alg_enabled_t, TLoc, T, Tk...>;

///
/// Vector reduction kernel evaluators.
///
template<typename TEnable, alg_loc_t TLoc, typename T, typename... Tk>
struct alg_eval_reduce_kernel_impl;
template<alg_loc_t TLoc, typename T, typename... Tk>
using alg_eval_reduce_kernel_impl_t = alg_eval_reduce_kernel_impl<alg_enabled_t, TLoc, T, Tk...>;

// ------------------------------------------------------------------------------------ //
// -----                       Host Vector kernel evaluators.                     ----- //
// ------------------------------------------------------------------------------------ //

#if ALG_EVAL_ON_HOST
///
/// Host Vector kernel evaluator for types without acceleration support.
///
template<typename T, typename... Tk>
struct alg_eval_kernel_impl<
	alg_enable_if_t<!ALG_HAS_SIMD_SUPPORT(T)>,
	alg_loc_host, T, Tk...>
{
public:
	template<typename TKernel>
	ALG_INLINE static void eval(const TKernel& kernel,
        alg_vector_host_t<T>& x0, const alg_vector_host_t<Tk>&... xk)
	{
		size_t size;
		alg_verify_prop(size, x0, xk...);
		alg_nonparallel_for(size, [&](const size_t i) { kernel(x0(i), xk(i)...); });
	}
};	// struct alg_eval_kernel_impl</*no SIMD for T*/, alg_loc_host, T, Tk...>
#endif	// if ALG_EVAL_ON_HOST

#if ALG_EVAL_ON_HOST
///
/// Host Vector kernel evaluator for types with SIMD acceleration support.
///
template<typename T, typename... Tk>
struct alg_eval_kernel_impl<
	alg_enable_if_t<alg_has_simd_support<T>>, 
	alg_loc_host, T, Tk...>
{
public:
	template<typename TKernel>
	ALG_INLINE static void eval(const TKernel& kernel,
        alg_vector_host_t<T>& x0, const alg_vector_host_t<Tk>&... xk)
	{
		size_t size;
		alg_verify_prop(size, x0, xk...);
		eval_block(size, kernel, x0.data(), xk.data()...);
	}

private:
	template<typename TKernel, size_t TBlockSize = alg_simd_block_max_supported_size_v<T>>
	ALG_INLINE static void eval_block(
		const alg_enable_if_t<TBlockSize != alg_simd_block_min_supported_size_v<T> / 2, size_t> n,
		const TKernel& kernel, T* const x0_mem, const Tk* const... xk_mem)
	{
		const size_t nb = n / TBlockSize;
		const size_t nl = n % TBlockSize;
		if (nb != 0)
		{
			// Evaluating for current block size..
			alg_nonparallel_for(nb, [&](const size_t i)
			{
				alg_store_simd_block<T, TBlockSize> x0_i(x0_mem + i * TBlockSize);
				kernel(x0_i, alg_simd_block<Tk, TBlockSize>(xk_mem + i * TBlockSize)...);
			});
			// Evaluating for half-sized blocks..
			eval_block<TKernel, TBlockSize / 2>(nl, kernel,
				x0_mem + nb * TBlockSize, xk_mem + nb * TBlockSize...);
			return;
		}
		// This memory is smaller than a block size. Lets try with half of the block size..
		eval_block<TKernel, TBlockSize / 2>(nl, kernel,
			x0_mem, xk_mem...);
	}
	template<typename TKernel, size_t TBlockSize>
	ALG_INLINE static void eval_block(
		const alg_enable_if_t<TBlockSize == alg_simd_block_min_supported_size_v<T> / 2, size_t> n,
		const TKernel& kernel, T* const x0_mem, const Tk* const... xk_mem)
	{
		// ...And we've finally reached the leftover.
		alg_nonparallel_for(n, [&](const size_t i) { kernel(x0_mem[i], xk_mem[i]...); });
	}
};	// struct alg_eval_kernel_impl</*has SIMD for T*/, alg_loc_host, T, Tk...>
#endif	// if ALG_EVAL_ON_HOST

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#if ALG_EVAL_ON_HOST
///
/// Vector reduction kernel evaluator for types without acceleration support.
///
template<typename T, typename... Tk>
struct alg_eval_reduce_kernel_impl<
	alg_enable_if_t<!ALG_HAS_SIMD_SUPPORT(T) && !ALG_HAS_SIMD_SUPPORT(alg_related_scalar_t<alg_vector_host_t<T>>)>,
	alg_loc_host, T, Tk...>
{
	using alg_scalar_t = alg_related_scalar_t<alg_vector_host_t<T>>;
public:
	template<typename TKernel, typename TRUKernel>
	ALG_INLINE static alg_scalar_t eval(const TKernel& kernel, const TRUKernel& ru_kernel,
        const alg_vector_host_t<T>& x0, const alg_vector_host_t<Tk>&... xk)
	{
		size_t size;
		alg_verify_prop(size, x0, xk...);
		alg_scalar_t r(kernel(x0(0), xk(0)...));
		alg_nonparallel_for(size, [&](const size_t i) { ru_kernel(r, kernel(x0(i), xk(i)...)); }, 1);
		return r;
	}
};	// struct alg_eval_reduce_kernel_impl</*no SIMD for both T and related scalar*/, alg_loc_host, T, Tk...>
#endif	// if ALG_EVAL_ON_HOST

#if ALG_EVAL_ON_HOST
///
/// Vector reduction kernel evaluator for types, related scalars of whose have SIMD acceleration support.
///
template<typename T, typename... Tk>
struct alg_eval_reduce_kernel_impl<
	alg_enable_if_t<!ALG_HAS_SIMD_SUPPORT(T) && ALG_HAS_SIMD_SUPPORT(alg_related_scalar_t<alg_vector_host_t<T>>)>,
	alg_loc_host, T, Tk...>
{
	using alg_scalar_t = alg_related_scalar_t<alg_vector_host_t<T>>;
public:
	template<typename TKernel, typename TRUKernel>
	ALG_INLINE static alg_scalar_t eval(const TKernel& kernel, const TRUKernel& ru_kernel,
        const alg_vector_host_t<T>& x0, const alg_vector_host_t<Tk>&... xk)
	{
		size_t size;
		alg_verify_prop(size, x0, xk...);
		return eval_block(size, kernel, ru_kernel, x0, xk...);
	}

private:
	template<typename TKernel, typename TRUKernel, size_t TBlockSize = alg_simd_block_max_supported_size_v<T>>
	ALG_INLINE static alg_scalar_t eval_block(
		const alg_enable_if_t<TBlockSize != alg_simd_block_min_supported_size_v<T>, size_t> n,
		const TKernel& kernel, const TRUKernel& ru_kernel,
		const alg_vector_host_t<T>& x0, const alg_vector_host_t<Tk>&... xk)
	{
		if (n >= TBlockSize)
		{
			// This block size is enough to start reduction..
			const alg_simd_block<T, TBlockSize> partial([&](const size_t j) { return kernel(x0(j), xk(j)...); });
			return eval_block_end<TKernel, TRUKernel, TBlockSize>(n - TBlockSize, kernel, ru_kernel, partial,
				TBlockSize, x0, xk...);
		}
		// Trying to start reduction with half block size..
		return eval_block<TKernel, TRUKernel, TBlockSize / 2>(n, kernel, ru_kernel,
			x0, xk...);
	}
	template<typename TKernel, typename TRUKernel, size_t TBlockSize>
	ALG_INLINE static alg_scalar_t eval_block_end(
		const alg_enable_if_t<TBlockSize != alg_simd_block_min_supported_size_v<T>, size_t> n,
		const TKernel& kernel, const TRUKernel& ru_kernel, alg_simd_block<T, TBlockSize> partial,
		const size_t j0, const alg_vector_host_t<T>& x0, const alg_vector_host_t<Tk>&... xk)
	{
		const size_t nb = n / TBlockSize;
		const size_t nl = n % TBlockSize;
		if (nb != 0)
		{
			// Evaluating for current block size..
			alg_nonparallel_for(nb, [&](const size_t i)
			{
				partial = ru_kernel(partial, kernel(
					alg_simd_block<T,  TBlockSize>([&](const size_t j) { return kernel(x0(j0 + j), xk(j0 + j)...); })));
			});
			// Evaluating for half-sized blocks..
			return eval_block_end<TKernel, TRUKernel, TBlockSize / 2>(nl,
				kernel, ru_kernel, ru_kernel(partial, alg_reduce()),
				j0 + nb * TBlockSize, x0, xk...);
		}
		// This memory is smaller than a block size. Lets try with half of the block size..
		return eval_block_end<TKernel, TRUKernel, TBlockSize / 2>(nl,
			kernel, ru_kernel, ru_kernel(partial, alg_reduce()),
			j0, x0, xk...);
	}
	template<typename TKernel, typename TRUKernel, size_t TBlockSize>
	ALG_INLINE static alg_scalar_t eval_block_end(
		const alg_enable_if_t<TBlockSize == alg_simd_block_min_supported_size_v<T>, size_t> n,
		const TKernel& kernel, const TRUKernel& ru_kernel, alg_simd_block<T, TBlockSize> partial,
		const size_t j0, const alg_vector_host_t<T>& x0, const alg_vector_host_t<Tk>&... xk)
	{
		// ...And we've finally reached the minimum block size.
		const size_t nb = n / TBlockSize;
		const size_t nl = n % TBlockSize;
		if (nb != 0)
		{
			// Evaluating for current block size..
			alg_nonparallel_for(nb, [&](const size_t i)
			{
				partial = ru_kernel(partial, kernel(
					alg_simd_block<T, TBlockSize>([&](const size_t j) { return kernel(x0(j0 + j), xk(j0 + j)...); })));
			});
			// Evaluating for half-sized blocks..
			return eval_block_end<TKernel, TRUKernel, TBlockSize / 2>(nl,
				kernel, ru_kernel, ru_kernel(partial, alg_reduce()),
				j0 + nb * TBlockSize, x0, xk...);
		}
		// Evaluating for the leftover..
		const size_t i0 = j0 + nb * TBlockSize;
		alg_scalar_t r = ru_kernel(partial, alg_reduce());
		alg_nonparallel_for(nl, [&](const size_t i)
		{
			ru_kernel(r, kernel(x0(i0 + i), xk(i0 + i)...));
		});
		return r;
	}
};	// struct alg_eval_reduce_kernel_impl</*no SIMD for T, but not for related scalar*/, alg_loc_host, T, Tk...>
#endif	// if ALG_EVAL_ON_HOST

#if ALG_EVAL_ON_HOST
///
/// Vector reduction kernel evaluator for types with SIMD acceleration support.
///
template<typename T, typename... Tk>
struct alg_eval_reduce_kernel_impl<
	alg_enable_if_t<alg_has_simd_support<T>>,
	alg_loc_host, T, Tk...>
{
	using alg_scalar_t = alg_related_scalar_t<alg_vector_host_t<T>>;
public:
	template<typename TKernel, typename TRUKernel>
	ALG_INLINE static alg_scalar_t eval(const TKernel& kernel, const TRUKernel& ru_kernel,
        const alg_vector_host_t<T>& x0, const alg_vector_host_t<Tk>&... xk)
	{
		size_t size;
		alg_verify_prop(size, x0, xk...);
		return eval_block(size, kernel, ru_kernel, x0.data(), xk.data()...);
	}

private:
	template<typename TKernel, typename TRUKernel, size_t TBlockSize = alg_simd_block_max_supported_size_v<T>>
	ALG_INLINE static alg_scalar_t eval_block(
		const alg_enable_if_t<TBlockSize != alg_simd_block_min_supported_size_v<T>, size_t> n,
		const TKernel& kernel, const TRUKernel& ru_kernel,
		T* const x0_mem, const Tk* const... xk_mem)
	{
		if (n >= TBlockSize)
		{
			// This block size is enough to start reduction..
			const alg_simd_block<T, TBlockSize> partial = kernel(
				alg_simd_block<T,  TBlockSize>(x0_mem),
				alg_simd_block<Tk, TBlockSize>(xk_mem)...);
			return eval_block_end<TKernel, TRUKernel, TBlockSize>(n - TBlockSize, kernel, ru_kernel, partial,
			    x0_mem + TBlockSize, (xk_mem + TBlockSize)...);
		}
		// Trying to start reduction with half block size..
		return eval_block<TKernel, TRUKernel, TBlockSize / 2>(n, kernel, ru_kernel,
            x0_mem, xk_mem...);
	}
	template<typename TKernel, typename TRUKernel, size_t TBlockSize>
	ALG_INLINE static alg_scalar_t eval_block_end(
		const alg_enable_if_t<TBlockSize != alg_simd_block_min_supported_size_v<T>, size_t> n,
		const TKernel& kernel, const TRUKernel& ru_kernel, alg_simd_block<T, TBlockSize> partial,
		T* const x0_mem, const Tk* const... xk_mem)
	{
		const size_t nb = n / TBlockSize;
		const size_t nl = n % TBlockSize;
		if (nb != 0)
		{
			// Evaluating for current block size..
			alg_nonparallel_for(nb, [&](const size_t i)
			{
				partial = ru_kernel(partial, kernel(
					alg_simd_block<T,  TBlockSize>(x0_mem + i * TBlockSize),
	                alg_simd_block<Tk, TBlockSize>(xk_mem + i * TBlockSize)...));
			});
			// Evaluating for half-sized blocks..
			return eval_block_end<TKernel, TRUKernel, TBlockSize / 2>(nl,
			    kernel, ru_kernel, ru_kernel(partial, alg_reduce()),
                x0_mem + nb * TBlockSize, xk_mem + nb * TBlockSize...);
		}
		// This memory is smaller than a block size. Lets try with half of the block size..
		return eval_block_end<TKernel, TRUKernel, TBlockSize / 2>(nl,
			kernel, ru_kernel, ru_kernel(partial, alg_reduce()),
			x0_mem, xk_mem...);
	}
	template<typename TKernel, typename TRUKernel, size_t TBlockSize>
	ALG_INLINE static alg_scalar_t eval_block_end(
		const alg_enable_if_t<TBlockSize == alg_simd_block_min_supported_size_v<T>, size_t> n,
		const TKernel& kernel, const TRUKernel& ru_kernel, alg_simd_block<T, TBlockSize> partial,
		T* const x0_mem, const Tk* const... xk_mem)
	{
		// ...And we've finally reached the minimum block size.
		const size_t nb = n / TBlockSize;
		const size_t nl = n % TBlockSize;
		if (nb != 0)
		{
			// Evaluating for minimum block size..
			alg_nonparallel_for(nb, [&](const size_t i)
			{
				partial = ru_kernel(partial, kernel(
					alg_simd_block<T,  TBlockSize>(x0_mem + i * TBlockSize),
		            alg_simd_block<Tk, TBlockSize>(xk_mem + i * TBlockSize)...));
			});
		}
		// Evaluating for the leftover..
		alg_scalar_t r = ru_kernel(partial, alg_reduce());
		alg_nonparallel_for(nl, [&](const size_t i)
		{
			ru_kernel(r, kernel(x0_mem[nb * TBlockSize + i], xk_mem[nb * TBlockSize + i]...));
		});
		return r;
	}
};	// struct alg_eval_reduce_kernel_impl</*has SIMD for T*/, alg_loc_host, T, Tk...>
#endif	// if ALG_EVAL_ON_HOST

// ------------------------------------------------------------------------------------ //
// -----                      Device Vector kernel evaluators.                    ----- //
// ------------------------------------------------------------------------------------ //

#if ALG_EVAL_ON_DEVICE
///
/// Device Vector kernel evaluator.
///
template<typename T, typename... Tk>
struct alg_eval_kernel_impl<void, alg_loc_device, T, Tk...> final
{
public:
	template<typename TKernel>
	ALG_INLINE static void eval(const TKernel& kernel,
        alg_vector_device_t<T>& x0, const alg_vector_host_t<Tk>&... xk)
	{
		size_t size;
		alg_verify_prop(size, x0, xk...);
		alg_gpu_sized_buffer<T> x0_device_data(x0.device_data(), size);
		kernel(x0_device_data, alg_gpu_buffer<T>(xk.device_data())...);
	}
};	// struct alg_eval_kernel_impl<alg_loc_device, T, Tk...>
#endif	// if ALG_EVAL_ON_DEVICE

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#if ALG_EVAL_ON_DEVICE
///
/// Device Vector reduction kernel evaluator.
///
template<typename T, typename... Tk>
struct alg_eval_reduce_kernel_impl<void, alg_loc_device, T, Tk...> final
{
	using alg_scalar_t = alg_related_scalar_t<alg_vector_device_t<T>>;
public:
	template<typename TKernel, typename TRUKernel>
	ALG_INLINE static alg_scalar_t eval(const TKernel& kernel, const TRUKernel&,
        const alg_vector_host_t<T>& x0, const alg_vector_host_t<Tk>&... xk)
	{
		size_t size;
		alg_verify_prop(size, x0, xk...);
		alg_gpu_sized_buffer<T> x0_device_data(x0.device_data(), size);
		return kernel(x0_device_data, alg_gpu_buffer<T>(xk.device_data())...);
	}
};  // struct alg_eval_reduce_kernel_impl<alg_loc_host, T, Tk...>
#endif	// if ALG_EVAL_ON_DEVICE

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/// @f$ kernel(x_i^k), i = 0..N-1 @f$
template<alg_loc_t TLoc, typename T, typename... Tk>
ALG_INLINE static auto alg_eval_kernel(
	alg_vector<T, TLoc>& x0, const alg_vector<Tk, TLoc>&... xk)
{
	return [&](const auto& kernel)
	{
		alg_eval_kernel_impl_t<TLoc, T, Tk...>::eval(kernel, x0, xk...);
	};
}

/// @f$ r_0 := kernel(x_0); r_i := kernelUpdate(r_{i-1}, kernel(x_i)), i = 1..N-1; r := r_{N-1} @f$
template<alg_loc_t TLoc, typename T, typename... Tk>
ALG_INLINE static auto alg_reduce_kernel(
	const alg_vector<T, TLoc>& x0, const alg_vector<Tk, TLoc>&... xk)
{
	return [&](const auto& kernel, const auto& ru_kernel) -> alg_related_scalar_t<alg_vector<T>>
	{
		return alg_eval_reduce_kernel_impl_t<TLoc, T, Tk...>::eval(kernel, ru_kernel,
			x0, xk...);
	};
}

// ------------------------------------------------------------------------------------ //
// -----                        Vector arithmetic functions.                      ----- //
// ------------------------------------------------------------------------------------ //

/// @defgroup alg_vector_arithmetic_functions Arithmetic functions
/// @{

/// @f$ \vec x := \alpha \vec y + \beta \vec z @f$
template<typename T, alg_loc_t TLoc>
static void alg_add(alg_vector<T, TLoc>& x,
    const alg_related_scalar_t<alg_vector<T, TLoc>> alpha, const alg_vector<T, TLoc>& y,
    const alg_related_scalar_t<alg_vector<T, TLoc>> beta,  const alg_vector<T, TLoc>& z)
{
	alg_verify_coef(alpha, beta);
	alg_eval_kernel(x, y, z)([&](auto& xi, const auto& yi, const auto& zi)
	{
		alg_add(xi, alpha, yi, beta, zi);
	});
}

/// @f$ \vec x += \alpha \vec y @f$
template<typename T>
static void alg_add_assign(alg_vector<T>& x,
    const alg_related_scalar_t<alg_vector<T>> alpha, const alg_vector<T>& y)
{
	alg_verify_coef(alpha);
	alg_eval_kernel(x, y)([&](auto& xi, const auto& yi)
	{
		alg_add_assign(xi, alpha, yi);
	});
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/// @f$ \vec x := \alpha \vec y @f$
template<typename T>
static void alg_mul(alg_vector<T>& x,
	const alg_related_scalar_t<alg_vector<T>> alpha, const alg_vector<T>& y)
{
	alg_verify_coef(alpha);
	alg_eval_kernel(x, y)([&](auto& xi, const auto& yi)
	{
		alg_mul(xi, alpha, yi);
	});
}

/// @f$ \vec x := \alpha \vec x @f$
template<typename T>
static void alg_mul_assign(alg_vector<T>& x,
   const alg_related_scalar_t<alg_vector<T>> alpha)
{
	alg_verify_coef(alpha);
	alg_eval_kernel(x)([&](auto& xi)
	{
		alg_mul_assign(xi, alpha);
	});
}

/// @}

// ------------------------------------------------------------------------------------ //
// -----                        Vector reduction operations.                      ----- //
// ------------------------------------------------------------------------------------ //

/// @defgroup alg_vector_reduction_operations Reduction operations
/// @{

/// @f$ r := \vec x \dot \vec y @f$
template<typename T>
static alg_related_scalar_t<alg_vector<T>> alg_dot(const alg_vector<T>& x, const alg_vector<T>& y)
{
	size_t size;
	alg_verify_prop(size, x, y);
	alg_related_scalar_t<alg_vector<T>> dot(0.0);
	alg_eval_kernel(const_cast<alg_vector<T>&>(x), y)([&](const auto& xi, const auto& yi)
	{
		dot += alg_dot(xi, yi);
	});

	return dot;
}

/// @f$ r := \vec x \dot \vec y @f$
template<typename T>
static alg_related_scalar_t<alg_vector<T>> alg_dot2(const alg_vector<T>& x, const alg_vector<T>& y)
{
	size_t size;
	alg_verify_prop(size, x, y);
	alg_related_scalar_t<alg_vector<T>> dot(0.0);
	alg_reduce_kernel(const_cast<alg_vector<T>&>(x), y)(
		[&](const auto& xi, const auto& yi) { dot += alg_dot(xi, yi); });

	return dot;
}

// ------------------------------------------------------------------------------------ //

/// @f$ r := ||\vec x||_{l_1} @f$
template<typename T>
static alg_related_scalar_t<alg_vector<T>> alg_norm_l1(const alg_vector<T>& x)
{
	size_t size;
	alg_verify_prop(size, x);
	alg_related_scalar_t<alg_vector<T>> norm(0.0);
	alg_eval_kernel(const_cast<alg_vector<T>&>(x))([&](const auto& xi)
	{
		norm += alg_norm_l1(xi);
	});

	return norm;
}

/// @f$ r := ||\vec x||_{l_2}^2 @f$
template<typename T>
static alg_related_scalar_t<alg_vector<T>> alg_norm_l2_sq(const alg_vector<T>& x)
{
	const alg_related_scalar_t<alg_vector<T>> norm_sq = alg_dot(x, x);
	return norm_sq;
}

/// @f$ r := ||\vec x||_{l_2} @f$
template<typename T>
static alg_related_scalar_t<alg_vector<T>> alg_norm_l2(const alg_vector<T>& x)
{
	const alg_related_scalar_t<alg_vector<T>> norm_sq = alg_norm_l2_sq(x);
	return std::sqrt(norm_sq);
}

/// @f$ r := ||\vec x||_{l_\infty} @f$
template<typename T>
static alg_related_scalar_t<alg_vector<T>> alg_norm_linf(const alg_vector<T>& x)
{
	size_t size;
	alg_verify_prop(size, x);
	alg_related_scalar_t<alg_vector<T>> norm(-1.0);
	alg_eval_kernel(const_cast<alg_vector<T>&>(x))([&](const auto& xi)
	{
		norm = std::max(norm, alg_norm_linf(xi));
	});

	return norm;
}

// ------------------------------------------------------------------------------------ //

/// @f$ r := max|x_i - y_i| \le \varepsilon @f$
template<typename T>
static bool alg_near_eql(
	const alg_vector<T>& x, const alg_vector<T>& y,
	const alg_related_scalar_t<alg_vector<T>> eps = alg_related_scalar_t<alg_vector<T>>(ALG_EPS))
{
	/// @todo Rewrite as kernel
	size_t size;
	alg_verify_prop(size, x, y);
	for (size_t i = 0; i < size; ++i)
	{
		if (!alg_near_eql(x(i), y(i), eps))
		{
			return false;
		}
	}
	return true;
}

/// @}

/// @}
