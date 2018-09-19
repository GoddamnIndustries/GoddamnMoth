// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once
#include "SimdBlock.h"
#include "SimdSSE.h"

/// @addtogroup alg_simd_block
/// @{

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++          SSE-accelerated operation for double-precision vectors.         +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

///
/// Generic SIMD-accelerated column vector.
///
/// @f$ \vec x = [ x_0, x_1, ..., x_{N-1} ]^T @f$
/// @see https://en.wikipedia.org/wiki/Vector_space
/// @{
///
template<>
class alg_simd_block<alg_float64_t, 2>
{
public:
    __m128d data;

public:
	ALG_INLINE explicit alg_simd_block(const alg_float64_t* const ptr)
		: data(_mm_loadu_pd(ptr)) {}
public:
    ALG_INLINE void store(alg_float64_t* const ptr) const
    {
        _mm_storeu_pd(ptr, data);
    }
};  // class alg_simd_block<alg_float64_t, 2>

#if ALG_AVX256_SUPPORT
template<>
class alg_simd_block<alg_float64_t, 4>
{
public:
    __m256d data;

public:
	ALG_INLINE explicit alg_simd_block(const alg_float64_t* const ptr)
		: data(_mm256_loadu_pd(ptr)) {}
public:
    ALG_INLINE void store(alg_float64_t* const ptr) const
    {
        _mm256_storeu_pd(ptr, data);
    }
};  // class alg_simd_block<alg_float64_t, 4>
#endif  // if ALG_AVX256_SUPPORT

#if ALG_AVX512_SUPPORT
template<>
class alg_simd_block<alg_float64_t, 8> final
{
public:
    __m512d data;

public:
	ALG_INLINE explicit alg_simd_block(const alg_float64_t* const ptr)
		: data(_mm512_loadu_pd(ptr)) {}
public:
    ALG_INLINE void store(alg_float64_t* const ptr) const
    {
        _mm512_storeu_pd(ptr, data);
    }
};  // class alg_simd_block<alg_float64_t, 8>
#endif  // if ALG_AVX512_SUPPORT
/// @}

///
/// Minimum machine size for specific types of the SIMD vector.
///
template<>
struct alg_simd_block_min_supported_size<alg_float64_t>
{
	enum { value = 2 };
};	// struct alg_simd_block_min_supported_size<alg_float64_t>

///
/// Maximum machine size for specific types of the SIMD vector.
///
template<>
struct alg_simd_block_max_supported_size<alg_float64_t>
{
#if ALG_AVX512_SUPPORT
	enum { value = 8 };
#endif	//if ALG_AVX512_SUPPORT
#if ALG_AVX256_SUPPORT && !ALG_AVX512_SUPPORT
	enum { value = 4 };
#endif	//if ALG_AVX256_SUPPORT && !ALG_AVX512_SUPPORT
#if !ALG_AVX256_SUPPORT && !ALG_AVX512_SUPPORT
	enum { value = 2 };
#endif	// if !ALG_AVX256_SUPPORT && !ALG_AVX512_SUPPORT
};	// struct alg_simd_block_max_size<alg_float64_t>

// ------------------------------------------------------------------------------------ //
// -----                         Vector-Vector operations.                        ----- //
// ------------------------------------------------------------------------------------ //

/// @f$ \vec x := \alpha \vec y + \beta \vec z @f$
/// @todo https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=_mm_fma&expand=2381,2409,2389
/// @{
ALG_INLINE static void ALG_VECTORCALL alg_add(alg_simd_block<alg_float64_t, 2>& x,
    const alg_float64_t alpha, const alg_simd_block<alg_float64_t, 2>& y,
    const alg_float64_t beta,  const alg_simd_block<alg_float64_t, 2>& z)
{
    x.data = _mm_add_pd(
        _mm_mul_pd(_mm_set1_pd(alpha), y.data),
        _mm_mul_pd(_mm_set1_pd(beta),  z.data));
}

#if ALG_AVX256_SUPPORT
ALG_INLINE static void ALG_VECTORCALL alg_add(alg_simd_block<alg_float64_t, 4>& x,
    const alg_float64_t alpha, const alg_simd_block<alg_float64_t, 4>& y,
    const alg_float64_t beta,  const alg_simd_block<alg_float64_t, 4>& z)
{
    x.data = _mm256_add_pd(
        _mm256_mul_pd(_mm256_set1_pd(alpha), y.data),
        _mm256_mul_pd(_mm256_set1_pd(beta),  z.data));
}
#endif  // if ALG_AVX256_SUPPORT

#if ALG_AVX512_SUPPORT
ALG_INLINE static void ALG_VECTORCALL alg_add(alg_simd_block<alg_float64_t, 8>& x,
    const alg_float64_t alpha, const alg_simd_block<alg_float64_t, 8>& y,
    const alg_float64_t beta,  const alg_simd_block<alg_float64_t, 8>& z)
{
    x.data = _mm512_add_pd(
        _mm512_mul_pd(_mm512_set1_pd(alpha),  y.data),
        _mm512_mul_pd(_mm512_set1_pd(beta), z.data));
}
#endif  // if ALG_AVX512_SUPPORT
/// @}

// ------------------------------------------------------------------------------------ //

/// @f$ \vec x += \alpha \vec y @f$
/// @todo https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=_mm_fma&expand=2381,2409,2389
/// @{
ALG_INLINE static void ALG_VECTORCALL alg_add_assign(alg_simd_block<alg_float64_t, 2>& x,
    const alg_float64_t alpha, const alg_simd_block<alg_float64_t, 2>& y)
{
    x.data = _mm_add_pd(x.data,
        _mm_mul_pd(_mm_set1_pd(alpha), y.data));
}

#if ALG_AVX256_SUPPORT
ALG_INLINE static void ALG_VECTORCALL alg_add_assign(alg_simd_block<alg_float64_t, 4>& x,
    const alg_float64_t alpha, const alg_simd_block<alg_float64_t, 4>& y)
{
    x.data = _mm256_add_pd(x.data,
        _mm256_mul_pd(_mm256_set1_pd(alpha), y.data));
}
#endif  // if ALG_AVX256_SUPPORT

#if ALG_AVX512_SUPPORT
ALG_INLINE static void ALG_VECTORCALL alg_add_assign(alg_simd_block<alg_float64_t, 8>& x,
    const alg_float64_t alpha, const alg_simd_block<alg_float64_t, 8>& y)
{
    x.data = _mm512_add_pd(x.data,
        _mm512_mul_pd(_mm512_set1_pd(alpha), y.data));
}
#endif  // if ALG_AVX512_SUPPORT
/// @}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/// @f$ \vec x := \alpha \vec y @f$
/// @{
ALG_INLINE static void ALG_VECTORCALL alg_mul(alg_simd_block<alg_float64_t, 2>& x,
    const alg_float64_t alpha, const alg_simd_block<alg_float64_t, 2>& y)
{
    x.data = _mm_mul_pd(_mm_set1_pd(alpha), y.data);
}

#if ALG_AVX256_SUPPORT
ALG_INLINE static void ALG_VECTORCALL alg_mul(alg_simd_block<alg_float64_t, 4>& x,
    const alg_float64_t alpha, const alg_simd_block<alg_float64_t, 4>& y)
{
    x.data = _mm256_mul_pd(_mm256_set1_pd(alpha), y.data);
}
#endif  // if ALG_AVX256_SUPPORT

#if ALG_AVX512_SUPPORT
ALG_INLINE static void ALG_VECTORCALL alg_mul(alg_simd_block<alg_float64_t, 8>& x,
    const alg_float64_t alpha, const alg_simd_block<alg_float64_t, 8>& y)
{
    x.data = _mm512_mul_pd(_mm512_set1_pd(alpha), y.data);
}
#endif  // if ALG_AVX512_SUPPORT
/// @}

// ------------------------------------------------------------------------------------ //

/// @f$ \vec x := \alpha \vec x @f$
/// @{
ALG_INLINE static void ALG_VECTORCALL alg_mul_assign(alg_simd_block<alg_float64_t, 2>& x,
    const alg_float64_t alpha)
{
    // This seems to be optimal.
    alg_mul(x, alpha, x);
}

#if ALG_AVX256_SUPPORT
ALG_INLINE static void ALG_VECTORCALL alg_mul_assign(alg_simd_block<alg_float64_t, 4>& x,
	const alg_float64_t alpha)
{
	// This seems to be optimal.
	alg_mul(x, alpha, x);
}
#endif  // if ALG_AVX256_SUPPORT

#if ALG_AVX512_SUPPORT
ALG_INLINE static void ALG_VECTORCALL alg_mul_assign(alg_simd_block<alg_float64_t, 8>& x,
	const alg_float64_t alpha)
{
	// This seems to be optimal.
	alg_mul(x, alpha, x);
}
#endif  // if ALG_AVX512_SUPPORT
/// @}

// ------------------------------------------------------------------------------------ //
// -----                        Vector reduction operations.                      ----- //
// ------------------------------------------------------------------------------------ //

/// @f$ r := \vec x \dot \vec y @f$
/// @{
ALG_INLINE static alg_float64_t ALG_VECTORCALL alg_dot(const alg_simd_block<alg_float64_t, 2>& x,
	const alg_simd_block<alg_float64_t, 2>& y)
{
	return _mm_cvtsd_f64(_mm_dp_pd(x.data, y.data, 0x31 /* 0011 0001 */));
}

#if ALG_AVX256_SUPPORT
ALG_INLINE static alg_float64_t ALG_VECTORCALL alg_dot(const alg_simd_block<alg_float64_t, 4>& x,
	const alg_simd_block<alg_float64_t, 4>& y)
{
	return _mm256_reduce_add_pd(_mm256_mul_pd(x.data, y.data));
}
#endif  // if ALG_AVX256_SUPPORT

#if ALG_AVX512_SUPPORT
ALG_INLINE static alg_float64_t ALG_VECTORCALL alg_dot(const alg_simd_block<alg_float64_t, 8>& x,
	const alg_simd_block<alg_float64_t, 8>& y)
{
	return _mm512_reduce_add_pd(_mm512_mul_pd(x.data, y.data));
}
#endif  // if ALG_AVX512_SUPPORT
/// @}

// ------------------------------------------------------------------------------------ //

/// @f$ r := ||\vec x||_{l_1} @f$
/// @{
ALG_INLINE static alg_float64_t ALG_VECTORCALL alg_norm_l1(const alg_simd_block<alg_float64_t, 2>& x)
{
	return _mm_reduce_add_pd(_mm_abs_pd(x.data));
}

#if ALG_AVX256_SUPPORT
ALG_INLINE static alg_float64_t ALG_VECTORCALL alg_norm_l1(const alg_simd_block<alg_float64_t, 4>& x)
{
	return _mm256_reduce_add_pd(_mm256_abs_pd(x.data));
}
#endif  // if ALG_AVX256_SUPPORT

#if ALG_AVX512_SUPPORT
ALG_INLINE static alg_float64_t ALG_VECTORCALL alg_norm_l1(const alg_simd_block<alg_float64_t, 8>& x)
{
	return _mm512_reduce_add_pd(_mm512_abs_pd(x.data));
}
#endif  // if ALG_AVX512_SUPPORT
/// @}

// ------------------------------------------------------------------------------------ //

/// @f$ r := ||\vec x||_{l_\infty} @f$
/// @{
ALG_INLINE static alg_float64_t ALG_VECTORCALL alg_norm_linf(const alg_simd_block<alg_float64_t, 2>& x)
{
	return _mm_reduce_add_pd(_mm_abs_pd(x.data));
}

#if ALG_AVX256_SUPPORT
ALG_INLINE static alg_float64_t ALG_VECTORCALL alg_norm_linf(const alg_simd_block<alg_float64_t, 4>& x)
{
	return _mm256_reduce_add_pd(_mm256_abs_pd(x.data));
}
#endif  // if ALG_AVX256_SUPPORT

#if ALG_AVX512_SUPPORT
ALG_INLINE static alg_float64_t ALG_VECTORCALL alg_norm_linf(const alg_simd_block<alg_float64_t, 8>& x)
{
	return _mm512_reduce_add_pd(_mm512_abs_pd(x.data));
}
#endif  // if ALG_AVX512_SUPPORT
/// @}

/// @}
