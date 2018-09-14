// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Utils/Assert.h"

#include <cmath>
#include <cassert>
#include <type_traits>
//#include <atomic>

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++                             Generic scalars.                             +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

enum alg_loc_t
{
	alg_loc_host,
	alg_loc_device,
};  // enum alg_loc_t

#if ALG_EVAL_ON_DEVICE 
template<typename T>
class alg_device_scalar final
{
public:
	alg_gpu_buffer<T> m_buffer;
	alg_gpu_event m_event;
public:
	/// @f$ \hat x := x @f$
	ALG_INLINE alg_device_scalar(const T scalar = -122.0)
			: m_buffer(1, &scalar)
	{
	}
public:
	ALG_INLINE operator T() const
	{
		T r;
		m_buffer.read_sync(&r, 1);
		return r;
	}
};  // class alg_device_scalar
#endif  // if ALG_EVAL_ON_DEVICE

///
/// Basic scalar floating-point types.
///
using alg_float32_t = float;
using alg_float64_t = double;

template<typename, alg_loc_t>
struct alg_scalar;
template<typename T>
struct alg_scalar<T, alg_loc_host> { using type = T; };
//template<typename T>
//struct alg_scalar<T, alg_loc_device> { using type = alg_device_scalar<T>; };
template<typename T, alg_loc_t TLoc>
using alg_scalar_t = typename alg_scalar<T, TLoc>::type;

///
/// Related scalar type.
///
template<typename>
struct alg_related_scalar;
template<>
struct alg_related_scalar<alg_float32_t> { using type = alg_float32_t; };
template<>
struct alg_related_scalar<alg_float64_t> { using type = alg_float64_t; };

template<typename U>
using alg_related_scalar_t = typename alg_related_scalar<U>::type;

///
/// Related inverse type.
///
template<typename T>
struct alg_related_inverse { using type = T; };
template<>
struct alg_related_inverse<alg_float32_t> { using type = alg_float32_t; };
template<>
struct alg_related_inverse<alg_float64_t> { using type = alg_float64_t; };

template<typename U>
using alg_related_inverse_t = typename alg_related_inverse<U>::type;

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

//template<alg_loc_t TLoc, typename T, typename... Tk>
//ALG_INLINE static auto alg_eval_kernel(
//	alg_scalar_t<T, TLoc>& x0, const alg_scalar_t<Tk, TLoc>&... xk)
//{
//	return [&](const auto& kernel)
//	{
//		alg_gpu_sized_buffer<T> x0_device_data(x0.m_buffer(), 1);
//		return kernel(x0_device_data, alg_gpu_buffer<T>(xk.m_buffer)...);
//	};
//}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/// @f$ \hat x := \alpha \hat y + \beta \hat z @f$
//template<typename T, alg_loc_t TLoc>
//static void alg_add(alg_scalar_t<T, TLoc>& x,
//    const T alpha, const alg_scalar_t<T, TLoc> y,
//    const T beta,  const alg_scalar_t<T, TLoc> z)
//{
//	alg_eval_kernel(x, y, z)([&]()
//	{
//		alg_add(x, alpha);
//	});
//}
static void alg_add(alg_float64_t& x, const alg_float64_t alpha, const alg_float64_t y, const alg_float64_t beta, const alg_float64_t z)
{
	alg_verify_coef(alpha, beta);
	x = fma(alpha, y, beta * z);
	alg_verify_coef(x);
}
static void alg_add(alg_float32_t& x, const alg_float32_t alpha, const alg_float32_t y, const alg_float32_t beta, const alg_float32_t z)
{
	alg_verify_coef(alpha, beta);
	x = fmaf(alpha, y, beta * z);
	alg_verify_coef(x);
}

// ------------------------------------------------------------------------------------ //

/// x += alpha * y
static void alg_add_assign(alg_float64_t& x, const alg_float64_t alpha, const alg_float64_t y)
{
	alg_verify_coef(alpha);
	x = fma(alpha, y, x);
	alg_verify_coef(x);
}
static void alg_add_assign(alg_float32_t& x, const alg_float32_t alpha, const alg_float32_t y)
{
	alg_verify_coef(alpha);
	x = fmaf(alpha, y, x);
	alg_verify_coef(x);
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/// x := alpha * y
static void alg_mul(alg_float64_t& x, const alg_float64_t alpha, const alg_float64_t y)
{
	alg_verify_coef(alpha);
	x = alpha * y;
	alg_verify_coef(x);
}
static void alg_mul(alg_float32_t& x, const alg_float32_t alpha, const alg_float32_t y)
{
	alg_verify_coef(alpha);
	x = alpha * y;
	alg_verify_coef(x);
}

/// x := alpha * y * z
static void alg_mul(alg_float64_t& x, const alg_float64_t alpha, const alg_float64_t y, const alg_float64_t z)
{
	alg_verify_coef(alpha);
	x = alpha * y * z;
	alg_verify_coef(x);
}
static void alg_mul(alg_float32_t& x, const alg_float32_t alpha, const alg_float32_t y, const alg_float32_t z)
{
	alg_verify_coef(alpha);
	x = alpha * y * z;
	alg_verify_coef(x);
}

/// x *= alpha
static void alg_mul_assign(alg_float64_t& x, const alg_float64_t alpha)
{
	alg_verify_coef(alpha);
	x *= alpha;
	alg_verify_coef(x);
}
static void alg_mul_assign(alg_float32_t& x, const alg_float32_t alpha)
{
	alg_verify_coef(alpha);
	x *= alpha;
	alg_verify_coef(x);
}

/// x *= alpha * y
static void alg_mul_assign(alg_float64_t& x, const alg_float64_t alpha, const alg_float64_t y)
{
	alg_verify_coef(alpha);
	x *= alpha * y;
	alg_verify_coef(x);
}
static void alg_mul_assign(alg_float32_t& x, const alg_float32_t alpha, const alg_float32_t y)
{
	alg_verify_coef(alpha);
	x *= alpha * y;
	alg_verify_coef(x);
}

// ------------------------------------------------------------------------------------ //

/// x := alpha * y * z + beta * d
static void alg_mul_add(alg_float64_t& x, const alg_float64_t alpha, const alg_float64_t y, const alg_float64_t z, const alg_float64_t beta, const alg_float64_t d)
{
	alg_verify_coef(alpha, beta);
	x = fma(alpha, y * z, beta * d);
	alg_verify_coef(x);
}
static void alg_mul_add(alg_float32_t& x, const alg_float32_t alpha, const alg_float32_t y, const alg_float32_t z, const alg_float32_t beta, const alg_float32_t d)
{
	alg_verify_coef(alpha, beta);
	x = fmaf(alpha, y * z, beta * d);
	alg_verify_coef(x);
}

/// x += alpha * y * z
static void alg_mul_add_assign(alg_float64_t& x, const alg_float64_t alpha, const alg_float64_t y, const alg_float64_t z)
{
	alg_verify_coef(alpha);
	x = fma(alpha, y * z, x);
	alg_verify_coef(x);
}
static void alg_mul_add_assign(alg_float32_t& x, const alg_float32_t alpha, const alg_float32_t y, const alg_float32_t z)
{
	alg_verify_coef(alpha);
	x = fmaf(alpha, y * z, x);
	alg_verify_coef(x);
}

// ------------------------------------------------------------------------------------ //

/// r := (x, y).
static alg_float64_t alg_dot(const alg_float64_t x, const alg_float64_t y)
{
	return x * y;
}
static alg_float32_t alg_dot(const alg_float32_t x, const alg_float32_t y)
{
	return x * y;
}

/// r := ||x||_L1
static alg_float64_t alg_norm_l1(const alg_float64_t x)
{
	return std::fabs(x);
}
static alg_float32_t alg_norm_l1(const alg_float32_t x)
{
	return std::fabs(x);
}

/// r := ||x||_L1
static alg_float64_t alg_norm_linf(const alg_float64_t x)
{
	return std::fabs(x);
}
static alg_float32_t alg_norm_linf(const alg_float32_t x)
{
	return std::fabs(x);
}

#define ALG_EPS 0.0000000001

/// r := |x - y| <= eps
static bool alg_near_eql(const alg_float64_t x, const alg_float64_t y, const alg_float64_t eps = alg_float64_t(ALG_EPS))
{
	return fabs(x - y) <= eps;
}

// ------------------------------------------------------------------------------------ //

/// y := x^-1
static void alg_inverse(alg_float64_t& y, const alg_float64_t x)
{
    y = 1.0 / x;
	ALG_DEBUG_ASSERT(isfinite(y));
}
static void alg_inverse(alg_float32_t& y, const alg_float32_t x)
{
	y = 1.0f / x;
	ALG_DEBUG_ASSERT(isfinite(y));
}
