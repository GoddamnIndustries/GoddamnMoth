// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "../Scalar.h"

#include <immintrin.h>
#define ALG_VECTORCALL //__vectorcall

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++                              SSE extensions.                             +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

// ------------------------------------------------------------------------------------ //
// -----                           Reduce add operations.                         ----- //
// ------------------------------------------------------------------------------------ //

/// @f$ r := \sum_{i=1}^{n} x_i @f$
/// @{
ALG_INLINE static alg_float32_t ALG_VECTORCALL _mm_reduce_add_ps(const __m128 x)
{
	// http://stackoverflow.com/a/35270026
	const __m128 x_shuf = _mm_movehdup_ps(x);
	const __m128 x_shuf__add__x = _mm_add_ps(x, x_shuf);
	const __m128 x_shuf___movehl___x_shuf__add__x = _mm_movehl_ps(x_shuf, x_shuf__add__x);
	return _mm_cvtss_f32(_mm_add_ss(x_shuf__add__x, x_shuf___movehl___x_shuf__add__x));
}

#if ALG_AVX256_SUPPORT
ALG_INLINE static alg_float32_t ALG_VECTORCALL _mm256_reduce_add_ps(const __m256 x)
{
	// http://stackoverflow.com/a/35270026
	const __m256 x_hadd_x = _mm256_hadd_ps(x, x);
	const __m256 x_hadd_x__hadd__x_hadd_x = _mm256_hadd_ps(x_hadd_x, x_hadd_x);
	const __m128 x_hadd_x__hadd__x_hadd_x__lo = _mm256_castps256_ps128(x_hadd_x__hadd__x_hadd_x);
	const __m128 x_hadd_x__hadd__x_hadd_x__hi = _mm256_extractf128_ps(x_hadd_x__hadd__x_hadd_x, 1);
	return _mm_cvtss_f32(_mm_add_ss(x_hadd_x__hadd__x_hadd_x__lo, x_hadd_x__hadd__x_hadd_x__hi));
}
#endif  // if ALG_AVX256_SUPPORT

#if ALG_AVX512_SUPPORT
ALG_INLINE static float32_t ALG_VECTORCALL _mm512_reduce_add_ps(const __m512 x)
{
	// http://stackoverflow.com/a/35270026
	const __m256 x__lo = _mm512_castps512_ps256(x);
	const __m256 x__hi = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(x), 1));
	return _mm256_reduce_add_ps(_mm256_add_ps(x__lo, x__hi));
}
#endif  // if ALG_AVX512_SUPPORT
/// @}

// ------------------------------------------------------------------------------------ //

/// @f$ r := \sum_{i=1}^{n} x_i @f$
/// @{
ALG_INLINE static alg_float64_t ALG_VECTORCALL _mm_reduce_add_pd(const __m128d x)
{
	return _mm_cvtsd_f64(_mm_hadd_pd(x, x));
}

#if ALG_AVX256_SUPPORT
ALG_INLINE static alg_float64_t ALG_VECTORCALL _mm256_reduce_add_pd(const __m256d x)
{
	// http://stackoverflow.com/a/35270026
	const __m256d x_hadd_x = _mm256_hadd_pd(x, x);
	const __m128d x_hadd_x__lo = _mm256_castpd256_pd128(x_hadd_x);
	const __m128d x_hadd_x__hi = _mm256_extractf128_pd(x_hadd_x, 1);
	return _mm_cvtsd_f64(_mm_add_sd(x_hadd_x__lo, x_hadd_x__hi));
}
#endif  // if ALG_AVX256_SUPPORT

#if ALG_AVX512_SUPPORT
ALG_INLINE static float64_t ALG_VECTORCALL _mm512_reduce_add_pd(__m512d const x)
{
	// http://stackoverflow.com/a/35270026
	const __m256d x__lo = _mm512_castpd512_pd256(x);
	const __m256d x__hi = _mm512_extractf64x4_pd(x, 1);
	return _mm256_reduce_add_pd(_mm256_add_pd(x__lo, x__hi));
}
#endif  // if ALG_AVX512_SUPPORT
/// @}

// ------------------------------------------------------------------------------------ //
// -----                              Abs operations.                             ----- //
// ------------------------------------------------------------------------------------ //

/// @f$ \vec r := |\vec x| @f$
/// @{
ALG_INLINE static __m128 ALG_VECTORCALL _mm_abs_ps(__m128 const x)
{
	// http://fastcpp.blogspot.ru/2011/03/changing-sign-of-float-values-using-sse.html
	static const __m128 sign_mask_128 = _mm_castsi128_ps(_mm_set1_epi32(~(1 << 31)));
	return _mm_and_ps(sign_mask_128, x);
}

#if ALG_AVX256_SUPPORT
ALG_INLINE static __m256 ALG_VECTORCALL _mm256_abs_ps(__m256 const x)
{
	// http://fastcpp.blogspot.ru/2011/03/changing-sign-of-float-values-using-sse.html
	static const __m256 sign_mask_256 = _mm256_castsi256_ps(_mm256_set1_epi64x(~(1 << 31)));
	return _mm256_and_ps(sign_mask_256, x);
}
#endif  // if ALG_AVX256_SUPPORT
/// @}

// ------------------------------------------------------------------------------------ //

/// @f$ \vec r := |\vec x| @f$
/// @{
ALG_INLINE static __m128d ALG_VECTORCALL _mm_abs_pd(__m128d const x)
{
	// http://fastcpp.blogspot.ru/2011/03/changing-sign-of-float-values-using-sse.html
	static const __m128d sign_mask_128d = _mm_castsi128_pd(_mm_set1_epi64x(~(1ull << 63)));
	return _mm_and_pd(sign_mask_128d, x);
}

#if ALG_AVX256_SUPPORT
ALG_INLINE static __m256d ALG_VECTORCALL _mm256_abs_pd(__m256d const x)
{
	// http://fastcpp.blogspot.ru/2011/03/changing-sign-of-float-values-using-sse.html
	static const __m256d sign_mask_256d = _mm256_castsi256_pd(_mm256_set1_epi64x(~(1ull << 63)));
	return _mm256_and_pd(sign_mask_256d, x);
}
#endif  // if ALG_AVX256_SUPPORT
/// @}
