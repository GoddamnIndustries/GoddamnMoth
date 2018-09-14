// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++                        SIMD-accelerated blocks.                          +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

#pragma once
#include "../Scalar.h"
#include "../Utils/Parallel.h"

/// @defgroup alg_simd_block SIMD-accelerated blocks
/// @{

struct alg_reduce {};

///
/// Generic SIMD-accelerated column vector.
///
/// @f$ \vec x = [ x_0, x_1, ..., x_{N-1} ]^T @f$
/// @see https://en.wikipedia.org/wiki/Vector_space
/// @{
///
template<typename T, size_t N>
class alg_simd_block;
template<typename T, size_t N>
class alg_store_simd_block final : public alg_simd_block<T, N>
{
private:
	T* const m_ptr;
public:
	ALG_INLINE explicit alg_store_simd_block(T* const ptr)
		: alg_simd_block<T, N>(ptr), m_ptr(ptr) {}
	ALG_INLINE ~alg_store_simd_block() { this->store(m_ptr); }
};  // class alg_store_simd_block
/// @}

///
/// Minimum machine size for specific types of the SIMD vector.
///
template<typename>
struct alg_simd_block_min_supported_size;

template<typename V>
static constexpr size_t alg_simd_block_min_supported_size_v = alg_simd_block_min_supported_size<V>::value;

///
/// Maximum machine size for specific types of the SIMD vector.
///
template<typename>
struct alg_simd_block_max_supported_size;

template<typename V>
static constexpr size_t alg_simd_block_max_supported_size_v = alg_simd_block_max_supported_size<V>::value;

/// @}

#define ALG_EVAL_ON_HOST 1
#define ALG_EVAL_ON_DEVICE 0

#if ALG_EVAL_ON_DEVICE
#include "../Gpgpu/GpgpuOpenCL.h"
#endif

/// @todo Correctly determine platform and so on..
#define ALG_SIMD_SUPPORT 1
#define ALG_HAS_SIMD_SUPPORT(T) (ALG_SIMD_SUPPORT && algstd::is_floating_point_v<T>)

template<typename T, typename...>
static constexpr bool alg_has_simd_support = ALG_HAS_SIMD_SUPPORT(T);

#define ALG_SSE_SUPPORT 1 && ALG_SIMD_SUPPORT
#define ALG_AVX256_SUPPORT 1
#if defined(__INTELLISENSE__)
#	define ALG_AVX512_SUPPORT 1
#else
#	define ALG_AVX512_SUPPORT 0
#endif

#if ALG_SSE_SUPPORT
#   include "SimdBlockSSE32.h"
#   include "SimdBlockSSE64.h"
#endif  // if ALG_SSE_SUPPORT
