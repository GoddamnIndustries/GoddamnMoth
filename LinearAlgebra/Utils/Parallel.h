// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <cstdlib>
#include <cstdint>
#include <type_traits>

#ifndef _MSC_VER
#	define __forceinline

#endif

#define ALG_INLINE __forceinline
#define ALG_HOST
#pragma warning(disable : 4505)

namespace algstd
{
	template<typename T, typename U>
	static constexpr bool is_same_v = std::is_same<T, U>::value;

	template<typename T>
	static constexpr bool is_const_v = std::is_const<T>::value;

	template<typename T>
	static constexpr bool is_pod_v = std::is_pod<T>::value;

	template<bool B, typename T, typename F>
	using conditional_t = typename std::conditional<B, T, F>::type;

	template<bool B, typename T = void>
	using enable_if_t = typename std::enable_if<B, T>::type;

	template<bool B, typename T = void>
	using disable_if_t = typename std::enable_if<!B, T>::type;

	template<typename T>
	static constexpr bool is_floating_point_v = std::is_floating_point<T>::value;
}

using alg_enabled_t = void;
template<bool B, typename T = alg_enabled_t>
using alg_enable_if_t = typename std::enable_if<B, T>::type;

template<size_t N>
using alg_constexpr_size_t = std::integral_constant<size_t, N>;

template<typename T, typename...>
using alg_first_t = T;

#ifndef ALG_DOXYGEN_SHOULD_SKIP_THIS
#define ALG_DOXYGEN_SHOULD_SKIP_THIS 1
#endif  // ifndef ALG_DOXYGEN_SHOULD_SKIP_THIS

// ------------------------------------------------------------------------------------ //

///
/// Simple 1D for-each loop.
///
template<typename P>
static ALG_INLINE void alg_nonparallel_for(const size_t n, 
	const P& predicate, 
	const size_t i_start = 0)
{
	for (size_t i = i_start; i < n; ++i)
	{
		predicate(i);
	}
}

///
/// Simple 1D reverse for-each loop.
///
template<typename P>
static ALG_INLINE void alg_nonparallel_for_reverse(const size_t n, 
	const P& predicate)
{
	for (size_t i = n - 1; i != SIZE_MAX; --i)
	{
		predicate(i);
	}
}

///
/// Simple 2D for-each loop.
///
template<typename P>
static ALG_INLINE void alg_nonparallel_for(const size_t n, const size_t m, 
	const P& predicate, 
	const size_t i_start = 0, const size_t j_start = 0)
{
	for (size_t i = i_start; i < n; ++i)
	{
		for (size_t j = j_start; j < m; ++j)
		{
			predicate(i, j);
		}
	}
}

// ------------------------------------------------------------------------------------ //

///
/// OMP-accelerated parallel execution for two labmdas.
///
template<typename P1, typename P2>
static ALG_INLINE void alg_parallel(
	const P1& predicate1, const P2& predicate2)
{
#pragma omp parallel sections 
	{
#pragma omp section
		{
			predicate1();
		}

#pragma omp section
		{
			predicate2();
		}
	}
}

///
/// OMP-accelerated parallel execution for two labmdas.
///
template<typename P1, typename P2, typename P3>
static ALG_INLINE void alg_parallel(
	const P1& predicate1, const P2& predicate2, const P3& predicate3)
{
#pragma omp parallel sections 
	{
#pragma omp section
		{
			predicate1();
		}

#pragma omp section
		{
			predicate2();
		}

#pragma omp section
		{
			predicate3();
		}
	}
}

// ------------------------------------------------------------------------------------ //

///
/// OMP-accelerated parallel 1D for loop.
///
template<typename P>
static ALG_INLINE void alg_parallel_for(const size_t n,
	const P& predicate, 
	const size_t i_start = 0)
{
#pragma omp parallel for
    for (int64_t _i = static_cast<int64_t>(i_start); _i < n; ++_i)
    {
		const size_t i = static_cast<size_t>(_i);
        predicate(i);
    }
}

///
/// OMP-accelerated parallel 2D for-each loop.
///
template<typename P>
static ALG_INLINE void alg_parallel_for(const size_t n, const size_t m,
	const P& predicate, 
	const size_t i_start = 0, const size_t j_start = 0)
{
#pragma omp parallel for
    for (int64_t _i = static_cast<int64_t>(i_start); _i < n; ++_i)
    {
		const size_t i = static_cast<size_t>(_i);
		for (size_t j = j_start; j < m; ++j)
        {
            predicate(i, j);
        }
    }
}
