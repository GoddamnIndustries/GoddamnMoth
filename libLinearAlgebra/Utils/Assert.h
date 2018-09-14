// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#define ALG_TEMPLATES_TESTING 1

#include <cassert>
#include <cmath>
#include <math.h>

#define ALG_ASSERT(...) do { if(!(__VA_ARGS__)) { puts(#__VA_ARGS__); abort(); } } while(false)
#define ALG_STATIC_ASSERT(...) static_assert((__VA_ARGS__), #__VA_ARGS__)
#define ALG_DEBUG_ASSERT(...) //assert((__VA_ARGS__))
#define ALG_FATAL(...) do { ALG_ASSERT(__VA_ARGS__); abort(); } while (false)

#include "Parallel.h"

ALG_INLINE static bool alg_verify_prop_impl(const size_t)
{
	return true;
}
template<typename T, typename... TRest>
ALG_INLINE static bool alg_verify_prop_impl(const size_t size, const T& first, const TRest&... rest)
{
	return size == first.size() && alg_verify_prop_impl(size, rest...);
}
template<typename T, typename... TRest>
ALG_INLINE static void _alg_verify_size(size_t& size, const T& first, const TRest&... rest)
{
	size = first.size();
	//ALG_ASSERT(_alg_verify_size_impl(size, rest...));
}
#define alg_verify_prop _alg_verify_size

template<typename T>
ALG_INLINE static void _alg_verify_coefs(T const& alpha, T const& beta = 1.0, T const& gamma = 1.0)
{
	/*if (!isfinite(alpha)) *(T*)&alpha = 0;
	if (!isfinite(beta)) *(T*)&beta = 0;
	if (!isfinite(gamma)) *(T*)&gamma = 0;*/

	ALG_ASSERT(isfinite(alpha));
	ALG_ASSERT(isfinite(beta));
	ALG_ASSERT(isfinite(gamma));
}
#define alg_verify_coef _alg_verify_coefs

template<typename T>
ALG_INLINE static void _alg_verify_addr(...)
{
}
#define alg_verify_addr(...)
