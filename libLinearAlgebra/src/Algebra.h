// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg. 
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once
#include <array>
#define __forceinline

#define ALG_TEMPLATES_TESTING 1
#if ALG_TEMPLATES_TESTING
#	include <type_traits>
#endif	// if ALG_TEMPLATES_TESTING

namespace std
{
    template<typename T, typename U>
    static constexpr bool is_same_v = std::is_same<T, U>::value;
}

template<typename T, typename U>
static __forceinline T alg_union_cast(const U& u)
{
}
