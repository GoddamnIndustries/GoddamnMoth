// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Utils/Array.h"
#include "Utils/Parallel.h"

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++                               Permutations.                              +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

///
/// Defines a permutation.
/// @todo Determine minimal size type to store in permutation. 
///
template<size_t N>
class alg_permutation final
{
private:
	alg_array<size_t> m_data;

public:

	/// Initializes an identical permutation:
	/// sigma(i) := i
	alg_permutation() : m_data(N)
	{
		alg_parallel_for<N>([&](const size_t i)
		{
			m_data[i] = i;
		});
	}

	/// sigma(i) := j
	void swap(size_t const i, size_t const j)
	{
		std::swap(m_data[i], m_data[j]);
	}

	/// r := sigma(i)
	size_t operator() (size_t const i) const
	{
		return m_data[i];
	}

};	// class alg_permutation
