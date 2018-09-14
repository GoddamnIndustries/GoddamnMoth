// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++                          Arithmetic progression.                         +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

/// Computes n'th element of the arithmetic progression.
/// an := a1 + (n - 1) * d.
template<typename T>
static constexpr T arithm_progression_an(const T a1, const T d, const T n)
{
	return a1 + (n - 1) * d;
}

/// Computes sum of the elements of the arithmetic progression.
/// Sn := (a1 + d * (n - 1) / 2) * n.
template<typename T>
static constexpr T arithm_progression_Sn(const T a1, const T d, const T n)
{
	return (a1 + d * (n - 1) / 2) * n;
}
