// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Matrix.h"
#include "Utils/Progression.h"

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++                     Triangluar and symmetric matrices.                   +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

///
/// @todo
///
enum alg_half_storage_matrix_type
{
	alg_half_storage_matrix_type_upper_triangular,
	alg_half_storage_matrix_type_lower_triangular,
	alg_half_storage_matrix_type_symmetric,
	alg_half_storage_matrix_type_skew_symmetric,
};	// enum alg_half_storage_matrix_type

///
/// Defines a generic matrix with only half degrees of freedom: triangular matrix, symmetric and others.
///
template<typename T, size_t N>
class alg_half_storage_matrix_base_base : public alg_matrix_base<T, N>
{
private:
	alg_array<T, arithm_progression_Sn(size_t(1), size_t(1), N)> m_data;

};	// class alg_half_storage_matrix_base_base<T, N>

template<typename T, size_t N>
class alg_half_storage_matrix_base;

template<typename T, size_t N>
class alg_upper_triangular_matrix;
template<typename T, size_t N>
class alg_upper_triangular_matrix_base;
template<typename T, size_t N>
class alg_lower_triangular_matrix;
template<typename T, size_t N>
class alg_lower_triangular_matrix_base;
template<typename T, size_t N>
class alg_symmetric_matrix;
template<typename T, size_t N>
class alg_symmetric_matrix_base;

// ------------------------------------------------------------------------------------ //
// -----                        System solving operations.                        ----- //
// ------------------------------------------------------------------------------------ //

/// B * x = g via an algorithm for upper triangluar matrices.
template<typename T, size_t N>
static void alg_solve_system(alg_related_vector_t<alg_upper_triangular_matrix_base<T, N>>& x,
	const alg_upper_triangular_matrix_base<T, N>& B,
	const alg_related_vector_t<alg_upper_triangular_matrix_base<T, N>>& g)
{
	using scalar_t = alg_related_scalar_t<alg_matrix_base<T, N>>;

	// Backward substitution steps of the Gaussian elemenation:
	alg_nonparallel_for_reverse<N>([&](const size_t i)
	{
		alg_related_inverse_t<T> Bii_inv;
		alg_inverse(Bii_inv, B(i, i));

		x(i) = g(i);
		alg_foreach_np<N>([&](const size_t j) { alg_mul_add_assign(x(i), scalar_t(-1.0), B(i, j), x(j)); }, i + 1);
		alg_mul_assign(x(i), scalar_t(1.0), Bii_inv);
	});
}

/// B * x = g via an algorithm for lower triangluar matrices.
template<typename T, size_t N>
static void alg_solve_system(alg_related_vector_t<alg_lower_triangular_matrix<T, N>>& x,
	const alg_lower_triangular_matrix<T, N>& B,
	const alg_related_vector_t<alg_lower_triangular_matrix<T, N>>& g)
{
	using scalar_t = alg_related_scalar_t<alg_matrix_base<T, N>>;

	alg_parallel_for<N>([&](const size_t i)
	{
		alg_related_inverse_t<T> Aii_inv;
		alg_inverse(Aii_inv, B(i, i));

		x(i) = g(i);
		/// @todo
		alg_mul_assign(x(i), scalar_t(1.0), Aii_inv);
	});
}

// ------------------------------------------------------------------------------------ //
// -----                          Decomposing operations.                         ----- //
// ------------------------------------------------------------------------------------ //

/// A = LU, L is lower trianglular and U is upper triangular.
/// @see https://en.wikipedia.org/wiki/LU_decomposition
template<typename T, size_t N>
static void alg_matrix_decompose_lu(alg_lower_triangular_matrix<T, N>& L, alg_upper_triangular_matrix<T, N>& U,
	const alg_matrix_base<T, N>& A)
{
}

/// A = QR, Q is orthogonal and R is upper triangular via Householder reflections.
/// @see https://en.wikipedia.org/wiki/QR_decomposition
template<typename T, size_t N>
static void alg_matrix_decompose_qr_householder(alg_matrix<T, N>& L, alg_upper_triangular_matrix<T, N>& R,
	const alg_matrix_base<T, N>& A)
{
}

/// A = QR, Q is orthogonal and R is upper triangular via Givens rotations.
/// @see https://en.wikipedia.org/wiki/QR_decomposition
template<typename T, size_t N>
static void alg_matrix_decompose_qr_givens(alg_matrix<T, N>& L, alg_upper_triangular_matrix<T, N>& R,
	const alg_matrix_base<T, N>& A)
{
}

/// A = LL^T, L is lower trianglular and A is symmetric.
/// @see https://en.wikipedia.org/wiki/Cholesky_decomposition
template<typename T, size_t N>
static void alg_matrix_decompose_llt_cholesky(alg_lower_triangular_matrix<T, N>& L,
	const alg_symmetric_matrix_base<T, N>& A)
{
}

// ------------------------------------------------------------------------------------ //
// -----                       Matrix inversion operations.                       ----- //
// ------------------------------------------------------------------------------------ //

/// B := A^1 via matrix decomposition and N system solving (for non-block matrices).
template<typename T, size_t N, typename = algstd::enable_if_t<algstd::is_floating_point_v<T>>>
static void alg_inverse(alg_matrix<T, N>& B,
	const alg_matrix<T, N>& A)
{
	alg_lower_triangular_matrix<T, N> L;
	alg_upper_triangular_matrix<T, N> U;
	alg_matrix_decompose_lu(L, U, A);

	alg_parallel_for<N>([&](const size_t i)
	{
		alg_related_vector_t<alg_matrix_base<T, N>> g, x, y;
		g(i) = T(1.0);

		alg_solve_system(y, L, g);
		alg_solve_system(x, U, y);
		alg_parallel_for<N>([&](const size_t j)
		{
			B(i, j) = x(j);
		});
	});
}

/// B := A^1 via matrix decomposition and N system solving (for non-block matrices), A is symmetric.
template<typename T, size_t N, typename = algstd::enable_if_t<algstd::is_floating_point_v<T>>>
static void alg_inverse(alg_symmetric_matrix<T, N>& B,
	const alg_symmetric_matrix_base<T, N>& A)
{
	alg_lower_triangular_matrix<T, N> L;
	alg_matrix_decompose_llt_cholesky(L, A);
	/// @todo
}
