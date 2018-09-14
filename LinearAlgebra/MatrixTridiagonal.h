// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Matrix.h"
#include "utils/Array.h"

#include <array>
#include <initializer_list>

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++                           Tridiagonal matrices.                          +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

///
/// Defines a generic row-major square tridiagonal matrix:
///
/// [ c0 b0 0  0  ... 0      0      0      ]
/// [ a1 c1 b1 0  ... 0      0      0      ]
/// [ 0  a2 c2 b2 ... 0      0      0      ]
/// [ 0  0  a3 c3 ... 0      0      0      ]
/// [             ...                      ]
/// [ 0  0  0  0  ... c(N-3) b(N-3) 0      ]
/// [ 0  0  0  0  ... a(N-2) c(N-2) b(N-2) ]
/// [ 0  0  0  0  ... 0      a(N-1) c(N-1) ]
///
/// @see https://en.wikipedia.org/wiki/Tridiagonal_matrix
///
template<typename T, size_t N>
class alg_tridiagonal_matrix_base : public alg_matrix_base<T, N>
{
public:

	/// Initializes a scalar matrix:
	/// Aij := delta_ij * t
	ALG_INLINE explicit alg_tridiagonal_matrix_base(const T& t = T(0.0)) : alg_matrix_base<T, N>(t)
	{
	}

public:

	/// Returns value on the diagonal at the i'th row:
	/// r := Aii = ci
	/// @{
	ALG_INLINE virtual T& diagonal(size_t const i) = 0;
	ALG_INLINE virtual const T& diagonal(size_t const i) const
	{
		return const_cast<const T&>(const_cast<alg_tridiagonal_matrix_base<T, N>*>(this)->diagonal(i));
	}
	/// @}

	/// Returns value on the diagonal at the first row:
	/// r := A00 = c0
	/// @{
	ALG_INLINE virtual T& first_diagonal()
	{
		return this->diagonal(0);
	}
	ALG_INLINE virtual const T& first_diagonal() const
	{
		return const_cast<const T&>(const_cast<alg_tridiagonal_matrix_base<T, N>*>(this)->first_diagonal());
	}
	/// @}

	/// Returns value on the diagonal at the last row:
	/// r := A(N-1)(N-1) = c(N-1)
	/// @{
	ALG_INLINE virtual T& last_diagonal()
	{
		return this->diagonal(N - 1);
	}
	ALG_INLINE virtual const T& last_diagonal() const
	{
		return const_cast<const T&>(const_cast<alg_tridiagonal_matrix_base<T, N>*>(this)->last_diagonal());
	}
	/// @}

	/// Returns value on the subdiagonal at the i'th row:
	/// r := Ai(i-1) = ai, i=1..N-1
	/// @{
	ALG_INLINE virtual T& subdiagonal(size_t const i) = 0;
	ALG_INLINE virtual const T& subdiagonal(size_t const i) const
	{
		return const_cast<const T&>(const_cast<alg_tridiagonal_matrix_base<T, N>*>(this)->subdiagonal(i));
	}
	/// @}

	/// Returns value on the diagonal at the last row:
	/// r := A(N-1)(N-2) = a(N-1)
	/// @{
	ALG_INLINE virtual T& last_subdiagonal() { return this->subdiagonal(N - 1); }
	ALG_INLINE virtual const T& last_subdiagonal() const
	{
		return const_cast<const T&>(const_cast<alg_tridiagonal_matrix_base<T, N>*>(this)->last_subdiagonal());
	}

	/// Returns value on the superdiagonal at the i'th row:
	/// r := Ai(i+1) = bi, i=0..N-2
	/// @{
	ALG_INLINE virtual T& superdiagonal(size_t const i) = 0;
	ALG_INLINE virtual const T& superdiagonal(size_t const i) const
	{
		return const_cast<const T&>(const_cast<alg_tridiagonal_matrix_base<T, N>*>(this)->superdiagonal(i));
	}
	/// @}

	/// Returns value on the diagonal at the first row:
	/// r := A01 = b0
	/// @{
	ALG_INLINE virtual T& first_superdiagonal() { return this->superdiagonal(0); }
	ALG_INLINE virtual const T& first_superdiagonal() const
	{
		return const_cast<const T&>(const_cast<alg_tridiagonal_matrix_base<T, N>*>(this)->first_superdiagonal());
	}
	/// @}

private:

	/// r := &Aij if |i-j|<=1 else nullptr
	ALG_INLINE const T* pointer_at(const size_t i, const size_t j) const
	{
		if (i == j) return &this->diagonal(i);
		if (i == j + 1) return &this->subdiagonal(i);
		if (i == j - 1) return &this->superdiagonal(i);
		return nullptr;
	}

public:
	ALG_INLINE T& operator() (const size_t i, const size_t j) override final
	{
		const T* const pointer = this->pointer_at(i, j);
		if (pointer != nullptr)
		{
			return const_cast<T&>(*pointer);
		}
		ALG_FATAL(
			"Attempting to get a mutable reference on the out of the diagonal elements "
			"of a tridiagonal matrix.");
	}
	ALG_INLINE const T& operator() (const size_t i, const size_t j) const override final
	{
		const T* const pointer = this->pointer_at(i, j);
		if (pointer != nullptr)
		{
			return *pointer;
		}
		static const T zero = {};
		return zero;
	}
};	// class alg_tridiagonal_matrix_base

ALG_DEFINE_RELATED_SCALAR_VECTOR(alg_tridiagonal_matrix_base);
ALG_DEFINE_RELATED_INVERSE(alg_tridiagonal_matrix_base, alg_matrix);

///
/// Implements a generic row-major square tridiagonal matrix.
/// Supports linear operations, inverting and system solving.
///
/// [ c0 b0 0  0  ... 0      0      0      ]
/// [ a1 c1 b1 0  ... 0      0      0      ]
/// [ 0  a2 c2 b2 ... 0      0      0      ]
/// [ 0  0  a3 c3 ... 0      0      0      ]
/// [             ...                      ]
/// [ 0  0  0  0  ... c(N-3) b(N-3) 0      ]
/// [ 0  0  0  0  ... a(N-2) c(N-2) b(N-2) ]
/// [ 0  0  0  0  ... 0      a(N-1) c(N-1) ]
///
/// @see https://en.wikipedia.org/wiki/Tridiagonal_matrix
///
template<typename T, size_t N>
class alg_tridiagonal_matrix final : public alg_tridiagonal_matrix_base<T, N>
{
private:
	bool m_is_identity = false;
	alg_array<T, N> m_diagonal;				// ci, i = 0..N-1
	alg_array<T, N - 1, 1> m_subdiagonal;	// ai, i = 1..N
	alg_array<T, N - 1> m_superdiagonal;	// bi, i = 0..N-2

public:

	/// Initializes a scalar matrix:
	/// Aij := delta_ij * t
	ALG_INLINE explicit alg_tridiagonal_matrix(const T& t) : alg_tridiagonal_matrix_base<T, N>(t)
	{
		m_diagonal.fill(t);
		m_subdiagonal.fill(T(0.0));
		m_superdiagonal.fill(T(0.0));
	}

	/// Initializes an identity matrix:
	/// Aij := delta_ij
	ALG_INLINE alg_tridiagonal_matrix() : alg_tridiagonal_matrix(T(1.0))
	{
		m_is_identity = true;
	}

	/// Initializes a tridiagonal matrix with values of some other tridiagonal matrix implementation:
	/// Aij = Bij, |i-j|<=1
	ALG_INLINE explicit alg_tridiagonal_matrix(const alg_tridiagonal_matrix_base<T, N>& B)
	{
		m_is_identity = B.is_identity();
		this->first_diagonal() = B.first_diagonal();
		this->first_superdiagonal() = B.first_superdiagonal();
		alg_nonparallel_for<N - 1>([&](const size_t i)
		{
			this->subdiagonal(i) = B.subdiagonal(i);
			this->diagonal(i) = B.diagonal(i);
			this->superdiagonal(i) = B.superdiagonal(i);
		}, 1);
		this->last_subdiagonal() = B.last_subdiagonal();
		this->last_diagonal() = B.last_diagonal();
	}
	
	/// Initializes a tridiagonal matrix with initial data: 
	/// Aij := dij
	ALG_INLINE alg_tridiagonal_matrix(const std::initializer_list<T>& data)
	{
		ALG_ASSERT(data.size() == 3 * N - 2);
		
		/// @todo Parallelize this!
		auto iterator = data.begin();
		this->first_diagonal() = *(iterator++);
		this->first_superdiagonal() = *(iterator++);
		alg_nonparallel_for<N - 1>([&](const size_t i)
		{
			this->subdiagonal(i) = *(iterator++);
			this->diagonal(i) = *(iterator++);
			this->superdiagonal(i) = *(iterator++);
		}, 1);
		this->last_subdiagonal() = *(iterator++);
		this->last_diagonal() = *(iterator);
	}

public:
	ALG_INLINE bool is_identity() const override final
	{
		return m_is_identity;
	}

public:
	ALG_INLINE T& diagonal(size_t const i) override final
	{
		return m_diagonal[i];
	}
	ALG_INLINE T& subdiagonal(size_t const i) override final
	{
		return m_subdiagonal[i];
	}
	ALG_INLINE T& superdiagonal(size_t const i) override final
	{
		return m_superdiagonal[i];
	}
};	// class alg_tridiagonal_matrix<T, N>

ALG_DEFINE_RELATED_SCALAR_VECTOR(alg_tridiagonal_matrix);
ALG_DEFINE_RELATED_INVERSE(alg_tridiagonal_matrix, alg_matrix);

///
/// Implements a generic row-major square tridiagonal matrix with same elements in all diagonals
/// except first and last rows.
/// Supports linear operations, inverting and system solving.
///
/// [ c0 b0 0  0 ... 0  0      0      ]
/// [ a  c  b  0 ... 0  0      0      ]
/// [ 0  a  c  b ... 0  0      0      ]
/// [ 0  0  a  c ... 0  0      0      ]
/// [            ...                  ]
/// [ 0  0  0  0 ... c  b      0      ]
/// [ 0  0  0  0 ... a  c      b      ]
/// [ 0  0  0  0 ... 0  a(N-1) c(N-1) ]
///
/// @see https://en.wikipedia.org/wiki/Tridiagonal_matrix
///
template<typename T, size_t N>
class alg_compressed_tridiagonal_matrix final : public alg_tridiagonal_matrix_base<T, N>
{
private:
	bool m_is_identity = false;
	T m_first_diagonal, m_first_superdiagonal;	// c0, b0
	T m_subdiagonal, m_diagonal, m_superdiagonal;	// a, c, b
	T m_last_subdiagonal, m_last_diagonal;	// a(N-1), c(N-1)

public:

	/// Initializes a scalar matrix:
	/// Aij := delta_ij * t
	ALG_INLINE explicit alg_compressed_tridiagonal_matrix(const T& t) : alg_tridiagonal_matrix_base<T, N>(t)
	{
		m_first_diagonal = t; 
		m_first_superdiagonal = T(0.0);
		m_subdiagonal = T(0.0);
		m_diagonal = t;
		m_superdiagonal = T(0.0);
		m_last_subdiagonal = T(0.0);
		m_last_diagonal = t;
	}

	/// Initializes an identity matrix:
	/// Aij := delta_ij
	ALG_INLINE alg_compressed_tridiagonal_matrix() : alg_compressed_tridiagonal_matrix(T(1.0))
	{
		m_is_identity = true;
	}

	/// Initializes a tridiagonal matrix with initial data: 
	/// Aij := dij
	ALG_INLINE alg_compressed_tridiagonal_matrix(const std::initializer_list<T>& data)
	{
		ALG_ASSERT(data.size() == 7);
		m_first_diagonal = *(data.begin() + 0);
		m_first_superdiagonal = *(data.begin() + 1);
		m_subdiagonal = *(data.begin() + 2);
		m_diagonal = *(data.begin() + 3);
		m_superdiagonal = *(data.begin() + 4);
		m_last_subdiagonal = *(data.begin() + 5);
		m_last_diagonal = *(data.begin() + 6);
	}

public:
	ALG_INLINE bool is_identity() const override final
	{
		return m_is_identity;
	}

private:
	/// @todo

public:
	ALG_INLINE virtual T& diagonal(size_t const i) override
	{
		return m_diagonal;
	}
	ALG_INLINE virtual const T& diagonal(size_t const i) const override
	{
		return m_diagonal;
	}

	ALG_INLINE T& first_diagonal() override
	{
		return m_first_diagonal;
	}

	ALG_INLINE T& last_diagonal() override
	{
		return m_last_diagonal;
	}

	ALG_INLINE T& subdiagonal(size_t const i) override
	{
		return m_subdiagonal;
	}
	ALG_INLINE const T& subdiagonal(size_t const i) const override
	{
		return m_subdiagonal;
	}

	ALG_INLINE T& last_subdiagonal() override
	{
		return m_last_subdiagonal;
	}

	ALG_INLINE T& superdiagonal(size_t const i) override
	{
		return m_superdiagonal;
	}
	ALG_INLINE const T& superdiagonal(size_t const i) const override
	{
		return m_superdiagonal;
	}

	ALG_INLINE T& first_superdiagonal() override
	{
		return m_first_superdiagonal;
	}
};	// class alg_compressed_tridiagonal_matrix<T, N>

ALG_DEFINE_RELATED_SCALAR_VECTOR(alg_compressed_tridiagonal_matrix);
ALG_DEFINE_RELATED_INVERSE(alg_compressed_tridiagonal_matrix, alg_matrix);

// ------------------------------------------------------------------------------------ //
// -----                        System solving operations.                        ----- //
// ------------------------------------------------------------------------------------ //

/// A * x = f via Thomas algorithm.
/// @see https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
/// @todo Add stability and singularity tests.
template<typename T, size_t N>
static void alg_solve_system(alg_related_vector_t<alg_matrix_base<T, N>>& x, 
	const alg_tridiagonal_matrix_base<T, N>& A,
	const alg_related_vector_t<alg_matrix_base<T, N>>& f)
{
	using scalar_t = alg_related_scalar_t<alg_matrix_base<T, N>>;
	using vector_t = alg_component_t < alg_related_vector_t<alg_matrix_base<T, N>>>;

	alg_array<alg_related_inverse_t<T>, N - 1> alpha;
	alg_array<vector_t, N> beta;

	// Forward steps of the elimination:
	T delta0 = A.first_diagonal();
	alg_related_inverse_t<T> delta;		// Since delta stores inversion results of the previous steps.
	alg_related_inverse_t<T> delta_inv;
	alg_inverse(delta_inv, delta0);
	alg_mul(alpha[0], scalar_t(-1.0), delta_inv, A.first_superdiagonal());
	alg_mul(beta [0], scalar_t(+1.0), delta_inv, f(0));

	vector_t omega;
	alg_nonparallel_for<N - 1>([&](const size_t i)
	{
		alg_mul_add(delta, scalar_t(+1.0), A.subdiagonal(i), alpha[i - 1], scalar_t(+1.0), A.diagonal(i));
		alg_inverse(delta_inv, delta);

		alg_mul_add(omega, scalar_t(-1.0), A.subdiagonal(i), beta[i - 1], scalar_t(+1.0), f(i));
		alg_mul(alpha[i], scalar_t(-1.0), delta_inv, A.superdiagonal(i));
		alg_mul(beta [i], scalar_t(+1.0), delta_inv, omega);
	}, 1);

	alg_mul_add(delta, scalar_t(+1.0), A.last_subdiagonal(), alpha[N - 2], scalar_t(+1.0), A.last_diagonal());
	alg_inverse(delta_inv, delta);

	alg_mul_add(omega, scalar_t(-1.0), A.last_subdiagonal(), beta[N - 2], scalar_t(+1.0), f(N - 1));
	alg_mul(beta[N - 1], scalar_t(1.0), delta_inv, omega);

	// Backward substitution steps:
	x(N - 1) = beta[N - 1];
	alg_nonparallel_for_reverse<N - 1>([&](const size_t i)
	{
		alg_mul_add(x(i), scalar_t(1.0), alpha[i], x(i + 1), scalar_t(1.0), beta[i]);
	});
}

// ------------------------------------------------------------------------------------ //
// -----                       Matrix inversion operations.                       ----- //
// ------------------------------------------------------------------------------------ //

/// B := A^1 via special algorithm for tridiagonal matrices.
//template<typename T, size_t N>
//static void alg_inverse(alg_matrix<T, N>& B, const alg_tridiagonal_matrix_base<T, N>& A)
//{
//	/*std::array<T, N + 1> theta;
//	theta[0] = T(1.0);
//	theta[1] = A.diagonal(0);
//	alg_nonparallel_for<N>([&](const size_t i)
//	{
//
//	});
//
//	std::array<T, N> phi;*/
//
//}
