// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Scalar.h"
#include "Vector.h"
#include "Utils/Parallel.h"
#include "Utils/Array.h"

#include <initializer_list>
#include <iostream>

/// @defgroup alg_dense_matrices (Dense matrices)
/// @{

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++                          Generic dense matrices.                         +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

///
/// Defines a generic row-major square matrix.
///
/// [ a00     a01     a02     ... a0(N-1)     ]
/// [ a10     a11     a12     ... a1(N-1)     ]
/// [                         ...             ]
/// [ a(N-1)0 a(N-1)1 a(N-1)2 ... a(N-1)(N-1) ]
///
/// @see https://en.wikipedia.org/wiki/Matrix_(mathematics)
///
template<typename T>
class alg_matrix final : public alg_vector<T>
{
private:
	size_t m_size;
	bool m_is_identity = false;

public:

	/// Initializes an identity matrix:
	/// Aij := delta_ij
	ALG_INLINE explicit alg_matrix(const size_t n = 0)
		: alg_matrix(n, T(1.0))
	{
		m_is_identity = true;
	}

	/// Initializes a scalar matrix:
	/// Aij := delta_ij t
	ALG_INLINE explicit alg_matrix(const size_t n, const T& t)
		: alg_vector<T>(n * n), m_size(n)
	{
		alg_parallel_for(size(), size(), [&](const size_t i, const size_t j)
		{
			(*this)(i, j) = i != j ? T(0.0) : t;
		});
	}

	/// Initializes a dense matrix with values of some other matrix implementation:
	/// Aij := Bij
	template<template<typename> class TMatrix>
	ALG_INLINE explicit alg_matrix(const TMatrix<T>& B)
		: alg_vector<T>(B.size() * B.size()), m_size(B.size()), m_is_identity(B.is_identity())
	{
		alg_parallel_for(size(), size(), [&](const size_t i, const size_t j)
		{
			(*this)(i, j) = B(i, j);
		});
	}

	/// Initializes a dense matrix with initial data:
	/// Aij := dij
	ALG_INLINE alg_matrix(const std::initializer_list<T>& data)
		: alg_vector<T>(data), m_size(static_cast<size_t>(std::sqrt(data.size())))
	{}

public:

	/// Returns size of the matrix.
	ALG_INLINE size_t size() const
	{
		return m_size;
	}

	/// Resizes a vector.
	ALG_INLINE void resize(const size_t n)
	{
		m_size = n;
		alg_vector<T>::resize(n * n);
	}

	/// Returns true if this matrix is surely identity.
	ALG_INLINE bool is_identity() const
	{
		return m_is_identity;
	}

public:

	ALG_INLINE alg_vector<T> row(const size_t i) const
	{
		//return alg_vector<T>(this->data() + i * m_size, m_size);
		return alg_vector<T>(*this, i * m_size, m_size);
	}

	ALG_INLINE alg_vector<T> column(const size_t i) const
	{
		return alg_vector<T>(*this, i, m_size, m_size);
		//return alg_vector<T>(this->data() + i, m_size, m_size);
	}

public:
	ALG_INLINE T& operator() (const size_t i, const size_t j)
	{
		m_is_identity = false;
		return alg_vector<T>::operator()(i * m_size + j);
	}
	ALG_INLINE const T& operator() (const size_t i, const size_t j) const
	{
		return alg_vector<T>::operator()(i * m_size + j);
	}
};  // class alg_matrix<T>

///
/// Related scalar type.
///
template<>
struct alg_related_scalar<alg_matrix<alg_float32_t>> { using type = alg_float32_t; };
template<>
struct alg_related_scalar<alg_matrix<alg_float64_t>> { using type = alg_float64_t; };
template<typename T>
struct alg_related_scalar<alg_matrix<T>> { using type = typename alg_related_scalar<T>::type; };

#if ALG_TEMPLATES_TESTING
static_assert(algstd::is_same_v<alg_float32_t, alg_related_scalar_t<alg_matrix<alg_float32_t>>>, "");
static_assert(algstd::is_same_v<alg_float32_t, alg_related_scalar_t<alg_matrix<alg_vector<alg_float32_t>>>>, "");
#endif	// if ALG_TEMPLATES_TESTING

///
/// Related vector type.
///
template<typename>
struct alg_related_vector;
template<>
struct alg_related_vector<alg_matrix<alg_float32_t>> { using type = alg_vector<alg_float32_t>; };
template<>
struct alg_related_vector<alg_matrix<alg_float64_t>> { using type = alg_vector<alg_float64_t>; };
template<typename T>
struct alg_related_vector<alg_matrix<T>> { using type = alg_vector<typename alg_related_vector<T>::type>; };

template<typename M>
using alg_related_vector_t = typename alg_related_vector<M>::type;

#if ALG_TEMPLATES_TESTING
static_assert(algstd::is_same_v<alg_vector<alg_float32_t>, alg_related_vector_t<alg_matrix<alg_float32_t>>>, "");
static_assert(algstd::is_same_v<alg_vector<alg_vector<alg_float32_t>>, alg_related_vector_t<alg_matrix<alg_matrix<alg_float32_t>>>>, "");
#endif	// if ALG_TEMPLATES_TESTING

///
/// Wraps related scalar and vector types for specific matrices implementations.
///
#define ALG_DEFINE_RELATED_SCALAR_VECTOR(Matrix) \
	template<typename T> \
	struct alg_related_scalar<Matrix<T>> { using type = alg_related_scalar_t<alg_matrix<T>>; }; \
	template<typename T> \
	struct alg_related_vector<Matrix<T>> { using type = alg_related_vector_t<alg_matrix<T>>; }

///
/// Wraps related inverse type for matrices.
///
#define ALG_DEFINE_RELATED_INVERSE(BaseMatrix, InverseMatrix) \
	template<typename T> \
	struct alg_related_inverse<BaseMatrix<T>> { using type = InverseMatrix<T>; }

// ------------------------------------------------------------------------------------ //
// -----                          Matrix-Vector products.                         ----- //
// ------------------------------------------------------------------------------------ //

/// @f$ r := (A \vec x)_i @f$
template<typename T>
static void alg_mul_i(alg_component_t<alg_related_vector_t<alg_matrix<T>>>& r,
    const alg_matrix<T>& A, const alg_related_vector_t<alg_matrix<T>>& x,
	const size_t i)
{
    r = alg_dot(A.row(i), x);
}

/// @f$ \vec x := \alpha A \vec y @f$
template<typename T>
static void alg_mul(alg_related_vector_t<alg_matrix<T>>& x,
    const alg_related_scalar_t<alg_matrix<T>> alpha, const alg_matrix<T>& A, const alg_related_vector_t<alg_matrix<T>>& y)
{
	// x and y should not be same vectors.
	alg_verify_addr(x, y);

	size_t size;
	alg_verify_prop(size, x, y, A);
	alg_verify_coef(alpha);
    alg_parallel_for(size, [&](const size_t i)
    {
        alg_mul_i(x(i), A, y, i);
        alg_mul_assign(x(i), alpha);
    });
}

/// @f$ \vec x := \alpha A \vec x @f$
template<typename T>
static void alg_mul_assign(alg_related_vector_t<alg_matrix<T>>& x,
    const alg_related_scalar_t<alg_matrix<T>> alpha, const alg_matrix<T>& A)
{
	alg_vector<T> y;
	alg_mul(y, alpha, A, x);
	x = y;
}

// ------------------------------------------------------------------------------------ //

/// @f$ x := \alpha A \vec y + \beta \vec z @f$
template<typename T>
static void alg_mul_add(alg_related_vector_t<alg_matrix<T>>& x,
    const alg_related_scalar_t<alg_vector<T>> alpha,  const alg_matrix<T>& A, const alg_related_vector_t<alg_matrix<T>>& y,
    const alg_related_scalar_t<alg_vector<T>> beta,                           const alg_related_vector_t<alg_matrix<T>>& z)
{
	// x and y should not be same operands.
	alg_verify_addr(x, y, z);

	size_t size;
	alg_verify_prop(size, x, y, z, A);
	alg_verify_coef(alpha, beta);
    alg_parallel_for(size, [&](const size_t i)
    {
        alg_mul_i(x(i), A, y, i);
        alg_add(x(i), alpha, x(i), beta, z(i));
    });
}

/// @f$ \vec x += \alpha A \vec y @f$
template<typename T>
static void alg_mul_add_assign(alg_related_vector_t<alg_matrix<T>>& x,
    const alg_related_scalar_t<alg_vector<T>> alpha, const alg_matrix<T>& A, const alg_related_vector_t<alg_matrix<T>>& y)
{
	size_t size;
	alg_verify_prop(size, x, y, A);
	alg_verify_coef(alpha);
	alg_parallel_for(size, [&](const size_t i)
    {
		alg_component_t<alg_related_vector_t<alg_matrix<T>>> r;
		alg_mul_i(r, A, y, i);
        alg_add_assign(x(i), alpha, r);
    });
}

// ------------------------------------------------------------------------------------ //
// -----                          Matrix-Matrix products.                         ----- //
// ------------------------------------------------------------------------------------ //

/// @f$ r := (B C)_{i,j} @f$
template<typename T>
static void alg_mul_ij(T& r,
    const alg_matrix<T>& B, const alg_matrix<T>& C,
    const size_t i, const size_t j)
{
	r = alg_dot(B.row(i), C.column(j));
}

/// @f$ A := \alpha B C @f$
/// @note Result matrix should not be same with any operands.
template<typename T>
static void alg_mul(alg_matrix<T>& A,
    const alg_related_scalar_t<alg_matrix<T>> alpha, const alg_matrix<T>& B,
                                                     const alg_matrix<T>& C)
{
	// Neither B nor C should not be same operands with A.
	alg_verify_addr(A, B, C);

	size_t size;
	alg_verify_prop(size, A, B, C);
	alg_verify_coef(alpha);
	alg_parallel_for(size, [&](const size_t i, const size_t j)
	{
		alg_mul_ij(A(i, j), B, C, i, j);
		alg_mul_assign(A(i, j), alpha);
	});
}

/// @f$ A := \alpha B A @f$
template<typename T>
static void alg_mul_assign(alg_matrix<T>& A,
    const alg_related_scalar_t<alg_matrix<T>> alpha, const alg_matrix<T>& B)
{
	alg_matrix<T> D;
	alg_mul(D, alpha, B, A);
	A = D;
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/// @f$ A := \alpha B C + \beta D @f$
/// @note Result matrix should not be same with any of multiplication operands.
template<typename T>
static void alg_mul_add(alg_matrix<T>& A,
    const alg_related_scalar_t<alg_matrix<T>> alpha, const alg_matrix<T>& B, const alg_matrix<T>& C,
    const alg_related_scalar_t<alg_matrix<T>> beta,  const alg_matrix<T>& D)
{
	// Neither B nor C should not be same operands with A.
	alg_verify_addr(A, B, C);

	size_t size;
	alg_verify_prop(size, A, B, C, D);
	alg_verify_coef(alpha, beta);
    alg_parallel_for(size, [&](const size_t i, const size_t j)
    {
        alg_mul_ij(A(i, j), B, C, i, j);
        alg_add(A(i, j), alpha, A(i, j), beta, D(i, j));
    });
}

/// @f$ A += \alpha B C @f$
template<typename T>
static void alg_mul_add_assign(alg_matrix<T>& A,
    const alg_related_scalar_t<alg_matrix<T>> alpha, const alg_matrix<T>& B, const alg_matrix<T>& C)
{
	// Neither B nor C should not be same operands with A.
	alg_verify_addr(A, B, C);

	size_t size;
	alg_verify_prop(size, A, B, C);
	alg_verify_coef(alpha);
	alg_parallel_for(size, [&](const size_t i, const size_t j)
	{
		T r;
		alg_mul_ij(r, B, C, i, j);
		alg_add_assign(A(i, j), alpha, r);
	});
}

// ------------------------------------------------------------------------------------ //
// -----                        System solving operations.                        ----- //
// ------------------------------------------------------------------------------------ //

/// Solves @f$ B \vec x = \vec g @f$ via Gaussian elimination with leading element selection.
/// @todo Add stability and singularity tests.
template<typename T>
static void alg_solve_system(alg_related_vector_t<alg_matrix<T>>& x,
	const alg_matrix<T>& B,
	const alg_related_vector_t<alg_matrix<T>>& g)
{
    using scalar_t = alg_related_scalar_t<alg_matrix<T>>;

	size_t size;
	alg_verify_prop(size, B, x, g);

	// Now we are solving the C * x = f system. C and f are modifiable, C is dense.
	alg_matrix<T> A(B);
    alg_related_vector_t<alg_matrix<T>> f(g);

	// Column permutation for leading element selection.
	// Finally, we are solving the A x = f, A and f are modifiable, A is permutable.
	/*alg_permutation<N> row_permutation;
	const auto A = [&row_permutation, &C](const size_t i, const size_t j) -> T&
	{
		return C(i, row_permutation(j));
	};*/
	
    // Forward steps of the elimination:
    alg_nonparallel_for(size, [&](const size_t i)
    {
		/// @todo Leading element?!
		/// @todo Use fnear here..
        //assert(A(i, i) != T(0.0) && "Near singular system.");
		alg_related_inverse_t<T> Aii_inv;
        alg_inverse(Aii_inv, A(i, i));

        alg_nonparallel_for(size, [&](const size_t j) { alg_mul_assign(A(i, j), scalar_t(1.0), Aii_inv); }, i + 1);
        alg_mul_assign(f(i), scalar_t(1.0), Aii_inv);
        A(i, i) = T(1.0);

        alg_nonparallel_for(size, [&](const size_t p)
        {
            alg_nonparallel_for(size, [&](const size_t q) { alg_mul_add_assign(A(p, q), scalar_t(-1.0), A(p, i), A(i, q)); }, i + 1);
            alg_mul_add_assign(f(p), scalar_t(-1.0), A(p, i), f(i));
            A(p, i) = T(0.0);
        }, i + 1);
    });

    // Backward substitution steps:
    alg_nonparallel_for_reverse(size, [&](const size_t i)
    {
        x(i) = f(i);
        alg_nonparallel_for(size, [&](const size_t j) { alg_mul_add_assign(x(i), scalar_t(-1.0), A(i, j), x(j)); }, i + 1);
    });
}

/// Solves @f$ A \vec x = \vec f @f$ via Biconjugate gradient stabilized method.
/// @see https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
template<typename T, template<typename> class TMatrix>
static void alg_solve_system_iter_bicgstab(alg_related_vector_t<TMatrix<T>>& x,
	const TMatrix<T>& A,
	const alg_related_vector_t<TMatrix<T>>& f,
	const alg_related_scalar_t<TMatrix<T>> eps = alg_related_scalar_t<TMatrix<T>>(ALG_EPS))
{
	using scalar_t = alg_related_scalar_t<alg_matrix<T>>;

	size_t size;
	alg_verify_prop(size, x, f, A);

	//x = f;

	// Preparations:
	alg_vector<T>& x_0 = x, &x_km = x_0, x_k(x_km);
	//alg_vector<T> x_0(size, scalar_t(0.0)), &x_km = x_0, x_k(x_km);
	alg_vector<T> r_0(size, scalar_t(0.0)), &r_km = r_0, r_k(r_km);
	alg_vector<T> v_0(size, scalar_t(0.0)), &v_km = v_0, v_k(v_km);
	alg_vector<T> p_0(size, scalar_t(0.0)), &p_km = p_0, p_k(p_km);

	alg_vector<T> s_k(size, scalar_t(0.0));
	alg_vector<T> t_k(size, scalar_t(0.0));

	scalar_t rho_km = scalar_t(1.0);
	scalar_t alpha_km = scalar_t(1.0);
	scalar_t omega_km = scalar_t(1.0);

	alg_mul_add(r_0, scalar_t(-1.0), A, x_0, scalar_t(1.0), f);
	const scalar_t keke = alg_norm_l2(r_0);
	if (keke < eps*eps) {
		//std::cerr << k << std::endl;
		return;
	}

	const alg_vector<T> r_hat(r_0);

	// Iterations:
	for (size_t k = 1; k < size; ++k)
	{
		const scalar_t rho_k = alg_dot(r_hat, r_km);
		const scalar_t beta_k = rho_k / rho_km * alpha_km / omega_km;
		alg_add(p_k, beta_k, p_km, -beta_k * omega_km, v_km);
		alg_add_assign(p_k, scalar_t(1.0), r_km);
		alg_mul(v_k, scalar_t(1.0), A, p_k);

		const scalar_t alpha_k = rho_k / alg_dot(r_hat, v_k);
		alg_add(s_k, scalar_t(1.0), r_km, -alpha_k, v_k);
		alg_mul(t_k, scalar_t(1.0), A, s_k);

		const scalar_t omega_k = alg_dot(t_k, s_k) / alg_norm_l2_sq(t_k);
		alg_add(x_k, omega_k, s_k, alpha_k, p_k);
		alg_add_assign(x_k, scalar_t(1.0), x_km);
		alg_add(r_k, scalar_t(1.0), s_k, -omega_k, t_k);

		const scalar_t keke = alg_norm_l2(r_k) / alg_norm_l2(f);
		//std::cerr << keke << std::endl;
		if (keke < eps) {
			//std::cerr << k << std::endl;
			break;
		}

		x_km = x_k;
		r_km = r_k;
		v_km = v_k;
		p_km = p_k;

		rho_km = rho_k;
		alpha_km = alpha_k;
		omega_km = omega_k;
	}

	x = x_km;
}

// ------------------------------------------------------------------------------------ //
// -----                       Matrix inversion operations.                       ----- //
// ------------------------------------------------------------------------------------ //

/// @f$ B := A^1 @f$ via N system solving (for non-block matrices).
template<typename T, typename = algstd::enable_if_t<algstd::is_floating_point_v<T>>>
static void alg_inverse(alg_matrix<T>& B, 
	const alg_matrix<T>& A)
{
    alg_parallel_for(B.size(), [&](const size_t j)
    {
        alg_related_vector_t<alg_matrix<T>> g, x;
        g(j) = T(1.0);

        alg_solve_system(x, A, g);
        alg_parallel_for(B.size(), [&](const size_t i)
        {
            B(i, j) = x(i);
        });
    });
}

/// @}
