// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MatrixTridiagonal.h"
#include "Test.h"

using vector3 = alg_vector<alg_float64_t, 3>;
using bvector3 = alg_vector<alg_vector<alg_float64_t, 3>, 3>;
using matrix3 = alg_matrix<alg_float64_t, 3>;
using tmatrix3 = alg_tridiagonal_matrix<alg_float64_t, 3>;
using ctmatrix3 = alg_compressed_tridiagonal_matrix<alg_float64_t, 3>;
using bctmatrix3 = alg_compressed_tridiagonal_matrix<alg_compressed_tridiagonal_matrix<alg_float64_t, 3>, 3>;

// ------------------------------------------------------------------------------------ //
// -----                        System solving operations.                        ----- //
// ------------------------------------------------------------------------------------ //

testing_unit_test(tridiagonal, matrix_solve_system)
{
	vector3 x;
	const tmatrix3 A = {
		+2.0, +1.0, 
		+1.0, +2.0, +1.0,
		      +1.0, +2.0,
	};
	const vector3 f = {
		+1.0,
		+0.0,
		+1.0,
	};
	const vector3 r = {
		+1.0,
		-1.0,
		+1.0,
	};

	bool equals;
	alg_solve_system(x, A, f);
	equals = alg_near_eql(x, r);
	testing_verify(equals);
};

testing_unit_test(ctridiagonal, matrix_solve_system)
{
	vector3 x;
	const ctmatrix3 A = {
		+2.0, +1.0,
		+1.0, +2.0, +1.0,
		      +1.0, +2.0,
	};
	const vector3 f = {
		+1.0,
		+0.0,
		+1.0,
	};
	const vector3 r = {
		+1.0,
		-1.0,
		+1.0,
	};

	bool equals;
	alg_solve_system(x, A, f);
	equals = alg_near_eql(x, r);
	testing_verify(equals);
};

using bvector7 = alg_vector<alg_vector<alg_float64_t, 7>, 3>;
using ctmatrix7 = alg_compressed_tridiagonal_matrix<alg_float64_t, 7>;
using bctmatrix7 = alg_compressed_tridiagonal_matrix<alg_compressed_tridiagonal_matrix<alg_float64_t, 7>, 3>;

testing_unit_test(tridiagonal2, matrix_solve_system)
{
	bvector7 x;
	static const ctmatrix7 I, E(2.0);
	static const ctmatrix7 C = {
		-4, +1,
		+1, -4, +1,
	        +1, -4,
	};
	static const bctmatrix7 A = {
		C, E,
		I, C, I,
		   E, C,
	};
	static const bvector7 f = { {-3,0,0,0,0,0,-3},{ -3,0,0,0,0,0,-3 },{ -3,0,0,0,0,0,-3 } };
	static const bvector7 r = { {3,3,3,3,3,3,3},{3,3,3,3,3,3,3},{3,3,3,3,3,3,3} };

	bool equals;
	alg_solve_system(x, A, f);
	equals = alg_near_eql(x, r);
	testing_verify(equals);
};
