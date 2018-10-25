#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif  // ifndef _USE_MATH_DEFINES

#ifndef NOMINMAX
#define NOMINMAX
#endif  // #ifdef NOMINMAX

#include <cstdlib>
#include <cassert>
#include <cmath>

#include <stdexcept>
#include <algorithm>
#include <iostream>

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#if __NVCC__
#define MOTH_HOST __host__
#define MOTH_DEVICE __device__
#define MOTH_GLOBAL __global__
#else
#define MOTH_HOST
#define MOTH_DEVICE
#define MOTH_GLOBAL
#endif

#if __cplusplus >= 201703L || _MSVC_LANG >= 201703
#define MOTH_CPP17 1
#else
#define MOTH_CPP17 0
#endif

#define MOTH_CORE //__declspec(dllexport)

#define MOTH_NPOS SIZE_MAX

#ifdef M_PI
#define MOTH_PI M_PI
#else
#define MOTH_PI (4.0 * atan(1.0))
#endif

#ifdef M_PI_2
#define MOTH_PI_2 M_PI_2
#else
#define MOTH_PI_2 (2.0 * atan(1.0))
#endif

#ifdef M_PI_4
#define MOTH_PI_4 M_PI_4
#else
#define MOTH_PI_4 (1.0 * atan(1.0))
#endif

using moth_real_t = double;
using moth_size_t = std::size_t;
using moth_diff_t = std::ptrdiff_t;

using moth_radians_t = moth_real_t;
using moth_degrees_t = moth_real_t;

template<typename T>
int fsgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Determinant of the 2x2 matrix.
 */
MOTH_HOST MOTH_DEVICE
inline moth_real_t moth_det(moth_real_t a11, moth_real_t a12,
                            moth_real_t a21, moth_real_t a22)
{
    return a11 * a22 - a12 * a21;
}

/**
 * Determinant of the 3x3 matrix.
 */
MOTH_HOST MOTH_DEVICE
inline moth_real_t moth_det(moth_real_t a11, moth_real_t a12, moth_real_t a13,
                            moth_real_t a21, moth_real_t a22, moth_real_t a23,
                            moth_real_t a31, moth_real_t a32, moth_real_t a33)
{
    moth_real_t det{};
    det += a11 * moth_det(a22, a23,
                          a32, a33);
    det -= a12 * moth_det(a21, a23,
                          a31, a33);
    det += a13 * moth_det(a21, a22,
                          a31, a32);
    return det;
}
