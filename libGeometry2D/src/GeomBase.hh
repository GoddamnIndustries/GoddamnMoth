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
#include <iosfwd>

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

#define MOTH_CORE //__declspec(dllexport)

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
