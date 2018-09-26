#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif  // ifndef _USE_MATH_DEFINES

#ifndef NOMINMAX
#define NOMINMAX
#endif  // #ifdef NOMINMAX

#include <cstdlib>
#include <cmath>

#include <stdexcept>
#include <algorithm>
#include <iosfwd>

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#define GEOM_CORE __declspec(dllexport)

#ifdef M_PI
#define GEOM_PI M_PI
#else
#define GEOM_PI (4.0 * atan(1.0))
#endif

#ifdef M_PI_2
#define GEOM_PI_2 M_PI_2
#else
#define GEOM_PI_2 (2.0 * atan(1.0))
#endif

#ifdef M_PI_4
#define GEOM_PI_4 M_PI_4
#else
#define GEOM_PI_4 (1.0 * atan(1.0))
#endif

using geom_real_t = double;
using geom_size_t = std::size_t;

using geom_radians_t = geom_real_t;
using geom_degrees_t = geom_real_t;

template<typename T>
int fsgn(T val)
{
    return (T(0) < val) - (val < T(0));
}
