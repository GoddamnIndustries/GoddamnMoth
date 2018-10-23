#pragma once

#include "libGeometry2D/src/GeomPoint.hh"

/**
 * Triangle in 3D space.
 */
struct MOTH_CORE moth_tri2d
{
    moth_p2d p1{};
    moth_p2d p2{};
    moth_p2d p3{};

public:

};  // struct moth_tri2d

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Triangle in 3D space.
 */
struct MOTH_CORE moth_tri3d
{
    moth_p3d p1{};
    moth_p3d p2{};
    moth_p3d p3{};

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t area(const moth_tri3d& t1)
    {
        moth_p3d n{moth_p3d::cross(t1.p2 - t1.p1, t1.p3 - t1.p1)};
        return 0.5 * moth_p3d::len(n);
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_p3d normal(const moth_tri3d& t1)
    {
        moth_p3d n{moth_p3d::cross(t1.p2 - t1.p1, t1.p3 - t1.p1)};
        moth_real_t l{moth_p3d::len(n)};
        if (l != 0.0) {
            return n / l;
        } else {
            return {0.0, 0.0, 0.0};
        }
    }
};  // struct moth_tri3d
