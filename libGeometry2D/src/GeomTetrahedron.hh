#pragma once

#include "libGeometry2D/src/GeomTriangle.hh"

/**
 * Tetrahedron in 3D space.
 */
struct MOTH_CORE moth_tetrahedron
{
public:
    moth_p3d p1;
    moth_p3d p2;
    moth_p3d p3;
    moth_p3d p4;

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t volume(const moth_tetrahedron& T)
    {
        moth_real_t v{moth_p3d::mixed(T.p1 - T.p4, T.p2 - T.p4, T.p3 - T.p4) / 6.0};
        return std::fabs(v);
    }
};  // struct moth_tetrahedron
