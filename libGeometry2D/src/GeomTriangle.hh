#pragma once

#include "libGeometry2D/src/GeomPoint.hh"

/**
 * Triangle in 2D space.
 */
struct MOTH_CORE moth_tri2d
{
    moth_p2d p1{};
    moth_p2d p2{};
    moth_p2d p3{};

public:
    MOTH_HOST MOTH_DEVICE
    static moth_p2d circumcenter(const moth_tri2d& t1)
    {
        moth_p3d a = {t1.p1.x, t1.p1.y}, b = {t1.p2.x,t1.p2.y}, c = {t1.p3.x,t1.p3.y};
        moth_p3d ac = c - a;
        moth_p3d ab = b - a;
        moth_p3d ab_ac = moth_p3d::cross(ab, ac);
        moth_p3d to_cc = moth_p3d::cross(ab_ac, ab) * moth_p3d::dot(ac, ac) +
                         moth_p3d::cross(ac, ab_ac) * moth_p3d::dot(ab, ab);
        to_cc /= 2.0 * moth_p3d::dot(ab_ac, ab_ac);
        return a + to_cc;
    }

public:
    MOTH_HOST MOTH_DEVICE
    static bool circle(const moth_tri2d& t1, const moth_p2d& p1)
    {
        moth_real_t a13 = (t1.p1.x*t1.p1.x - p1.x*p1.x) + (t1.p1.y*t1.p1.y - p1.y*p1.y);
        moth_real_t a23 = (t1.p2.x*t1.p2.x - p1.x*p1.x) + (t1.p2.y*t1.p2.y - p1.y*p1.y);
        moth_real_t a33 = (t1.p3.x*t1.p3.x - p1.x*p1.x) + (t1.p3.y*t1.p3.y - p1.y*p1.y);
        moth_real_t det{moth_det(t1.p1.x - p1.x, t1.p1.y - p1.y, a13,
                                 t1.p2.x - p1.x, t1.p2.y - p1.y, a23,
                                 t1.p3.x - p1.x, t1.p3.y - p1.y, a33)};
        return det > 0;
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t len(const moth_tri2d& t1)
    {
        moth_real_t l{};
        l += moth_p2d::len(t1.p1 - t1.p2);
        l += moth_p2d::len(t1.p2 - t1.p3);
        l += moth_p2d::len(t1.p3 - t1.p1);
        return l;
    }
    MOTH_HOST MOTH_DEVICE
    static moth_real_t area(const moth_tri2d& t1)
    {
        moth_real_t n{moth_p2d::det(t1.p2 - t1.p1, t1.p3 - t1.p1)};
        return 0.5 * n;
    }
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
    static moth_p3d barycenter(const moth_tri3d& t1)
    {
        return (t1.p1 + t1.p2 + t1.p3)/3.0;
    }
    MOTH_HOST MOTH_DEVICE
    static moth_p3d circumcenter(const moth_tri3d& t1)
    {
        moth_p3d a = t1.p1, b = t1.p2, c = t1.p3;
        moth_p3d ac = c - a;
        moth_p3d ab = b - a;
        moth_p3d ab_ac = moth_p3d::cross(ab, ac);
        moth_p3d to_cc = moth_p3d::cross(ab_ac, ab) * moth_p3d::dot(ac, ac) +
                         moth_p3d::cross(ac, ab_ac) * moth_p3d::dot(ab, ab);
        to_cc /= 2.0 * moth_p3d::dot(ab_ac, ab_ac);
        return a + to_cc;
    }

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
            n /= l;
            return n;
        } else {
            return {0.0, 0.0, 0.0};
        }
    }
};  // struct moth_tri3d

