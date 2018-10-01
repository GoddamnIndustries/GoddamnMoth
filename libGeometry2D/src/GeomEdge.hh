#pragma once

#include "GeomBase.hh"
#include "GeomPoint.hh"

///
/// Edge (half-open interval) in 2D space:
/// @f$
/// \hat{e_1}(\vec{s}, \vec{t}) := {\vec{p}: \vec{p} = \alpha\vec{s} + (1-\alpha)\vec{t}, 0 < \alpha \le 1}.
/// @f$
///
struct GEOM_CORE geom_e2d final
{
	geom_p2d s{};
	geom_p2d t{};

public:
    GEOM_HOST GEOM_DEVICE
    bool operator==(const geom_e2d& e) const
    {
        return (s == e.s) && (t == e.t);
    }
    GEOM_HOST GEOM_DEVICE
    bool operator!=(const geom_e2d& e) const
    {
        return (s != e.s) || (t != e.t);
    }

public:
    GEOM_HOST GEOM_DEVICE
    geom_e2d operator+(const geom_p2d& p) const
    {
        return {s + p, t + p};
    }
    GEOM_HOST GEOM_DEVICE
    geom_e2d& operator+=(const geom_p2d& p)
    {
        return *this = *this + p;
    }

    GEOM_HOST GEOM_DEVICE
    geom_e2d operator-(const geom_p2d& p) const
    {
        return {s - p, t - p};
    }
    GEOM_HOST GEOM_DEVICE
    geom_e2d& operator-=(const geom_p2d& p)
    {
        return *this = *this - p;
    }

public:
    GEOM_HOST GEOM_DEVICE
    static geom_real_t dot(const geom_e2d& e1, const geom_e2d& e2)
    {
        return geom_p2d::dot(e1.t - e1.s, e2.t - e2.s);
    }
    GEOM_HOST GEOM_DEVICE
    static geom_real_t len(const geom_e2d& e1)
    {
        return geom_p2d::len(e1.t - e1.s);
    }

public:
    GEOM_HOST GEOM_DEVICE
    static geom_p2d normal(const geom_e2d& e1)
    {
        return geom_p2d::normal(e1.t - e1.s);
    }
    GEOM_HOST GEOM_DEVICE
    static geom_real_t angle(const geom_e2d& e1, const geom_e2d& e2)
    {
        return geom_p2d::angle(e1.t - e1.s, e2.t - e2.s);
    }

public:
    GEOM_HOST GEOM_DEVICE
    static bool intersect(const geom_e2d& e1, geom_e2d e2, geom_e2d& e)
    {
        geom_p2d v1 = e1.t - e1.s;
        geom_p2d v2 = e2.t - e2.s;
        geom_real_t det = geom_p2d::det(v1, v2);
        if(det != 0.0) {
            /* Edges are not collinear. */
            geom_p2d v = e2.s - e1.s;
            geom_real_t alpha = geom_p2d::det(v, v2) / det;
            geom_real_t gamma = geom_p2d::det(v, v1) / det;
            if ((0.0 < alpha && alpha <= 1.0) && (0.0 < gamma && gamma <= 1.0)) {
                /* Edges intersect on point. */
                e.s = e1.s + alpha * v1;
                e.t = e.s;
                return true;
            }
        } else {
            /* Edges are collinear. */
            if (geom_p2d::dot(v1, v2) < 0.0) {
                std::swap(e2.s, e2.t);
            }

            geom_p2d q1 = e2.s - e1.t;
            geom_p2d q2 = e2.t - e1.s;
            if (geom_p2d::det(q1, q2) == 0.0) {
                /* Edges are on one line. */
                if (geom_p2d::len(v1) + geom_p2d::len(v2) != fabs(geom_p2d::len(q2) - geom_p2d::len(q1))) {
                    /* Edges intersect on line. */
                    e.s = geom_p2d::max(geom_p2d::min(e2.s, e2.t), geom_p2d::min(e1.s, e1.t));
                    e.t = geom_p2d::min(geom_p2d::max(e2.s, e2.t), geom_p2d::max(e1.s, e1.t));
                    return true;
                }
            }
        }
        return false;
    }

public:
    GEOM_HOST GEOM_DEVICE
    static bool reflect(const geom_e2d& e1, const geom_e2d& e2, geom_e2d& e)
    {
        if (intersect(e1, e2, e) && e.s == e.t) {
            /* Edges intersect and are not collinear. */
            geom_p2d p = e.s;
            geom_real_t l1 = geom_p2d::len(e2.s - p);
            geom_real_t l2 = geom_p2d::len(e2.t - p);
            geom_p2d v = (e2.s - p) / l1;
            geom_p2d n = geom_p2d::normal(e1.s - p);
            geom_real_t sin = geom_p2d::det(v, n);
            geom_real_t cos = geom_p2d::dot(v, n);
            if (cos > 0.0) {
                /* Edge reflect from the opaque side of the diode wall. */
                geom_p2d q = geom_p2d::rotate(n, sin, cos);
                e = {e2.s, p + l2 * q};
            } else {
                /* Edge goes through the transparent side of the diode wall. */
                e = {e2.s, e2.t};
            }
            return true;
        } else {
            e = {e2.s, e2.t};
        }
        return false;
    }

public:
    GEOM_HOST
    static std::string str(const geom_e2d& e);
};	// struct geom_e2d

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

GEOM_HOST GEOM_CORE
std::ostream& operator<<(std::ostream& stream, const geom_e2d& e);
GEOM_HOST GEOM_CORE
std::istream& operator>>(std::istream& stream, geom_e2d& e);
