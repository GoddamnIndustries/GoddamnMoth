#pragma once

#include "GeomBase.hh"
#include "GeomPoint.hh"

/**
 * Edge (half-open interval) in 2D space.
 */
struct MOTH_CORE moth_e2d final
{
	moth_p2d s{};
	moth_p2d t{};

public:
    MOTH_HOST MOTH_DEVICE
    bool operator==(const moth_e2d& e) const
    {
        return (s == e.s) && (t == e.t);
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_e2d& e) const
    {
        return (s != e.s) || (t != e.t);
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_e2d operator+(const moth_p2d& p) const
    {
        return {s + p, t + p};
    }
    MOTH_HOST MOTH_DEVICE
    moth_e2d& operator+=(const moth_p2d& p)
    {
        return *this = *this + p;
    }

    MOTH_HOST MOTH_DEVICE
    moth_e2d operator-(const moth_p2d& p) const
    {
        return {s - p, t - p};
    }
    MOTH_HOST MOTH_DEVICE
    moth_e2d& operator-=(const moth_p2d& p)
    {
        return *this = *this - p;
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t dot(const moth_e2d& e1, const moth_e2d& e2)
    {
        return moth_p2d::dot(e1.t - e1.s, e2.t - e2.s);
    }
    MOTH_HOST MOTH_DEVICE
    static moth_real_t len(const moth_e2d& e1)
    {
        return moth_p2d::len(e1.t - e1.s);
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_p2d normal(const moth_e2d& e1)
    {
        return moth_p2d::normal(e1.t - e1.s);
    }
    MOTH_HOST MOTH_DEVICE
    static moth_real_t angle(const moth_e2d& e1, const moth_e2d& e2)
    {
        return moth_p2d::angle(e1.t - e1.s, e2.t - e2.s);
    }

public:
    /**
     * Calculate intersection of the two edges.
     */
    MOTH_HOST MOTH_DEVICE
    static bool intersect(const moth_e2d& e1, moth_e2d e2, moth_e2d& e)
    {
        moth_p2d v1 = e1.t - e1.s;
        moth_p2d v2 = e2.t - e2.s;
        moth_real_t det = moth_p2d::det(v1, v2);
        if(det != 0.0) {
            /* Edges are not collinear. */
            moth_p2d v = e2.s - e1.s;
            moth_real_t alpha = moth_p2d::det(v, v2) / det;
            moth_real_t gamma = moth_p2d::det(v, v1) / det;
            if ((0.0 < alpha && alpha <= 1.0) && (0.0 < gamma && gamma <= 1.0)) {
                /* Edges intersect on point. */
                e.s = e1.s + alpha * v1;
                e.t = e.s;
                return true;
            }
        } else {
            /* Edges are collinear. */
            if (moth_p2d::dot(v1, v2) < 0.0) {
                std::swap(e2.s, e2.t);
            }

            moth_p2d q1 = e2.s - e1.t;
            moth_p2d q2 = e2.t - e1.s;
            if (moth_p2d::det(q1, q2) == 0.0) {
                /* Edges are on one line. */
                if (moth_p2d::len(v1) + moth_p2d::len(v2) != fabs(moth_p2d::len(q2) - moth_p2d::len(q1))) {
                    /* Edges intersect on line. */
                    e.s = moth_p2d::max(moth_p2d::min(e2.s, e2.t), moth_p2d::min(e1.s, e1.t));
                    e.t = moth_p2d::min(moth_p2d::max(e2.s, e2.t), moth_p2d::max(e1.s, e1.t));
                    return true;
                }
            }
        }
        return false;
    }
    MOTH_HOST MOTH_DEVICE
    static bool intersect(const moth_e2d& e1, moth_e2d e2, moth_p2d& p)
    {
        moth_e2d e;
        if (intersect(e1, e2, e)) {
            if (e.s == e.t) {
                p = e.s;
                return true;
            }
        }
        return false;
    }

public:
    /**
     * Reflect trajectory on the diode wall.
     */
    MOTH_HOST MOTH_DEVICE [[deprecated]]
    static bool reflect(const moth_e2d& e1, const moth_e2d& e2, moth_e2d& e)
    {
        if (intersect(e1, e2, e) && e.s == e.t) {
            /* Edges intersect and are not collinear. */
            moth_p2d p = e.s;
            moth_real_t l1 = moth_p2d::len(e2.s - p);
            moth_real_t l2 = moth_p2d::len(e2.t - p);
            moth_p2d v = (e2.s - p) / l1;
            moth_p2d n = moth_p2d::normal(e1.s - p);
            moth_real_t sin = moth_p2d::det(v, n);
            moth_real_t cos = moth_p2d::dot(v, n);
            if (cos > 0.0) {
                /* Edge reflect from the opaque side of the diode wall. */
                moth_p2d q = moth_p2d::rotate(n, sin, cos);
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
    MOTH_HOST
    static std::string str(const moth_e2d& e);
};	// struct moth_e2d

MOTH_HOST MOTH_CORE
extern std::ostream& operator<<(std::ostream& stream, const moth_e2d& e);
MOTH_HOST MOTH_CORE
extern std::istream& operator>>(std::istream& stream, moth_e2d& e);

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Edge (half-open interval) in 2D space.
 */
struct MOTH_CORE moth_e3d final
{
    moth_p3d s{};
    moth_p3d t{};

public:
    MOTH_HOST MOTH_DEVICE
    bool operator==(const moth_e3d& e) const
    {
        return (s == e.s) && (t == e.t);
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_e3d& e) const
    {
        return (s != e.s) || (t != e.t);
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_e3d operator+(const moth_p3d& p) const
    {
        return {s + p, t + p};
    }
    MOTH_HOST MOTH_DEVICE
    moth_e3d& operator+=(const moth_p3d& p)
    {
        return *this = *this + p;
    }

    MOTH_HOST MOTH_DEVICE
    moth_e3d operator-(const moth_p3d& p) const
    {
        return {s - p, t - p};
    }
    MOTH_HOST MOTH_DEVICE
    moth_e3d& operator-=(const moth_p3d& p)
    {
        return *this = *this - p;
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t dot(const moth_e3d& e1, const moth_e3d& e2)
    {
        return moth_p3d::dot(e1.t - e1.s, e2.t - e2.s);
    }
    MOTH_HOST MOTH_DEVICE
    static moth_real_t len(const moth_e3d& e1)
    {
        return moth_p3d::len(e1.t - e1.s);
    }

public:
    /*@todo angle*/

public:
    /*@todo intersect*/

public:
    /*@todo reflect*/

public:
    MOTH_HOST
    static std::string str(const moth_e3d& e);
};  // struct moth_e3d

MOTH_HOST MOTH_CORE
extern std::ostream& operator<<(std::ostream& stream, const moth_e3d& e);
MOTH_HOST MOTH_CORE
extern std::istream& operator>>(std::istream& stream, moth_e3d& e);
