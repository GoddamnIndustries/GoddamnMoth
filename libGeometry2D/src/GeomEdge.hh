#pragma once

#include "GeomBase.hh"
#include "GeomPoint.hh"

/**
 * Edge (half-open interval) in 2D space.
 */
struct MOTH_CORE moth_e2d final
{
	moth_p2d p1{};
	moth_p2d p2{};

public:
    MOTH_HOST MOTH_DEVICE
    moth_p2d vec() const
    {
        return p2 - p1;
    }

public:
    MOTH_HOST MOTH_DEVICE
    bool operator==(const moth_e2d& e) const
    {
        return (p1 == e.p1) && (p2 == e.p2);
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_e2d& e) const
    {
        return (p1 != e.p1) || (p2 != e.p2);
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_e2d operator+(const moth_p2d& p) const
    {
        return {p1 + p, p2 + p};
    }
    MOTH_HOST MOTH_DEVICE
    moth_e2d& operator+=(const moth_p2d& p)
    {
        return *this = *this + p;
    }

    MOTH_HOST MOTH_DEVICE
    moth_e2d operator-(const moth_p2d& p) const
    {
        return {p1 - p, p2 - p};
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
        moth_p2d p1{e1.vec()};
        moth_p2d p2{e2.vec()};
        return moth_p2d::dot(p1, p2);
    }
    MOTH_HOST MOTH_DEVICE
    static moth_real_t len(const moth_e2d& e1)
    {
        moth_p2d p1{e1.vec()};
        return moth_p2d::len(p1);
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_p2d cross(const moth_e2d& e1, const moth_e2d& e2, const moth_e2d& e3)
    {
        moth_p2d p1{e1.vec()};
        moth_p2d p2{e2.vec()};
        moth_p2d p3{e3.vec()};
        return moth_p2d::cross(p1, p2, p3);
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t det(const moth_e2d& e1, const moth_e2d& e2)
    {
        moth_p2d p1{e1.vec()};
        moth_p2d p2{e2.vec()};
        return moth_p2d::det(p1, p2);
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t angle(const moth_e2d& e1, const moth_e2d& e2)
    {
        moth_p2d p1{e1.vec()};
        moth_p2d p2{e2.vec()};
        return moth_p2d::angle(p1, p2);
    }

    MOTH_HOST MOTH_DEVICE
    static moth_e2d rotate(const moth_e2d& e, const moth_p2d& p, moth_real_t sin, moth_real_t cos)
    {
        moth_p2d p1 = e.p1 - p;
        moth_p2d p2 = e.p2 - p;
        moth_e2d rot{moth_p2d::rotate(p1, sin, cos),
                     moth_p2d::rotate(p2, sin, cos)};
        rot += p;
        return rot;
    }
    MOTH_HOST MOTH_DEVICE
    static moth_e2d rotate(const moth_e2d& e, const moth_p2d& p, moth_radians_t theta)
    {
        return rotate(e, p, std::sin(theta), std::cos(theta));
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_p2d normal(const moth_e2d& e)
    {
        moth_p2d p = e.vec();
        moth_p2d n{-p.y, p.x};
        moth_real_t l{moth_p2d::len(n)};
        if (l != 0.0) {
            return n / l;
        } else {
            std::cerr << "Warning: normal to singular edge." << std::endl;
            return {0.0, 0.0};
        }
    }

public:
    MOTH_HOST MOTH_DEVICE
    static bool intersect(const moth_e2d& e1, moth_e2d e2, moth_e2d& e)
    {
        moth_p2d v1{e1.vec()};
        moth_p2d v2{e2.vec()};
        moth_real_t det = moth_p2d::det(v1, v2);
        if(det != 0.0) {
            /* Edges are not collinear. */
            moth_p2d v = e2.p1 - e1.p1;
            moth_real_t alpha = moth_p2d::det(v, v2) / det;
            moth_real_t gamma = moth_p2d::det(v, v1) / det;
            if ((0.0 < alpha && alpha <= 1.0) && (0.0 < gamma && gamma <= 1.0)) {
                /* Edges intersect on point. */
                e.p1 = e1.p1 + alpha * v1;
                e.p2 = e.p1;
                return true;
            }
        } else {
            /* Edges are collinear. */
            if (moth_p2d::dot(v1, v2) < 0.0) {
                std::swap(e2.p1, e2.p2);
            }

            moth_p2d q1 = e2.p1 - e1.p2;
            moth_p2d q2 = e2.p2 - e1.p1;
            if (moth_p2d::det(q1, q2) == 0.0) {
                /* Edges are on one line. */
                if (moth_p2d::len(v1) + moth_p2d::len(v2) != std::fabs(moth_p2d::len(q2) - moth_p2d::len(q1))) {
                    /* Edges intersect on line. */
                    e.p1 = moth_p2d::max(moth_p2d::min(e2.p1, e2.p2), moth_p2d::min(e1.p1, e1.p2));
                    e.p2 = moth_p2d::min(moth_p2d::max(e2.p1, e2.p2), moth_p2d::max(e1.p1, e1.p2));
                    return true;
                }
            }
        }
        return false;
    }
    MOTH_HOST MOTH_DEVICE
    static bool intersect(const moth_e2d& e1, const moth_e2d& e2, moth_p2d& p)
    {
        moth_e2d e;
        if (intersect(e1, e2, e)) {
            if (e.p1 == e.p2) {
                p = e.p1;
                return true;
            }
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
    moth_p3d p1{};
    moth_p3d p2{};

public:
    MOTH_HOST MOTH_DEVICE
    moth_p3d vec() const
    {
        return p2 - p1;
    }

public:
    MOTH_HOST MOTH_DEVICE
    bool operator==(const moth_e3d& e) const
    {
        return (p1 == e.p1) && (p2 == e.p2);
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_e3d& e) const
    {
        return (p1 != e.p1) || (p2 != e.p2);
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_e3d operator+(const moth_p3d& p) const
    {
        return {p1 + p, p2 + p};
    }
    MOTH_HOST MOTH_DEVICE
    moth_e3d& operator+=(const moth_p3d& p)
    {
        return *this = *this + p;
    }

    MOTH_HOST MOTH_DEVICE
    moth_e3d operator-(const moth_p3d& p) const
    {
        return {p1 - p, p2 - p};
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
        moth_p3d p1{e1.vec()};
        moth_p3d p2{e2.vec()};
        return moth_p3d::dot(p1, p2);
    }
    MOTH_HOST MOTH_DEVICE
    static moth_real_t len(const moth_e3d& e1)
    {
        moth_p3d p1{e1.vec()};
        return moth_p3d::len(p1);
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_p3d cross(const moth_e3d& e1, const moth_e3d& e2)
    {
        moth_p3d p1{e1.vec()};
        moth_p3d p2{e2.vec()};
        return moth_p3d::cross(p1, p2);
    }
    MOTH_HOST MOTH_DEVICE
    static moth_p3d cross(const moth_e3d& e1, const moth_e3d& e2, const moth_e3d& e3)
    {
        moth_p3d p1{e1.vec()};
        moth_p3d p2{e2.vec()};
        moth_p3d p3{e3.vec()};
        return moth_p3d::cross(p1, p2, p3);
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t mixed(const moth_e3d& e1, const moth_e3d& e2, const moth_e3d& e3)
    {
        moth_p3d p1{e1.vec()};
        moth_p3d p2{e2.vec()};
        moth_p3d p3{e3.vec()};
        return moth_p3d::mixed(p1, p2, p3);
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t det(const moth_e3d& e1, const moth_e3d& e2)
    {
        moth_p3d p1{e1.vec()};
        moth_p3d p2{e2.vec()};
        return moth_p3d::det(p1, p2);
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t angle(const moth_e3d& e1, const moth_e3d& e2)
    {
        moth_p3d p1{e1.vec()};
        moth_p3d p2{e2.vec()};
        return moth_p3d::angle(p1, p2);
    }

    MOTH_HOST MOTH_DEVICE
    static moth_e3d rotate(const moth_e3d& e, const moth_p3d& p, const moth_p3d& u, moth_real_t sin, moth_real_t cos)
    {
        moth_p3d p1 = e.p1 - p;
        moth_p3d p2 = e.p2 - p;
        moth_e3d rot{moth_p3d::rotate(p1, u, sin, cos),
                     moth_p3d::rotate(p2, u, sin, cos)};
        rot += p;
        return rot;
    }
    MOTH_HOST MOTH_DEVICE
    static moth_e3d rotate(const moth_e3d& e, const moth_p3d& p, const moth_p3d& u, moth_radians_t theta)
    {
        return rotate(e, p, u, std::sin(theta), std::cos(theta));
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_p3d normal(const moth_e3d& e1, const moth_e3d& e2)
    {
        moth_p3d n{cross(e1, e2)};
        moth_real_t l{moth_p2d::len(n)};
        if (l != 0.0) {
            return n / l;
        } else {
            std::cerr << "Warning: normal to singular edges." << std::endl;
            return {0.0, 0.0};
        }
    }

public:
    MOTH_HOST MOTH_DEVICE
    static bool intersect(const moth_e3d& e1, moth_e3d e2, moth_e3d& e)
    {
        moth_p3d v1{e1.vec()};
        moth_p3d v2{e2.vec()};
        moth_real_t det = moth_p2d::det(v1, v2);
        if(det != 0.0) {
            /* Edges are not collinear. */
            moth_p3d v = e2.p1 - e1.p1;
            moth_real_t alpha = moth_p3d::det(v, v2) / det;
            moth_real_t gamma = moth_p3d::det(v, v1) / det;
            if ((0.0 < alpha && alpha <= 1.0) && (0.0 < gamma && gamma <= 1.0)) {
                /* Edges intersect on point. */
                e.p1 = e1.p1 + alpha * v1;
                e.p2 = e.p1;
                return true;
            }
        } else {
            /* Edges are collinear. */
            if (moth_p3d::dot(v1, v2) < 0.0) {
                std::swap(e2.p1, e2.p2);
            }

            moth_p3d q1 = e2.p1 - e1.p2;
            moth_p3d q2 = e2.p2 - e1.p1;
            if (moth_p3d::det(q1, q2) == 0.0) {
                /* Edges are on one line. */
                if (moth_p3d::len(v1) + moth_p2d::len(v2) != std::fabs(moth_p2d::len(q2) - moth_p2d::len(q1))) {
                    /* Edges intersect on line. */
                    e.p1 = moth_p3d::max(moth_p3d::min(e2.p1, e2.p2), moth_p3d::min(e1.p1, e1.p2));
                    e.p2 = moth_p3d::min(moth_p3d::max(e2.p1, e2.p2), moth_p3d::max(e1.p1, e1.p2));
                    return true;
                }
            }
        }
        return false;
    }
    MOTH_HOST MOTH_DEVICE
    static bool intersect(const moth_e3d& e1, const moth_e3d& e2, moth_p3d& p)
    {
        moth_e3d e;
        if (intersect(e1, e2, e)) {
            if (e.p1 == e.p2) {
                p = e.p1;
                return true;
            }
        }
        return false;
    }

public:
    MOTH_HOST
    static std::string str(const moth_e3d& e);
};  // struct moth_e3d

MOTH_HOST MOTH_CORE
extern std::ostream& operator<<(std::ostream& stream, const moth_e3d& e);
MOTH_HOST MOTH_CORE
extern std::istream& operator>>(std::istream& stream, moth_e3d& e);
