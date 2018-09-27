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
    static bool intersect(const geom_e2d& e1, geom_e2d e2, geom_e2d& e);

    GEOM_HOST GEOM_DEVICE
    static bool reflect(const geom_e2d& e1, const geom_e2d& e2, geom_e2d& e);

public:
    /// String representation of the edge
    /// acceptable by the scripting engine.
    GEOM_HOST
    static std::string str(const geom_e2d& e);
};	// struct geom_e2d

GEOM_HOST GEOM_CORE
std::ostream& operator<<(std::ostream& stream, const geom_e2d& e);
GEOM_HOST GEOM_CORE
std::istream& operator>>(std::istream& stream, geom_e2d& e);
