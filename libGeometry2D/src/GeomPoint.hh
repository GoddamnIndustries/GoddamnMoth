#pragma once
#include "GeomBase.hh"

///
/// Point in 2D space.
///
struct GEOM_CORE geom_p2d
{
    geom_real_t x{}, y{};
    geom_real_t u{}, v{};

public:
    GEOM_HOST GEOM_DEVICE
    bool operator==(const geom_p2d& p1) const
    {
        return (x == p1.x) && (y == p1.y);
    }
    GEOM_HOST GEOM_DEVICE
    bool operator!=(const geom_p2d& p1) const
    {
        return (x != p1.x) || (y != p1.y);
    }

public:
    GEOM_HOST GEOM_DEVICE
    geom_p2d operator+() const
    {
        return *this;
    }
    GEOM_HOST GEOM_DEVICE
    geom_p2d operator+(const geom_p2d& p1) const
    {
        return {x + p1.x, y + p1.y};
    }
    GEOM_HOST GEOM_DEVICE
    geom_p2d& operator+=(const geom_p2d& p1)
    {
        return *this = *this + p1;
    }

    GEOM_HOST GEOM_DEVICE
    geom_p2d operator-() const
    {
        return {-x, -y};
    }
    GEOM_HOST GEOM_DEVICE
    geom_p2d operator-(const geom_p2d& p1) const
    {
        return {x - p1.x, y - p1.y};
    }
    GEOM_HOST GEOM_DEVICE
    geom_p2d& operator-=(const geom_p2d& p1)
    {
        return *this = *this - p1;
    }

public:
    GEOM_HOST GEOM_DEVICE
    geom_p2d operator*(geom_real_t a) const
    {
        return {a * x, a * y};
    }
    GEOM_HOST GEOM_DEVICE
    geom_p2d& operator*=(geom_real_t a)
    {
        return *this = a * *this;
    }
    GEOM_HOST GEOM_DEVICE
    friend geom_p2d operator*(geom_real_t a, const geom_p2d& p1)
    {
        return p1 * a;
    }

    GEOM_HOST GEOM_DEVICE
    geom_p2d operator/(geom_real_t a) const
    {
        return {x / a, y / a};
    }
    GEOM_HOST GEOM_DEVICE
    geom_p2d& operator/=(geom_real_t a)
    {
        return *this = *this / a;
    }
    GEOM_HOST GEOM_DEVICE
    friend geom_p2d operator/(geom_real_t a, const geom_p2d& p1)
    {
        return p1 / a;
    }

public:
    GEOM_HOST GEOM_DEVICE
    static geom_real_t dot(const geom_p2d& p1, const geom_p2d& p2)
    {
        return p1.x * p2.x + p1.y * p2.y;
    }
    GEOM_HOST GEOM_DEVICE
    static geom_real_t len(const geom_p2d& p1)
    {
        return std::hypot(p1.x, p1.y);
    }

public:
    GEOM_HOST GEOM_DEVICE
    static geom_real_t det(const geom_p2d& p1, const geom_p2d& p2)
    {
        return p1.x * p2.y - p2.x * p1.y;
    }

public:
    GEOM_HOST GEOM_DEVICE
    static geom_p2d rotate(const geom_p2d& n, geom_real_t sin, geom_real_t cos)
    {
        return {n.x * cos - n.y * sin, n.x * sin + n.y * cos};
    }
    GEOM_HOST GEOM_DEVICE
    static geom_p2d rotate(const geom_p2d& n, geom_radians_t theta)
    {
        return rotate(n, std::sin(theta), std::cos(theta));
    }

public:
    GEOM_HOST GEOM_DEVICE
    static geom_p2d normal(const geom_p2d& p1)
    {
        geom_p2d n{-p1.y, p1.x};
        geom_real_t l = len(n);
        if (l != 0.0) {
            return n / l;
        } else {
            return {0.0, 0.0};
        }
    }
    GEOM_HOST GEOM_DEVICE
    static geom_radians_t angle(const geom_p2d& p1, const geom_p2d& p2)
    {
        return std::atan2(det(p1, p2), dot(p1, p2));
    }

public:
    GEOM_HOST GEOM_DEVICE
    static geom_p2d max(const geom_p2d& p1, const geom_p2d& p2)
    {
        return {std::max(p1.x, p2.x), std::max(p1.y, p2.y)};
    }
    GEOM_HOST GEOM_DEVICE
    static geom_p2d min(const geom_p2d& p1, const geom_p2d& p2)
    {
        return {std::min(p1.x, p2.x), std::min(p1.y, p2.y)};
    }

public:
    /// String representation of the edge
    /// acceptable by the scripting engine.
    GEOM_HOST
    static std::string str(const geom_p2d& p1);
};	// struct geom_p2d

GEOM_HOST GEOM_CORE
std::ostream& operator<<(std::ostream& stream, const geom_p2d& p);
GEOM_HOST GEOM_CORE
std::istream& operator>>(std::istream& stream, geom_p2d& p);
