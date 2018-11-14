#pragma once

#include "libGeometry2D/src/GeomBase.hh"

/**
 * Point in 2D space.
 */
struct MOTH_CORE moth_p2d
{
    moth_real_t x{};
    moth_real_t y{};

public:
    MOTH_HOST MOTH_DEVICE
    bool operator==(const moth_p2d& p) const
    {
        return (x == p.x) && (y == p.y);
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_p2d& p) const
    {
        return (x != p.x) || (y != p.y);
    }

    MOTH_HOST MOTH_DEVICE
    bool operator>(const moth_p2d& p) const
    {
        if (x == p.x) {
            return y > p.y;
        } else {
            return x > p.x;
        }
    }
    MOTH_HOST MOTH_DEVICE
    bool operator>=(const moth_p2d& p) const
    {
        return (*this == p) || (*this > p);
    }

    MOTH_HOST MOTH_DEVICE
    bool operator<(const moth_p2d& p) const
    {
        if (x == p.x) {
            return y < p.y;
        } else {
            return x < p.x;
        }
    }
    MOTH_HOST MOTH_DEVICE
    bool operator<=(const moth_p2d& p) const
    {
        return (*this == p) || (*this < p);
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_p2d operator+() const
    {
        return {+x, +y};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p2d operator+(const moth_p2d& p) const
    {
        return {x + p.x, y + p.y};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p2d& operator+=(const moth_p2d& p)
    {
        return *this = *this + p;
    }

    MOTH_HOST MOTH_DEVICE
    moth_p2d operator-() const
    {
        return {-x, -y};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p2d operator-(const moth_p2d& p) const
    {
        return {x - p.x, y - p.y};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p2d& operator-=(const moth_p2d& p)
    {
        return *this = *this - p;
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_p2d operator*(moth_real_t a) const
    {
        return {a * x, a * y};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p2d& operator*=(moth_real_t a)
    {
        return *this = *this * a;
    }
    MOTH_HOST MOTH_DEVICE
    friend moth_p2d operator*(moth_real_t a, const moth_p2d& p)
    {
        return p * a;
    }

    MOTH_HOST MOTH_DEVICE
    moth_p2d operator/(moth_real_t a) const
    {
        return {x / a, y / a};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p2d& operator/=(moth_real_t a)
    {
        return *this = *this / a;
    }
    MOTH_HOST MOTH_DEVICE
    friend moth_p2d operator/(moth_real_t a, const moth_p2d& p)
    {
        return p / a;
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_p2d max(const moth_p2d& p1, const moth_p2d& p2)
    {
        return {std::max(p1.x, p2.x), std::max(p1.y, p2.y)};
    }
    MOTH_HOST MOTH_DEVICE
    static moth_p2d min(const moth_p2d& p1, const moth_p2d& p2)
    {
        return {std::min(p1.x, p2.x), std::min(p1.y, p2.y)};
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t dot(const moth_p2d& p1, const moth_p2d& p2)
    {
        return p1.x * p2.x + p1.y * p2.y;
    }
    MOTH_HOST MOTH_DEVICE
    static moth_real_t len(const moth_p2d& p1)
    {
        return std::hypot(p1.x, p1.y);
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_p2d cross(const moth_p2d& p1, const moth_p2d& p2, const moth_p2d& p3)
    {
        moth_p2d n{-p3.y, p3.x};
        moth_real_t l{len(n)};
        if (l != 0.0) {
            n /= l;
            n *= det(p1, p2);
            return n / l;
        } else {
            return {0.0, 0.0};
        }
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t det(const moth_p2d& p1, const moth_p2d& p2)
    {
        return p1.x * p2.y - p2.x * p1.y;
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_radians_t angle(const moth_p2d& p1, const moth_p2d& p2)
    {
        //moth_real_t sin{det(p1, p2)};
        moth_real_t cos{dot(p1, p2)/ len(p1) / len(p2)};
        return std::acos(cos);
        //return std::atan2(sin, cos);
    }

    MOTH_HOST MOTH_DEVICE
    static moth_p2d rotate(const moth_p2d& p, moth_real_t sin, moth_real_t cos)
    {
        moth_p2d rot{p.x * cos - p.y * sin,
                     p.x * sin + p.y * cos};
        return rot;
    }
    MOTH_HOST MOTH_DEVICE
    static moth_p2d rotate(const moth_p2d& p, moth_radians_t theta)
    {
        return rotate(p, std::sin(theta), std::cos(theta));
    }

public:
    MOTH_HOST
    static std::string str(const moth_p2d& p1);
};	// struct moth_p2d

MOTH_HOST MOTH_CORE
extern std::ostream& operator<<(std::ostream& stream, const moth_p2d& p);
MOTH_HOST MOTH_CORE
extern std::istream& operator>>(std::istream& stream, moth_p2d& p);

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Point in 3D space.
 */
struct MOTH_CORE moth_p3d final
{
    moth_real_t x{};
    moth_real_t y{};
    moth_real_t z{};

public:
    MOTH_HOST MOTH_DEVICE
    operator moth_p2d() const
    {
        return {x, y};
    }

public:
    MOTH_HOST MOTH_DEVICE
    bool operator==(const moth_p3d& p) const
    {
        return (x == p.x) && (y == p.y) && (z == p.z);
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_p3d& p) const
    {
        return (x != p.x) || (y != p.y) || (z != p.z);
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_p3d operator+() const
    {
        return {+x, +y, +z};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p3d operator+(const moth_p3d& p) const
    {
        return {x + p.x, y + p.y, z + p.z};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p3d& operator+=(const moth_p3d& p)
    {
        return *this = *this + p;
    }

    MOTH_HOST MOTH_DEVICE
    moth_p3d operator-() const
    {
        return {-x, -y, -z};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p3d operator-(const moth_p3d& p) const
    {
        return {x - p.x, y - p.y, z - p.z};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p3d& operator-=(const moth_p3d& p)
    {
        return *this = *this - p;
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_p3d operator*(moth_real_t a) const
    {
        return {a * x, a * y, a * z};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p3d& operator*=(moth_real_t a)
    {
        return *this = *this * a;
    }
    MOTH_HOST MOTH_DEVICE
    friend moth_p3d operator*(moth_real_t a, const moth_p3d& p)
    {
        return p * a;
    }

    MOTH_HOST MOTH_DEVICE
    moth_p3d operator/(moth_real_t a) const
    {
        return {x / a, y / a, z / a};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p3d& operator/=(moth_real_t a)
    {
        return *this = *this / a;
    }
    MOTH_HOST MOTH_DEVICE
    friend moth_p3d operator/(moth_real_t a, const moth_p3d& p)
    {
        return p / a;
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_p3d max(const moth_p3d& p1, const moth_p3d& p2)
    {
        return {std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z)};
    }
    MOTH_HOST MOTH_DEVICE
    static moth_p3d min(const moth_p3d& p1, const moth_p3d& p2)
    {
        return {std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z)};
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t dot(const moth_p3d& p1, const moth_p3d& p2)
    {
        return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
    }
    MOTH_HOST MOTH_DEVICE
    static moth_real_t len(const moth_p3d& p1)
    {
#if MOTH_CPP17
        return std::hypot(p1.x, p1.y, p1.z);
#else   // if MOTH_CPP17
        return std::sqrt(dot(p1, p1));
#endif  // if MOTH_CPP17
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_p3d cross(const moth_p3d& p1, const moth_p3d& p2)
    {
        return {p1.y * p2.z - p2.y * p1.z,
                p1.z * p2.x - p2.z * p1.x,
                p1.x * p2.y - p2.x * p1.y};
    }
    MOTH_HOST MOTH_DEVICE
    static moth_p3d cross(const moth_p3d& p1, const moth_p3d& p2, const moth_p3d& p3)
    {
        return cross(cross(p1, p2), p3);
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t mixed(const moth_p3d& p1, const moth_p3d& p2, const moth_p3d& p3)
    {
        return dot(p1, cross(p2, p3));
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t det(const moth_p3d& p1, const moth_p3d& p2)
    {
        moth_p3d n{cross(p1, p2)};
        moth_real_t l{len(n)};
        if (l != 0.0) {
            n /= l;
            return mixed(n, p1, p2);
        } else {
            return 0.0;
        }
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_radians_t angle(const moth_p3d& p1, const moth_p3d& p2)
    {
        moth_real_t sin{det(p1, p2)};
        moth_real_t cos{dot(p1, p2)};
        return std::atan2(sin, cos);
    }

    MOTH_HOST MOTH_DEVICE
    static moth_p3d rotate(const moth_p3d& p, const moth_p3d& u, moth_real_t sin, moth_real_t cos)
    {
        moth_real_t ver = 1.0 - cos;
        moth_p3d a{u.x*u.x * ver + cos, u.x*u.y * ver - u.z*sin, u.x*u.z * ver + u.y*sin};
        moth_p3d b{u.y*u.x * ver + u.z*sin, u.y*u.y * ver + cos, u.y*u.z * ver - u.x*sin};
        moth_p3d c{u.z*u.x * ver - u.y*sin, u.z*u.y * ver + u.x*sin, u.z*u.z * ver + cos};
        moth_p3d rot{dot(a, p), dot(b, p), dot(c, p)};
        return rot;
    }
    MOTH_HOST MOTH_DEVICE
    static moth_p3d rotate(const moth_p3d& p, const moth_p3d& u, moth_radians_t theta)
    {
        return rotate(p, u, std::sin(theta), std::cos(theta));
    }

public:
    MOTH_HOST
    static std::string str(const moth_p3d& p1);
};  // struct moth_p3d

MOTH_HOST MOTH_CORE
extern std::ostream& operator<<(std::ostream& stream, const moth_p3d& p);
MOTH_HOST MOTH_CORE
extern std::istream& operator>>(std::istream& stream, moth_p3d& p);

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

extern moth_p3d geom_p_inf;
extern moth_p3d geom_m_inf;
