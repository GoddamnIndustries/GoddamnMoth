#pragma once

#include "libGeometry2D/src/GeomBase.hh"

/**
 * Point in 2D space.
 */
struct MOTH_CORE moth_p2d final
{
    moth_real_t x{};
    moth_real_t y{};

public:
    MOTH_HOST MOTH_DEVICE
    bool operator==(const moth_p2d& p1) const
    {
        return (x == p1.x) && (y == p1.y);
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_p2d& p1) const
    {
        return (x != p1.x) || (y != p1.y);
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_p2d operator+() const
    {
        return {+x, +y};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p2d operator+(const moth_p2d& p1) const
    {
        return {x + p1.x, y + p1.y};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p2d& operator+=(const moth_p2d& p1)
    {
        return *this = *this + p1;
    }

    MOTH_HOST MOTH_DEVICE
    moth_p2d operator-() const
    {
        return {-x, -y};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p2d operator-(const moth_p2d& p1) const
    {
        return {x - p1.x, y - p1.y};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p2d& operator-=(const moth_p2d& p1)
    {
        return *this = *this - p1;
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
    friend moth_p2d operator*(moth_real_t a, const moth_p2d& p1)
    {
        return p1 * a;
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
    friend moth_p2d operator/(moth_real_t a, const moth_p2d& p1)
    {
        return p1 / a;
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
    static moth_real_t det(const moth_p2d& p1, const moth_p2d& p2)
    {
        return p1.x * p2.y - p2.x * p1.y;
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_radians_t angle(const moth_p2d& p1, const moth_p2d& p2)
    {
        moth_real_t sin{det(p1, p2)};
        moth_real_t cos{dot(p1, p2)};
        return std::atan2(sin, cos);
    }

    MOTH_HOST MOTH_DEVICE
    static moth_p2d rotate(const moth_p2d& n, moth_real_t sin, moth_real_t cos)
    {
        moth_p2d rot{n.x * cos - n.y * sin,
                     n.x * sin + n.y * cos};
        return rot;
    }
    MOTH_HOST MOTH_DEVICE
    static moth_p2d rotate(const moth_p2d& n, moth_radians_t theta)
    {
        return rotate(n, std::sin(theta), std::cos(theta));
    }

public:
    MOTH_HOST MOTH_DEVICE [[deprecated]]
    static moth_p2d normal(const moth_p2d& p1)
    {
        moth_p2d n{-p1.y, p1.x};
        moth_real_t l = len(n);
        if (l != 0.0) {
            return n / l;
        } else {
            std::cerr << "Warning: normal to null vector." << std::endl;
            return {0.0, 0.0};
        }
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
    bool operator==(const moth_p3d& p1) const
    {
        return (x == p1.x) && (y == p1.y) && (z == p1.z);
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_p3d& p1) const
    {
        return (x != p1.x) || (y != p1.y) || (z != p1.z);
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_p3d operator+() const
    {
        return {+x, +y, +z};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p3d operator+(const moth_p3d& p1) const
    {
        return {x + p1.x, y + p1.y, z + p1.z};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p3d& operator+=(const moth_p3d& p1)
    {
        return *this = *this + p1;
    }

    MOTH_HOST MOTH_DEVICE
    moth_p3d operator-() const
    {
        return {-x, -y, -z};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p3d operator-(const moth_p3d& p1) const
    {
        return {x - p1.x, y - p1.y, z - p1.z};
    }
    MOTH_HOST MOTH_DEVICE
    moth_p3d& operator-=(const moth_p3d& p1)
    {
        return *this = *this - p1;
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
    friend moth_p3d operator*(moth_real_t a, const moth_p3d& p1)
    {
        return p1 * a;
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
    friend moth_p3d operator/(moth_real_t a, const moth_p3d& p1)
    {
        return p1 / a;
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
    static moth_p3d rotate(const moth_p3d& n, const moth_p3d& u, moth_real_t sin, moth_real_t cos)
    {
        moth_real_t ver = 1.0 - cos;
        moth_p3d rot{n.x * (u.x*u.x * ver + cos) + n.y * (u.x*u.y * ver - u.z*sin) + n.z * (u.x*u.z * ver + u.y*sin),
                     n.x * (u.y*u.x * ver + u.z*sin) + n.y * (u.y*u.y * ver + cos) + n.z * (u.y*u.z * ver - u.x*sin),
                     n.x * (u.z*u.x * ver - u.y*sin) + n.y * (u.z*u.y * ver + u.x*sin) + n.z * (u.z*u.z * ver + cos)};
        return rot;
    }
    MOTH_HOST MOTH_DEVICE
    static moth_p3d rotate(const moth_p3d& n, const moth_p3d& u, moth_radians_t theta)
    {
        return rotate(n, u, std::sin(theta), std::cos(theta));
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

static const moth_p3d geom_p_inf{+500000.0, +500000.0, +500000.0};
static const moth_p3d geom_m_inf{-500000.0, -500000.0, -500000.0};
