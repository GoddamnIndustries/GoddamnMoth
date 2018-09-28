#pragma once

#include "libGeometry2D/src/GeomPoint.hh"
#include "libGeometry2D/src/GeomEdge.hh"

template<typename geom_p2d_vector = std::vector<geom_p2d>>
struct geom_poly2d;
template<typename geom_p2d_vector = std::vector<geom_p2d>>
struct geom_poly2d_iter;

///
/// 2D poly (list of points in 2D space) iterator.
///
template<typename geom_p2d_vector>
struct GEOM_CORE geom_poly2d_iter final
{
    using geom_poly2d = geom_poly2d<geom_p2d_vector>;
    const geom_poly2d* ptr = nullptr;
    geom_diff_t offset = 1;
    geom_size_t index = 0;

public:
    GEOM_HOST GEOM_DEVICE
    bool operator==(const geom_poly2d_iter& poly) const
    {
        return (ptr == poly.ptr) && (index == poly.index);
    }
    GEOM_HOST GEOM_DEVICE
    bool operator!=(const geom_poly2d_iter& poly) const
    {
        return (ptr != poly.ptr) || (index == poly.index);
    }

public:
    GEOM_HOST GEOM_DEVICE
    geom_poly2d_iter operator+(const geom_diff_t delta) const
    {
        return {ptr, offset, (index + offset * delta) % ptr->points.size()};
    }
    GEOM_HOST GEOM_DEVICE
    geom_poly2d_iter& operator+=(const geom_diff_t delta)
    {
        return *this = *this + delta;
    }

    GEOM_HOST GEOM_DEVICE
    geom_poly2d_iter operator-(const geom_diff_t delta) const
    {
        return {ptr, offset, (index - offset * delta) % ptr->points.size()};
    }
    GEOM_HOST GEOM_DEVICE
    geom_poly2d_iter& operator-=(const geom_diff_t delta)
    {
        return *this = *this - delta;
    }

public:
    GEOM_HOST GEOM_DEVICE
    geom_poly2d_iter& operator++()
    {
        return *this += 1;
    }
    GEOM_HOST GEOM_DEVICE
    const geom_poly2d_iter operator++(int)
    {
        return (*this += 1) - 1;
    }

    GEOM_HOST GEOM_DEVICE
    geom_poly2d_iter& operator--()
    {
        return *this -= 1;
    }
    GEOM_HOST GEOM_DEVICE
    const geom_poly2d_iter operator--(int)
    {
        return (*this -= 1) + 1;
    }

public:
    GEOM_HOST GEOM_DEVICE
    static bool move(geom_poly2d_iter& iter, const geom_poly2d& poly)
    {
        return (++poly).ptr != &poly;
    }

public:
    GEOM_HOST GEOM_DEVICE
    geom_p2d point() const
    {
        return (*ptr)[index];
    }
    GEOM_HOST GEOM_DEVICE
    geom_e2d edge() const
    {
        return {(*ptr)[index], (*ptr)[index + 1]};
    }
};  // struct geom_poly2d_iter

///
/// 2D poly (list of points in 2D space).
///
template<typename geom_p2d_vector>
struct GEOM_CORE geom_poly2d final
{
    using geom_poly2d_iter = geom_poly2d_iter<geom_p2d_vector>;
    geom_p2d_vector points;

public:
    GEOM_HOST GEOM_DEVICE
    const geom_p2d& operator[](geom_diff_t i) const
    {
        return points[i % points.size()];
    }

public:
    GEOM_HOST GEOM_DEVICE
    geom_poly2d_iter iter() const
    {
        return {this, +1};
    }
    GEOM_HOST GEOM_DEVICE
    geom_poly2d_iter iter_rev() const
    {
        return {this, -1};
    }

public:
    GEOM_HOST GEOM_DEVICE
    static geom_real_t len(const geom_poly2d& poly);

    GEOM_HOST GEOM_DEVICE
    static geom_real_t area(const geom_poly2d& poly);

public:
    GEOM_HOST GEOM_DEVICE
    static bool contains(const geom_poly2d& poly, const geom_p2d& p);

public:
    GEOM_HOST
    static std::string str(const geom_poly2d& poly);
    GEOM_HOST
    static std::string plt(const geom_poly2d& poly);
};  // struct geom_poly2d

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

GEOM_HOST GEOM_CORE
std::ostream& operator<<(std::ostream& stream, const geom_poly2d<>& poly);
GEOM_HOST GEOM_CORE
std::istream& operator>>(std::istream& stream, geom_poly2d<>& poly);

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

GEOM_HOST GEOM_DEVICE
template<typename geom_p2d_vector>
inline geom_real_t geom_poly2d<geom_p2d_vector>::len(const geom_poly2d& poly)
{
    geom_real_t length = 0.0;
    geom_poly2d_iter iter = poly.iter();
    do {
        length += geom_e2d::len(iter.edge());
    } while ((++iter) != poly.iter());
    return length;
}

GEOM_HOST GEOM_DEVICE
template<typename geom_p2d_vector>
inline geom_real_t geom_poly2d<geom_p2d_vector>::area(const geom_poly2d& poly)
{
    geom_real_t area = 0.0;
    geom_poly2d_iter iter = poly.iter();
    do {
        area += geom_p2d::det(iter.edge().t, iter.edge().s) * 0.5;
    } while ((++iter) != poly.iter());
    return area;
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

GEOM_HOST GEOM_DEVICE
template<typename geom_p2d_vector>
inline bool geom_poly2d<geom_p2d_vector>::contains(const geom_poly2d& poly, const geom_p2d& p)
{
    geom_size_t intersections = 0;
    geom_poly2d_iter iter = poly.iter();
    do {
        geom_e2d e;
        geom_e2d e1 = iter.edge();
        geom_e2d e2 = {p, p + geom_p_inf};
        if (geom_e2d::intersect(e1, e2, e)) {
            if (e.s == e.t) {
                intersections += 1;
            } else {
                intersections += 2;
            }
        }
    } while ((++iter) != poly.iter());
    return intersections % 2 == 1;
}
