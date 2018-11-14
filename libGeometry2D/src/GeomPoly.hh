#pragma once

#include "libGeometry2D/src/GeomPoint.hh"
#include "libGeometry2D/src/GeomEdge.hh"

#include <vector>
#include <cassert>

using geom_p2d_array = std::vector<moth_p2d>;

/**
 * 2D polygon base container.
 */
struct MOTH_CORE moth_poly2d_base : public geom_p2d_array
{
    friend struct geom_poly2d_iter;

public:
    MOTH_HOST MOTH_DEVICE
    explicit moth_poly2d_base(const moth_p2d& p): geom_p2d_array{{p}} {}
    MOTH_HOST MOTH_DEVICE
    explicit moth_poly2d_base() = default;
    MOTH_HOST MOTH_DEVICE
    virtual ~moth_poly2d_base() = default;
};  // struct moth_poly2d_base

/**
 * 2D polygon vertex/edge iterator.
 */
struct MOTH_CORE moth_poly2d_iter final
{
    const moth_poly2d_base* poly = nullptr;
    moth_diff_t offset = 1;
    moth_size_t index = 0;

public:
    MOTH_HOST MOTH_DEVICE
    bool operator==(const moth_poly2d_iter& iter) const
    {
        return (poly == iter.poly) && (index == iter.index);
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_poly2d_iter& iter) const
    {
        return (poly != iter.poly) || (index != iter.index);
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_poly2d_iter operator+(const moth_diff_t delta) const
    {
        return {poly, offset, (index + offset * delta) % poly->size()};
    }
    MOTH_HOST MOTH_DEVICE
    moth_poly2d_iter& operator+=(const moth_diff_t delta)
    {
        return *this = *this + delta;
    }

    MOTH_HOST MOTH_DEVICE
    moth_poly2d_iter operator-(const moth_diff_t delta) const
    {
        return {poly, offset, (index - offset * delta) % poly->size()};
    }
    MOTH_HOST MOTH_DEVICE
    moth_poly2d_iter& operator-=(const moth_diff_t delta)
    {
        return *this = *this - delta;
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_poly2d_iter& operator++()
    {
        return *this += 1;
    }
    MOTH_HOST MOTH_DEVICE
    const moth_poly2d_iter operator++(int)
    {
        return (*this += 1) - 1;
    }

    MOTH_HOST MOTH_DEVICE
    moth_poly2d_iter& operator--()
    {
        return *this -= 1;
    }
    MOTH_HOST MOTH_DEVICE
    const moth_poly2d_iter operator--(int)
    {
        return (*this -= 1) + 1;
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_poly2d_iter next() const
    {
        return *this + 1;
    }
    MOTH_HOST MOTH_DEVICE
    moth_poly2d_iter prev() const
    {
        return *this - 1;
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_p2d point() const
    {
        return (*poly)[index];
    }
    MOTH_HOST MOTH_DEVICE
    moth_e2d edge() const
    {
        return {point(), next().point()};
    }
};  // struct moth_poly2d_iter

/**
 * 2D polygon.
 */
struct MOTH_CORE moth_poly2d final : public moth_poly2d_base
{
public:
    MOTH_HOST MOTH_DEVICE
    explicit moth_poly2d(const moth_p2d& p): moth_poly2d_base{p} {}
    MOTH_HOST MOTH_DEVICE
    explicit moth_poly2d() = default;

public:
    MOTH_HOST MOTH_DEVICE
    moth_poly2d_iter iter() const
    {
        return {this, +1};
    }
    MOTH_HOST MOTH_DEVICE
    moth_poly2d_iter iter_rev() const
    {
        return {this, -1};
    }

public:
    MOTH_HOST
    static void push(moth_poly2d& poly, const moth_p2d& p)
    {
        poly.push_back(p);
    }
    MOTH_HOST
    static void pop(moth_poly2d& poly)
    {
        poly.pop_back();
    }

public:
    MOTH_HOST MOTH_DEVICE
    static moth_real_t len(const moth_poly2d& poly)
    {
        moth_real_t length = 0.0;
        moth_poly2d_iter iter = poly.iter();
        do {
            length += moth_e2d::len(iter.edge());
        } while ((++iter) != poly.iter());
        return length;
    }

    MOTH_HOST MOTH_DEVICE
    static moth_real_t area(const moth_poly2d& poly)
    {
        moth_real_t area = 0.0;
        moth_poly2d_iter iter = poly.iter();
        do {
            area += moth_p2d::det(iter.edge().p1, iter.edge().p2) * 0.5;
        } while ((++iter) != poly.iter());
        return area;
    }

public:
    /**
     * Check if point is inside polygon.
     */
    MOTH_HOST MOTH_DEVICE
    static bool contains(const moth_poly2d& poly, const moth_p2d& p)
    {
        size_t intersections = 0;
        moth_poly2d_iter iter = poly.iter();
        do {
            moth_e2d e;
            moth_e2d e1 = iter.edge();
            moth_e2d e2 = {p, p + geom_p_inf};
            if (moth_e2d::intersect(e1, e2, e)) {
                if (e.p1 == e.p2) {
                    intersections += 1;
                } else {
                    intersections += 2;
                }
            }
        } while ((++iter) != poly.iter());
        return intersections % 2 == 1;
    }

public:
    MOTH_HOST
    static bool simple(const moth_poly2d& poly);

public:
    /**
     * Intersection of two polygons.
     */
    MOTH_HOST
    friend moth_poly2d operator*(const moth_poly2d& poly1, const moth_poly2d& poly2)
    {
        std::abort();
    }

    /**
     * Union of two polygons.
     */
    MOTH_HOST
    friend moth_poly2d operator+(const moth_poly2d& poly1, const moth_poly2d& poly2)
    {
        std::abort();
    }

    /**
     * Difference of two polygons.
     */
    MOTH_HOST
    friend moth_poly2d operator-(const moth_poly2d& poly1, const moth_poly2d& poly2)
    {
        std::abort();
    }

public:
    MOTH_HOST
    static std::string str(const moth_poly2d& poly);
    MOTH_HOST
    static std::string plt(const moth_poly2d& poly);
};  // struct moth_poly2d

MOTH_HOST MOTH_CORE
std::ostream& operator<<(std::ostream& stream, const moth_poly2d& poly);
MOTH_HOST MOTH_CORE
std::istream& operator>>(std::istream& stream, moth_poly2d& poly);

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * 2D polygon utilities.
 */
struct geom_poly2d_primitives
{
public:
    /**
     * Returns rectangle with given south-west and north-east points.
     */
    MOTH_HOST
    static moth_poly2d rect(const moth_p2d& p_sw, const moth_p2d& p_ne)
    {
        moth_p2d p_se{p_ne.x, p_sw.y};
        moth_p2d p_nw{p_sw.x, p_ne.y};
        moth_poly2d poly{};
        moth_poly2d::push(poly, p_sw);
        moth_poly2d::push(poly, p_se);
        moth_poly2d::push(poly, p_ne);
        moth_poly2d::push(poly, p_nw);
        return poly;
    }

    /**
     * Circle of the given radius with CCW orientation.
     */
    MOTH_HOST
    static moth_poly2d circle(const moth_p2d& c, moth_real_t r, moth_size_t n = 10)
    {
        assert(r > 0.0);
        assert(n > 1);
        moth_poly2d poly{};
        for (moth_size_t i = 0; i < n; ++i) {
            moth_real_t phi = 2.0 * MOTH_PI * i / n;
            moth_p2d p = c + moth_p2d{r * std::cos(phi), r * std::sin(phi)};
            moth_poly2d::push(poly, p);
        }
        return poly;
    }

    /**
     * Star of given outer and inner radius with CCW orientation.
     */
    MOTH_HOST
    static moth_poly2d star(const moth_p2d &c, moth_real_t r1, moth_real_t r2, moth_size_t n = 10)
    {
        assert(r1 > 0.0 && r2 > 0.0);
        assert(n > 1);
        moth_poly2d poly{};
        for (moth_size_t i = 0; i < n; ++i) {
            moth_real_t phi1 = 2.0 * MOTH_PI * (i + 0.0) / n;
            moth_p2d p1 = c + moth_p2d{r1 * std::cos(phi1), r1 * std::sin(phi1)};
            moth_poly2d::push(poly, p1);

            moth_real_t phi2 = 2.0 * MOTH_PI * (i + 0.5) / n;
            moth_p2d p2 = c + moth_p2d{r2 * std::cos(phi2), r2 * std::sin(phi2)};
            moth_poly2d::push(poly, p2);
        }
        return poly;
    }

public:
    /**
     * Returns convex hull of the given set of points.
     */
    MOTH_HOST
    static moth_poly2d convex_hull(geom_p2d_array& p)
    {
        std::iter_swap(p.begin(), std::min_element(p.begin(), p.end(),
            [&](const moth_p2d& p1, const moth_p2d& p2) {
                return p1.x < p2.x;
            }));
        std::sort(p.begin() + 1, p.end(),
            [&](const moth_p2d& p1, const moth_p2d& p2) {
                return moth_p2d::det(p1 - p[0], p2 - p[0]) < 0;
            });

        moth_poly2d hull;
        moth_poly2d::push(hull, p[0]);
        moth_poly2d::push(hull, p[1]);
        moth_poly2d::push(hull, p[2]);
        for (moth_size_t i = 3; i < p.size(); ++i) {
            moth_p2d p_i = p[i];
            moth_e2d e = hull.iter_rev().edge();
            while (moth_p2d::det(e.p2 - p_i, e.p1 - p_i) >= 0.0) {
                moth_poly2d::pop(hull);
            }
            moth_poly2d::push(hull, p_i);
        }
        return hull;
    }
};  // struct geom_poly2d_primitives
