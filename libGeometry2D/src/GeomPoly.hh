#pragma once

#include "libGeometry2D/src/GeomPoint.hh"
#include "libGeometry2D/src/GeomEdge.hh"

#include <vector>
#include <cassert>

using geom_p2d_array = std::vector<geom_p2d>;

/**
 * 2D polygon base container.
 */
struct GEOM_CORE geom_poly2d_base : protected geom_p2d_array
{
    friend struct geom_poly2d_iter;

public:
    GEOM_HOST GEOM_DEVICE
    explicit geom_poly2d_base(const geom_p2d& p): geom_p2d_array{{p}} {}
    GEOM_HOST GEOM_DEVICE
    explicit geom_poly2d_base() = default;
    GEOM_HOST GEOM_DEVICE
    virtual ~geom_poly2d_base() = default;
};  // struct geom_poly2d_base

/**
 * 2D polygon vertex/edge iterator.
 */
struct GEOM_CORE geom_poly2d_iter final
{
    const geom_poly2d_base* poly = nullptr;
    geom_diff_t offset = 1;
    geom_size_t index = 0;

public:
    GEOM_HOST GEOM_DEVICE
    bool operator==(const geom_poly2d_iter& iter) const
    {
        return (poly == iter.poly) && (index == iter.index);
    }
    GEOM_HOST GEOM_DEVICE
    bool operator!=(const geom_poly2d_iter& iter) const
    {
        return (poly != iter.poly) || (index != iter.index);
    }

public:
    GEOM_HOST GEOM_DEVICE
    geom_poly2d_iter operator+(const geom_diff_t delta) const
    {
        return {poly, offset, (index + offset * delta) % poly->size()};
    }
    GEOM_HOST GEOM_DEVICE
    geom_poly2d_iter& operator+=(const geom_diff_t delta)
    {
        return *this = *this + delta;
    }

    GEOM_HOST GEOM_DEVICE
    geom_poly2d_iter operator-(const geom_diff_t delta) const
    {
        return {poly, offset, (index - offset * delta) % poly->size()};
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
    geom_poly2d_iter next() const
    {
        return *this + 1;
    }
    GEOM_HOST GEOM_DEVICE
    geom_poly2d_iter prev() const
    {
        return *this - 1;
    }

public:
    GEOM_HOST GEOM_DEVICE
    geom_p2d point() const
    {
        return (*poly)[index];
    }
    GEOM_HOST GEOM_DEVICE
    geom_e2d edge() const
    {
        return {point(), next().point()};
    }
};  // struct geom_poly2d_iter

/**
 * 2D polygon.
 */
struct GEOM_CORE geom_poly2d final : public geom_poly2d_base
{
public:
    GEOM_HOST GEOM_DEVICE
    explicit geom_poly2d(const geom_p2d& p): geom_poly2d_base{p} {}
    GEOM_HOST GEOM_DEVICE
    explicit geom_poly2d() = default;

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
    GEOM_HOST
    static geom_e2d_list* lst(const geom_poly2d& poly)
    {
        geom_e2d_list* list = nullptr;
        geom_poly2d_iter iter = poly.iter();
        do {
            geom_e2d_list::pull(list, iter.point());
        } while ((++iter) != poly.iter());
        return list->next;
    }
    GEOM_HOST
    static geom_e2d_list* lst_rev(const geom_poly2d& poly)
    {
        geom_e2d_list* list = nullptr;
        geom_poly2d_iter iter = poly.iter_rev();
        do {
            geom_e2d_list::pull(list, iter.point());
        } while ((++iter) != poly.iter());
        return list;
    }

public:
    GEOM_HOST
    static void push(geom_poly2d& poly, const geom_p2d& p)
    {
        poly.push_back(p);
    }
    GEOM_HOST
    static void pop(geom_poly2d& poly)
    {
        poly.pop_back();
    }

public:
    GEOM_HOST GEOM_DEVICE
    static geom_real_t len(const geom_poly2d& poly)
    {
        geom_real_t length = 0.0;
        geom_poly2d_iter iter = poly.iter();
        do {
            length += geom_e2d::len(iter.edge());
        } while ((++iter) != poly.iter());
        return length;
    }

    GEOM_HOST GEOM_DEVICE
    static geom_real_t area(const geom_poly2d& poly)
    {
        geom_real_t area = 0.0;
        geom_poly2d_iter iter = poly.iter();
        do {
            area += geom_p2d::det(iter.edge().s, iter.edge().t) * 0.5;
        } while ((++iter) != poly.iter());
        return area;
    }

public:
    /**
     * Check if point is inside polygon.
     */
    GEOM_HOST GEOM_DEVICE
    static bool contains(const geom_poly2d& poly, const geom_p2d& p)
    {
        size_t intersections = 0;
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

public:
    /**
     * Intersection of two polygons.
     */
    GEOM_HOST
    friend geom_poly2d operator*(const geom_poly2d& poly1, const geom_poly2d& poly2)
    {
        std::abort();
    }

    /**
     * Union of two polygons.
     */
    GEOM_HOST
    friend geom_poly2d operator+(const geom_poly2d& poly1, const geom_poly2d& poly2)
    {
        /* Make CCW list for polygons. */
        geom_e2d_list* list1;
        if (area(poly1) > 0.0) {
            list1 = lst(poly1);
        } else {
            list1 = lst_rev(poly1);
        }
        geom_e2d_list* list2;
        if (area(poly2) > 0.0) {
            list2 = lst(poly2);
        } else {
            list2 = lst_rev(poly2);
        }

        geom_poly2d poly;
        geom_e2d_list* head1 = list1;
        geom_e2d_list::clip(list1, list2);
        do {
            geom_poly2d::push(poly, list1->point);
            if (list1->next_ext != nullptr) {
                list1 = list1->next_ext;
            } else {
                list1 = list1->next;
            }
        } while (list1 != head1);

        delete list1;
        delete list2;
        return poly;
    }

    /**
     * Difference of two polygons.
     */
    GEOM_HOST
    friend geom_poly2d operator-(const geom_poly2d& poly1, const geom_poly2d& poly2)
    {
        std::abort();
    }

public:
    GEOM_HOST
    static std::string str(const geom_poly2d& poly);
    GEOM_HOST
    static std::string plt(const geom_poly2d& poly);
};  // struct geom_poly2d

GEOM_HOST GEOM_CORE
std::ostream& operator<<(std::ostream& stream, const geom_poly2d& poly);
GEOM_HOST GEOM_CORE
std::istream& operator>>(std::istream& stream, geom_poly2d& poly);

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
    GEOM_HOST
    static geom_poly2d rect(const geom_p2d& p_sw, const geom_p2d& p_ne)
    {
        geom_p2d p_se{p_ne.x, p_sw.y};
        geom_p2d p_nw{p_sw.x, p_ne.y};
        geom_poly2d poly{};
        geom_poly2d::push(poly, p_sw);
        geom_poly2d::push(poly, p_se);
        geom_poly2d::push(poly, p_ne);
        geom_poly2d::push(poly, p_nw);
        return poly;
    }

    /**
     * Circle of the given radius with CCW orientation.
     */
    GEOM_HOST
    static geom_poly2d circle(const geom_p2d& c, geom_real_t r, geom_size_t n = 10)
    {
        assert(r > 0.0);
        assert(n > 1);
        geom_poly2d poly{};
        for (geom_size_t i = 0; i < n; ++i) {
            geom_real_t phi = 2.0 * GEOM_PI * i / n;
            geom_p2d p = c + geom_p2d{r * std::cos(phi), r * std::sin(phi)};
            geom_poly2d::push(poly, p);
        }
        return poly;
    }

    /**
     * Star of given outer and inner radius with CCW orientation.
     */
    GEOM_HOST
    static geom_poly2d star(const geom_p2d &c, geom_real_t r1, geom_real_t r2, geom_size_t n = 10)
    {
        assert(r1 > 0.0 && r2 > 0.0);
        assert(n > 1);
        geom_poly2d poly{};
        for (geom_size_t i = 0; i < n; ++i) {
            geom_real_t phi1 = 2.0 * GEOM_PI * (i + 0.0) / n;
            geom_p2d p1 = c + geom_p2d{r1 * std::cos(phi1), r1 * std::sin(phi1)};
            geom_poly2d::push(poly, p1);

            geom_real_t phi2 = 2.0 * GEOM_PI * (i + 0.5) / n;
            geom_p2d p2 = c + geom_p2d{r1 * std::cos(phi1), r2 * std::sin(phi1)};
            geom_poly2d::push(poly, p2);
        }
        return poly;
    }

public:
    /**
     * Returns convex hull of the given set of points.
     */
    GEOM_HOST
    static geom_poly2d convex_hull(geom_p2d_array& p)
    {
        std::iter_swap(p.begin(), std::min_element(p.begin(), p.end(),
            [&](const geom_p2d& p1, const geom_p2d& p2) {
                return p1.x < p2.x;
            }));
        std::sort(p.begin() + 1, p.end(),
            [&](const geom_p2d& p1, const geom_p2d& p2) {
                return geom_p2d::det(p1 - p[0], p2 - p[0]) < 0;
            });

        geom_poly2d hull;
        geom_poly2d::push(hull, p[0]);
        geom_poly2d::push(hull, p[1]);
        geom_poly2d::push(hull, p[2]);
        for (geom_size_t i = 3; i < p.size(); ++i) {
            geom_p2d p_i = p[i];
            geom_e2d e = hull.iter_rev().edge();
            while (geom_p2d::det(e.t - p_i, e.s - p_i) >= 0.0) {
                geom_poly2d::pop(hull);
            }
            geom_poly2d::push(hull, p_i);
        }
        return hull;
    }
};  // struct geom_poly2d_primitives
