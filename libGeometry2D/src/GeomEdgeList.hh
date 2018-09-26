#pragma once
#include "GeomEdge.hh"

#include <algorithm>
#include <iostream>
#include <istream>
#include <vector>
#include <memory>

#include <cassert>
#include <cmath>

#undef min
#undef max

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

///
/// 2D Polygon (list of points in 2D space).
///
struct geom_e2d_list final
{
public:
    geom_p2d point{};
    geom_e2d_list* next{this};

public:
    geom_e2d_list(geom_e2d_list&&) = delete;
    geom_e2d_list(const geom_e2d_list&) = delete;
    geom_e2d_list& operator= (geom_e2d_list&&) = delete;
    geom_e2d_list& operator= (const geom_e2d_list&) = delete;

public:
    explicit geom_e2d_list(const geom_p2d& p): point(p) {}
    virtual ~geom_e2d_list()
    {
        if (next != this) {
            do {
                remove();
            } while (next != this);
        }
    }

public:
    geom_e2d edge() const
    {
        return {point, next->point};
    }

public:
    static bool move(geom_e2d_list*& list, const geom_e2d_list* head = nullptr)
    {
        list = list->next;
        return list != head;
    }
    static bool move(const geom_e2d_list*& list, const geom_e2d_list* head = nullptr)
    {
        return move(const_cast<geom_e2d_list*&>(list), head);
    }

public:
    void insert(const geom_p2d& p)
    {
        geom_e2d_list* list = this;
        geom_e2d_list* node = new geom_e2d_list{p};
        node->next = list->next;
        list->next = node;
    }
    void remove()
    {
        geom_e2d_list* list = this;
        geom_e2d_list* node = list->next;
        list->next = node->next;
        node->next = node;
        delete node;
    }

public:
    static void push(geom_e2d_list*& poly, const geom_p2d& p)
    {
        if (poly == nullptr) {
            poly = new geom_e2d_list{p};
        } else {
            geom_e2d_list* node = new geom_e2d_list{p};
            node->next = poly->next;
            poly->next = node;
            std::swap(poly->point, poly->next->point);
        }
    }
    static void pull(geom_e2d_list*& poly, const geom_p2d& p)
    {
        if (poly == nullptr) {
            poly = new geom_e2d_list{p};
        } else {
            geom_e2d_list* node = new geom_e2d_list{p};
            node->next = poly->next;
            poly->next = node;
            move(poly);
        }
    }
    static void pop(geom_e2d_list*& poly)
    {
        if (poly->next == poly) {
            delete poly;
            poly = nullptr;
        } else {
            std::swap(poly->point, poly->next->point);
            poly->remove();
        }
    }

public:
    static geom_real_t len(const geom_e2d_list* poly)
    {
        geom_real_t length = 0.0;
        const geom_e2d_list* head = poly;
        do {
            length += geom_e2d::len(poly->edge());
        } while (move(poly, head));
        return length;
    }

    static geom_real_t area(const geom_e2d_list* poly)
    {
        geom_real_t area = 0.0;
        const geom_e2d_list* head = poly;
        do {
            geom_e2d e = poly->edge();
            area += geom_p2d::det(e.t, e.s) * 0.5;
        } while (move(poly, head));
        return area;
    }

public:
    static bool contains(const geom_e2d_list* poly, const geom_p2d& p)
    {
        geom_size_t intersections = 0;
        const geom_e2d_list* head = poly;
        do {
            geom_e2d e;
            geom_e2d e1 = poly->edge();
            if (e1.s == p || e1.t == p) {
                /* Boundary points are explicitly treated as internal. */
                return true;
            }

            geom_e2d e2 = {p, {std::max(p.x + 50000.0, std::max(e1.s.x, e1.t.x)), p.y + 1.0}};
            if (geom_e2d::intersect(e1, e2, e)) {
                if (e.s == e.t) {
                    intersections += 1;
                } else {
                    intersections += 2;
                }
            }
        } while (move(poly, head));
        return intersections % 2 == 1;
    }

    static bool reflect(const geom_e2d_list* poly, const geom_e2d& e2, geom_e2d e)
    {
        const geom_e2d_list* head = poly;
        do {
            geom_e2d e1 = poly->edge();
            if (geom_e2d::reflect(e1, e2, e)) {
                return true;
            }
        } while (move(poly, head));
        return false;
    }

public:
    static std::string str(const geom_e2d_list* poly);
    static std::string plt(const geom_e2d_list* poly, ...);
};  // struct geom_e2d_list

std::ostream& operator<< (std::ostream& stream, const geom_e2d_list* poly);
std::istream& operator>> (std::istream& stream, geom_e2d_list*& poly);

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

///
/// 2D polygon factory.
///
struct geom_e2d_list_factory final
{
public:
    /// Returns copy of the polygon.
    static geom_e2d_list* copy(const geom_e2d_list* poly)
    {
        assert(poly != nullptr);
        geom_e2d_list* copy = nullptr;
        const geom_e2d_list* head = poly;
        do {
            geom_e2d_list::pull(copy, poly->point);
        } while (geom_e2d_list::move(poly, head));
        geom_e2d_list::move(copy);
        return copy;
    }

    /// Returns reverse copy of the polygon.
    static geom_e2d_list* copy_rev(const geom_e2d_list* poly)
    {
        assert(poly != nullptr);
        geom_e2d_list* copy = nullptr;
        const geom_e2d_list* head = poly;
        do {
            geom_e2d_list::push(copy, poly->point);
        } while (geom_e2d_list::move(poly, head));
        return copy;
    }

public:
    /// Returns rectangle with given south-west and north-east points.
    static geom_e2d_list* new_rect_ccw(const geom_p2d& p_sw, const geom_p2d& p_ne)
    {
        assert(p_sw.x < p_ne.x);
        assert(p_sw.y < p_ne.y);
        geom_p2d p_se{p_sw.x, p_ne.y};
        geom_p2d p_nw{p_ne.x, p_sw.y};
        geom_e2d_list* poly = new geom_e2d_list(p_sw);
        poly->insert(p_nw);
        poly->insert(p_ne);
        poly->insert(p_se);
        return poly;
    }

public:
    /// Returns circle of the given radius with CCW orientation.
    static geom_e2d_list* new_circle_ccw(const geom_p2d& c, geom_real_t r, geom_size_t n = 10)
    {
        assert(r > 0.0);
        assert(n > 1);
        geom_e2d_list* poly = nullptr;
        for (geom_size_t i = 0; i < n; ++i) {
            geom_real_t phi = 2.0 * GEOM_PI * i / n;
            geom_p2d p = c + geom_p2d{r * cos(phi), r * sin(phi)};
            if (poly == nullptr) {
                poly = new geom_e2d_list(p);
            } else {
                poly->insert(p);
            }
        }
        return poly;
    }

    /// Returns start of given outer and inner radius with CCW orientation.
    static geom_e2d_list* new_star_ccw(const geom_p2d &c, geom_real_t r1, geom_real_t r2, geom_size_t n = 10)
    {
        assert(r1 > 0.0 && r2 > 0.0);
        assert(n > 1);
        geom_e2d_list* poly = nullptr;
        for (geom_size_t i = 0; i < n; ++i) {
            geom_real_t phi1 = 2.0 * GEOM_PI * i / n;
            geom_p2d p1 = c + geom_p2d{r1 * cos(phi1), r1 * sin(phi1)};
            if (poly == nullptr) {
                poly = new geom_e2d_list(p1);
            } else {
                poly->insert(p1);
            }

            geom_real_t phi2 = 2.0 * GEOM_PI * (i + 0.5) / n;
            geom_p2d p2 = c + geom_p2d{r2 * cos(phi2), r2 * sin(phi2)};
            poly->insert(p2);
        }
        return poly;
    }

public:
    /// Returns convex hull of the given set of points.
    static geom_e2d_list* convex_hull(const geom_e2d_list* poly)
    {
        std::vector<geom_p2d> p{};
        const geom_e2d_list* head = poly;
        do {
            p.push_back(poly->point);
        } while (geom_e2d_list::move(poly, head));

        /* Find the left-most point and move it to index 0. */
        std::iter_swap(p.begin(), std::min_element(p.begin(), p.end(),
                                                   [&](const geom_p2d& p1, const geom_p2d& p2) {
                                                       return p1.x < p2.x;
                                                   }));
        /* Sort points by the polar angles with p[0]. */
        std::sort(p.begin() + 1, p.end(),
                  [&](const geom_p2d& p1, const geom_p2d& p2) {
                      return geom_p2d::det(p1 - p[0], p2 - p[0]) < 0;
                  });

        geom_e2d_list* hull = nullptr;
        geom_e2d_list::push(hull, p[0]);
        geom_e2d_list::push(hull, p[1]);
        geom_e2d_list::push(hull, p[2]);
        for (geom_size_t i = 3; i < p.size(); ++i) {
            geom_p2d p_i = p[i];
            geom_e2d e = hull->edge();
            while (geom_p2d::det(e.t - p_i, e.s - p_i) >= 0.0) {
                geom_e2d_list::pop(hull);
            }
            geom_e2d_list::push(hull, p_i);
        }
        return hull;
    }

public:
    /// Returns result of the generic order-dependant boolean operation.
    static geom_e2d_list* boolean_generic(const geom_e2d_list* poly1, const geom_e2d_list* poly2)
    {
        geom_e2d_list* un = nullptr;
        const geom_e2d_list* start = poly1;
        const geom_e2d_list* head1 = poly1;
        const geom_e2d_list* head2 = poly2;
        do {
            do {
                geom_e2d e1 = poly1->edge();
                geom_e2d_list::push(un, poly1->point);
                do {
                    geom_e2d e0{};
                    geom_e2d e2 = poly2->edge();
                    if (geom_e2d::intersect(e1, e2, e0)) {
                        auto a = geom_e2d::angle(e1, e2);
                        std::cerr << a << std::endl;
                        geom_e2d_list::push(un, e0.s);
                        std::swap(poly1, poly2);
                        std::swap(head1, head2);
                        break;
                    }
                } while (geom_e2d_list::move(poly2, head2));
            } while (geom_e2d_list::move(poly1, head1));
        } while (start != head1);
        return un;
    }

    static geom_e2d_list* boolean_union(const geom_e2d_list* poly1, const geom_e2d_list* poly2)
    {
        std::unique_ptr<geom_e2d_list> poly1_ptr{};
        geom_real_t area1 = geom_e2d_list::area(poly1);
        geom_real_t area2 = geom_e2d_list::area(poly2);
        if (area1 * area2 < 0.0) {
            poly1_ptr.reset(geom_e2d_list_factory::copy_rev(poly1));
            poly1 = poly1_ptr.get();
        }

        const geom_e2d_list* head1 = poly1;
        do {
            if (!geom_e2d_list::contains(poly2, poly1->point)) {
                return boolean_generic(poly1, poly2);
            }
        } while (geom_e2d_list::move(poly1, head1));
        abort();
    }

    static geom_e2d_list* boolean_inter(const geom_e2d_list* poly1, const geom_e2d_list* poly2)
    {
        std::unique_ptr<geom_e2d_list> poly1_ptr{};
        geom_real_t area1 = geom_e2d_list::area(poly1);
        geom_real_t area2 = geom_e2d_list::area(poly2);
        if (area1 * area2 < 0.0) {
            poly1_ptr.reset(geom_e2d_list_factory::copy_rev(poly1));
            poly1 = poly1_ptr.get();
        }

        const geom_e2d_list* head1 = poly1;
        do {
            if (geom_e2d_list::contains(poly2, poly1->point)) {
                return boolean_generic(poly1, poly2);
            }
        } while (geom_e2d_list::move(poly1, head1));
        abort();
    }

    static geom_e2d_list* boolean_minus(const geom_e2d_list* poly1, const geom_e2d_list* poly2)
    {
        std::unique_ptr<geom_e2d_list> poly1_ptr{};
        geom_real_t area1 = geom_e2d_list::area(poly1);
        geom_real_t area2 = geom_e2d_list::area(poly2);
        if (area1 * area2 > 0.0) {
            poly1_ptr.reset(geom_e2d_list_factory::copy_rev(poly1));
            poly1 = poly1_ptr.get();
        }

        const geom_e2d_list* head1 = poly1;
        do {
            if (!geom_e2d_list::contains(poly2, poly1->point)) {
                return boolean_generic(poly1, poly2);
            }
        } while (geom_e2d_list::move(poly1, head1));
        abort();
    }
};  // struct geom_e2d_list_factor

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

enum geom_c2d_type
{
    GEOM_C2D_NONE,
    GEOM_C2D_TOUCH,
    GEOM_C2D_INTERSECT,
};	// enum geom_c2d_type

enum geom_c2d_dirc_type
{
    GEOM_C2D_DIRECTION_NONE,
    GEOM_C2D_DIRECTION_IN,
    GEOM_C2D_DIRECTION_OUT,
};  // enum geom_c2d_dirc_type

struct geom_c2d_list
{
    geom_c2d_type type{};
    geom_c2d_dirc_type direction{};
    geom_e2d e{};
    struct geom_e2d_list* c1 = nullptr;
    struct geom_e2d_list* c2 = nullptr;
    geom_c2d_list* next = nullptr;
    ~geom_c2d_list()
    {
        delete next;
    }
};  // struct geom_c2d_list

inline
geom_c2d_type
geom_collide(const geom_e2d &e1, const geom_e2d &e2, geom_e2d& e_int)
{
    extern void
    geom_collide(geom_e2d, geom_e2d, geom_c2d_list*);

    geom_c2d_list collision{};
    geom_collide(e1, e2, &collision);
    e_int = collision.e;
    return collision.type;
}
