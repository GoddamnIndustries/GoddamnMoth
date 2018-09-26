#pragma once

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

#define GEOM_PI (4 * atan(1.0))

using geom_real_t = double;
using geom_size_t = size_t;

template <typename T>
int fsgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

///
/// Point in 2D space.
///
struct geom_p2d
{
	geom_real_t x{};
	geom_real_t y{};
    geom_real_t u{};
    geom_real_t v{};

public:
    geom_p2d operator+() const
    {
        return { x, y };
    }
    geom_p2d operator+(const geom_p2d& p) const
    {
        return { x + p.x, y + p.y };
    }

    geom_p2d operator-() const
    {
        return { -x, -y };
    }
	geom_p2d operator-(const geom_p2d& p) const
	{
		return { x - p.x, y - p.y };
	}

public:
    geom_p2d operator*(geom_real_t a) const
    {
        return { a * x, a * y };
    }
	friend geom_p2d operator*(geom_real_t a, const geom_p2d& p)
	{
		return { a * p.x, a * p.y };
	}

    geom_p2d operator/(geom_real_t a) const
    {
        return { x / a, y / a };
    }
    friend geom_p2d operator/(geom_real_t a, const geom_p2d& p)
    {
        return { p.x / a, p.y / a };
    }

public:
    bool operator==(const geom_p2d& p) const
    {
        return (x == p.x) && (y == p.y);
    }
    bool operator!=(const geom_p2d& p) const
    {
        return (x != p.x) || (y != p.y);
    }

public:
    static geom_real_t dot(const geom_p2d& p1, const geom_p2d& p2)
    {
        return p1.x * p2.x + p1.y * p2.y;
    }
    static geom_real_t len(const geom_p2d& p1)
    {
        return sqrt(dot(p1, p1));
    }

public:
    static geom_real_t det(const geom_p2d& v1, const geom_p2d& v2)
    {
        return v1.x * v2.y - v2.x * v1.y;
    }
    static geom_real_t ang(const geom_p2d& v1, const geom_p2d& v2)
    {
        return atan2(det(v1, v2), dot(v1, v2));
    }

public:
	static geom_p2d max(const geom_p2d& p1, const geom_p2d& p2)
    {
        return { std::max(p1.x, p2.x), std::max(p1.y, p2.y) };
    }
    static geom_p2d min(const geom_p2d& p1, const geom_p2d& p2)
    {
        return { std::min(p1.x, p2.x), std::min(p1.y, p2.y) };
    }

public:
    static geom_p2d normal(const geom_p2d& p1)
    {
        geom_p2d n{p1.y, -p1.x};
        return n / len(n);
    }

public:
    static geom_p2d rot(const geom_p2d& n, geom_real_t a)
    {
        geom_real_t c = cos(a);
        geom_real_t s = sin(a);
        return{ n.x * c - n.y * s, n.x * s + n.y * c };
    }

public:
    static std::string str(const geom_p2d& p);
};	// struct geom_p2d

std::ostream& operator<<(std::ostream& stream, const geom_p2d& p);
std::istream& operator>>(std::istream& stream, geom_p2d& p);

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

///
/// Edge in 2D space.
///
struct geom_e2d final
{
	geom_p2d s{};
	geom_p2d t{};

public:
    geom_e2d operator+(const geom_p2d& p) const
    {
        return { s + p, t + p };
    }
    geom_e2d operator-(const geom_p2d& p) const
    {
        return { s - p, t - p };
    }

public:
    bool operator==(const geom_e2d& e) const
    {
        return (s == e.s) && (t == e.t);
    }
    bool operator!=(const geom_e2d& e) const
    {
        return (s != e.s) || (t != e.t);
    }

public:
    static geom_real_t dot(const geom_e2d& e)
    {
        return geom_p2d::dot(e.t, e.s);
    }
    static geom_real_t dot(const geom_e2d& e1, const geom_e2d& e2)
    {
        return geom_p2d::dot(e1.t - e1.s, e2.t - e2.s);
    }

    static geom_real_t len(const geom_e2d& e)
    {
        return geom_p2d::len(e.t - e.s);
    }

public:
    static geom_real_t det(const geom_e2d& e, const geom_p2d& p = {})
    {
        return geom_p2d::det(e.t - p, e.s - p);
    }
    static geom_real_t det(const geom_e2d& e1, const geom_e2d& e2)
    {
        return geom_p2d::det(e1.t - e1.s, e2.t - e2.s);
    }
    static geom_real_t ang(const geom_e2d& e1, const geom_e2d& e2)
    {
        return geom_p2d::ang(e1.t - e1.s, e2.t - e2.s);
    }

public:
    /// Calculates intersection of two edges.
    static bool intersect(const geom_e2d& e1, geom_e2d e2, geom_e2d& e)
    {
        geom_p2d v1 = e1.t - e1.s;
        geom_p2d v2 = e2.t - e2.s;
        geom_real_t det = geom_p2d::det(v1, v2);
        if(det != 0.0) {
            // Edges are not collinear.
            geom_p2d b = e2.s - e1.s;
            geom_real_t alpha = geom_p2d::det(b, v2) / det;
            geom_real_t gamma = geom_p2d::det(b, v1) / det;
            if (alpha >= 0 && alpha <= 1 && gamma >= 0 && gamma <= 1) {
                // Edges intersect on point.
                e.s = e1.s + alpha * v1;
                e.t = e.s;
                return true;
            }
        } else {
            // Edges are collinear.
            if (geom_p2d::dot(v1, v2) < 0.0) {
                std::swap(e2.s, e2.t);
                v2 = -v2;
            }

            geom_p2d q1 = e2.s - e1.t;
            geom_p2d q2 = e2.t - e1.s;
            if (geom_p2d::det(q1, q2) == 0.0) {
                // Edges are on one line.
                if (geom_p2d::len(v1) + geom_p2d::len(v2) != fabs(geom_p2d::len(q2) - geom_p2d::len(q1))) {
                    // Edges intersect on line.
                    e.s = geom_p2d::max(geom_p2d::min(e2.s, e2.t), geom_p2d::min(e1.s, e1.t));
                    e.t = geom_p2d::min(geom_p2d::max(e2.s, e2.t), geom_p2d::max(e1.s, e1.t));
                    return true;
                }
            }
        }
        return false;
    }

    /// Calculates reflection of two edges.
    static bool reflect(const geom_e2d& e1, const geom_e2d& e2, geom_e2d& e)
    {
        if (intersect(e1, e2, e) && e.s == e.t) {
            geom_p2d p = e.s;
            geom_real_t l1 = geom_p2d::len(e2.s - p);
            geom_real_t l2 = geom_p2d::len(e2.t - p);
            geom_p2d v = (e2.s - p) / l1;
            geom_p2d n = geom_p2d::normal(e1.s - p);

            geom_real_t cos = geom_p2d::dot(v, n);
            geom_real_t sin = geom_p2d::det(v, n);
            geom_p2d q{n.x * cos - n.y * sin, n.x * sin + n.y * cos};
            geom_p2d w{n.x * cos - n.y * sin, n.x * sin + n.y * cos};
            e.s = e2.s;
            e.t = p + l2 * q;
            return true;
        }
        e = e2;
        return false;
    }

public:
    static std::string str(const geom_e2d& e);
};	// struct geom_e2d

std::ostream& operator<<(std::ostream& stream, const geom_e2d& e);
std::istream& operator>>(std::istream& stream, geom_e2d& e);

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
            area += geom_e2d::det(poly->edge()) * 0.5;
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
            while (geom_e2d::det(hull->edge(), p_i) >= 0.0) {
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
					    auto a = geom_e2d::ang(e1, e2);
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
