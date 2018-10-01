#pragma once

#include "GeomBase.hh"
#include "GeomPoint.hh"

/**
 * Edge (half-open interval) in 2D space:
 * @f$
 * \hat{e_1}(\vec{s}, \vec{t}) := {\vec{p}: \vec{p} = \alpha\vec{s} + (1-\alpha)\vec{t}, 0 < \alpha \le 1}.
 * @f$
 */
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
    /**
     * Calculate intersection of the two edges.
     */
    GEOM_HOST GEOM_DEVICE
    static bool intersect(const geom_e2d& e1, geom_e2d e2, geom_e2d& e)
    {
        geom_p2d v1 = e1.t - e1.s;
        geom_p2d v2 = e2.t - e2.s;
        geom_real_t det = geom_p2d::det(v1, v2);
        if(det != 0.0) {
            /* Edges are not collinear. */
            geom_p2d v = e2.s - e1.s;
            geom_real_t alpha = geom_p2d::det(v, v2) / det;
            geom_real_t gamma = geom_p2d::det(v, v1) / det;
            if ((0.0 < alpha && alpha <= 1.0) && (0.0 < gamma && gamma <= 1.0)) {
                /* Edges intersect on point. */
                e.s = e1.s + alpha * v1;
                e.t = e.s;
                return true;
            }
        } else {
            /* Edges are collinear. */
            if (geom_p2d::dot(v1, v2) < 0.0) {
                std::swap(e2.s, e2.t);
            }

            geom_p2d q1 = e2.s - e1.t;
            geom_p2d q2 = e2.t - e1.s;
            if (geom_p2d::det(q1, q2) == 0.0) {
                /* Edges are on one line. */
                if (geom_p2d::len(v1) + geom_p2d::len(v2) != fabs(geom_p2d::len(q2) - geom_p2d::len(q1))) {
                    /* Edges intersect on line. */
                    e.s = geom_p2d::max(geom_p2d::min(e2.s, e2.t), geom_p2d::min(e1.s, e1.t));
                    e.t = geom_p2d::min(geom_p2d::max(e2.s, e2.t), geom_p2d::max(e1.s, e1.t));
                    return true;
                }
            }
        }
        return false;
    }
    GEOM_HOST GEOM_DEVICE
    static bool intersect(const geom_e2d& e1, geom_e2d e2, geom_p2d& p)
    {
        geom_e2d e;
        if (intersect(e1, e2, e)) {
            if (e.s == e.t) {
                p = e.s;
                return true;
            }
        }
        return false;
    }

public:
    /**
     * Reflect trajectory on the diode wall.
     */
    GEOM_HOST GEOM_DEVICE
    static bool reflect(const geom_e2d& e1, const geom_e2d& e2, geom_e2d& e)
    {
        if (intersect(e1, e2, e) && e.s == e.t) {
            /* Edges intersect and are not collinear. */
            geom_p2d p = e.s;
            geom_real_t l1 = geom_p2d::len(e2.s - p);
            geom_real_t l2 = geom_p2d::len(e2.t - p);
            geom_p2d v = (e2.s - p) / l1;
            geom_p2d n = geom_p2d::normal(e1.s - p);
            geom_real_t sin = geom_p2d::det(v, n);
            geom_real_t cos = geom_p2d::dot(v, n);
            if (cos > 0.0) {
                /* Edge reflect from the opaque side of the diode wall. */
                geom_p2d q = geom_p2d::rotate(n, sin, cos);
                e = {e2.s, p + l2 * q};
            } else {
                /* Edge goes through the transparent side of the diode wall. */
                e = {e2.s, e2.t};
            }
            return true;
        } else {
            e = {e2.s, e2.t};
        }
        return false;
    }

public:
    GEOM_HOST
    static std::string str(const geom_e2d& e);
};	// struct geom_e2d

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

GEOM_HOST GEOM_CORE
std::ostream& operator<<(std::ostream& stream, const geom_e2d& e);
GEOM_HOST GEOM_CORE
std::istream& operator>>(std::istream& stream, geom_e2d& e);

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * List of 2D edges.
 */
struct GEOM_CORE geom_e2d_list final
{
public:
    geom_p2d point{};
    geom_e2d_list* next{this};
    geom_e2d_list* next_int{};
    geom_e2d_list* next_ext{};

public:
    GEOM_HOST
    geom_e2d_list(geom_e2d_list&&) = delete;
    GEOM_HOST
    geom_e2d_list(const geom_e2d_list&) = delete;
    GEOM_HOST
    geom_e2d_list& operator= (geom_e2d_list&&) = delete;
    GEOM_HOST
    geom_e2d_list& operator= (const geom_e2d_list&) = delete;

public:
    GEOM_HOST
    explicit geom_e2d_list(const geom_p2d& p): point(p) {}
    GEOM_HOST
    ~geom_e2d_list()
    {
        if (next != this) {
            do {
                remove();
            } while (next != this);
        }
    }

public:
    GEOM_HOST
    geom_e2d edge() const
    {
        return {point, next->point};
    }

public:
    GEOM_HOST
    static bool move(geom_e2d_list*& list, const geom_e2d_list* head = nullptr)
    {
        list = list->next;
        return list != head;
    }
    GEOM_HOST
    static bool move(const geom_e2d_list*& list, const geom_e2d_list* head = nullptr)
    {
        return move(const_cast<geom_e2d_list*&>(list), head);
    }

public:
    GEOM_HOST
    void insert(const geom_p2d& p)
    {
        geom_e2d_list* list = this;
        geom_e2d_list* node = new geom_e2d_list{p};
        node->next = list->next;
        list->next = node;
    }
    GEOM_HOST
    void remove()
    {
        geom_e2d_list* list = this;
        geom_e2d_list* node = list->next;
        list->next = node->next;
        node->next = node;
        delete node;
    }

public:
    GEOM_HOST
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
    GEOM_HOST
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
    GEOM_HOST
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
    GEOM_HOST [[deprecated]]
    static bool contains(const geom_e2d_list* poly, const geom_p2d& p)
    {
        geom_size_t intersections = 0;
        const geom_e2d_list* head = poly;
        do {
            geom_e2d e;
            geom_e2d e1 = poly->edge();
            geom_e2d e2 = {p, p + geom_p_inf};
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

public:
    GEOM_HOST
    static void clip(geom_e2d_list* list1, geom_e2d_list* list2)
    {
        geom_e2d_list* head1 = list1;
        do {
            geom_e2d_list* head2 = list2;
            do {
                geom_p2d p;
                geom_e2d e1 = list1->edge();
                geom_e2d e2 = list2->edge();
                if (geom_e2d::intersect(e1, e2, p)) {
                    bool n0 = contains(list1, e2.s);
                    bool n1 = contains(list1, e2.t);
                    //assert(n0 ^ n1);
                    list1->insert(p);
                    list2->insert(p);
                    if (n0 && !n1) {
                        list1->next->next_ext = list2->next->next;
                        list2->next->next_int = list1->next->next;
                    } else {
                        list1->next->next_int = list2->next->next;
                        list2->next->next_ext = list1->next->next;
                    }
                    break;
                }
            } while (move(list2, head2));
        } while (move(list1, head1));
    }

public:
    GEOM_HOST [[deprecated]]
    static bool reflect(const geom_e2d_list* poly, const geom_e2d& e2, geom_e2d e)
    {
        /// @todo Incorrect!
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
    GEOM_HOST
    static std::string str(const geom_e2d_list* poly);
};  // struct geom_e2d_list

GEOM_HOST GEOM_CORE
std::ostream& operator<< (std::ostream& stream, const geom_e2d_list* poly);
GEOM_HOST GEOM_CORE
std::istream& operator>> (std::istream& stream, geom_e2d_list*& poly);

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

///
/// 2D polygon factory.
///
struct GEOM_CORE geom_e2d_list_factory final
{
public:
    /**
     * Reverse copy the polygon.
     */
    GEOM_HOST [[deprecated]]
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
    /**
     * Rectangle with given south-west and north-east points.
     */
    GEOM_HOST [[deprecated]]
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

    /**
     * Circle of the given radius with CCW orientation.
     */
    GEOM_HOST [[deprecated]]
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
};  // struct geom_e2d_list_factory

