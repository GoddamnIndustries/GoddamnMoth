#pragma once

#include <algorithm>
#include <ostream>
#include <istream>
#include <memory>

#include <cassert>
#include <cmath>

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

using geom_real_t = double;

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
struct geom_p2d final
{
	geom_real_t x{};
	geom_real_t y{};

public:
    geom_p2d operator+(const geom_p2d& p) const
    {
        return { x + p.x, y + p.y };
    }
	geom_p2d operator-(const geom_p2d& p) const
	{
		return { x - p.x, y - p.y };
	}

    geom_p2d operator*(geom_real_t a) const
    {
        return {a * x, a * y};
    }
	friend geom_p2d operator*(geom_real_t a, const geom_p2d& p)
	{
		return {a * p.x, a * p.y};
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
    friend std::ostream& operator<<(std::ostream& stream, const geom_p2d& p);
    friend std::istream& operator>>(std::istream& stream, geom_p2d& p);

    static std::string str(const geom_p2d& p);

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
	    return {p1.y, -p1.x};
    }
};	// struct geom_p2d

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
    friend std::ostream& operator<<(std::ostream& stream, const geom_e2d& e);
    friend std::istream& operator>>(std::istream& stream, geom_e2d& e);

    static std::string str(const geom_e2d& e);

public:
    static geom_real_t det(const geom_e2d& e)
    {
        return geom_p2d::det(e.t, e.s);
    }

public:
    static geom_real_t len(const geom_e2d& e)
    {
        return geom_p2d::len(e.t - e.s);
    }
};	// struct geom_e2d

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

///
/// List of edges in 2D space.
///
struct geom_e2d_list final
{
private:
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
        geom_e2d_list* poly = this;
        do {
            remove_nxt(poly);
        } while (poly->next != this);
    }

public:
    friend std::ostream& operator<< (std::ostream& stream, const geom_e2d_list* poly);
    friend std::istream& operator>> (std::istream& stream, geom_e2d_list* poly);

    static std::string str(const geom_e2d_list* poly);

public:
    /// Returns copy of the polygon.
    static geom_e2d_list* cpy(const geom_e2d_list* poly)
    {
        geom_e2d_list* copy = nullptr;
        const geom_e2d_list* head = poly;
        do {
            if (copy == nullptr) {
                copy = new geom_e2d_list{poly->point};
            } else {
                copy->insert_next(poly->point);
            }
        } while (move_next(copy), move_next(poly, head));
        move_next(copy);
        return copy;
    }

    /// Returns reverse copy of the polygon.
    static geom_e2d_list* rev(const geom_e2d_list* poly)
    {
        geom_e2d_list* copy = nullptr;
        const geom_e2d_list* head = poly;
        do {
            if (copy == nullptr) {
                copy = new geom_e2d_list{poly->point};
            } else {
                copy->insert_next(poly->point);
            }
        } while (move_next(poly, head));
        return copy;
    }

public:
    geom_e2d edge() const
    {
        return {point, next->point};
    }

public:
    /// Returns perimeter of the polygon.
    static geom_real_t len(const geom_e2d_list* poly)
    {
        geom_real_t length = 0.0;
        const geom_e2d_list* head = poly;
        do {
            length += geom_e2d::len(poly->edge());
        } while (move_next(poly, head));
        return length;
    }

    /// Returns signed area of the polygon.
    static geom_real_t sqr(const geom_e2d_list* poly)
    {
        geom_real_t square = 0.0;
        const geom_e2d_list* head = poly;
        do {
            square += geom_e2d::det(poly->edge()) * 0.5;
        } while (move_next(poly, head));
        return square;
    }

public:
    static bool move_next(geom_e2d_list*& list, const geom_e2d_list* head = nullptr)
    {
        list = list->next;
        return list != head;
    }
    static bool move_prev(geom_e2d_list*& list, const geom_e2d_list* head = nullptr)
    {
        for (const geom_e2d_list* last = list; list->next != last; move_next(list));
        return list != head;
    }

    static bool move_next(const geom_e2d_list*& list, const geom_e2d_list* head = nullptr)
    {
        return move_next(const_cast<geom_e2d_list*&>(list), head);
    }
    static bool move_prev(const geom_e2d_list*& list, const geom_e2d_list* head = nullptr)
    {
        return move_prev(const_cast<geom_e2d_list*&>(list), head);
    }

public:
    void insert_next(const geom_p2d& p)
    {
        geom_e2d_list* list = this;
        auto node = new geom_e2d_list{p};
        node->next = list->next;
        list->next = node;
    }
    void insert_prev(const geom_p2d& p)
    {
        geom_e2d_list* list = this;
        move_prev(list);
        list->insert_next(p);
    }

    static void remove_nxt(geom_e2d_list* list)
    {
        auto node = list->next;
        list->next = node->next;
        node->next = node;
        delete node;
    }
    static void remove_prv(geom_e2d_list* list)
    {
        move_prev(list);
        move_prev(list);
        remove_nxt(list);
    }

private:
    static bool contains(const geom_e2d_list* E, const geom_p2d& p);
    static void clip(geom_e2d_list* E1, geom_e2d_list* E2);
};  // struct geom_e2d_list

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
    geom_collide(const geom_e2d&, const geom_e2d&, geom_c2d_list*);

    geom_c2d_list collision{};
    geom_collide(e1, e2, &collision);
    e_int = collision.e;
    return collision.type;
}
