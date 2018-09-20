#pragma once
#include <ostream>
#include <array>
#include <fstream>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <memory>


using geom_real_t = double;

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

struct geom_p2d
{
	geom_real_t x = 0.0;
	geom_real_t y = 0.0;

    bool operator == (const geom_p2d& p) const
    {
        return (x == p.x) && (y == p.y);
    }
    bool operator != (const geom_p2d& p) const
	{
		return (x != p.x) || (y != p.y);
	}

    geom_p2d operator+ (const geom_p2d& p) const
    {
        return { x + p.x, y + p.y };
    }
	geom_p2d operator- (const geom_p2d& p) const
	{
		return { x - p.x, y - p.y };
	}

    geom_p2d operator* (geom_real_t a) const
    {
        return {a * x, a * y};
    }
	friend geom_p2d operator* (geom_real_t a, const geom_p2d& p)
	{
		return {a * p.x, a * p.y};
	}

    static geom_real_t dot(const geom_p2d& p1, const geom_p2d& p2)
    {
        return p1.x * p2.x + p1.y * p2.y;
    }
    static geom_real_t len(const geom_p2d& p1)
    {
        return sqrt(dot(p1, p1));
    }

    static geom_real_t det(const geom_p2d& v1, const geom_p2d& v2)
    {
        return v1.x * v2.y - v2.x * v1.y ;
    }

	static geom_p2d max(const geom_p2d& p1, const geom_p2d& p2)
    {
        return { std::max(p1.x, p2.x), std::max(p1.y, p2.y) };
    }
    static geom_p2d min(const geom_p2d& p1, const geom_p2d& p2)
    {
        return { std::min(p1.x, p2.x), std::min(p1.y, p2.y) };
    }

    static geom_p2d normal(const geom_p2d& p1)
    {
	    return {p1.y, -p1.x};
    }
};	// struct geom_p2d

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

struct geom_e2d
{
	geom_p2d s{};
	geom_p2d t{};
};	// struct geom_e2d

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
