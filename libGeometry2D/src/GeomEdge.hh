#pragma once
#include <ostream>
#include <array>
#include <fstream>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <memory>


using geom_real_t = double;
template <typename T>
int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

///
/// Point in 2D space:
/// @f$ (x, y)^T @f$.
///
struct geom_p2d
{
	geom_real_t x;
	geom_real_t y;

	geom_p2d(geom_real_t X = 0.0, geom_real_t Y = 0.0): x(X), y(Y)
	{}

	/// @todo To be removed.
	bool operator != (const geom_p2d& p) const
	{
		return (x != p.x) || (y != p.y);
	}
	geom_p2d operator- (const geom_p2d& p) const
	{
		return { x - p.x, y - p.y };
	}
	geom_p2d operator+ (const geom_p2d& p) const
	{
		return { x + p.x, y + p.y };
	}
	friend geom_p2d operator*(geom_real_t a, const geom_p2d& p)
	{
		return {a * p.x, a * p.y};
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
};	// struct geom_p2d

//static geom_real_t geom_dot(const geom_p2d&, const geom_p2d&);

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

///
/// Edge in 2D space:
/// @f$ [\vec s, \vec t] := \alpha \vec s + (1 - \alpha) \vec t @f$.
///
struct geom_e2d
{
	geom_p2d s;
	geom_p2d t;

	geom_e2d(): s(), t()
	{}
	geom_e2d(geom_p2d const& point1, geom_p2d const& point2): s(point1), t(point2)
	{}
	/// @todo To be removed.
	bool if_point_on_segment(geom_p2d const& p) const
	{
		auto const det = (t.x - s.x) * (p.y - s.y) - (p.x - s.x) * (t.y - s.y);
		if(det != 0) // p is not on the line [s,t]
			return false;

        return (p.x >= std::min(s.x, t.x)) && (p.x <= std::max(s.x, t.x)) &&
               (p.y >= std::min(s.y, t.y)) && (p.y <= std::max(s.y, t.y));
    }

};	// struct geom_e2d

struct geom_e2d_list : public geom_p2d
{
    geom_e2d_list* next = this;
    geom_e2d_list* next_in = nullptr;
    geom_e2d_list* next_out = nullptr;
    geom_e2d edge() const { return { *this, *next }; }
    void insert(const geom_p2d& p)
    {
        auto* n = new geom_e2d_list{};
        (geom_p2d&)*n = p;
        n->next = next;
        next = n;
    }


    static int orientation(const geom_e2d_list* E);

    static bool contains(const geom_e2d_list* E, const geom_p2d& p);

    static void clip(geom_e2d_list* E1, geom_e2d_list* E2);
};  // struct geom_e2d_list

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
    geom_c2d_list* next = nullptr;
    ~geom_c2d_list()
    {
        delete next;
    }
};  // struct geom_c2d_list

void
geom_collide(const geom_e2d &e1, const geom_e2d &e2, geom_c2d_list* collision);

inline
geom_c2d_type
geom_collide(const geom_e2d &e1, const geom_e2d &e2, geom_e2d& e_int)
{
    geom_c2d_list collision{};
    geom_collide(e1, e2, &collision);
    e_int = collision.e;
    return collision.type;
}
