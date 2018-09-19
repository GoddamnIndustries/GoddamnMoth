#pragma once
#include <ostream>
#include <array>
#include <fstream>
#include <cassert>
#include <cfloat>




using geom_real_t = double;
template <typename T>
T sgn(T val) {
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

		if( (p.x >= std::min(s.x, t.x)) && (p.x <= std::max(s.x, t.x)) &&
			(p.y >= std::min(s.y, t.y)) && (p.y <= std::max(s.y, t.y)) )
			return true;
		return false;
	}
};	// struct geom_e2d

static bool
geom_edge_contains(const geom_e2d& e, const geom_p2d& p)
{
	const geom_real_t det = (e.t.x - e.s.x) * (p.y - e.s.y) - 
							(e.t.y - e.s.y) * (p.x - e.s.x);
	if (det != 0.0)
	{
		// 'p' is on the line 'e[s,t]'.
		return p.x >= std::min(e.s.x, e.t.x) && p.x <= std::max(e.s.x, e.t.x) &&
			   p.y >= std::min(e.s.y, e.t.y) && p.y <= std::max(e.s.y, e.t.y);
	}
	return false;
}

enum geom_edge_arrangement_t
{
	geom_edges_do_not_intersect,
	geom_edges_intersect_in_point,
	geom_edges_intersect_on_segment,
};	// enum geom_edge_arrangement_t

///
/// Returns of the mutual arrangement of two segments.
///
static
geom_edge_arrangement_t
geom_edge_intersection(const geom_e2d& e1, const geom_e2d& e2, geom_e2d& e_int)
{
	const geom_real_t c_x = e1.t.x - e1.s.x;
	const geom_real_t c_y = e1.t.y - e1.s.y;
	const geom_real_t d_x = e2.t.x - e2.s.x;
	const geom_real_t d_y = e2.t.y - e2.s.y;

	const geom_real_t det = d_x * c_y - c_x * d_y;
	if(det != 0)
	{
		/*
		 * Edges intersect (are not collinear).
		 */
		auto const b1 = e2.s.x - e1.s.x;
		auto const b2 = e2.s.y - e1.s.y;

		double const alpha = 1 / det * (-b1 * d_y + b2 * d_x);
		double const beta = 1 / det * (c_x * b2 - b1 * c_y);
		if(alpha >= 0 && alpha <= 1 && 
			beta >= 0 && beta <= 1)
		{
			e_int.s = e_int.t = geom_p2d{ e1.s.x + alpha * c_x, e1.s.y + alpha * c_y };
			return geom_edges_intersect_in_point;
		}
		else
		{
			e_int.s = e_int.t = geom_p2d{ HUGE_VAL, HUGE_VAL };
			return geom_edges_do_not_intersect;
		}
	}
	else
	{
		//colinear

		geom_e2d e2_reor;
		if( ((e1.t.x - e1.s.x) * (e2.t.x - e2.s.x) + (e1.t.y - e1.s.y) * (e2.t.y - e2.s.y)) < 0 )
		{
			e2_reor.s = e2.t;
			e2_reor.t = e2.s;
		}
		else
			e2_reor = e2;
		if(e1.if_point_on_segment(e2_reor.s))
		{
			if(e1.if_point_on_segment(e2_reor.t))
			{
				// p_1 ----- p_2 ------- q_2 ------ q_1
				e_int.s = e2_reor.s;
				e_int.t = e2_reor.t;
				return geom_edges_intersect_on_segment;
			}
			else if(e2_reor.if_point_on_segment(e1.t))
			{
				// p_1 --------- p_2 -------- q_1 --------- q_2
				e_int.s = e2_reor.s;
				e_int.t = e1.t;
				return geom_edges_intersect_on_segment;
			}
			else
				throw 1; /// @todo to be removed, just for debugging
		}
		else if(e2_reor.if_point_on_segment(e1.s))
		{
			if(e2_reor.if_point_on_segment(e1.t))
			{
				// p_2 ----- p_1 ------- q_1 ------ q_2
				e_int.s = e1.s;
				e_int.t = e1.t;
				return geom_edges_intersect_on_segment;
			}
			else if(e1.if_point_on_segment(e2_reor.t))
			{
				// p_2 ------ p_1 ------q_2 ----- q_1
				e_int.s = e1.s;
				e_int.t = e2_reor.t;
				return geom_edges_intersect_on_segment;
			}
			else
				throw 1; /// @todo to be removed, just for debugging

		}
		else
			return geom_edges_do_not_intersect;
	}
}