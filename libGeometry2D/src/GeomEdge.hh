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
/// @f$ [\vec p, \vec q] := \alpha \vec p + (1 - \alpha) \vec q @f$.
///
struct geom_e2d
{
	geom_p2d p;
	geom_p2d q;

	geom_e2d(): p(), q()
	{}
	geom_e2d(geom_p2d const& point1, geom_p2d const& point2): p(point1), q(point2)
	{}
	/// @todo To be removed.
	bool if_point_on_segment(geom_p2d const& point) const
	{
		auto const det = (q.x - p.x) * (point.y - p.y) - (point.x - p.x) * (q.y - p.y);
		if(det != 0) // point is not on the line [p,q]
			return false;

		if( (point.x >= std::min(p.x, q.x)) && (point.x <= std::max(p.x, q.x)) &&
			(point.y >= std::min(p.y, q.y)) && (point.y <= std::max(p.y, q.y)) )
			return true;
		return false;
	}
};	// struct geom_e2d




enum geom_edge_mutual_arrangement_t
{
	geom_edges_do_not_intersect,
	geom_edges_intersect_in_point,
	geom_edges_intersect_on_segment,
};	// enum geom_edge_mutual_arrangement_t

///
/// Returns of the mutual arrangement of two segments.
///
static
geom_edge_mutual_arrangement_t
geom_mutual_arrangement(const geom_e2d& e1, const geom_e2d& e2, geom_e2d& intersection_segment)
{


	double const c_vec_0 = e1.q.x - e1.p.x;
	double const c_vec_1 = e1.q.y - e1.p.y;
	double const d_vec_0 = e2.q.x - e2.p.x;
	double const d_vec_1 = e2.q.y - e2.p.y;

	auto const det = -c_vec_0 * d_vec_1 + d_vec_0 * c_vec_1;

	if(det != 0)
	{
		//not colinear
		auto const b1 = e2.p.x - e1.p.x;
		auto const b2 = e2.p.y - e1.p.y;

		double const alpha = 1 / det * (-b1 * d_vec_1 + b2 * d_vec_0);
		double const beta = 1 / det * (c_vec_0 * b2 - b1 * c_vec_1);

		if( (alpha >= 0) && (alpha <= 1) && (beta >= 0) && (beta <= 1) )
		{
			geom_p2d const int_p(e1.p.x + alpha * c_vec_0, e1.p.y + alpha * c_vec_1);
			intersection_segment.p = int_p;
			intersection_segment.q = int_p;
			return geom_edges_intersect_in_point;
		}
		else
			return geom_edges_do_not_intersect;
	}
	else
	{
		//colinear

		geom_e2d e2_reor;
		if( ((e1.q.x - e1.p.x) * (e2.q.x - e2.p.x) + (e1.q.y - e1.p.y) * (e2.q.y - e2.p.y)) < 0 )
		{
			e2_reor.p = e2.q;
			e2_reor.q = e2.p;
		}
		else
			e2_reor = e2;
		if(e1.if_point_on_segment(e2_reor.p))
		{
			if(e1.if_point_on_segment(e2_reor.q))
			{
				// p_1 ----- p_2 ------- q_2 ------ q_1
				intersection_segment.p = e2_reor.p;
				intersection_segment.q = e2_reor.q;
				return geom_edges_intersect_on_segment;
			}
			else if(e2_reor.if_point_on_segment(e1.q))
			{
				// p_1 --------- p_2 -------- q_1 --------- q_2
				intersection_segment.p = e2_reor.p;
				intersection_segment.q = e1.q;
				return geom_edges_intersect_on_segment;
			}
			else
				throw 1; /// @todo to be removed, just for debugging
		}
		else if(e2_reor.if_point_on_segment(e1.p))
		{
			if(e2_reor.if_point_on_segment(e1.q))
			{
				// p_2 ----- p_1 ------- q_1 ------ q_2
				intersection_segment.p = e1.p;
				intersection_segment.q = e1.q;
				return geom_edges_intersect_on_segment;
			}
			else if(e1.if_point_on_segment(e2_reor.q))
			{
				// p_2 ------ p_1 ------q_2 ----- q_1
				intersection_segment.p = e1.p;
				intersection_segment.q = e2_reor.q;
				return geom_edges_intersect_on_segment;
			}
			else
				throw 1; /// @todo to be removed, just for debugging

		}
		else
			return geom_edges_do_not_intersect;
	}
}