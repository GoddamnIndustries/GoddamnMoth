#pragma once
#include <ostream>
#include <array>
#include <fstream>
#include <cassert>

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

	/// @todo To be removed.
	bool operator!= (const geom_p2d& p) const
	{
		return x != p.x && y != p.y;
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

	/// @todo To be removed.
	bool if_point_on_segment(geom_p2d const& point) const
	{
		// Fixed.
		if (q.x - p.x != 0.0) {
			double const alpha = (point.x - p.x) / (q.x - p.x);
			if ((alpha < 0) || (alpha > 1))
				return false;
			else
				return true;
		}
		if (q.y - p.y != 0.0) {
			double const alpha = (point.y - p.y) / (q.y - p.y);
			if ((alpha < 0) || (alpha > 1))
				return false;
			else
				return true;
		}
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
geom_mutual_arrangement(const geom_e2d& e1, const geom_e2d& e2)
{
	/// This classification is needed for correct clipping.
	/// @todo Implements this based on the article https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
#if 0
	inline
		geom_real_t
		geom_orientation(const geom_p2d& p, const geom_p2d& q, const geom_p2d& r)
	{
		const auto prod = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
		return sgn(prod);
	}
	inline
		geom_real_t
		geom_orientation(const geom_e2d& pq, const geom_p2d& r)
	{
		return geom_orientation(pq.p, pq.q, r);
	}
	inline
		geom_real_t
		geom_orientation(const geom_p2d& p, const geom_e2d& qr)
	{
		return geom_orientation(p, qr.p, qr.q);
	}
	inline
		bool
		geom_intersects(const geom_e2d& e1, const geom_e2d& e2)
	{
		const auto orient1 = geom_orientation(e1, e2.p);
		const auto orient2 = geom_orientation(e1, e2.q);
		const auto orient3 = geom_orientation(e2, e1.p);
		const auto orient4 = geom_orientation(e2, e1.q);
		if (orient1 != orient2 &&
			orient3 != orient4) {
			return true;
		}

		return false;
	}
#endif
	return geom_edges_do_not_intersect;
}

static
geom_p2d 
geom_intersects_p(const geom_e2d& e1, const geom_e2d& e2)
{
	std::array<double, 2> c_vec = { e1.q.x - e1.p.x, e1.q.y - e1.p.y };
	std::array<double, 2> d_vec = { e2.q.x - e2.p.x, e2.q.y - e2.p.y };

	auto const det = -c_vec[0] * d_vec[1] + d_vec[0] * c_vec[1];
	auto const b1 = e2.p.x - e1.p.x;
	auto const b2 = e2.p.y - e1.p.y;

	double const alpha = 1 / det * (-b1 * d_vec[1] + b2 * d_vec[0]);

	return { e1.p.x + alpha * c_vec[0], e1.p.y + alpha * c_vec[1] };
}

#if 1
class geom_polygon2d //: public std::enable_shared_from_this<geom_polygon2d>
{
public:
	geom_p2d point;
	geom_polygon2d* next;
	geom_polygon2d* other_in = nullptr;
	geom_polygon2d* other_out = nullptr;


public:

	geom_polygon2d(geom_p2d const& p) : point(p), next(this)
	{}

	void insert(geom_p2d const& p)
	{
		auto new_pol_list = new geom_polygon2d(p);
		new_pol_list->next = next;
		next = new_pol_list;
	}

	geom_e2d make_line() const
	{
		return { point, next->point };
	}

	bool is_internal(const geom_p2d& point) const
	{
		geom_e2d line_from_point{ point, {point.x + 1.0, point.y} };
		unsigned inters_counter = 0;
		auto list_iter = this;
		do
		{
			auto line = list_iter->make_line();
			auto const inters_point = geom_intersects_p(line_from_point, line);
			if ((inters_point.x > point.x) && (line.if_point_on_segment(inters_point)))
				++inters_counter;

			list_iter = list_iter->next;

		} while (list_iter != this);

		if ((inters_counter % 2) == 0)
			return false;
		return true;
	}

	geom_p2d get_point() const
	{
		return point;
	}

	void print(std::string const& filepath) const
	{
		std::ofstream file;
		file.open(filepath);

		auto list_iter = this;

		do
		{
			auto edge = list_iter->make_line();
			file << edge.p.x << " " << edge.p.y << std::endl;
			file << edge.q.x << " " << edge.q.y << std::endl;
			file << std::endl;

			list_iter = list_iter->next;
		} while (list_iter != this);

	}

	void plot() const
	{
		print("A.txt");
		system("gnuplot -e \"plot 'A.txt' with lines; pause -1;\"");
		remove("A.txt");
	}

	void set_other_in(geom_polygon2d* const list)
	{
		other_in = list;
	}
	void set_other_out(geom_polygon2d* const list)
	{
		other_out = list;
	}

	geom_polygon2d* get_other_in() const
	{
		return other_in;
	}

	geom_polygon2d* get_other_out() const
	{
		return other_out;
	}
};

void geom_clip(geom_polygon2d* A, geom_polygon2d* B)
{
	auto A_iter = A;
	do
	{
		auto B_iter = B;
		do
		{
			auto A_edge_seg = A_iter->make_line();
			auto B_edge_seg = B_iter->make_line();

			auto p = geom_intersects_p(A_edge_seg, B_edge_seg);

			bool intersects = A_edge_seg.if_point_on_segment(p);
			intersects &= B_edge_seg.if_point_on_segment(p);
			intersects &= p != A_edge_seg.p;
			intersects &= p != A_edge_seg.q;
			intersects &= p != B_edge_seg.p;
			intersects &= p != B_edge_seg.q;
			if (intersects)
			{
				bool n_0 = A->is_internal(B_edge_seg.p);
				bool n_1 = A->is_internal(B_edge_seg.q);
				//assert(n_0 ^ n_1);

				A_iter->insert(p);
				B_iter->insert(p);

				if (n_0 && !n_1)
				{
					A_iter->next->set_other_out(B_iter->next);
					B_iter->next->set_other_in(A_iter->next);
				}
				else
				{
					A_iter->next->set_other_in(B_iter->next);
					B_iter->next->set_other_out(A_iter->next);
				}

				break;
			}

			B_iter = B_iter->next;
		} while (B_iter != B);


		A_iter = A_iter->next;
	} while (A_iter != A);
}
#endif

void geom_minus(geom_polygon2d* A, geom_polygon2d* B)
{
	auto A_iter = A;
	do
	{
		if (A_iter->other_out != nullptr)
			A_iter->next = A_iter->other_out;

		A_iter = A_iter->next;
	} while (A_iter != A);
}
