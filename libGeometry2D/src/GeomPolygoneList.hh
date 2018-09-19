#pragma once
#include "GeomEdge.hh"



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

	void insert_back(geom_p2d const& p)
	{
		auto list_iter = this;
		while(list_iter->next != this)
		{
			list_iter = list_iter->next;
		}
		list_iter->insert(p);
	}

	geom_e2d make_line() const
	{
		return { point, next->point };
	}

	bool is_internal(const geom_p2d& point) const
	{
		double d = -1e6;
		auto list_iter = this;
		do
		{
			if(list_iter->point.x > d)
				d = list_iter->point.x;
			list_iter = list_iter->next;
		} while (list_iter != this);
		geom_e2d line_from_point{ point, {point.x + d + 1.0, point.y} };
		unsigned inters_counter = 0;

		list_iter = this;
		do
		{
			auto line = list_iter->make_line();
//			auto const inters_point = geom_intersects_p(line_from_point, line);
//			if ((inters_point.x > point.x) && (line.if_point_on_segment(inters_point)))
//				++inters_counter;

			geom_e2d int_e;
			auto const mutual_arrangement_type = geom_edge_intersection(line_from_point, line, int_e);
			if(mutual_arrangement_type == geom_edges_intersect_in_point)
				++inters_counter;

			list_iter = list_iter->next;

		} while (list_iter != this);

		if ((inters_counter % 2) == 0)
			return false;
		return true;
	}


	void print(std::string const& filepath) const
	{
		std::ofstream file;
		file.open(filepath);

		auto list_iter = this;

		do
		{
			auto edge = list_iter->make_line();
			file << edge.s.x << " " << edge.s.y << std::endl;
			file << edge.t.x << " " << edge.t.y << std::endl;
			file << std::endl;

			list_iter = list_iter->next;
		} while (list_iter != this);

	}

	void plot() const
	{
		print("A.txt");
		system("gnuplot -e \"plot 'A.txt' with lines; pause -1;\"");
//		remove("A.txt");
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
			auto const A_edge_seg = A_iter->make_line();
			auto const B_edge_seg = B_iter->make_line();

//			auto s = geom_intersects_p(A_edge_seg, B_edge_seg);

			geom_e2d int_edge;

			auto const mutual_arr_type = geom_edge_intersection(A_edge_seg, B_edge_seg, int_edge);

//			bool intersects = A_edge_seg.if_point_on_segment(s);
//			intersects &= B_edge_seg.if_point_on_segment(s);

			bool intersects = false;
			geom_p2d p;
			if(mutual_arr_type == geom_edges_intersect_in_point)
			{
				p = int_edge.s;
				intersects = true;
				intersects &= p != A_edge_seg.s;
				intersects &= p != A_edge_seg.t;
				intersects &= p != B_edge_seg.s;
				intersects &= p != B_edge_seg.t;
			}
			if (intersects)
			{
				bool const n_0 = A->is_internal(B_edge_seg.s);
				if((B_edge_seg.t.x == 3) && (B_edge_seg.t.y == 1))
				{
					double a = 0.0;
				}
				bool const n_1 = A->is_internal(B_edge_seg.t);
				assert(n_0 ^ n_1);

				A_iter->insert(p);
				B_iter->insert(p);

				if (n_0 && !n_1)
				{
					A_iter->next->other_out = B_iter->next;
					B_iter->next->other_in = A_iter->next;
				}
				else
				{
					A_iter->next->other_in = B_iter->next;
					B_iter->next->other_out = A_iter->next;
				}

				break;
			}

			B_iter = B_iter->next;
		} while (B_iter != B);


		A_iter = A_iter->next;
	} while (A_iter != A);
}
#endif


void geom_union(geom_polygon2d *A, geom_polygon2d *B)
{
	auto A_iter = A;
	do
	{
		if (A_iter->other_out != nullptr)
			A_iter->next = A_iter->other_out;

		A_iter = A_iter->next;
	} while (A_iter != A);
}
