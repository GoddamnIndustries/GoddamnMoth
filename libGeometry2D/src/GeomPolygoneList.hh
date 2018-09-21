#pragma once
#include "GeomEdge.hh"



#if 1
class geom_vertex2d //: public std::enable_shared_from_this<geom_vertex2d>
{
public:
	geom_p2d point;
	geom_vertex2d* next;
	geom_vertex2d* other_in = nullptr;
	geom_vertex2d* other_out = nullptr;


public:

	geom_vertex2d(geom_p2d const& p) : point(p), next(this)
	{}


	void insert_next(geom_p2d const &p)
	{
		auto new_pol_list = new geom_vertex2d(p);
		new_pol_list->next = next;
		next = new_pol_list;
	}

	void remove_next()
	{
		assert(next != this);

		auto node = next;
		next = next->next;

		node->next = node;
		delete node;
	}


	geom_e2d make_line() const
	{
		return { point, next->point };
	}

};

class geom_polygon2d
{

public:
	geom_vertex2d* head;


public:

	geom_polygon2d(): head(nullptr)
	{}

	geom_polygon2d(geom_p2d const& p): head(new geom_vertex2d(p))
	{
		head->next = head;
	}

	geom_polygon2d(geom_polygon2d const& other): head(nullptr)
	{
		auto iter = other.head;
		do
		{
			// Only non-clipped polygons can be copied
			assert(iter->other_in == nullptr);
			assert(iter->other_out == nullptr);
			insert_back(iter->point);
			iter = iter->next;

		} while(iter != other.head);
	}

	geom_polygon2d& operator=(geom_polygon2d const& other)
	{
		if(head != nullptr)
		{
			while (head->next != head)
				head->remove_next();
			delete head;
			head = nullptr;
		}
		auto iter = other.head;
		do
		{
			// Only non-clipped polygons can be copied
			assert(iter->other_in == nullptr);
			assert(iter->other_out == nullptr);
			insert_back(iter->point);
			iter = iter->next;

		} while(iter != other.head);

		return (*this);

	}

	~geom_polygon2d()
	{
		if(head != nullptr)
		{
			while(head->next != head)
			{
				head->remove_next();
			}
			delete head;
		}
	}


	void insert_back(geom_p2d const& p)
	{

		if(head != nullptr)
		{
			auto list_iter = head;
			while (list_iter->next != head)
			{
				list_iter = list_iter->next;

			}
			list_iter->insert_next(p);
		}
		else
		{
			head = new geom_vertex2d(p);
			head->next = head;
		}
	}

	bool is_internal(const geom_p2d& point) const
	{
		double d = -1e6;
		auto list_iter = head;
		do
		{
			if(list_iter->point.x > d)
				d = list_iter->point.x;
			list_iter = list_iter->next;
		} while (list_iter != head);
		geom_e2d line_from_point{ point, {point.x + d + 1.0, point.y} };
		unsigned inters_counter = 0;

		list_iter = head;
		do
		{
			auto line = list_iter->make_line();
//			auto const inters_point = geom_intersects_p(line_from_point, line);
//			if ((inters_point.x > point.x) && (line.if_point_on_segment(inters_point)))
//				++inters_counter;

			geom_e2d int_e;
			auto const mutual_arrangement_type = geom_collide(line_from_point, line, int_e);
			if(mutual_arrangement_type == GEOM_C2D_INTERSECT)
				++inters_counter;

			list_iter = list_iter->next;

		} while (list_iter != head);

		if ((inters_counter % 2) == 0)
			return false;
		return true;
	}



	void print(std::string const& filepath) const
	{
		std::ofstream file;
		file.open(filepath);

		auto list_iter = head;

		do
		{
			auto edge = list_iter->make_line();
			file << edge.s.x << " " << edge.s.y << std::endl;
			file << edge.t.x << " " << edge.t.y << std::endl;
			file << std::endl;

			list_iter = list_iter->next;
		} while (list_iter != head);

	}

	void plot() const
	{
		print("A.txt");
		system("gnuplot -e \"plot 'A.txt' with lines; pause -1;\"");
		remove("A.txt");
	}


};

void geom_copy_polygon(geom_polygon2d const& source, geom_polygon2d& dest)
{
	auto iter = source.head;
	assert(dest.head == nullptr);
	do
	{
		// Only non-clipped polygons can be copied
		assert(iter->other_in == nullptr);
		assert(iter->other_out == nullptr);
		dest.insert_back(iter->point);
		iter = iter->next;

	} while(iter != source.head);

}

template <typename T>
int sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}

int
geom_orientation(geom_polygon2d const& E)
{
	geom_real_t E_square = 0.0;
	const geom_vertex2d* E_i = E.head;
	do
	{
		const geom_vertex2d* E_ip = E.head->next;
		E_square += geom_p2d::det(E_i->point, E_ip->point);
	} while ((E_i = E_i->next) != E.head);
	return sgn(E_square);
}

void geom_clip(geom_polygon2d const& A, geom_polygon2d const& B, geom_polygon2d& A_res, geom_polygon2d& B_res)
{
	geom_copy_polygon(A, A_res);
	geom_copy_polygon(B, B_res);

	auto A_iter = A_res.head;
	do
	{
		auto B_iter = B_res.head;
		do
		{
			auto const A_edge_seg = A_iter->make_line();
			auto const B_edge_seg = B_iter->make_line();


			geom_e2d int_edge;

			auto const mutual_arr_type = geom_collide(A_edge_seg, B_edge_seg, int_edge);


			bool intersects = false;
			geom_p2d p;
			if(mutual_arr_type == GEOM_C2D_INTERSECT)
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
				bool const n_0 = A_res.is_internal(B_edge_seg.s);
				bool const n_1 = A_res.is_internal(B_edge_seg.t);
				assert(n_0 ^ n_1);

				A_iter->insert_next(p);
				B_iter->insert_next(p);

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
		} while (B_iter != B_res.head);


		A_iter = A_iter->next;
	} while (A_iter != A_res.head);
}
#endif



geom_polygon2d geom_union(geom_polygon2d const&A, geom_polygon2d const &B)
{
/*	auto A_iter = A.head;
	do
	{
		if (A_iter->other_out != nullptr)
			A_iter->next = A_iter->other_out;

		A_iter = A_iter->next;
	} while (A_iter != A.head);
*/
	geom_polygon2d result;
	auto A_iter = A.head;
	do
	{
		result.insert_back(A_iter->point);
		if(A_iter->other_out == nullptr)
			A_iter = A_iter->next;
		else
			A_iter = A_iter->other_out;

	} while (A_iter != A.head);

	return result;
}

void geom_polygon_orientation_change(geom_polygon2d& A)
{
	geom_vertex2d* curr = A.head;
	geom_vertex2d* next = nullptr;
	geom_vertex2d* prev = nullptr;

	while(curr->next != A.head)
	{
		//Saving next
		next = curr->next;
		//Reverting
		curr->next = prev;

		//Iteration
		prev = curr;
		curr = next;
	}

	curr->next = prev;
	A.head->next = curr;
	A.head = curr;
}

geom_polygon2d geom_minus(geom_polygon2d const& A, geom_polygon2d& B)
{
	geom_polygon2d result;
	auto B_iter = B.head;
	bool B_inside_A = true;
	do
	{
		B_inside_A &= A.is_internal(B_iter->point);
		B_iter = B_iter->next;
	} while(B_iter != B.head);

	geom_polygon_orientation_change(B);
	if(!B_inside_A)
	{
		result = geom_union(A, B);
	}
	else
	{
		///@todo implement
	}
	return result;
}