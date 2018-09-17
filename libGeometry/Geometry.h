#ifndef GODDAMNOPUWENIJSOLVER_GEOMETRY_H
#define GODDAMNOPUWENIJSOLVER_GEOMETRY_H

#include <Fade_2D.h>
#include <vector>
#include <functional>

using geom_real_t = double;

// 2D point in space.
class geom_point_2d
{
	geom_real_t x, y;

public:
	geom_point_2d() = default;

	geom_point_2d(geom_real_t X, geom_real_t Y): x(X), y(Y)
	{}

	geom_point_2d(GEOM_FADE2D::Point2 const& point): x(point.x()), y(point.y())
	{}

	geom_real_t X() const
	{
		return x;
	}
	geom_real_t Y() const
	{
		return y;
	}
};

// Line in 2D space (segment)
class geom_line_2d
{
	geom_point_2d point0, point1;

public:
	geom_line_2d() = default;

	geom_line_2d(geom_point_2d const& first, geom_point_2d const& second): point0(first), point1(second)
	{}

	geom_point_2d P0() const
	{
		return point0;
	}

	geom_point_2d P1() const
	{
		return point1;
	}

	bool if_point_on_segment(geom_point_2d const& point) const
	{
		double const alpha = (point.X() - point0.X()) / (point1.X() - point0.X());
		if((alpha < 0) || (alpha > 1))
			return false;

		double const beta = (point.Y() - point0.Y()) / (point1.Y() - point0.Y());
		if(alpha != beta)
			return false;
		return true;
	}

};

// Returns intersection point of two lines
geom_point_2d geom_intersects(geom_line_2d const& c, geom_line_2d const& d);

// Convex shape class
class geom_convex_shape
{
	std::vector<geom_line_2d> lines;

	std::vector<geom_convex_shape*> neighbours;

public:

	geom_convex_shape(): lines(0), neighbours(0)
	{}

	explicit geom_convex_shape(std::vector<geom_line_2d> const& lines_v): lines(lines_v), neighbours(lines_v.size(), nullptr)
	{}

	bool is_internal(const geom_point_2d& point) const
	{
		geom_line_2d line_from_point(point, geom_point_2d(point.X() + 1.0, point.Y()));
		unsigned inters_counter = 0;
		for(auto const& line : lines)
		{
			auto const inters_point = geom_intersects(line_from_point, line);
			if((inters_point.X() > point.X()) && (line.if_point_on_segment(inters_point)))
				++inters_counter;
		}
		if((inters_counter % 2) == 0)
			return false;
		return true;
	}

	void set_neighbour(int neigh_num, geom_convex_shape* const p_neigh)
	{
		neighbours[neigh_num] = p_neigh;
	}

	geom_convex_shape* get_neighbour(int neigh_num) const
	{
		return neighbours[neigh_num];
	}
};

class geom_convex_shape_vector
{
	std::vector<geom_convex_shape> parts;
// and other useful stuff.
public:

	geom_convex_shape& operator[](unsigned const i)
	{
		return parts[i];
	}

	geom_convex_shape const& operator[](unsigned const i) const
	{
		return parts[i];
	}

	void resize(size_t const size)
	{
		parts.resize(size);
	}

	bool is_internal(const geom_point_2d& point) const
	{
		for(const geom_convex_shape& shape : parts)
			if(shape.is_internal(point))
				return true;
			return false;
	}
};

geom_convex_shape geom_triangle2_to_shape(GEOM_FADE2D::Triangle2* triangle);

class geom_constraint_function: public std::function<GEOM_FADE2D::Point2(double)> // function of argument t from 0 to 1
{
	unsigned m_discr_point_number;

public:

	explicit geom_constraint_function(std::function<GEOM_FADE2D::Point2(double)> const& func, unsigned discrPointNumber):
			std::function<GEOM_FADE2D::Point2(double)>(func), m_discr_point_number(discrPointNumber)
	{}


	std::vector<GEOM_FADE2D::Point2> Discretize() const
	{
		auto const h = 1.0 / m_discr_point_number;
		std::vector<GEOM_FADE2D::Point2> vPoints;

		for(unsigned i = 0; i < m_discr_point_number; ++i)
		{
			auto point = (*this)(i*h);
			vPoints.push_back(point);
		}

		return vPoints;
	}


};


//!@brief Zone class describes a zone bounded by some functions.
class geom_zone
{
	const geom_constraint_function m_function;
	const bool m_is_inside;


//		const std::function<std::array<double, 4>(GEOM_FADE2D::Point2)> m_boundaryConditions; //returns {ro, u, v, P}

public:
	geom_zone(geom_constraint_function const& func, bool isInside): m_function(func),
														 m_is_inside(isInside)
	{}


	std::vector<GEOM_FADE2D::Point2> discretize() const
	{
		return m_function.Discretize();
	}

	bool is_inside() const
	{
		return m_is_inside;
	}

};

class geom_triangle_params
{
public:
	double min_angle_degree = 20.0;
	double min_edge_length = 0.001;
	double max_edge_length = DBL_MAX;
	GEOM_FADE2D::Vector2 grid_vector = {1.0, 0.0};
	double grid_length = 0;
	double grow_factor = DBL_MAX;
	double grow_factor_min_area = 0.001;
	double cap_aspect_limit = 10.0;

};


class geom_area
{
	std::vector<geom_zone> m_zones;




public:

	geom_area(std::vector<geom_zone> const &zones) : m_zones(zones)
	{}


	/**@brief trianglulates area consisted of zones
	**/
	geom_convex_shape_vector triangulate(geom_triangle_params const &triangleProperties);

private:

	/**@brief Finds a position of the Triangle in the triangle vector
	 *
	 * @param vpTriangle2 Triangle vector
	 * @param pTr triangle pointer
	 * @return an index
	 */
	int find_triangle_index(std::vector<GEOM_FADE2D::Triangle2 *> const& vpTriangle2,
							GEOM_FADE2D::Triangle2* const pTr) const
	{

		auto const si = vpTriangle2.size();
		if(pTr == nullptr)
			return -1;

		int i = 0;
		bool found = false;

		while ((!found) && (i < vpTriangle2.size()))
		{
			if (vpTriangle2[i] == pTr)
				found = true;
			else
				++i;
		}
		if (found)
			return i;
		return -1;
	}

	void make_convex_shape_vector(geom_convex_shape_vector & destination,
			std::vector<GEOM_FADE2D::Triangle2*>& source) const
	{
		destination.resize(source.size());
		for(int trngl_cntr = 0; trngl_cntr < source.size(); ++trngl_cntr)
		{
			destination[trngl_cntr] = geom_triangle2_to_shape(source[trngl_cntr]);
		}

		for(int trngl_cntr = 0; trngl_cntr < source.size(); ++trngl_cntr)
		{
			for(int ith = 0; ith < 3; ++ith)
			{
				auto const opp_triangle = source[trngl_cntr]->getOppositeTriangle(ith);
				auto const index = find_triangle_index(source, opp_triangle);

				if(index != -1)
					destination[trngl_cntr].set_neighbour((ith + 1) % 3, &destination[index]);
				else
					destination[trngl_cntr].set_neighbour((ith + 1) % 3, nullptr);
			}
		}
	}


};


// Polygon class
class geom_polygon
{
//	std::vector<geom_line_2d> lines;
	std::vector<geom_point_2d> vertices;


// and other useful stuff.
// with Hertel-Mehlhorn algorithm
	geom_convex_shape_vector split_convex() const;
};



#endif //GODDAMNOPUWENIJSOLVER_GEOMETRY_H
