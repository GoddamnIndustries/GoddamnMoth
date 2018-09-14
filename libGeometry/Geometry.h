#ifndef GODDAMNOPUWENIJSOLVER_GEOMETRY_H
#define GODDAMNOPUWENIJSOLVER_GEOMETRY_H

#include <vector>

using geom_real_t = double;

// 2D point in space.
class geom_point_2d
{
	geom_real_t x, y;

public:
	geom_point_2d(geom_real_t X, geom_real_t Y): x(X), y(Y)
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

public:

	geom_convex_shape(std::vector<geom_line_2d> const& lines_v): lines(lines_v)
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
};

class geom_convex_shape_vector
{
	std::vector<geom_convex_shape> parts;
// and other useful stuff.
public:
	bool is_internal(const geom_point_2d& point) const
	{
		for(const geom_convex_shape& shape : parts)
			if(shape.is_internal(point))
				return true;
			return false;
	}
};


// Polygon class
class geom_polygon
{
	std::vector<geom_line_2d> lines;
// and other useful stuff.
// with Hertel-Mehlhorn algorithm
	geom_convex_shape_vector split_convex() const;
};



#endif //GODDAMNOPUWENIJSOLVER_GEOMETRY_H
