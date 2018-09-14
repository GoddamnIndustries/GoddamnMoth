#include <iostream>

#include "libGeometry/Geometry.h"
#include "libLinearAlgebra/Matrix.h"

int main()
{
	geom_line_2d first(geom_point_2d(0,0), geom_point_2d(1, 1));
	geom_line_2d second(geom_point_2d(1, 1), geom_point_2d(2, 0));
	geom_line_2d third(geom_point_2d(2, 0), geom_point_2d(0,0));

	std::vector<geom_line_2d> line_v;
	line_v.emplace_back(first);
	line_v.emplace_back(second);
	line_v.emplace_back(third);

	geom_convex_shape shape(line_v);

	geom_point_2d point(1.0, 0.5);

	std::cout << shape.is_internal(point) << std::endl;



	return 0;
}