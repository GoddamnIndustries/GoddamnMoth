#include <iostream>

#include "Geometry/Geometry.h"

int main()
{
	geom_line_2d first(geom_point_2d(-2,-2), geom_point_2d(1, 1));
	geom_line_2d second(geom_point_2d(0, 0), geom_point_2d(1, 1));

	auto const int_p = geom_intersects(first, second);

	std::cout << int_p.X() << " " << int_p.Y() << std::endl;

	return 0;
}