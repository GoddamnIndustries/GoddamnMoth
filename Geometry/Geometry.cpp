#include "Geometry.h"

#include <array>


geom_point_2d geom_intersects(geom_line_2d const& c, geom_line_2d const& d)
{
	std::array<double, 2> c_vec = {c.P1().X() - c.P0().X(), c.P1().Y() - c.P0().Y()};
	std::array<double, 2> d_vec = {d.P1().X() - d.P0().X(), d.P1().Y() - d.P0().Y()};

	auto const det = -c_vec[0] * d_vec[1] + d_vec[0] * c_vec[1];


	if(det == 0) // c == d or has now intersection points
		throw geom_error();

	auto const b1 = d.P0().X() - c.P0().X();
	auto const b2 = d.P0().Y() - c.P0().Y();

	double const alpha = 1 / det * (-b1 * d_vec[1] + b2 * d_vec[0]);

	return geom_point_2d(c.P0().X() + alpha * c_vec[0], c.P0().Y() + alpha * c_vec[1]);

}



