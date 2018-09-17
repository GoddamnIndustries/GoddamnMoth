#include <iostream>

#include "libGeometry/Geometry.h"
//#include "libLinearAlgebra/Matrix.h"

/*
#define X_LEFT (-2.0)
#define X_RIGHT 8.0
#define Y_BOT (-5.0)
#define Y_TOP 5.0 */

#define X_LEFT (0.0)
#define X_RIGHT (3.0)
#define Y_TOP (2.0)
#define Y_BOT (0.0)

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


	geom_constraint_function circle1([](double t)
									  {
										  auto x = 0.2 * cos(2.0 * M_PI * t);
										  auto y = 0.2 * sin(2.0 * M_PI * t);
										  return GEOM_FADE2D::Point2(x, y);
									  }, 40); //0,2 50

	geom_constraint_function square([](double t)
									 {
										 double x, y;
										 if((t >= 0) && (t < 0.25))
										 {
											 x = X_LEFT;
											 //y = -4 * t + 1;
											 y = Y_TOP + (Y_BOT - Y_TOP) * 4 * t;
											 // y = -Y_TOP * (t - 0.25) / 0.25 + Y_BOT * t / 0.25;
										 }
										 else if((t >= 0.25) && (t < 0.5))
										 {
											 y = Y_BOT;
											 //x = 4*(t - 0.25) - 1;
											 x = X_LEFT + (X_RIGHT - X_LEFT) * 4 * (t - 0.25);
										 }
										 else if((t >= 0.5) && (t < 0.75))
										 {
											 x = X_RIGHT;
											 // y = 4*(t - 0.5) - 1;
											 y = Y_BOT + (Y_TOP - Y_BOT) * 4 * (t - 0.5);
										 }
										 else
										 {
											 y = Y_TOP;
											 //x = -4*(t - 0.75) + 1;
											 x = X_RIGHT + (X_LEFT - X_RIGHT) * 4 * (t - 0.75);
										 }

										 return GEOM_FADE2D::Point2(x, y);


									 }, 4);

	geom_zone circle_zone(circle1, false);
	geom_zone square_zone(square, true);

	std::vector<geom_zone> zones_v;
//	zones_v.push_back(circle_zone);
	zones_v.push_back(square_zone);

	geom_triangle_params params;

	geom_area area(zones_v);


	auto const convex_mesh = area.triangulate(params);



	return 0;
}