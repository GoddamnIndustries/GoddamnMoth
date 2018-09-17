#include "Geometry.h"

#include <array>


geom_point_2d geom_intersects(geom_line_2d const& c, geom_line_2d const& d)
{
	std::array<double, 2> c_vec = {c.P1().X() - c.P0().X(), c.P1().Y() - c.P0().Y()};
	std::array<double, 2> d_vec = {d.P1().X() - d.P0().X(), d.P1().Y() - d.P0().Y()};

	auto const det = -c_vec[0] * d_vec[1] + d_vec[0] * c_vec[1];

/*
	if(det == 0) // c == d or has now intersection points
		throw geom_error();
*/


	auto const b1 = d.P0().X() - c.P0().X();
	auto const b2 = d.P0().Y() - c.P0().Y();

	double const alpha = 1 / det * (-b1 * d_vec[1] + b2 * d_vec[0]);

	return geom_point_2d(c.P0().X() + alpha * c_vec[0], c.P0().Y() + alpha * c_vec[1]);

}

geom_convex_shape geom_triangle2_to_shape(GEOM_FADE2D::Triangle2* triangle)
{
	geom_point_2d const first_vertex(*(triangle->getCorner(0)));
	geom_point_2d const second_vertex(*(triangle->getCorner(1)));
	geom_point_2d const third_vertex(*(triangle->getCorner(2)));


	std::vector<geom_line_2d> lines_v;
	lines_v.emplace_back(geom_line_2d(first_vertex, second_vertex));
	lines_v.emplace_back(geom_line_2d(second_vertex, third_vertex));
	lines_v.emplace_back(geom_line_2d(third_vertex, first_vertex));

	geom_convex_shape shape(lines_v);

	return shape;

}

geom_convex_shape_vector geom_area::triangulate(geom_triangle_params const& triangleProperties)
{

	GEOM_FADE2D::Fade_2D global_area;

	//making bounding rectangle
	GEOM_FADE2D::Point2 p1(-1, -1), p2(-1, 1), p3(1, -1), p4(1, 1);
	global_area.insert(p1);
	global_area.insert(p2);
	global_area.insert(p3);
	global_area.insert(p4);


	std::vector<std::vector<GEOM_FADE2D::Segment2>> vvSegment(m_zones.size());
	std::vector<GEOM_FADE2D::ConstraintGraph2 *> vCSG;
	std::vector<GEOM_FADE2D::Zone2 *> vpZonesDealunay;

	for (int i = 0; i < m_zones.size(); ++i)
	{
		//discretizing zone's constraints
		auto vPoints = m_zones[i].discretize();

		//creating segments
		for (int j = 0; j < vPoints.size(); ++j)
		{
			vvSegment[i].push_back(GEOM_FADE2D::Segment2(vPoints[j], vPoints[(j + 1) % vPoints.size()]));
		}

		//creating constraint graphs
		vCSG.push_back(global_area.createConstraint(vvSegment[i], GEOM_FADE2D::CIS_CONFORMING_DELAUNAY));


		//creating fade2D delaunay zones
		m_zones[i].is_inside() ? vpZonesDealunay.push_back(global_area.createZone(vCSG[i], GEOM_FADE2D::ZL_INSIDE)) :
		vpZonesDealunay.push_back(global_area.createZone(vCSG[i], GEOM_FADE2D::ZL_OUTSIDE));


	}

	//auto pGrowZone = global_area.createZone(vCSG, GEOM_FADE2D::ZL_GROW, p1);

	GEOM_FADE2D::Zone2* pGrowZone = global_area.createZone(nullptr, GEOM_FADE2D::ZL_GLOBAL);

	//calculating final zone
	auto const size = static_cast<int>(vpZonesDealunay.size());
	for (int i = 0; i < size; ++i)
	{
		//pGrowZone = GEOM_FADE2D::zoneUnion(pGrowZone,
		//GEOM_FADE2D::zoneIntersection(vpZonesDealunay[i], vpZonesDealunay[i + 1]));
		pGrowZone = GEOM_FADE2D::zoneIntersection(pGrowZone, vpZonesDealunay[i]);
	}


	global_area.applyConstraintsAndZones();



//	pGrowZone->show("kekas.ps", false, true);


	//refining final zone
	auto pBoundedZone(pGrowZone->convertToBoundedZone());

	GEOM_FADE2D::MeshGenParams params(pBoundedZone); //changed
	params.minAngleDegree = triangleProperties.min_angle_degree;
	params.minEdgeLength = triangleProperties.min_edge_length;
	params.maxEdgeLength = triangleProperties.max_edge_length;
	params.gridVector = triangleProperties.grid_vector;
	params.gridLength = triangleProperties.grid_length;
	params.growFactor = triangleProperties.grow_factor;
	params.growFactorMinArea = triangleProperties.grow_factor_min_area;
	params.capAspectLimit = triangleProperties.cap_aspect_limit;

//	global_area.refine(pBoundedZone, triangleProperties[0], triangleProperties[1], triangleProperties[2], true);
	global_area.refineAdvanced(&params);

	pBoundedZone->show("lul.ps", false, true);


//	global_area.show("kek.ps");


	std::vector<GEOM_FADE2D::Triangle2*> vTriangles2;

	pBoundedZone->getTriangles(vTriangles2);


	geom_convex_shape_vector shape_v;
	make_convex_shape_vector(shape_v, vTriangles2);


	return shape_v;

}





