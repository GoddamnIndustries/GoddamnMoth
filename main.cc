
#include "libPosix/src/PosixFork.hh"
#include "libGeometry2D/src/GeomPolygoneList.hh"

int main()
{
	//system("gnuplot -e \"plot 'H:\\GoddamnOpuwenijSolver\\quad1.txt' with lines; pause -1;\"");

	/*geom_e2d e1{ {0.0, 0.0}, {1.0, 0.0} };
	geom_e2d e2{ {2.5, -.5}, {2.5, 0.5} };
	geom_intersects(e1, e2);*/

	geom_polygon2d p({ 0.0, 0.0 });
	p.insert({ 0.0, 4.0 });
	p.insert({ 4.0, 4.0 }); 
	p.insert({ 4.0, 0.0 });

	geom_polygon2d b({ 3.0, 1.0 });
	b.insert({ 3.0, 3.0 });
	b.insert({ 5.0, 3.0 });
	b.insert({ 5.0, 1.0 });


	geom_e2d e1({0.0, 0.0}, {1.0, 1.0});
	geom_e2d e2({0.5, 0.5}, {0.7, 0.7});

	geom_e2d e3({0.0, 1.0}, {1.0, 0.0});

	geom_e2d e4({0.0, -1.0}, {1.0, 0.0});


	geom_e2d int_seg_1, int_seg_2, int_seg_3;
	auto const m_a_1 = geom_mutual_arrangement(e1, e2, int_seg_1); // geom_edges_intersect_on_segment
	auto const m_a_2 = geom_mutual_arrangement(e1, e3, int_seg_2); //geom_edges_intersect_in_point
	auto const m_a_3 = geom_mutual_arrangement(e1, e4, int_seg_3); // do not intesect

	assert(m_a_1 == geom_edges_intersect_on_segment);
	assert(m_a_2 == geom_edges_intersect_in_point);
	assert(m_a_3 == geom_edges_do_not_intersect);

	assert((int_seg_1.p.x == 0.5) && (int_seg_1.p.y == 0.5) && (int_seg_1.q.x == 0.7) && (int_seg_1.q.y == 0.7));
	assert((int_seg_2.p.x == 0.5) && (int_seg_1.p.y == 0.5) && (int_seg_2.q.x == 0.5) && (int_seg_2.q.y == 0.5));


//	geom_clip(&p, &b);
//	geom_minus(&p, &b);
//	p.plot();

	p.print("A.txt");

	//int a = fork();
	return 0;
}
