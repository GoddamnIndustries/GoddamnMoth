
#include "libPosix/src/PosixFork.hh"
#include "libGeometry2D/src/GeomPolygoneList.hh"

int main()
{
	//system("gnuplot -e \"plot 'H:\\GoddamnOpuwenijSolver\\quad1.txt' with lines; pause -1;\"");

	/*geom_e2d e1{ {0.0, 0.0}, {1.0, 0.0} };
	geom_e2d e2{ {2.5, -.5}, {2.5, 0.5} };
	geom_intersects(e1, e2);*/

	auto* poly = new geom_e2d_list{{0.0, 0.0}};
    poly->insert_next({0.0, 1.0});
    poly->insert_next({1.0, 1.0});
    poly->insert_next({1.0, 0.0});
    auto poly_p = geom_e2d_list::len(poly);
    auto poly_s = geom_e2d_list::sqr(poly);
    auto poly_c = geom_e2d_list::str(poly);
    assert(poly_c == "((0, 0), (1, 0), (1, 1), (0, 1))");
    assert(poly_p == 4.0 && poly_s == -1.0);

    auto poly_cpy = geom_e2d_list::cpy(poly);
    poly_p = geom_e2d_list::len(poly_cpy);
    poly_s = geom_e2d_list::sqr(poly_cpy);
    poly_c = geom_e2d_list::str(poly_cpy);
    assert(poly_c == "((0, 0), (1, 0), (1, 1), (0, 1))");
    assert(poly_p == 4.0 && poly_s == -1.0);
    delete poly_cpy;

    auto poly_rev = geom_e2d_list::rev(poly);
    poly_p = geom_e2d_list::len(poly_rev);
    poly_s = geom_e2d_list::sqr(poly_rev);
    poly_c = geom_e2d_list::str(poly_rev);
    assert(poly_c == "((0, 0), (0, 1), (1, 1), (1, 0))");
    assert(poly_p == 4.0 && poly_s == +1.0);
    delete poly_rev;

	geom_e2d e1{{0.0, 0.0}, {1.0, 1.0}};
	geom_e2d e2{{0.5, 0.5}, {0.7, 0.7}};

	geom_e2d e3{{0.0, 1.0}, {1.0, 0.0}};

	geom_e2d e4{{0.0, -1.0}, {1.0, 0.0}};

	geom_e2d e5{{0.7, 0.7}, {0.5, 0.5}};
	geom_e2d e6{{0.5, 0.5}, {5.0, 5.0}};
	geom_e2d e7{{-5.0, -5.0}, {0.5, 0.5}};
	geom_e2d e8{{0.5, 0.5}, {-5.0, -5.0}};


	geom_e2d e9{{4, 4}, {4, 3}};
	geom_e2d e10{{3,1}, {8, 1}};

	geom_e2d int_seg_1, int_seg_2, int_seg_3, int_seg_4, int_seg_5, int_seg_6, int_seg_7, int_seg_8;
	auto const m_a_1 = geom_collide(e1, e2, int_seg_1); // GEOM_C2D_TOUCH
	auto const m_a_2 = geom_collide(e1, e3, int_seg_2); //GEOM_C2D_INTERSECT
	auto const m_a_3 = geom_collide(e1, e4, int_seg_3); // do not intesect

	auto const m_a_4 = geom_collide(e1, e5, int_seg_4);
	auto const m_a_5 = geom_collide(e1, e6, int_seg_5);
	auto const m_a_6 = geom_collide(e1, e7, int_seg_6);
	auto const m_a_7 = geom_collide(e1, e8, int_seg_7);

	auto const m_a_8 = geom_collide(e10, e9, int_seg_8);



	assert(m_a_1 == GEOM_C2D_TOUCH);
	assert(m_a_2 == GEOM_C2D_INTERSECT);
	assert(m_a_3 == GEOM_C2D_NONE);

	assert(m_a_4 == GEOM_C2D_TOUCH);
	assert(m_a_5 == GEOM_C2D_TOUCH);
	assert(m_a_6 == GEOM_C2D_TOUCH);
	assert(m_a_7 == GEOM_C2D_TOUCH);


	assert(m_a_8 == GEOM_C2D_NONE);


	assert((int_seg_1.s.x == 0.5) && (int_seg_1.s.y == 0.5) && (int_seg_1.t.x == 0.7) && (int_seg_1.t.y == 0.7));
	assert((int_seg_2.s.x == 0.5) && (int_seg_1.s.y == 0.5) && (int_seg_2.t.x == 0.5) && (int_seg_2.t.y == 0.5));


	geom_polygon2d p({ 0.0, 0.0 });
	p.insert_back({ 4.0, 0.0 });
	p.insert_back({ 4.0, 4.0 });
	p.insert_back({ 0.0, 4.0 });

	geom_p2d p1 = {3.0, 1.0};
	geom_p2d p2 = {5.0, 3.0};

	assert(p.is_internal(p1));
	assert(!p.is_internal(p2));



	geom_polygon2d b({ 3.0, 1.0 });
	b.insert_back({ 5.0, 1.0 });
	b.insert_back({ 5.0, 3.0 });
	b.insert_back({ 3.0, 3.0 });

	geom_polygon2d b_buf, p_buf;

	geom_clip(p, b, p_buf, b_buf);
	auto res = geom_minus(p_buf, b_buf);
	res.print("A.txt");

/*	geom_vertex2d c({2.0, 1.0});
	b.insert_back({ 2.0, 2.0 });
	b.insert_back({ 3.0, 2.0 });
	b.insert_back({ 3.0, 1.0 });

	geom_clip(&p, &c);
	geom_minus(&p, &c);
	p.print("A.txt");
*/
//	p.plot();


	//int a = fork();
	return 0;
}
