
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

	geom_clip(&p, &b);
	geom_minus(&p, &b);
	p.plot();

	//int a = fork();
	return 0;
}
