cmake_minimum_required(VERSION 3.10)
project(GoddamnOpuwenijSolver)

include_directories(.)

if(MSVC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /std:c++latest")
endif(MSVC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
add_executable(GoddamnOpuwenijSolver main.cc
        libCommon/src/CommTest.hh

        libGeometry2D/src/GeomBase.hh
        libGeometry2D/src/GeomPoint.hh libGeometry2D/src/GeomPoint.cc
        libGeometry2D/src/GeomBox.hh libGeometry2D/src/GeomBox.cc
        libGeometry2D/src/GeomEdge.hh libGeometry2D/src/GeomEdge.cc
        libGeometry2D/src/GeomPoly.cc libGeometry2D/src/GeomPoly.hh
        libGeometry2D/src/GeomTriangle.cc libGeometry2D/src/GeomTriangle.hh libGeometry2D/src/GeomTetrahedron.cc libGeometry2D/src/GeomTetrahedron.hh libGeometry2D/src/GeomDelaunay.cc libGeometry2D/src/GeomDelaunay.hh libGeometry2D/src/GeomTriangulation.cc libGeometry2D/src/GeomTriangulation.hh)
