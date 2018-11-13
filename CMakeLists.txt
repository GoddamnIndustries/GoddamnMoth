cmake_minimum_required(VERSION 3.10)
project(GoddamnOpuwenijSolver)

include_directories(.)


set(CMAKE_CXX_STANDARD 17)

if(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /openmp")
else(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -ffast-math -mavx -g")
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -ffast-math -mavx -DNDEBUG=1")
	if(APPLE)
	else(APPLE)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
	endif()
endif(MSVC)

add_executable(GoddamnOpuwenijSolver main.cc
        libCommon/src/CommTest.hh

        libGeometry2D/src/GeomBase.hh
        libGeometry2D/src/GeomPoint.hh libGeometry2D/src/GeomPoint.cc
		libGeometry2D/src/GeomEdge.hh libGeometry2D/src/GeomEdge.cc
        libGeometry2D/src/GeomPoly.cc libGeometry2D/src/GeomPoly.hh
        libGeometry2D/src/GeomTriangle.cc libGeometry2D/src/GeomTriangle.hh libGeometry2D/src/GeomTetrahedron.cc libGeometry2D/src/GeomTetrahedron.hh libGeometry2D/src/GeomMesh.cc libGeometry2D/src/GeomMesh.hh libGeometry2D/src/GeomMesh2.hh libGeometry2D/src/GeomMesh2.cc libGeometry2D/src/GeomSort.cc libGeometry2D/src/GeomSort.hh)
