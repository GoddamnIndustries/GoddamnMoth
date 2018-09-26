cmake_minimum_required(VERSION 3.10)
project(GoddamnOpuwenijSolver)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
add_executable(GoddamnOpuwenijSolver main.cc 
			libGeometry2D/src/GeomPolygoneList.hh libGeometry2D/src/GeomEdge.hh libGeometry2D/src/GeomEdge.cc)
