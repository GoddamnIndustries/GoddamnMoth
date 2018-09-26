cmake_minimum_required(VERSION 3.10)
project(GoddamnOpuwenijSolver)

include_directories(.)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
add_executable(GoddamnOpuwenijSolver main.cc 
        libGeometry2D/src/GeomPolygoneList.hh
        libGeometry2D/src/GeomBase.hh
        libGeometry2D/src/GeomEdge.hh libGeometry2D/src/GeomEdge.cc
        libGeometry2D/src/GeomPoint.hh libGeometry2D/src/GeomPoint.cc
        libGeometry2D/src/GeomEdgeList.cc libGeometry2D/src/GeomEdgeList.hh libCommon/src/CommTest.hh)
