cmake_minimum_required(VERSION 3.10)
project(GoddamnOpuwenijSolver)

set(CMAKE_CXX_STANDARD 14)

add_executable(GoddamnOpuwenijSolver main.cpp 
				libGeometry/Geometry.h 
				libGeometry/Geometry.cpp)