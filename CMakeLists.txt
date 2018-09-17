cmake_minimum_required(VERSION 3.12)
project(GoddamnOpuwenijSolver)

set(CMAKE_CXX_STANDARD 14)

IF(APPLE)

    set(CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE FILEPATH /usr/local/bin/clang-omp CACHE FILEPATH "C compiler" FORCE)
    set(CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE FILEPATH "C++ compiler" FORCE)

    FIND_PACKAGE(OpenMP REQUIRED)
    if(OPENMP_FOUND)
        message("OPENMP FOUND")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        include_directories(${OpenMP_INCLUDE_DIRS})
    endif()

    FIND_LIBRARY(FADE2D_LIBRARY fade2d ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/lib)
    IF(FADE2D_LIBRARY)
        message("kek")
    ENDIF(FADE2D_LIBRARY)

    set(FADE2D_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/fade2d/include)

    include_directories(${FADE2D_INCLUDE_DIRS})

ELSEIF(UNIX)

    set(CMAKE_CXX_COMPILER "clang")
    set(CMAKE_CXX_COMPILER "clang++")

    FIND_PACKAGE(OpenMP REQUIRED)
    if(OPENMP_FOUND)
        message("OPENMP FOUND")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        include_directories(${OpenMP_INCLUDE_DIRS})
    endif()

    FIND_LIBRARY(FADE2D_LIBRARY fade2d ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/lib/fade2d_ubuntu)
    IF(FADE2D_LIBRARY)
        message("kek")
    ENDIF(FADE2D_LIBRARY)

    set(FADE2D_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/fade2d/include)
    include_directories(${FADE2D_INCLUDE_DIRS})


    FIND_PACKAGE(MKL REQUIRED)


ELSE()


    FIND_LIBRARY(FADE2D_LIBRARY fade2D_x64_v141_Release ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/lib/fade2d_win64)
    IF(FADE2D_LIBRARY)
        message("kek")
    ENDIF()


    #FIND_LIBRARY(GMP_LIBRARY libgmp-10 ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/lib/fade2d_win64)


    set(FADE2D_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lib_source/fade2d/include)
    include_directories(${FADE2D_INCLUDE_DIRS})



    FIND_PACKAGE(OpenMP REQUIRED)
    if(OPENMP_FOUND)
        message("OPENMP FOUND")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        include_directories(${OpenMP_INCLUDE_DIRS})
    endif()



ENDIF()


add_executable(GoddamnOpuwenijSolver main.cpp libGeometry/Geometry.h libGeometry/Geometry.cpp)

IF(APPLE)
    target_link_libraries(GoddamnOpuwenijSolver ${FADE2D_LIBRARY})
ELSEIF(UNIX)
    target_link_libraries(GoddamnOpuwenijSolver ${FADE2D_LIBRARY} ${MKL_LIBRARIES})
ELSE()
    #target_link_libraries(Triangular_Solvers ${FADE2D_LIBRARY} ${LAPACK_LIBRARY} ${GMP_LIBRARY} ${BLASS_LIBRARY})
    target_link_libraries(GoddamnOpuwenijSolver ${FADE2D_LIBRARY})
ENDIF()
