#pragma once

#include "libGeometry2D/src/GeomTriangle.hh"

#include <vector>

struct MOTH_CORE moth_tri_cell2d
{
public:
    moth_size_t p1{MOTH_NPOS};
    moth_size_t p2{MOTH_NPOS};
    moth_size_t p3{MOTH_NPOS};
};  // struct moth_tri_cell2d

/**
 * Abstract Delaunay triangulation builder.
 */
struct MOTH_CORE moth_tri_grid2d
{
public:
    /**
     * Adds given set of points onto the set.
     */
    MOTH_HOST
    virtual void add_points(std::vector<moth_tri_cell2d>& cells,
                            std::vector<moth_p2d> p_i) = 0;
};  // struct moth_tri_grid2d
