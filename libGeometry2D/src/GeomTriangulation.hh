#pragma once

#include "libGeometry2D/src/GeomTriangle.hh"

#include <vector>

struct MOTH_CORE moth_cell2d_tri
{
public:
    moth_size_t k1{MOTH_NPOS};
    moth_size_t k2{MOTH_NPOS};
    moth_size_t k3{MOTH_NPOS};
    bool bad = false;
};  // struct moth_tri_cell2d

struct MOTH_CORE moth_cell2d_edge
{
    moth_size_t k1{MOTH_NPOS};
    moth_size_t k2{MOTH_NPOS};
    bool bad = false;

public:
    bool operator==(const moth_cell2d_edge& e) const
    {
        return (k1 == e.k1 && k2 == e.k2) || (k2 == e.k1 && k1 == e.k2);
    }
};  // struct moth_cell2d_edge

/**
 * Delaunay triangulation builder.
 */
struct MOTH_CORE moth_tri_grid2d_builder
{
private:
    std::vector<moth_cell2d_tri> tri;
    std::vector<moth_p2d> tri_points;
public:
    std::vector<moth_e2d> tri_edges;

    MOTH_HOST
    bool _triangle(moth_triangle2d& tr, const moth_cell2d_tri& cell) const;
    MOTH_HOST
    moth_triangle2d triangle(const moth_cell2d_tri& C, bool* boundary = nullptr) const;

public:
    MOTH_HOST
    moth_tri_grid2d_builder();

public:
    MOTH_HOST
    void insert(const moth_p2d& p);

public:
    MOTH_HOST
    void refine();

public:
    MOTH_HOST
    void print();
};  // struct moth_tri_grid2d_builder
