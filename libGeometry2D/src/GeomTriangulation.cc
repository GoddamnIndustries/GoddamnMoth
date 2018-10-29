#include "libGeometry2D/src/GeomTriangulation.hh"

#include <unordered_set>
#include <unordered_map>

#include <vector>
#include <tuple>

#include <iostream>
#include <fstream>

static const moth_tri2d Sx{{-100, -100}, {+100, -100}, {0, 100}};

MOTH_HOST
moth_tri_grid2d_builder::moth_tri_grid2d_builder()
{
    /* Create and add points of super cell. */
    moth_cell2d_tri S{~1u, ~2u, ~3u};
    tri.push_back(S);
}

MOTH_HOST
moth_tri2d moth_tri_grid2d_builder::triangle(const moth_cell2d_tri& tri_cell,  bool* boundary) const
{
    if (boundary != nullptr) {
        *boundary = false;
    }

    moth_tri2d T{};
    switch (tri_cell.k1) {
        case ~1u:
            T.p1 = Sx.p1;
            if (boundary != nullptr) {
                *boundary = true;
            }
            break;
        case ~2u:
            T.p1 = Sx.p2;
            if (boundary != nullptr) {
                *boundary = true;
            }
            break;
        case ~3u:
            T.p1 = Sx.p3;
            if (boundary != nullptr) {
                *boundary = true;
            }
            break;
        default:
            T.p1 = tri_points[tri_cell.k1];
            break;
    }
    switch (tri_cell.k2) {
        case ~1u:
            T.p2 = Sx.p1;
            if (boundary != nullptr) {
                *boundary = true;
            }
            break;
        case ~2u:
            T.p2 = Sx.p2;
            if (boundary != nullptr) {
                *boundary = true;
            }
            break;
        case ~3u:
            T.p2 = Sx.p3;
            if (boundary != nullptr) {
                *boundary = true;
            }
            break;
        default:
            T.p2 = tri_points[tri_cell.k2];
            break;
    }
    switch (tri_cell.k3) {
        case ~1u:
            T.p3 = Sx.p1;
            if (boundary != nullptr) {
                *boundary = true;
            }
            break;
        case ~2u:
            T.p3 = Sx.p2;
            if (boundary != nullptr) {
                *boundary = true;
            }
            break;
        case ~3u:
            T.p3 = Sx.p3;
            if (boundary != nullptr) {
                *boundary = true;
            }
            break;
        default:
            T.p3 = tri_points[tri_cell.k3];
            break;
    }
    return T;
}
MOTH_HOST
bool moth_tri_grid2d_builder::_triangle(moth_tri2d& tr, const moth_cell2d_tri& cell) const
{
    bool boundary;
    tr = triangle(cell, &boundary);
    return !boundary;
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST
void moth_tri_grid2d_builder::insert(const moth_p2d& p)
{
    /* Insert new point to the end. */
    moth_size_t i = tri_points.size();
    tri_points.push_back(p);

    /* Find bad triangles. */
    std::vector<moth_cell2d_edge> tr_bad_poly;
    for (moth_cell2d_tri& tri_cell : tri) {
        tri_cell.bad = moth_tri2d::circle(triangle(tri_cell), p);
        if (tri_cell.bad) {
            tr_bad_poly.push_back({tri_cell.k1, tri_cell.k2});
            tr_bad_poly.push_back({tri_cell.k2, tri_cell.k3});
            tr_bad_poly.push_back({tri_cell.k3, tri_cell.k1});
        }
    }
    /* Remove bad triangles. */
    tri.erase(std::remove_if(
              tri.begin(), tri.end(), [](const auto& tri_cell){ return tri_cell.bad; }),
              tri.end());

    /* Find shared edges of the polygon triangulation. */
    for (moth_size_t k = 0; k < tr_bad_poly.size(); ++k) {
        for (moth_size_t m = k + 1; m < tr_bad_poly.size(); ++m) {
            moth_cell2d_edge& e1 = tr_bad_poly[k];
            moth_cell2d_edge& e2 = tr_bad_poly[m];
            if (e1 == e2) {
                e1.bad = true;
                e2.bad = true;
            }
        }
    }
    /* Remove shared edges. */
    tr_bad_poly.erase(std::remove_if(
                      tr_bad_poly.begin(), tr_bad_poly.end(), [](const moth_cell2d_edge& e) { return e.bad; }),
                      tr_bad_poly.end());

    /* Add new triangles. */
    for (const auto& e : tr_bad_poly) {
        moth_cell2d_tri C{e.k1, e.k2, i};
        moth_tri2d T{triangle(C)};

        /* Make sure that the triangle is CCW. */
        if (moth_tri2d::area(T) < 0.0) {
            std::swap(C.k1, C.k2);
        }
        tri.push_back(C);
    }
}

/**
 * @see https://en.wikipedia.org/wiki/Ruppert%27s_algorithm
 * @see https://people.eecs.berkeley.edu/~jrs/papers/2dj.pdf
 */
MOTH_HOST
void moth_tri_grid2d_builder::refine()
{
    std::unordered_set<moth_size_t> tri_bad_edges;
    std::unordered_set<moth_size_t> tri_bad_triangles;
    for (moth_size_t n = 0; n < 20; ++n) {
        std::cerr << n << " " << tri.size() << std::endl;

        /* Phase 1: Fix all encroached edges. */
        while (true) {
            for (moth_size_t m = 0; m < tri_edges.size(); ++m) {
                const moth_e2d& e = tri_edges[m];

                /* Find if this edge is encroached by any point. */
                moth_p2d c = 0.5 * (e.p1 + e.p2);
                moth_real_t r = 0.5 * moth_p2d::len(e.vec());
                for (const moth_p2d& p : tri_points) {
                    if (p != e.p1 && p != e.p2 &&
                        moth_p2d::len(p - c) < r) {
                        tri_bad_edges.insert(m);
                        break;
                    }
                }
            }
            if (!tri_bad_edges.empty()) {
                tri_edges.reserve(tri_edges.capacity() + tri_bad_edges.size());
                for (moth_size_t m : tri_bad_edges) {
                    moth_e2d& e = tri_edges[m];

                    /* Calculate center point and add it to triangulation. */
                    moth_p2d c = 0.5 * (e.p1 + e.p2);
                    insert(c);

                    moth_p2d p = e.p2;
                    e.p2 = c;
                    tri_edges.push_back({c, p});
                }
                tri_bad_edges.clear();
            } else {
                break;
            }
        }

        /* Phase 2: Fix all bad triangles. */
        for (moth_size_t k = 0; k < tri.size(); ++k) {
            const moth_cell2d_tri& C = tri[k];
            moth_tri2d T{};
            if (_triangle(T, C)) {
                /* Assume all triangles are poor. */
                moth_e2d e12 = T.edge(12);
                moth_e2d e23 = T.edge(23);
                moth_e2d e31 = T.edge(31);
                moth_real_t theta1{std::fabs(moth_e2d::angle(e31, e12)) * 180.0 / MOTH_PI};
                moth_real_t theta2{std::fabs(moth_e2d::angle(e12, e23)) * 180.0 / MOTH_PI};
                moth_real_t theta3{std::fabs(moth_e2d::angle(e23, e31)) * 180.0 / MOTH_PI};
                moth_real_t theta{std::min(theta1, std::min(theta2, theta3)) / 2};
                if (theta < 35.0) {
                    tri_bad_triangles.insert(k);
                }
            }
        }
        if (!tri_bad_triangles.empty()) {
            for (moth_size_t k : tri_bad_triangles) {
                const moth_cell2d_tri& C = tri[k];
                moth_tri2d T{};
                if (_triangle(T, C)) {
                    bool cc_bad = false;
                    moth_p2d cc{moth_tri2d::circumcenter(T)};
                    for (moth_size_t m = 0; m < tri_edges.size(); ++m) {
                        const moth_e2d& e = tri_edges[m];

                        /* Find if this edge is encroached by the circumcenter. */
                        moth_p2d c = 0.5 * (e.p1 + e.p2);
                        moth_real_t r = 0.5 * moth_p2d::len(e.vec());
                        if (moth_p2d::len(cc - c) < r) {
                            tri_bad_edges.insert(m);
                            cc_bad = true;
                            break;
                        }
                    }
                    if (!cc_bad) {
                        insert(cc);
                    }
                }
            }
            tri_bad_triangles.clear();
        }
    }
}

MOTH_HOST
void moth_tri_grid2d_builder::print()
{
    std::string TrFilePath("res/tr-" + std::to_string(99999) + ".txt");
    std::ofstream TrFile(TrFilePath);

    for (moth_size_t k = 0; k < tri.size(); ++k) {
        bool super_boundary;
        const moth_cell2d_tri& C = tri[k];
        moth_tri2d T{triangle(C, &super_boundary)};
        if (!super_boundary) {
            TrFile << T.p1 << std::endl
                   << T.p2 << std::endl
                   << T.p3 << std::endl
                   << T.p1 << std::endl
                   << std::endl;
        }
    }
}
