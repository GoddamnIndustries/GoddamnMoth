#include "libGeometry2D/src/GeomTriangulation.hh"

#include <tuple>
#include <iostream>

/**
 * Delaunay triangulation builder based on Bowyer–Watson algorithm.
 * @see https://en.wikipedia.org/wiki/Bowyer–Watson_algorithm
 */
struct MOTH_CORE moth_tri_grid2d_bowyer_watson : public moth_tri_grid2d
{
    std::vector<moth_tri_cell2d> s_cells;
public:
    MOTH_HOST
    virtual void add_points(std::vector<moth_tri_cell2d>& cells,
                            std::vector<moth_p2d> p_i) override;

    void add_point(std::vector<moth_tri_cell2d>& tr_cells,
                   std::vector<moth_p2d>& points,
                   const moth_p2d& p) const;
};  // struct moth_tri_grid2d


MOTH_HOST
void moth_tri_grid2d_bowyer_watson::add_point(std::vector<moth_tri_cell2d>& tr_cells,
                                              std::vector<moth_p2d>& points,
                                              const moth_p2d& p) const
{
    points.push_back(p);
    if (points.size() >= 3) {
        moth_size_t i = points.size() - 1;

        /* Create and add the super triangle. */
        moth_tri2d Ts{{-100, -100}, {+100, -100}, {0, 200}};
        points.reserve(points.size() + 3);
        points.push_back(Ts.p1);
        points.push_back(Ts.p2);
        points.push_back(Ts.p3);

        /* Create and add the super cell. */
        moth_tri_cell2d Cs{points.size() - 3, points.size() - 2, points.size() - 1};
        tr_cells.push_back(Cs);

        /* Find bad triangles and create polygon of them. */
        std::vector<bool> tr_cells_bad(tr_cells.size());
        std::vector<std::tuple<moth_size_t, moth_size_t>> tr_poly;
        for (moth_size_t k = 0; k < tr_cells.size(); ++k) {
            const moth_tri_cell2d& C = tr_cells[k];
            if (tr_cells_bad[k] = moth_tri2d::circle({points[C.p1], points[C.p2], points[C.p3]}, p)) {
                tr_poly.emplace_back(C.p1, C.p2);
                tr_poly.emplace_back(C.p2, C.p3);
                tr_poly.emplace_back(C.p3, C.p1);
            }
        }
        /* Remove bad triangles. */
        std::vector<moth_tri_cell2d> tr_cells1;
        for (moth_size_t k = 0; k < tr_cells.size(); ++k) {
            if (!tr_cells_bad[k]) {
                tr_cells1.push_back(tr_cells[k]);
            }
        }
        tr_cells = std::move(tr_cells1);

        /* Find shared edges of the polygon triangulation. */
        std::vector<bool> tr_poly_shared(tr_poly.size());
        for (moth_size_t k = 0; k < tr_poly.size(); ++k) {
            for (moth_size_t m = k + 1; m < tr_poly.size(); ++m) {
                const std::tuple<moth_size_t, moth_size_t>& e1 = tr_poly[k];
                const std::tuple<moth_size_t, moth_size_t>& e2 = tr_poly[m];
                if (e1 == e2 ||
                    e1 == std::make_tuple(std::get<1>(e2), std::get<0>(e2))) {
                    tr_poly_shared[k] = true;
                    tr_poly_shared[m] = true;
                }
            }
        }
        /* Remove shared edges. */
        std::vector<std::tuple<moth_size_t, moth_size_t>> tr_poly1;
        for (moth_size_t k = 0; k < tr_poly.size(); ++k) {
            if (!tr_poly_shared[k]) {
                tr_poly1.push_back(tr_poly[k]);
            }
        }
        tr_poly = std::move(tr_poly1);

        /* Add new triangles. */
        for (const std::tuple<moth_size_t, moth_size_t>& e : tr_poly) {
            moth_tri_cell2d C{std::get<0>(e), std::get<1>(e), i};
            moth_tri2d T{points[C.p1], points[C.p2], points[C.p3]};
            switch (fsgn(moth_tri2d::area(T))) {
                case +1:
                    /* Nonsingular CCW triangle. */
                    tr_cells.push_back(C);
                    break;
                case -1:
                    /* Nonsingular CW triangle, reorder to CCW. */
                    std::swap(C.p1, C.p2);
                    tr_cells.push_back(C);
                    break;
                default:
                    /* Singular triangle. */
                    std::cerr << "Singular triangle" << std::endl;
                    std::abort();
            }
        }

        /* Remove the super triangle. */
        std::vector<moth_tri_cell2d> tr_cells2;
        for (const moth_tri_cell2d& C : tr_cells) {
            if (C.p1 != Cs.p1 && C.p1 != Cs.p2 && C.p1 != Cs.p3 &&
                C.p2 != Cs.p1 && C.p2 != Cs.p2 && C.p2 != Cs.p3 &&
                C.p3 != Cs.p1 && C.p3 != Cs.p2 && C.p3 != Cs.p3) {
                tr_cells2.push_back(C);
            }
        }
        tr_cells = std::move(tr_cells2);

        /* Remove the vertices of super cell. */
        points.pop_back();
        points.pop_back();
        points.pop_back();
    }
}


















MOTH_HOST
void moth_tri_grid2d_bowyer_watson::add_points(std::vector<moth_tri_cell2d>& tr,
                                               std::vector<moth_p2d> P)
{
    /* Create and add the super triangle. */
    moth_tri2d Ts{{-100, -100}, {+100, -100}, {0, 200}};
    P.push_back(Ts.p1);
    P.push_back(Ts.p2);
    P.push_back(Ts.p3);

    /* Create and add the super cell. */
    moth_tri_cell2d Cs{P.size() - 3, P.size() - 2, P.size() - 1};
    tr.push_back(Cs);

    for (moth_size_t i = 0; i < P.size() - 3; ++i) {
        const moth_p2d& p = P[i];

        /* Find bad triangles and create polygon of them. */
        std::vector<bool> tr_bad(tr.size());
        std::vector<std::tuple<moth_size_t, moth_size_t>> tr_bad_poly;
        for (moth_size_t k = 0; k < tr.size(); ++k) {
            const moth_tri_cell2d& C = tr[k];
            if (tr_bad[k] = moth_tri2d::circle({P[C.p1], P[C.p2], P[C.p3]}, p)) {
                tr_bad_poly.emplace_back(C.p1, C.p2);
                tr_bad_poly.emplace_back(C.p2, C.p3);
                tr_bad_poly.emplace_back(C.p3, C.p1);
            }
        }
        /* Remove bad triangles. */
        std::vector<moth_tri_cell2d> tr1;
        for (moth_size_t k = 0; k < tr.size(); ++k) {
            if (!tr_bad[k]) {
                tr1.push_back(tr[k]);
            }
        }
        tr = std::move(tr1);

        /* Find shared edges of the polygon triangulation. */
        std::vector<bool> tr_bad_poly_shared(tr_bad_poly.size());
        for (moth_size_t k = 0; k < tr_bad_poly.size(); ++k) {
            for (moth_size_t m = k + 1; m < tr_bad_poly.size(); ++m) {
                const std::tuple<moth_size_t, moth_size_t>& e1 = tr_bad_poly[k];
                const std::tuple<moth_size_t, moth_size_t>& e2 = tr_bad_poly[m];
                if (e1 == e2 ||
                    e1 == std::make_tuple(std::get<1>(e2), std::get<0>(e2))) {
                    tr_bad_poly_shared[k] = true;
                    tr_bad_poly_shared[m] = true;
                }
            }
        }
        /* Remove shared edges. */
        std::vector<std::tuple<moth_size_t, moth_size_t>> tr_bad_poly1;
        for (moth_size_t k = 0; k < tr_bad_poly.size(); ++k) {
            if (!tr_bad_poly_shared[k]) {
                tr_bad_poly1.push_back(tr_bad_poly[k]);
            }
        }
        tr_bad_poly = std::move(tr_bad_poly1);

        /* Add new triangles. */
        for (const std::tuple<moth_size_t, moth_size_t>& e : tr_bad_poly) {
            moth_tri_cell2d C{std::get<0>(e), std::get<1>(e), i};

            /* Make sure that the triangle is CCW. */
            moth_tri2d T{P[C.p1], P[C.p2], P[C.p3]};
            if (moth_tri2d::area(T) < 0.0) {
                std::swap(C.p1, C.p2);
            }

            tr.push_back(C);
        }
    }

    /* Remove the super triangle. */
    std::vector<moth_tri_cell2d> tr1;
    for (const moth_tri_cell2d& C : tr) {
        if (C.p1 != Cs.p1 && C.p1 != Cs.p2 && C.p1 != Cs.p3 &&
            C.p2 != Cs.p1 && C.p2 != Cs.p2 && C.p2 != Cs.p3 &&
            C.p3 != Cs.p1 && C.p3 != Cs.p2 && C.p3 != Cs.p3) {
            tr1.push_back(C);
        }
    }
    tr = std::move(tr1);
}
