#pragma once

#include "libGeometry2D/src/GeomTriangle.hh"
#include "libGeometry2D/src/GeomEdge.hh"

#include <vector>
#include <tuple>

#include <fstream>
#include <iostream>

#include <string>
#include <algorithm>
#include <unordered_set>

struct MOTH_CORE moth_cell2d
{
    moth_size_t p1{MOTH_NPOS}, p2{MOTH_NPOS}, p3{MOTH_NPOS};
    moth_size_t n12{MOTH_NPOS}, n23{MOTH_NPOS}, n31{MOTH_NPOS};
};  // moth_cell2d

/**
 * Build Delaunay triangulation for the given set of points.
 * @see https://en.wikipedia.org/wiki/Bowyer–Watson_algorithm
 */
MOTH_HOST MOTH_CORE
static std::vector<moth_cell2d> moth_triangulate_bowyer_watson(std::vector<moth_p2d> P)
{
    std::vector<moth_cell2d> tr{};

    /* Add the super */
    moth_tri2d Sx{{-100, -100}, {+100, -100}, {0, 200}};
    P.push_back(Sx.p1);
    P.push_back(Sx.p2);
    P.push_back(Sx.p3);

    moth_cell2d S{P.size() - 3, P.size() - 2, P.size() - 1};
    tr.push_back(S);

    for (moth_size_t i = 0; i < P.size() - 3; ++i) {
        const moth_p2d& p = P[i];

        /* Find bad triangles and create polygon of them. */
        std::vector<bool> tr_bad(tr.size());
        std::vector<std::tuple<moth_size_t, moth_size_t>> tr_bad_poly;
        for (moth_size_t k = 0; k < tr.size(); ++k) {
            const moth_cell2d& C = tr[k];
            if (tr_bad[k] = moth_tri2d::circle({P[C.p1], P[C.p2], P[C.p3]}, p)) {
                tr_bad_poly.emplace_back(C.p1, C.p2);
                tr_bad_poly.emplace_back(C.p2, C.p3);
                tr_bad_poly.emplace_back(C.p3, C.p1);
            }
        }
        /* Remove bad triangles. */
        std::vector<moth_cell2d> tr1;
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
            moth_cell2d C{std::get<0>(e), std::get<1>(e), i};

            /* Make sure that the triangle is CCW. */
            moth_tri2d T{P[C.p1], P[C.p2], P[C.p3]};
            if (moth_tri2d::area(T) < 0.0) {
                std::swap(C.p1, C.p2);
            }
            tr.push_back(C);
        }
    }

    /* Remove the super triangle. */
#if 0
    std::vector<moth_cell2d> tr1;
    for (const moth_cell2d& C : tr) {
        if (C.p1 != S.p1 && C.p1 != S.p2 && C.p1 != S.p3 &&
            C.p2 != S.p1 && C.p2 != S.p2 && C.p2 != S.p3 &&
            C.p3 != S.p1 && C.p3 != S.p2 && C.p3 != S.p3) {
            tr1.push_back(C);
        }
    }
    tr = std::move(tr1);
#endif

#if 1
    std::string TrFilePath("res/tr-" + std::to_string(99999) + ".txt");
    std::ofstream TrFile(TrFilePath);
    for (const moth_cell2d& T : tr) {
        moth_tri2d Tx{P[T.p1], P[T.p2], P[T.p3]};
        TrFile << Tx.p1 << std::endl
               << Tx.p2 << std::endl
               << Tx.p3 << std::endl
               << Tx.p1 << std::endl
               << std::endl;
    }
#endif

    return tr;
}

/**
 *
 * @see https://en.wikipedia.org/wiki/Ruppert%27s_algorithm
 */
MOTH_HOST MOTH_CORE
static std::vector<moth_cell2d> moth_triangulate_ruppert(std::vector<moth_e2d> E)
{
    std::unordered_set<moth_p2d> Ps{};
    for (const moth_e2d& e : E) {
        Ps.insert(e.p1);
        Ps.insert(e.p2);
    }
    std::vector<moth_p2d> P{Ps.begin(), Ps.end()};


    std::vector<moth_cell2d> tr = moth_triangulate_bowyer_watson(P);
    for (moth_size_t m = 0; m < 3; ++m) {

        /* Find the bad-quality triangles. */
        std::vector<moth_size_t> tr_bad{};
        for (moth_size_t k = 0; k < tr.size(); ++k) {
            const moth_cell2d& C = tr[k];

            moth_tri2d T = {P[C.p1], P[C.p2], P[C.p3]};
            moth_real_t T_area = std::fabs(moth_tri2d::area(T));
            moth_real_t T_len = moth_tri2d::len(T);
            moth_real_t T_quality = T_area / T_len;
            if (std::fabs(1.0 - T_quality) < 1.0 - 0.0) {
                tr_bad.push_back(k);
            }
        }
        if (tr_bad.empty()) {
            break;
        }

        if (m == 0)
            P.push_back({{},{}}); else
        for (moth_size_t k : tr_bad) {
            const moth_cell2d& C = tr[k];
            moth_tri2d T = {P[C.p1], P[C.p2], P[C.p3]};
            moth_p2d CC{moth_tri2d::circumcenter(T)};
            P.push_back(CC);
        }

        tr = moth_triangulate_bowyer_watson(P);
    }

    std::string TrFilePath("res/tr-" + std::to_string(99999) + ".txt");
    std::ofstream TrFile(TrFilePath);
    for (const moth_cell2d& T : tr) {
        moth_tri2d Tx{P[T.p1], P[T.p2], P[T.p3]};
        TrFile << Tx.p1 << std::endl
               << Tx.p2 << std::endl
               << Tx.p3 << std::endl
               << Tx.p1 << std::endl
               << std::endl;
    }

    return tr;
}
