#include "libGeometry2D/src/GeomMesh2.hh"
#include "libGeometry2D/src/GeomSort.hh"

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST MOTH_CORE
moth_mesh2d::moth_mesh2d(moth_size_t capacity)
{
    /* Add super triangle points. */
    pPoints.reserve(2 * capacity);
    pPoints.push_back({-4.0, -4.0});
    pPoints.push_back({+4.0, -4.0});
    pPoints.push_back({ 0.0, +4.0});

    /* Add super triangle. */
    pTriangles.reserve(capacity);
    pTriangles.push_back({0, 1, 2});
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST MOTH_CORE
void moth_mesh2d::insert_unconstrained(const moth_p2d& p1, moth_real_t eps)
{
    /* Round-off the point. */
    moth_p2d p{p1};
    if (eps > 0.0) {
        p.x = std::round(p.x / eps) * eps;
        p.y = std::round(p.y / eps) * eps;
    }

    /* Insert the point. */
    moth_mesh2d_point_iter pP{*this, pPoints.size()};
    pPoints.push_back(p);

    /* Find first bad triangle, reach first bad border triangle and
     * proceed to walk-around it by using the connectivity -- ~O(log(n)). */
    moth_mesh2d_triangle_iter pT_bad{triangle_end()};
    for (moth_mesh2d_triangle_iter pT_cur{pT_bad - 1};
                                   pT_cur != pT_bad; --pT_cur) {
        pT_cur.set_bad(moth_triangle2d::circle(*pT_cur, *pP));
        pT_cur.not_visited();
        if (pT_cur.good()) {
            continue;
        }

        /* Move current triangle triangle to the end,
         * walk-through neighbors of the current triangle
         * until the border is reached -- ~O(log(n)). */
        moth_mesh2d_triangle_iter::swap_assign(pT_cur, --pT_bad);
        for (;;) {
            /* Check neighbors of the bad triangle,
             * also moving bad ones to the end. */
            if (pT_cur.triangle(1).not_visited()) {
                pT_cur.triangle(1).set_bad(moth_triangle2d::circle(*pT_cur.triangle(1), *pP));
                if (pT_cur.triangle(1).bad()) {
                    moth_mesh2d_triangle_iter::swap(pT_cur.triangle(1), --pT_bad);
                }
            }
            if (pT_cur.triangle(1).good()) {
                break;
            }
            if (pT_cur.triangle(2).not_visited()) {
                pT_cur.triangle(2).set_bad(moth_triangle2d::circle(*pT_cur.triangle(2), *pP));
                if (pT_cur.triangle(2).bad()) {
                    moth_mesh2d_triangle_iter::swap(pT_cur.triangle(2), --pT_bad);
                }
            }
            if (pT_cur.triangle(2).good()) {
                /* Presort the triangle to start the walk-around:
                 * make a border oppose the first point. */
                moth_mesh2d_triangle_iter::lshift(pT_cur);
                break;
            }
            if (pT_cur.triangle(3).not_visited()) {
                pT_cur.triangle(3).set_bad(moth_triangle2d::circle(*pT_cur.triangle(3), *pP));
                if (pT_cur.triangle(3).bad()) {
                    moth_mesh2d_triangle_iter::swap(pT_cur.triangle(3), --pT_bad);
                }
            }
            if (pT_cur.triangle(3).good()) {
                /* Presort the triangle to start the walk-around:
                 * make a border oppose the first point. */
                moth_mesh2d_triangle_iter::rshift(pT_cur);
                break;
            }

            /* The border if not reached. Check whether
             * we aren't going for a loop and move on to the next
             * neighbor. */
            if (pT_cur.triangle(1).triangle(1) == pT_cur) {
                moth_mesh2d_triangle_iter::lshift(pT_cur.triangle(1));
            }
            pT_cur = pT_cur.triangle(1);
        }

        /* Walk around a border of the re-triangulation
         * area in CCW orientation inserting the triangles to the back -- ~O(log(n)). */
        moth_mesh2d_point_iter pP_frs{pT_cur.point(2)};
        for (moth_mesh2d_triangle_iter pT_new{*this}, pT_prv{*this}, pT_frs{*this};;) {

            /* Add the new triangle. */
            pTriangles.push_back({pP.nP, pT_cur.point(2).nP,
                                         pT_cur.point(3).nP});
            pT_prv = pT_new;
            pT_new = triangle_end() - 1;
            if (pT_frs.invalid()) {
                pT_frs = pT_new;
            }

            /* Carefully set the neighbors. */
            if (pT_cur.triangle(1).valid()) {
                /* Link the triangle with the outer neighbor.
                 * ( The outer triangle may not be presorted. ) */
                while (pT_cur.point(2) != pT_cur.triangle(1).point(3)) {
                    moth_mesh2d_triangle_iter::lshift(pT_cur.triangle(1));
                }
                pT_cur.triangle(1).set_triangle(1, pT_new);
                pT_new.set_triangle(1, pT_cur.triangle(1));
            }
            if (pT_prv.valid()) {
                /* Link the triangle with the previous one. */
                pT_new.set_triangle(3, pT_prv);
                pT_prv.set_triangle(2, pT_new);
            }
            if (pT_cur.point(3) == pP_frs) {
                /* End point of the last segment is the starting point.
                 * Link the last triangle with the first one and break. */
                pT_frs.set_triangle(3, pT_new);
                pT_new.set_triangle(2, pT_frs);
                break;
            }

            /* Carefully proceed to the new edge:
             * if possible, rotate the current triangle, otherwise select the
             * new triangle with edge CCW continuation. */
            if (pT_cur.triangle(2).not_visited()) {
                pT_cur.triangle(2).set_bad(moth_triangle2d::circle(*pT_cur.triangle(2), *pP));
                if (pT_cur.triangle(2).bad()) {
                    moth_mesh2d_triangle_iter::swap(pT_cur.triangle(2), --pT_bad);
                }
            }
            if (pT_cur.triangle(2).good()) {
                moth_mesh2d_triangle_iter::lshift(pT_cur);
            } else {
                while (pT_cur.point(1) != pT_cur.triangle(2).point(1)) {
                    moth_mesh2d_triangle_iter::lshift(pT_cur.triangle(2));
                }
                pT_cur = pT_cur.triangle(2);
                while (true) {
                    if (pT_cur.triangle(1).not_visited()) {
                        pT_cur.triangle(1).set_bad(moth_triangle2d::circle(*pT_cur.triangle(1), *pP));
                        if (pT_cur.triangle(1).bad()) {
                            moth_mesh2d_triangle_iter::swap(pT_cur.triangle(1), --pT_bad);
                        }
                    }
                    if (pT_cur.triangle(1).bad()) {
                        while (pT_cur.point(3) != pT_cur.triangle(1).point(1)) {
                            moth_mesh2d_triangle_iter::lshift(pT_cur.triangle(1));
                        }
                        pT_cur = pT_cur.triangle(1);
                    } else {
                        break;
                    }
                }
            }
        }
        break;
    }

    /* Walk through bad triangles, move them to the end
     * and delete -- ~O(log(n)). */
    for (moth_mesh2d_triangle_iter pT_cur{pT_bad}, pT_end{triangle_end()};
                                   pT_cur != pT_end; ++pT_cur) {
        if (pT_cur.bad()) {
            moth_mesh2d_triangle_iter::swap(pT_cur, --pT_end);
            pTriangles.pop_back();
        } else {
            break;
        }
    }
}   // void moth_mesh2d::insert_unconstrained(moth_p2d)

MOTH_HOST
void moth_mesh2d::insert_unconstrained(moth_p2d* pP_beg, moth_p2d* pP_end)
{
    /* Presort the triangle so that insertion
     * will take constant time -- ~O(n). */
    const static moth_size_t nSort_cutoff{100};
    if (pP_end - pP_beg > nSort_cutoff + 1) {
        moth_sort(pP_beg, pP_end);
    }

    /* Insert the presorted triangles -- ~O(n). */
    for (moth_p2d* pP_cur{pP_beg};
                   pP_cur != pP_end; ++pP_cur) {
        insert_unconstrained(*pP_cur);
    }
}   // void moth_mesh2d::insert_unconstrained(moth_p2d*, moth_p2d*)

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#if 0
/**
 * Apply constrains to the pre-built Delaunay triangulation.
 */
MOTH_HOST
void moth_mesh2d::apply_constraints()
{
    for (moth_mesh2d_constraint_iter pC_cur{constraint_begin()},
                                     pC_end{constraint_end()};
                                     pC_end != pC_cur; ++pC_cur) {

        /* Walk around the triangles around constraint edge point
         * and check if the constraint is met or
         * subdivide the edge. */
        for (moth_mesh2d_triangle_around_point_iter pT_beg{pC_cur.point(1)},
                                                    pT_cur{pT_beg};;) {




            if (pT_beg == ++pT_cur) {
                break;
            }
        }
    }
}   // void moth_mesh2d::apply_constraints()
#endif
