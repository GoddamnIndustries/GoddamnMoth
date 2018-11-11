#include "libGeometry2D/src/GeomMesh2.hh"

int DT::cnt = 0;

MOTH_HOST MOTH_CORE
void DT::moth_mesh2d::insert(const moth_p2d& p1, moth_real_t eps)
{
    ++cnt;

    /* Round-off the point. */
    moth_p2d p{p1};
    if (eps > 0.0) {
        p.x = std::round(p.x / eps) * eps;
        p.y = std::round(p.y / eps) * eps;
    }

    /* Insert the point. */
    moth_mesh2d_point_iter pP{pPoints.data(), pTriangles.data(), pPoints.size()};
    pPoints.push_back(p);

    /* Find first bad triangle, reach first bad border triangle and
     * proceed to walk-around it by using the connectivity -- ~O(log(n)). */
    moth_mesh2d_triangle_iter pT_bad = triangle_end();
    for (moth_mesh2d_triangle_iter pT_cur{pT_bad - 1};
                                   pT_cur != pT_bad; --pT_cur) {
        pT_cur.bbad() = moth_triangle2d::circle(*pT_cur, *pP);
        pT_cur.unvisited();
        if (pT_cur.good()) {
            continue;
        }

        /* Move current triangle triangle to the end,
         * walk-through neighbors of the current triangle
         * until the border is reached -- ~O(1). */
        moth_mesh2d_triangle_iter::swap_assign(pT_cur, --pT_bad);
        for (;;) {
            /* Check neighbors of the bad triangle,
             * also moving bad ones to the end. */
            if (pT_cur.triangle(1).unvisited()) {
                pT_cur.triangle(1).bbad() = moth_triangle2d::circle(*pT_cur.triangle(1), *pP);
                if (pT_cur.triangle(1).bad()) {
                    moth_mesh2d_triangle_iter::swap(pT_cur.triangle(1), --pT_bad);
                }
            }
            if (pT_cur.triangle(1).good()) {
                break;
            }

            if (pT_cur.triangle(2).unvisited()) {
                pT_cur.triangle(2).bbad() = moth_triangle2d::circle(*pT_cur.triangle(2), *pP);
                if (pT_cur.triangle(2).bad()) {
                    moth_mesh2d_triangle_iter::swap(pT_cur.triangle(2), --pT_bad);
                }
            }
            if (pT_cur.triangle(2).good()) {
                /* Presort the triangle to start the walk-around:
                 * make a border oppose the first point. */
                moth_mesh2d_utils::lshift(pT_cur);
                break;
            }

            if (pT_cur.triangle(3).unvisited()) {
                pT_cur.triangle(3).bbad() = moth_triangle2d::circle(*pT_cur.triangle(3), *pP);
                if (pT_cur.triangle(3).bad()) {
                    moth_mesh2d_triangle_iter::swap(pT_cur.triangle(3), --pT_bad);
                }
            }
            if (pT_cur.triangle(3).good()) {
                /* Presort the triangle to start the walk-around:
                 * make a border oppose the first point. */
                moth_mesh2d_utils::rshift(pT_cur);
                break;
            }

            /* The border if not reached. Check whether
             * we aren't going for a loop and move on to the next
             * neighbor. */
            if (pT_cur.triangle(1).triangle(1) == pT_cur) {
                moth_mesh2d_utils::lshift(pT_cur.triangle(1));
            }
            pT_cur = pT_cur.triangle(1);
        }

        /* Walk around a border of the re-triangulation
         * area in CCW orientation inserting the triangles to the back. */
        moth_mesh2d_point_iter pP_frs{pT_cur.point(2)};
        for (moth_mesh2d_triangle_iter pT_new{}, pT_prv{}, pT_frs{};;) {

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
                    moth_mesh2d_utils::lshift(pT_cur.triangle(1));
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

            /* Check the second neighbor of the bad triangle,
             * also moving bad it to the end, if it is bad. */
            if (pT_cur.triangle(2).unvisited()) {
                pT_cur.triangle(2).bbad() = moth_triangle2d::circle(*pT_cur.triangle(2), *pP);
                if (pT_cur.triangle(2).bad()) {
                    moth_mesh2d_triangle_iter::swap(pT_cur.triangle(2), --pT_bad);
                }
            }

            /* Carefully proceed to the new edge:
             * if possible, rotate the current triangle, otherwise select the
             * new triangle with edge CCW continuation. */
            if (pT_cur.triangle(2).good()) {
                moth_mesh2d_utils::lshift(pT_cur);
            } else {
                while (pT_cur.point(1) != pT_cur.triangle(2).point(1)) {
                    moth_mesh2d_utils::lshift(pT_cur.triangle(2));
                }
                pT_cur = pT_cur.triangle(2);
                while (true) {
                    if (pT_cur.triangle(1).unvisited()) {
                        pT_cur.triangle(1).bbad() = moth_triangle2d::circle(*pT_cur.triangle(1), *pP);
                        if (pT_cur.triangle(1).bad()) {
                            moth_mesh2d_triangle_iter::swap(pT_cur.triangle(1), --pT_bad);
                        }
                    }
                    if (pT_cur.triangle(1).bad()) {
                        while (pT_cur.point(3) != pT_cur.triangle(1).point(1)) {
                            moth_mesh2d_utils::lshift(pT_cur.triangle(1));
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
}
