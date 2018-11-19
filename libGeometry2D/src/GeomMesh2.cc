#include "libGeometry2D/src/GeomMesh2.hh"
#include "libGeometry2D/src/GeomSort.hh"

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST MOTH_CORE
moth_mesh2d::moth_mesh2d(moth_size_t capacity)
{
    /* Add super triangle points. */
    pPoints.reserve(2 * capacity);
    pPoints.push_back({-400.0, -400.0});
    pPoints.push_back({+400.0, -400.0});
    pPoints.push_back({   0.0, +400.0});

    /* Add super triangle. */
    pTriangles.reserve(capacity);
    pTriangles.push_back({0, 1, 2});
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST MOTH_CORE
moth_mesh2d_point_iter
moth_mesh2d::insert_unconstrained(const moth_p2d& p1, moth_real_t eps)
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
     * proceed to walk-around by using the connectivity -- ~O(log(n)). */
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
         * until the border is reached -- ~O(1). */
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
         * area in CCW orientation inserting the triangles to the back -- ~O(1). */
        moth_mesh2d_point_iter pP_beg{pT_cur.point(2)};
        for (moth_mesh2d_triangle_iter pT_new{*this}, pT_prv{*this}, pT_beg{*this};;) {
            /* Add the new triangle. */
            pTriangles.push_back({pP.nP, pT_cur.point(2).nP,
                                         pT_cur.point(3).nP});
            pT_prv = pT_new;
            pT_new = triangle_end() - 1;
            if (pT_beg.invalid()) {
                pT_beg = pT_new;
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
            if (pT_cur.point(3) == pP_beg) {
                /* End point of the last segment is the starting point.
                 * Link the last triangle with the first one and break. */
                pT_beg.set_triangle(3, pT_new);
                pT_new.set_triangle(2, pT_beg);
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
     * and delete -- ~O(1). */
    for (moth_mesh2d_triangle_iter pT_cur{pT_bad}, pT_end{triangle_end()};
                                   pT_cur != pT_end; ++pT_cur) {
        if (pT_cur.bad()) {
            moth_mesh2d_triangle_iter::swap(pT_cur, --pT_end);
            pTriangles.pop_back();
        } else {
            break;
        }
    }

    return pP;
}   // moth_mesh2d::insert_unconstrained

MOTH_HOST
moth_mesh2d_point_iter
moth_mesh2d::insert_unconstrained(moth_p2d* pP_beg,
                                  moth_p2d* pP_end)
{
    /* Presort the triangle so that insertion
     * will take constant time -- ~O(n). */
    const static moth_size_t nSort_cutoff{100};
    if (pP_end - pP_beg > nSort_cutoff + 1) {
        moth_sort(pP_beg, pP_end);
    }

    /* Insert the presorted triangles -- ~O(n). */
    moth_mesh2d_point_iter pP{insert_unconstrained(*pP_beg)};
    for (moth_p2d* pP_cur{pP_beg + 1};
                   pP_cur != pP_end; ++pP_cur) {
        insert_unconstrained(*pP_cur);
    }

    return pP;
}   // moth_mesh2d::insert_unconstrained

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST
void moth_mesh2d::apply_constrain_conforming(const moth_mesh2d_cedge_iter& pE)
{
    /* Find first bad intersected triangle, proceed to walk-around
     * the intersected area by using the connectivity -- ~O(1). */
    moth_mesh2d_cedge_iter pE_cur{pE}, pE_int{*this};
    for (moth_mesh2d_point_iter pP_cur{pE.point(1)}, pP_end{pE.point(2)};
                                pP_cur != pP_end;) {
        for (moth_mesh2d_triangle_iter pT_cur{pP_cur.triangle()};;
                                       pT_cur = pT_cur.triangle(2)) {
            /* Presort the triangle to start the
             * forward walk-around. */
            while (pP_cur != pT_cur.point(1)) {
                moth_mesh2d_triangle_iter::lshift(pT_cur);
            }

            /* Check if the constraint was met,
             * or the triangle intersects the constraint edge,
             * otherwise select the next triangle around the starting point. */
            if (pT_cur.point(2) == pP_end ||
                pT_cur.point(3) == pP_end) {
                pP_cur = pP_end;
                break;
            } else {
                moth_p2d p_int{};
                if (moth_e2d::intersect((*pT_cur).edge(1), *pE_cur, p_int)) {
                    /* Tessellate the constraint edge
                     * by inserting the intersection point. */
                    moth_mesh2d_point_iter pP_int{insert_unconstrained(p_int)};
                    pConstraints.push_back({pP_int.nP,
                                                   pP_end.nP});

                    pE_int = --constraint_end();
                    pE_cur.set_point(2, pP_int);
                    if (pE_cur.edge(2).valid()) {
                        pE_cur.edge(2).set_edge(1, pE_int);
                        pE_int.set_edge(2, pE_cur.edge(2));
                    }
                    pE_int.set_edge(1, pE_cur);
                    pE_cur.set_edge(2, pE_int);

                    /* Proceed to the next triangle. */
                    pE_cur = pE_int;
                    pP_cur = pP_int;
                    break;
                }
            }
        }
    }
}   // moth_mesh2d::apply_constrain_conforming

MOTH_HOST
void moth_mesh2d::apply_constrains_conforming()
{
    /* Apply all conforming constraints. */
    for (moth_mesh2d_cedge_iter pE_cur{constraint_begin()}, pE_end{constraint_end()};
                                pE_cur != pE_end; ++pE_cur) {
        apply_constrain_conforming(pE_cur);
    }
}   // moth_mesh2d::apply_constrains_conforming

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST
void moth_mesh2d::apply_constrain_ignoring(const moth_mesh2d_cedge_iter& pE)
{
    moth_mesh2d_point_iter pP_cur = pE.point(1);
    moth_mesh2d_point_iter pP_end = pE.point(2);

    /* Find first bad intersected triangle, proceed to walk-around
     * the intersected area by using the connectivity -- ~O(1). */
    moth_mesh2d_triangle_iter pT_bad{triangle_end()};
    for (moth_mesh2d_triangle_iter pT_cur{pP_cur.triangle()};;
                                   pT_cur = pT_cur.triangle(2)) {
        /* Presort the triangle to start the
         * forward walk-around. */
        while (pP_cur != pT_cur.point(1)) {
            moth_mesh2d_triangle_iter::lshift(pT_cur);
        }

        /* Check if the constraint was met,
         * or the triangle intersects the constraint edge,
         * otherwise select the next triangle around the starting point. */
        if (pT_cur.point(2) == pP_end ||
            pT_cur.point(3) == pP_end) {
            break;
        } else {
            if (!moth_e2d::intersect((*pT_cur).edge(1), *pE)) {
                continue;
            }
        }

        /* Move the current triangle to the end,
         * walk around edges until the end point is reached,
         * also moving the intersected triangles the end. */
        pT_cur.set_bad();
        moth_mesh2d_triangle_iter::swap_assign(pT_cur, --pT_bad);
        for (;;) {
            /* Carefully proceed to the new edge:
             * select the new triangle with edge CCW continuation. */
            moth_mesh2d_triangle_iter::swap(pT_cur.triangle(1), --pT_bad);
            while (pT_cur.point(3) != pT_cur.triangle(1).point(3)) {
                moth_mesh2d_triangle_iter::lshift(pT_cur.triangle(1));
            }
            pT_cur = pT_cur.triangle(1);
            pT_cur.set_bad();
            while (!moth_e2d::intersect((*pT_cur).edge(1), *pE)) {
                moth_mesh2d_triangle_iter::swap(pT_cur.triangle(3), --pT_bad);
                while (pT_cur.point(1) != pT_cur.triangle(3).point(1)) {
                    moth_mesh2d_triangle_iter::lshift(pT_cur.triangle(3));
                }
                pT_cur = pT_cur.triangle(3);
                pT_cur.set_bad();
            }

            if (pT_cur.point(2) == pP_end) {
                /* End point is reached,
                 * recursively construct triangles in the intersected area,
                 * using the Growing Circles method and break. */
                moth_mesh2d_triangle_iter pT_lwr{apply_constrain_ignoring_impl(pP_cur, pP_end)};
                moth_mesh2d_triangle_iter pT_upr{apply_constrain_ignoring_impl(pP_end, pP_cur)};
                pT_lwr.set_triangle(1, pT_upr);
                pT_upr.set_triangle(1, pT_lwr);
                break;
            }
        }
        break;
    }

    /* Walk through intersected triangles, move them to the end
     * and delete -- ~O(1). */
    for (moth_mesh2d_triangle_iter pT_cur{pT_bad}, pT_end{triangle_end()};
                                   pT_cur != pT_end; ++pT_cur) {
        if (pT_cur.bad()) {
            moth_mesh2d_triangle_iter::swap(pT_cur, --pT_end);
            pTriangles.pop_back();
        } else {
            break;
        }
    }
}   // moth_mesh2d::apply_constrain_ignoring

MOTH_HOST
moth_mesh2d_triangle_iter
moth_mesh2d::apply_constrain_ignoring_impl(moth_mesh2d_point_iter pP_beg,
                                           moth_mesh2d_point_iter pP_end)
{
    moth_mesh2d_triangle_iter pT_new{*this}, pT_beg{*this};
    for (moth_mesh2d_triangle_iter pT_cur{pP_beg.triangle()};;
                                   pT_cur = pT_cur.triangle(2)) {
        /* Presort the triangle to start the walk-around. */
        while (pP_beg != pT_cur.point(1)) {
            moth_mesh2d_triangle_iter::lshift(pT_cur);
        }

        /* Check if the constraint was met,
         * or the triangle intersects the constraint edge,
         * otherwise select the next triangle around the starting point. */
        if (pT_cur.point(2) == pP_end ||
            pT_cur.point(3) == pP_end) {
            break;
        } else {
            if (!moth_e2d::intersect((*pT_cur).edge(1), {*pP_beg, *pP_end})) {
                continue;
            }
        }

        moth_mesh2d_point_iter pP_mid{pT_cur.point(2)};
        for (pT_beg = pT_cur;;) {
            /* Carefully proceed to the new edge:
             * select the new triangle with edge CCW continuation. */
            while (pT_cur.point(3) != pT_cur.triangle(1).point(3)) {
                moth_mesh2d_triangle_iter::lshift(pT_cur.triangle(1));
            }
            pT_cur = pT_cur.triangle(1);
            while (!moth_e2d::intersect((*pT_cur).edge(1), {*pP_beg, *pP_end})) {
                while (pT_cur.point(1) != pT_cur.triangle(3).point(1)) {
                    moth_mesh2d_triangle_iter::lshift(pT_cur.triangle(3));
                }
                pT_cur = pT_cur.triangle(3);
                pT_cur.set_bad();
            }

            /* Find the closest point to the edge. */
            if (pT_cur.point(2) != pP_end) {
                if (moth_e2d::dist({*pP_beg, *pP_end}, *pT_cur.point(2)) <
                    moth_e2d::dist({*pP_beg, *pP_end}, *pP_mid)) {
                    pP_mid = pT_cur.point(2);
                }
            } else {
                break;
            }
        }

        /* Add the new triangle. */
        pTriangles.push_back({pP_mid.nP, pP_end.nP,
                                         pP_beg.nP});
        pT_new = --triangle_end();

        /* Carefully set the neighbors. */
        if (pP_mid == pT_beg.point(1) ||
            pP_mid == pT_beg.point(2) ||
            pP_mid == pT_beg.point(3)) {
            /* Link the triangle with the outer neighbor.
             * ( The outer triangle may not be presorted. ) */
            while (pT_beg.point(2) != pT_beg.triangle(3).point(1)) {
                moth_mesh2d_triangle_iter::lshift(pT_beg.triangle(3));
            }
            pT_new.set_triangle(2, pT_beg.triangle(3));
            pT_new.triangle(2).set_triangle(3, pT_new);
        } else {
            /* Recursively re-triangulating the inner edge. */
            pT_new.set_triangle(2, apply_constrain_ignoring_impl(pP_beg, pP_mid));
            pT_new.triangle(2).set_triangle(1, pT_new);
        }
        if (pP_mid == pT_cur.point(1) ||
            pP_mid == pT_cur.point(2) ||
            pP_mid == pT_cur.point(3)) {
            /* Link the triangle with the outer neighbor.
             * ( The outer triangle may not be presorted. ) */
            while (pT_cur.point(2) != pT_cur.triangle(3).point(1)) {
                moth_mesh2d_triangle_iter::lshift(pT_cur.triangle(3));
            }
            pT_new.set_triangle(2, pT_beg.triangle(3));
            pT_new.triangle(2).set_triangle(3, pT_new);
        } else {
            /* Recursively re-triangulating the inner edge. */
            pT_new.set_triangle(3, apply_constrain_ignoring_impl(pP_mid, pP_end));
            pT_new.triangle(3).set_triangle(1, pT_new);
        }
        break;
    }
    return pT_new;
}   // moth_mesh2d::apply_constrain_ignoring_impl

MOTH_HOST
void moth_mesh2d::apply_constrains_ignoring()
{
    /* Apply all conforming constraints. */
    for (moth_mesh2d_cedge_iter pE_cur{constraint_begin()}, pE_end{constraint_end()};
                                pE_cur != pE_end; ++pE_cur) {
        apply_constrain_ignoring(pE_cur);
    }
}   // moth_mesh2d::apply_constrains_ignoring()

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST
moth_mesh2d_point_iter moth_mesh2d::insert_constrain(const moth_poly2d& poly)
{
    moth_poly2d_iter pE_cur{poly.iter()};
    moth_mesh2d_point_iter pP_frs{insert_unconstrained(pE_cur.point())};
    while (++pE_cur != poly.iter()) {
        insert_unconstrained(pE_cur.point());
        pConstraints.push_back({pPoints.size() - 2, pPoints.size() - 1});
    }
    pConstraints.push_back({pPoints.size() - 1, pP_frs.nP});
    apply_constrains_conforming();
    return {*this};
}   // moth_mesh2d_point_iter moth_mesh2d::insert_constrain(moth_poly2d)
