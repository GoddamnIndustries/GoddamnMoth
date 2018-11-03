#pragma once

#include <vector>
#include <algorithm>
#include <fstream>

#include <cassert>
#include <type_traits>

#include "libGeometry2D/src/GeomBase.hh"
#include "libGeometry2D/src/GeomTriangle.hh"

#undef assert
#define assert(...)

namespace DT {

struct moth_mesh2d_point;
struct moth_mesh2d_point_iter;
struct moth_mesh2d_triangle;
struct moth_mesh2d_triangle_iter;

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Triangle cell.
 */
struct MOTH_CORE moth_mesh2d_triangle
{
public:
    moth_size_t nP1{MOTH_NPOS}, nP2{MOTH_NPOS}, nP3{MOTH_NPOS};
    moth_size_t nT1{MOTH_NPOS}, nT2{MOTH_NPOS}, nT3{MOTH_NPOS};
    bool bad{};
    bool deleted{};

public:
    MOTH_HOST MOTH_DEVICE
    moth_size_t& nnP(size_t k)
    {
        switch (k) {
            case 1:
                return nP1;
            case 2:
                return nP2;
            default:
                return nP3;
        }
        //return *(&nP1 + (k - 1) % 3);
    }
    MOTH_HOST MOTH_DEVICE
    moth_size_t nnP(size_t k) const
    {
        return *(&nP1 + (k - 1) % 3);
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_size_t& nnT(size_t k)
    {
        return *(&nT1 + (k - 1) % 3);
    }
    MOTH_HOST MOTH_DEVICE
    moth_size_t nnT(size_t k) const
    {
        return *(&nT1 + (k - 1) % 3);
    }
};  // struct moth_mesh2d_triangle

    struct MOTH_CORE moth_mesh2d_utils
    {
    public:
        MOTH_HOST
        static void shift(moth_mesh2d_triangle& T)
        {
            std::swap(T.nP1, T.nP2);
            std::swap(T.nP2, T.nP3);
            std::swap(T.nT1, T.nT2);
            std::swap(T.nT2, T.nT3);
        }
        MOTH_HOST template<typename T = moth_mesh2d_triangle_iter>
        static void shift(const T& pT)
        {
            shift(pT.pTriangles[pT.nT]);
        }
    };  // struct moth_mesh2d_utils

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Triangulation points iterator.
 */
struct MOTH_CORE moth_mesh2d_point_iter
{
    moth_p2d* pPoints{};
    moth_mesh2d_triangle* pTriangles{};
    moth_size_t nP{MOTH_NPOS};

public:
    MOTH_HOST MOTH_DEVICE
    moth_p2d operator* () const
    {
        assert(nP != MOTH_NPOS);
        return pPoints[nP];
    }

public:
    MOTH_HOST MOTH_DEVICE
    bool operator==(const moth_mesh2d_point_iter& pP) const
    {
        assert(pP.pPoints == pPoints &&
               pP.pTriangles == pTriangles);
        return nP == pP.nP;
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_mesh2d_point_iter& pP) const
    {
        assert(pP.pPoints == pPoints &&
               pP.pTriangles == pTriangles);
        return nP != pP.nP;
    }


public:
    MOTH_HOST MOTH_DEVICE
    bool valid() const
    {
        return nP != MOTH_NPOS;
    }
    MOTH_HOST MOTH_DEVICE
    bool invalid() const
    {
        return nP == MOTH_NPOS;
    }
};  // struct MOTH_CORE moth_mesh2d_point_iter

/**
 * Triangle cells iterator.
 */
struct MOTH_CORE moth_mesh2d_triangle_iter
{
    moth_p2d* pPoints{};
    moth_mesh2d_triangle* pTriangles{};
    moth_size_t nT{MOTH_NPOS};

public:
    MOTH_HOST MOTH_DEVICE
    moth_triangle2d operator* () const
    {
        return {*point(1), *point(2), *point(3)};
    }

public:
    MOTH_HOST MOTH_DEVICE
    bool operator==(const moth_mesh2d_triangle_iter& pT) const
    {
        assert(pT.pPoints == pPoints &&
               pT.pTriangles == pTriangles);
        return nT == pT.nT;
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_mesh2d_triangle_iter& pT) const
    {
        assert(pT.pPoints == pPoints &&
               pT.pTriangles == pTriangles);
        return nT != pT.nT;
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_triangle_iter& operator++()
    {
        assert(nT != MOTH_NPOS);
        nT += 1;
        return *this;
    }
    MOTH_HOST MOTH_DEVICE
    const moth_mesh2d_triangle_iter operator++(int)
    {
        moth_mesh2d_triangle_iter copy(*this);
        ++*this;
        return copy;
    }

    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_triangle_iter& operator--()
    {
        assert(nT != MOTH_NPOS);
        nT -= 1;
        return *this;
    }
    MOTH_HOST MOTH_DEVICE
    const moth_mesh2d_triangle_iter operator--(int)
    {
        moth_mesh2d_triangle_iter copy(*this);
        --*this;
        return copy;
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_triangle_iter& operator+=(moth_diff_t d)
    {
        assert(nT != MOTH_NPOS);
        nT += d;
        return *this;
    }
    MOTH_HOST MOTH_DEVICE
    const moth_mesh2d_triangle_iter operator+(moth_diff_t d) const
    {
        moth_mesh2d_triangle_iter copy(*this);
        return copy += d;
    }

    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_triangle_iter& operator-=(moth_diff_t d)
    {
        assert(nT != MOTH_NPOS);
        nT -= d;
        return *this;
    }
    MOTH_HOST MOTH_DEVICE
    const moth_mesh2d_triangle_iter operator-(moth_diff_t d) const
    {
        moth_mesh2d_triangle_iter copy(*this);
        return copy -= d;
    }

public:
    MOTH_HOST MOTH_DEVICE
    bool valid() const
    {
        return nT != MOTH_NPOS;
    }
    MOTH_HOST MOTH_DEVICE
    bool invalid() const
    {
        return nT == MOTH_NPOS;
    }

    MOTH_HOST
    bool& bbad()
    {
        return pTriangles[nT].bad;
    }

    MOTH_HOST
    bool& ddel()
    {
        return pTriangles[nT].deleted;
    }

    MOTH_HOST
    bool bad() const
    {
        return valid() && pTriangles[nT].bad;
    }
    MOTH_HOST
    bool good()
    {
        return invalid() || !pTriangles[nT].bad;
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_point_iter point(size_t k) const
    {
        assert(nT != MOTH_NPOS);
        return {pPoints, pTriangles, pTriangles[nT].nnP(k)};
    }
    MOTH_HOST MOTH_DEVICE
    void set_point(size_t k, const moth_mesh2d_point_iter& pP)
    {
        assert(nT != MOTH_NPOS);
        assert(pP.pPoints == pPoints &&
               pP.pTriangles == pTriangles);
        pTriangles[nT].nnP(k) = MOTH_NPOS;
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_triangle_iter triangle(size_t k) const
    {
        assert(nT != MOTH_NPOS);
        return {pPoints, pTriangles, pTriangles[nT].nnT(k)};
    }
    MOTH_HOST MOTH_DEVICE
    void set_triangle(size_t k, const moth_mesh2d_triangle_iter& pT)
    {
        assert(nT != MOTH_NPOS);
        assert(pT.pPoints == pPoints &&
               pT.pTriangles == pTriangles);
        pTriangles[nT].nnT(k) = pT.nT;
    }

public:
    MOTH_HOST MOTH_DEVICE
    static void swap(moth_mesh2d_triangle_iter& pT1, moth_mesh2d_triangle_iter& pT2)
    {
        assert(pT1.nT != MOTH_NPOS && pT2.nT != MOTH_NPOS);
        assert(pT1.pPoints == pT2.pPoints &&
               pT1.pTriangles == pT2.pTriangles);
        if (pT1 != pT2) {
            if (pT1 != pT2.triangle(1) &&
                pT1 != pT2.triangle(2) &&
                pT1 != pT2.triangle(3)) {

                /* Re-link neighbors of the first triangle to the second one. */
                for (moth_size_t k = 1; k <= 3; ++k) {
                    moth_mesh2d_triangle_iter pT1_k = pT1.triangle(k);
                    if (pT1_k.valid()) {
                        for (moth_size_t m = 1; m <= 3; ++m) {
                            if (pT1_k.triangle(m) == pT1) {
                                pT1_k.set_triangle(m, pT2);
                                break;
                            }
                        }
                    }
                }
                /* Re-link neighbors of the second triangle to the first one. */
                for (moth_size_t k = 1; k <= 3; ++k) {
                    moth_mesh2d_triangle_iter pT2_k = pT2.triangle(k);
                    if (pT2_k.valid()) {
                        for (moth_size_t m = 1; m <= 3; ++m) {
                            if (pT2_k.triangle(m) == pT2) {
                                pT2_k.set_triangle(m, pT1);
                                break;
                            }
                        }
                    }
                }
            } else {

                /* Relink the triangles with common edges first. */
                while (pT1 != pT2.triangle(1)) {
                    moth_mesh2d_utils::shift(pT2);
                }
                while (pT2 != pT1.triangle(1)) {
                    moth_mesh2d_utils::shift(pT1);
                }
                std::swap(pT1.pTriangles[pT1.nT].nT1,
                          pT2.pTriangles[pT2.nT].nT1);

                /* Re-link other neighbors of the first triangle to the second one. */
                for (moth_size_t k = 2; k <= 3; ++k) {
                    moth_mesh2d_triangle_iter pT1_k = pT1.triangle(k);
                    if (pT1_k.valid()) {
                        for (moth_size_t m = 1; m <= 3; ++m) {
                            if (pT1_k.triangle(m) == pT1) {
                                pT1_k.set_triangle(m, pT2);
                                break;
                            }
                        }
                    }
                }
                /* Re-link other neighbors of the second triangle to the first one. */
                for (moth_size_t k = 2; k <= 3; ++k) {
                    moth_mesh2d_triangle_iter pT2_k = pT2.triangle(k);
                    if (pT2_k.valid()) {
                        for (moth_size_t m = 1; m <= 3; ++m) {
                            if (pT2_k.triangle(m) == pT2) {
                                pT2_k.set_triangle(m, pT1);
                                break;
                            }
                        }
                    }
                }
            }

            /* And finally swap the memory. */
            std::swap(pT1.pTriangles[pT1.nT], pT2.pTriangles[pT2.nT]);
        }
    }
};  // struct moth_mesh2d_triangle_iter

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * 2D Delaunay triangular mesh.
 */
struct MOTH_CORE moth_mesh2d
{
    std::vector<moth_p2d> pPoints;
    std::vector<moth_mesh2d_triangle> pTriangles;

public:
    MOTH_HOST
    moth_mesh2d(moth_size_t capacity = 10000000)
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

public:
    MOTH_HOST
    moth_mesh2d_triangle_iter triangle_begin()
    {
        return {pPoints.data(), pTriangles.data(), 0};
    }
    MOTH_HOST
    moth_mesh2d_triangle_iter triangle_end()
    {
        return {pPoints.data(), pTriangles.data(), pTriangles.size()};
    }

public:
    MOTH_HOST
    void insert(const moth_p2d& p1)
    {
        /* Round-off the point. */
        moth_p2d p{p1};
        //p.x = std::round(p.x / 0.0001) * 0.0001;
        //p.y = std::round(p.y / 0.0001) * 0.0001;

        /* Insert the point -- =O(1). */
        moth_mesh2d_point_iter pP{pPoints.data(), pTriangles.data(), pPoints.size()};
        pPoints.push_back(p);

        /* Find bad triangles and move them to the end -- O(n).
         * @todo Optimize this by using the connectivity -- O(log(n)). */
        for (moth_mesh2d_triangle_iter cur = triangle_begin(),
                                       end = triangle_end(); cur != end;) {
            cur.bbad() = moth_triangle2d::circle(*cur, *pP);
            if (cur.bad()) {
                moth_mesh2d_triangle_iter::swap(cur, --end);
            } else {
                ++cur;
            }
        }

        /* Find the first bad border triangle and proceed to walk-around it
         * by using the connectivity -- ~O(log(n)). */
        for (moth_mesh2d_triangle_iter cur = triangle_begin(),
                                       end = triangle_end(); cur != end; ++cur) {
            if (cur.good() || (cur.triangle(1).bad() &&
                               cur.triangle(2).bad() &&
                               cur.triangle(3).bad())) {
                continue;
            }

            /* Presort the triangle to start the walk-around:
             * make a border oppose the first point. */
            if (cur.triangle(2).good()) {
                moth_mesh2d_utils::shift(cur);  /* 1,2,3 -> 2,1,3 */
            } else if (cur.triangle(3).good()) {
                moth_mesh2d_utils::shift(cur);  /* 1,2,3 -> 2,1,3 */
                moth_mesh2d_utils::shift(cur);  /* 2,1,3 -> 3,2,1 */
            }

            /* Walk around a border of the re-triangulation
             * area in CCW orientation inserting the triangles to the back. */
            moth_mesh2d_point_iter pP_frs{cur.point(2)};
            for (moth_mesh2d_triangle_iter add{}, prv{}, frs{};;) {

                /* Add the new triangle. */
                pTriangles.push_back({pP.nP, cur.point(2).nP,
                                             cur.point(3).nP});
                prv = add;
                add = triangle_end() - 1;
                if (frs.invalid()) {
                    frs = add;
                }

                /* Carefully set the neighbors. */
                if (cur.triangle(1).valid()) {
                    /* Link the triangle with the outer neighbor.
                     * ( The outer triangle may not be presorted. )*/
                    while (cur.point(2) != cur.triangle(1).point(3)) {
                        moth_mesh2d_utils::shift(cur.triangle(1));
                    }
                    cur.triangle(1).set_triangle(1, add);
                    add.set_triangle(1, cur.triangle(1));
                }
                if (prv.valid()) {
                    /* Link the triangle with the previous one. */
                    add.set_triangle(3, prv);
                    prv.set_triangle(2, add);
                }
                if (cur.point(3) == pP_frs) {
                    /* End point of the last segment is the starting point.
                     * Link the last triangle with the first one and break. */
                    frs.set_triangle(3, add);
                    add.set_triangle(2, frs);
                    break;
                }

                /* Carefully proceed to the new edge:
                 * if possible, rotate the current triangle, otherwise select the
                 * new triangle with edge CCW continuation. */
                if (cur.triangle(2).good()) {
                    moth_mesh2d_utils::shift(cur);
                } else {
                    while (cur.point(1) != cur.triangle(2).point(1)) {
                        moth_mesh2d_utils::shift(cur.triangle(2));
                    }
                    cur = cur.triangle(2);
                    while (cur.triangle(1).bad()) {
                        while (cur.point(3) != cur.triangle(1).point(1)) {
                            moth_mesh2d_utils::shift(cur.triangle(1));
                        }
                        cur = cur.triangle(1);
                    }
                }
            }
            break;
        }

        /* Find bad triangles, move them to the end and delete -- O(n).
         * @todo Optimize this by using the connectivity -- O(log(n)). */
        for (moth_mesh2d_triangle_iter cur = triangle_begin(),
                                       end = triangle_end(); cur != end;) {
            if (cur.bad()) {
                moth_mesh2d_triangle_iter::swap(cur, --end);
                pTriangles.pop_back();
            } else {
                ++cur;
            }
        }
    }

    void print(bool center = false)
    {
        std::string file_path("res/tr-" + std::to_string(99999) + ".txt");
        std::ofstream file(file_path);

        for (moth_mesh2d_triangle_iter cur = triangle_begin(),
                                       end = triangle_end(); cur != end; ++cur) {
#if 0
            if (pT->pP1 == pPoints[0]) continue;
            if (pT->pP1 == pPoints[1]) continue;
            if (pT->pP1 == pPoints[2]) continue;
            if (pT->pP2 == pPoints[0]) continue;
            if (pT->pP2 == pPoints[1]) continue;
            if (pT->pP2 == pPoints[2]) continue;
            if (pT->pP3 == pPoints[0]) continue;
            if (pT->pP3 == pPoints[1]) continue;
            if (pT->pP3 == pPoints[2]) continue;
#endif
            if (cur.ddel()) continue;

            file << *cur.point(1) << std::endl
                 << *cur.point(2) << std::endl
                 << *cur.point(3) << std::endl
                 << *cur.point(1) << std::endl
                 << std::endl;
            if (center) {
                file << (*cur.point(1) + *cur.point(2) + *cur.point(3)) / 3.0 << std::endl
                     << (*cur.point(1) + *cur.point(2) + *cur.point(3)) / 3.0 << std::endl
                     << std::endl;
            }
        }
    }

};  // struct moth_mesh2d

}
