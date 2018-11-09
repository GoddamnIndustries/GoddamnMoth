#pragma once

#include <vector>
#include <algorithm>
#include <fstream>

#include <cassert>
#include <type_traits>

#include "libGeometry2D/src/GeomBase.hh"
#include "libGeometry2D/src/GeomTriangle.hh"
#include "libGeometry2D/src/GeomSort.hh"

#undef assert
#define assert(...)

namespace DT {

struct moth_mesh2d_point;
struct moth_mesh2d_point_iter;
struct moth_mesh2d_triangle;
struct moth_mesh2d_triangle_iter;

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

extern int cnt;

/**
 * Triangle cell.
 */
struct MOTH_CORE moth_mesh2d_triangle
{
public:
    moth_size_t nP1{MOTH_NPOS}, nP2{MOTH_NPOS}, nP3{MOTH_NPOS};
    moth_size_t nT1{MOTH_NPOS}, nT2{MOTH_NPOS}, nT3{MOTH_NPOS};
    bool bad{};
    int mcnt{};

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
    int& vvis()
    {
        return pTriangles[nT].mcnt;
    }

    MOTH_HOST
    bool visited() const
    {
        return invalid() && (pTriangles[nT].mcnt >= cnt);
    }
    MOTH_HOST
    bool unvisited() const
    {
        return valid() && (pTriangles[nT].mcnt < cnt);
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
    static void swap_assign(moth_mesh2d_triangle_iter& pT1, const moth_mesh2d_triangle_iter& pT2)
    {
        swap(pT1, pT2);
        pT1 = pT2;
    }
    MOTH_HOST MOTH_DEVICE
    static void swap(const moth_mesh2d_triangle_iter& pT1, const moth_mesh2d_triangle_iter& pT2)
    {
        assert(pT1.nT != MOTH_NPOS && pT2.nT != MOTH_NPOS);
        assert(pT1.pPoints == pT2.pPoints &&
               pT1.pTriangles == pT2.pTriangles);
        if (pT1 != pT2) {
            moth_mesh2d_triangle& T1{pT1.pTriangles[pT1.nT]};
            moth_mesh2d_triangle& T2{pT2.pTriangles[pT2.nT]};

            /* @todo Rewrite this without using additional vector. */
            std::vector<moth_size_t*> pnT_neighbors;
            pnT_neighbors.reserve(9);

            /* Find neighbors of the first triangle. */
            for (moth_size_t k = 1; k <= 3; ++k) {
                if (T1.nnT(k) != MOTH_NPOS) {
                    moth_mesh2d_triangle& T1_k{pT1.pTriangles[T1.nnT(k)]};
                    for (moth_size_t m = 1; m <= 3; ++m) {
                        if (T1_k.nnT(m) == pT1.nT) {
                            pnT_neighbors.push_back(&T1_k.nnT(m));
                            break;
                        }
                    }
                }
            }
            /* Find neighbors of the second triangle. */
            for (moth_size_t k = 1; k <= 3; ++k) {
                if (T2.nnT(k) != MOTH_NPOS) {
                    moth_mesh2d_triangle& T2_k{pT2.pTriangles[T2.nnT(k)]};
                    for (moth_size_t m = 1; m <= 3; ++m) {
                        if (T2_k.nnT(m) == pT2.nT) {
                            pnT_neighbors.push_back(&T2_k.nnT(m));
                            break;
                        }
                    }
                }
            }

            /* Find only unique neighbors and relink them. */
            std::sort(pnT_neighbors.begin(), pnT_neighbors.end());
            pnT_neighbors.erase(std::unique(pnT_neighbors.begin(), pnT_neighbors.end()), pnT_neighbors.end());
            for (moth_size_t* pnT : pnT_neighbors) {
                if (*pnT == pT1.nT) {
                    *pnT = pT2.nT;
                } else if (*pnT == pT2.nT) {
                    *pnT = pT1.nT;
                } else {
                    std::cerr << "Bug in triangles swapping." << std::endl;
                    std::abort();
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
    void insert(const moth_p2d& p1, moth_real_t eps = 0.0);
    MOTH_HOST
    void insert(moth_p2d* pP_beg, moth_p2d* pP_end)
    {
        moth_sort(pP_beg, pP_end);
        for (moth_p2d* pP_cur = pP_beg; pP_cur != pP_end; ++pP_cur) {
            static int lll=0;
            if (lll++%10000==0) std::cerr << lll << std::endl;
            insert(*pP_cur);
        }
    }

public:
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
            //if (cur.vvis()) continue;

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
