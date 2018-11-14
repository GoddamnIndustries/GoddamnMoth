#pragma once

#include <vector>
#include <algorithm>
#include <fstream>

#include <cassert>
#include <type_traits>

#include "libGeometry2D/src/GeomBase.hh"
#include "libGeometry2D/src/GeomTriangle.hh"
#include "libGeometry2D/src/GeomPoly.hh"

//#undef assert
//#define assert(...)

struct moth_mesh2d_point;
struct moth_mesh2d_point_iter;
struct moth_mesh2d_triangle_t;
struct moth_mesh2d_triangle_iter;
struct moth_mesh2d_constraint_t;
struct moth_mesh2d_cedge_iter;

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Mesh Triangle cell.
 */
struct MOTH_CORE moth_mesh2d_triangle_t
{
public:
    moth_size_t nP1{MOTH_NPOS}, nP2{MOTH_NPOS}, nP3{MOTH_NPOS};
    moth_size_t nT1{MOTH_NPOS}, nT2{MOTH_NPOS}, nT3{MOTH_NPOS};
    moth_size_t nCounter{};
    bool bad{};

public:
    MOTH_HOST
    moth_size_t& nP(size_t k)
    {
        return *(&nP1 + (k - 1) % 3);
    }
    MOTH_HOST
    moth_size_t& nT(size_t k)
    {
        return *(&nT1 + (k - 1) % 3);
    }
};  // struct moth_mesh2d_triangle_t

struct MOTH_CORE moth_mesh2d_constraint_t
{
    moth_size_t nP1{MOTH_NPOS}, nP2{MOTH_NPOS};
    moth_size_t nE1{MOTH_NPOS}, nE2{MOTH_NPOS};

public:
    MOTH_HOST
    moth_size_t& nP(size_t k)
    {
        return *(&nP1 + (k - 1) % 2);
    }
    MOTH_HOST
    moth_size_t& nE(size_t k)
    {
        return *(&nE1 + (k - 1) % 2);
    }
};  // struct moth_mesh2d_constraint_t

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * 2D Delaunay triangular mesh.
 */
struct MOTH_CORE moth_mesh2d
{
    std::vector<moth_p2d> pPoints;
    std::vector<moth_mesh2d_triangle_t> pTriangles;
    std::vector<moth_mesh2d_constraint_t> pConstraints;

public:
    MOTH_HOST
    moth_mesh2d(moth_size_t capacity = 100000000);

public:
    MOTH_HOST template<typename T = moth_mesh2d_triangle_iter>
    T triangle_begin()
    {
        return {*this, 0};
    }
    MOTH_HOST template<typename T = moth_mesh2d_triangle_iter>
    T triangle_end()
    {
        return {*this, pTriangles.size()};
    }

public:
    MOTH_HOST template<typename T = moth_mesh2d_cedge_iter>
    T constraint_begin()
    {
        return {*this, 0};
    }
    MOTH_HOST template<typename T = moth_mesh2d_cedge_iter>
    T constraint_end()
    {
        return {*this, pConstraints.size()};
    }

public:
    /**
     * @brief Insert a single point into the mesh ignoring \
     * previously applied constrains.
     *
     * Insertion is implemented using the Bowyer-Watson algorithm,
     * complexity in average is @f$ O(log(n)) @f$.
     * Performance can be improved significantly by inserting
     * multiple points at once.
     *
     * @returns Iterator to the inserted point.
     *
     * @see https://en.wikipedia.org/wiki/Bowyer–Watson_algorithm
     */
    MOTH_HOST
    moth_mesh2d_point_iter insert_unconstrained(const moth_p2d& p1, moth_real_t eps = 0.0001);

    /**
     * @brief Insert multiple points into the mesh ignoring \
     * previously applied constrains.
     *
     * Insertion is implemented using the Bowyer-Watson algorithm,
     * presorting points along a Hilbert curve,
     * complexity in average is @f$ O(n) @f$.
     *
     * @returns Iterator to the first inserted point.
     *
     * @see https://en.wikipedia.org/wiki/Bowyer–Watson_algorithm
     * @see https://en.wikipedia.org/wiki/Hilbert_curve
     */
    MOTH_HOST
    moth_mesh2d_point_iter insert_unconstrained(moth_p2d* pP_beg, moth_p2d* pP_end);

public:
    MOTH_HOST
    moth_mesh2d_cedge_iter apply_constrain_ignoring(const moth_mesh2d_cedge_iter& pE);

    MOTH_HOST
    void apply_constrains_ignoring();

public:
    MOTH_HOST
    void apply_constrain_conforming(const moth_mesh2d_cedge_iter& pE);

    MOTH_HOST
    void apply_constrains_conforming();

public:
    /**
     * @brief Insert points of a polygon into the mesh,
     *
     */
    MOTH_HOST
    moth_mesh2d_point_iter insert_constrain(const moth_poly2d& poly);

public:
    template<typename T = moth_mesh2d_triangle_iter>
    void print(bool center = false)
    {
        std::string file_path("res/tr-" + std::to_string(99999) + ".txt");
        std::ofstream file(file_path);

        static int k=0;

        for (T cur = triangle_begin<T>(),
               end = triangle_end<T>(); cur != end; ++cur) {
            if (k++%1000==0) std::cerr << k << std::endl;
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

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Triangulation points iterator.
 */
struct MOTH_CORE moth_mesh2d_point_iter
{
    moth_mesh2d& pMesh;
    moth_size_t nP{MOTH_NPOS};

public:
    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_point_iter& operator=(const moth_mesh2d_point_iter& pT)
    {
        nP = pT.nP;
        return *this;
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_p2d operator* () const
    {
        assert(nP != MOTH_NPOS);
        return pMesh.pPoints[nP];
    }

public:
    MOTH_HOST MOTH_DEVICE
    bool operator==(const moth_mesh2d_point_iter& pP) const
    {
        assert(&pP.pMesh == &pMesh);
        return nP == pP.nP;
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_mesh2d_point_iter& pP) const
    {
        assert(&pP.pMesh == &pMesh);
        return nP != pP.nP;
    }

public:
    template<typename T=moth_mesh2d_triangle_iter>
    T triangle() const
    {
        for (T pT_cur = pMesh.triangle_begin<T>(), pT_end = pMesh.triangle_end<T>();
               pT_cur != pT_end; ++pT_cur)
        {
            if (pT_cur.point(1) == *this) return pT_cur;
            if (pT_cur.point(2) == *this) return pT_cur;
            if (pT_cur.point(3) == *this) return pT_cur;
        }
        abort();
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
};  // struct moth_mesh2d_point_iter

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Triangle cells iterator.
 */
struct MOTH_CORE moth_mesh2d_triangle_iter
{
    moth_mesh2d& pMesh;
    moth_size_t nT{MOTH_NPOS};

public:
    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_point_iter point(moth_size_t k) const
    {
        assert(nT < pMesh.pTriangles.size());
        return {pMesh, pMesh.pTriangles[nT].nP(k)};
    }
    MOTH_HOST
    void set_point(moth_size_t k,
                   const moth_mesh2d_point_iter& pP)
    {
        assert(&pP.pMesh == &pMesh && "Incompatible iterators.");
        assert(nT < pMesh.pTriangles.size());
        pMesh.pTriangles[nT].nP(k) = pP.nP;
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_triangle_iter triangle(moth_size_t k) const
    {
        assert(nT < pMesh.pTriangles.size());
        return {pMesh, pMesh.pTriangles[nT].nT(k)};
    }
    MOTH_HOST
    void set_triangle(moth_size_t k,
                      const moth_mesh2d_triangle_iter& pT)
    {
        assert(&pT.pMesh == &pMesh && "Incompatible iterators.");
        assert(nT < pMesh.pTriangles.size());
        pMesh.pTriangles[nT].nT(k) = pT.nT;
    }

public:
    MOTH_HOST MOTH_DEVICE
    const moth_triangle2d operator*() const
    {
        return {*point(1),
                *point(2),
                *point(3)};
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_triangle_iter& operator=(const moth_mesh2d_triangle_iter& pT)
    {
        assert(&pT.pMesh == &pMesh);
        nT = pT.nT;
        return *this;
    }

public:
    MOTH_HOST MOTH_DEVICE
    bool operator==(const moth_mesh2d_triangle_iter& pT) const
    {
        assert(&pT.pMesh == &pMesh);
        return nT == pT.nT;
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(const moth_mesh2d_triangle_iter& pT) const
    {
        assert(&pT.pMesh == &pMesh);
        return nT != pT.nT;
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_triangle_iter& operator++()
    {
        assert(nT != MOTH_NPOS);
        ++nT;
        return *this;
    }
    MOTH_HOST MOTH_DEVICE
    const moth_mesh2d_triangle_iter operator++(int)
    {
        moth_mesh2d_triangle_iter copy{*this};
        ++*this;
        return copy;
    }

    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_triangle_iter& operator--()
    {
        assert(nT != MOTH_NPOS);
        --nT;
        return *this;
    }
    MOTH_HOST MOTH_DEVICE
    const moth_mesh2d_triangle_iter operator--(int)
    {
        moth_mesh2d_triangle_iter copy{*this};
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
        moth_mesh2d_triangle_iter copy{*this};
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
        moth_mesh2d_triangle_iter copy{*this};
        return copy -= d;
    }

public:
    MOTH_HOST
    static void lshift(const moth_mesh2d_triangle_iter& pT)
    {
        assert(pT.nT < pT.pMesh.pTriangles.size());
        moth_mesh2d_triangle_t& T{pT.pMesh.pTriangles[pT.nT]};
        std::swap(T.nP1, T.nP2);
        std::swap(T.nP2, T.nP3);
        std::swap(T.nT1, T.nT2);
        std::swap(T.nT2, T.nT3);
    }
    MOTH_HOST
    static void rshift(const moth_mesh2d_triangle_iter& pT)
    {
        assert(pT.nT < pT.pMesh.pTriangles.size());
        moth_mesh2d_triangle_t& T{pT.pMesh.pTriangles[pT.nT]};
        std::swap(T.nP1, T.nP3);
        std::swap(T.nP2, T.nP3);
        std::swap(T.nT1, T.nT3);
        std::swap(T.nT2, T.nT3);
    }

public:
    /**
     * Safely swap two triangles. */
    MOTH_HOST
    static void swap(const moth_mesh2d_triangle_iter& pT1,
                     const moth_mesh2d_triangle_iter& pT2)
    {
        assert(pT1.nT != MOTH_NPOS && pT2.nT != MOTH_NPOS);
        assert(&pT1.pMesh == &pT2.pMesh && "Incompatible iterators.");
        if (pT1 != pT2) {
            moth_size_t* pnT_neighbors[6]{};
            moth_size_t** pnT_neighbors_end = pnT_neighbors;

            /* Find neighbors of the first triangle. */
            moth_mesh2d_triangle_t& T1{pT1.pMesh.pTriangles[pT1.nT]};
            for (moth_size_t k = 1; k <= 3; ++k) {
                if (T1.nT(k) != MOTH_NPOS) {
                    moth_mesh2d_triangle_t& T1_k{pT1.pMesh.pTriangles[T1.nT(k)]};
                    for (moth_size_t m = 1; m <= 3; ++m) {
                        if (T1_k.nT(m) == pT1.nT) {
                            *(pnT_neighbors_end++) = &T1_k.nT(m);
                            break;
                        }
                    }
                }
            }
            /* Find neighbors of the second triangle. */
            moth_mesh2d_triangle_t& T2{pT2.pMesh.pTriangles[pT2.nT]};
            for (moth_size_t k = 1; k <= 3; ++k) {
                if (T2.nT(k) != MOTH_NPOS) {
                    moth_mesh2d_triangle_t& T2_k{pT2.pMesh.pTriangles[T2.nT(k)]};
                    for (moth_size_t m = 1; m <= 3; ++m) {
                        if (T2_k.nT(m) == pT2.nT) {
                            *(pnT_neighbors_end++) = &T2_k.nT(m);
                            break;
                        }
                    }
                }
            }

            /* Relink the neighbors. */
            for (moth_size_t** pnT_neighbors_cur{pnT_neighbors};
                               pnT_neighbors_cur != pnT_neighbors_end; ++pnT_neighbors_cur) {
                moth_size_t* pnT = *pnT_neighbors_cur;
                if (*pnT == pT1.nT) {
                    *pnT = pT2.nT;
                } else {
                    *pnT = pT1.nT;
                }
            }

            /* And finally swap the memory. */
            std::swap(pT1.pMesh.pTriangles[pT1.nT],
                      pT2.pMesh.pTriangles[pT2.nT]);
        }
    }

    MOTH_HOST MOTH_DEVICE
    static void swap_assign(moth_mesh2d_triangle_iter& pT1, const moth_mesh2d_triangle_iter& pT2)
    {
        swap(pT1, pT2);
        pT1 = pT2;
    }

public:
    /** Check if iterator is valid. */
    /** @{ */
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
    /** @} */

public:
    /**
     * Check if triangle is bad or good.
     * Invalid iterators implicitly are treated as good. */
    /** @{ */
    MOTH_HOST
    bool bad() const
    {
        return valid() && pMesh.pTriangles[nT].bad;
    }
    MOTH_HOST
    bool good() const
    {
        return invalid() || !pMesh.pTriangles[nT].bad;
    }
    /** @} */

    /** Mark the triangle bad. */
    MOTH_HOST
    void set_bad(bool b)
    {
        assert(nT < pMesh.pTriangles.size());
        pMesh.pTriangles[nT].bad = b;
    }

public:
    /**
     * Check if the triangle was not visited
     * during the current insertion and marks if visited.
     * Invalid iterators are implicitly treated as visited. */
    MOTH_HOST
    bool not_visited()
    {
        moth_size_t nCounter{pMesh.pPoints.size()};
        if (valid() && (pMesh.pTriangles[nT].nCounter < nCounter)) {
            pMesh.pTriangles[nT].nCounter = nCounter;
            return true;
        } else {
            return false;
        }
    }
};  // struct moth_mesh2d_triangle_iter

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Triangles around point iterator.
 */
struct MOTH_CORE moth_mesh2d_triangle_around_point_iter final
{
    moth_mesh2d_point_iter pP;
    moth_mesh2d_triangle_iter pT;

public:
    MOTH_HOST MOTH_DEVICE
    explicit moth_mesh2d_triangle_around_point_iter(
        const moth_mesh2d_point_iter& pP_): pP{pP_}, pT{pP_.triangle()}
    {
    }

public:
    MOTH_HOST MOTH_DEVICE
    operator moth_mesh2d_triangle_iter() const
    {
        return pT;
    }

public:
    /** Compare the iterators. */
    /** @{ */
    MOTH_HOST MOTH_DEVICE
    bool operator==(
        const moth_mesh2d_triangle_around_point_iter& pPT) const
    {
        return (pP == pPT.pP) && (pT == pPT.pT);
    }
    MOTH_HOST MOTH_DEVICE
    bool operator!=(
        const moth_mesh2d_triangle_around_point_iter& pPT) const
    {
        return (pP != pPT.pP) || (pT != pPT.pT);
    }
    /** @} */

public:
    /**
     * Increment the iterator: go to the next CCW triangle,
     * jumping over the borders. */
    /** @{ */
    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_triangle_around_point_iter& operator++()
    {
        /* @todo Implement without shift. */
        if (pT.triangle(2).valid()) {
            while (pT.triangle(2).point(1) != pP) {
                moth_mesh2d_triangle_iter::lshift(pT.triangle(2));
            }
            pT = pT.triangle(2);
        } else {
            while (pT.triangle(3).valid()) {
                --*this;
            }
        }
        return *this;
    }
    MOTH_HOST MOTH_DEVICE
    const moth_mesh2d_triangle_around_point_iter operator++(int)
    {
        moth_mesh2d_triangle_around_point_iter copy{*this};
        ++*this;
        return copy;
    }
    /** @} */

    /**
     * Decrement the iterator: go to the next CW triangle,
     * jumping over the borders. */
    /** @{ */
    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_triangle_around_point_iter& operator--()
    {
        /* @todo Implement without shift. */
        if (pT.triangle(3).valid()) {
            while (pT.triangle(3).point(1) != pP) {
                moth_mesh2d_triangle_iter::lshift(pT.triangle(3));
            }
            pT = pT.triangle(3);
        } else {
            while (pT.triangle(2).valid()) {
                ++*this;
            }
        }
        return *this;
    }
    MOTH_HOST MOTH_DEVICE
    const moth_mesh2d_triangle_around_point_iter operator--(int)
    {
        moth_mesh2d_triangle_around_point_iter copy{*this};
        ++*this;
        return copy;
    }
    /** @} */
};  // struct moth_mesh2d_triangle_around_point_iter

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

struct moth_mesh2d_cedge_iter
{
    moth_mesh2d& pMesh;
    moth_size_t nE{MOTH_NPOS};

    bool valid() const { return nE != MOTH_NPOS;}

public:
    moth_mesh2d_point_iter point(moth_size_t k) const
    {
        return {pMesh, pMesh.pConstraints[nE].nP(k)};
    }
    void set_point(moth_size_t k, moth_mesh2d_point_iter pP)
    {
        pMesh.pConstraints[nE].nP(k) = pP.nP;
    }

public:
    moth_mesh2d_cedge_iter edge(moth_size_t k) const
    {
        return {pMesh, pMesh.pConstraints[nE].nE(k)};
    }
    void set_edge(moth_size_t k, moth_mesh2d_cedge_iter pP)
    {
        pMesh.pConstraints[nE].nE(k) = pP.nE;
    }

public:
    moth_e2d operator*() const
    {
        return {*point(1), *point(2)};
    }

public:
    MOTH_HOST MOTH_DEVICE
    moth_mesh2d_cedge_iter& operator=(const moth_mesh2d_cedge_iter& pT)
    {
        nE = pT.nE;
        return *this;
    }

public:
    bool operator==(const moth_mesh2d_cedge_iter& pE) const
    {
        return nE == pE.nE;
    }
    bool operator!=(const moth_mesh2d_cedge_iter& pE) const
    {
        return nE != pE.nE;
    }

public:
    moth_mesh2d_cedge_iter& operator++()
    {
        ++nE;
        return *this;
    }
    moth_mesh2d_cedge_iter operator++(int);

    moth_mesh2d_cedge_iter& operator--()
    {
        --nE;
        return *this;
    }
};  // struct moth_mesh2d_constraint_iter
