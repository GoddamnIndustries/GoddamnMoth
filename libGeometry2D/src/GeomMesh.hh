#pragma once

#include <vector>
#include <algorithm>
#include <fstream>

#include "libGeometry2D/src/GeomBase.hh"
#include "libGeometry2D/src/GeomTriangle.hh"

struct moth_mesh2d;
struct moth_mesh2d_p;
struct moth_mesh2d_e;
struct moth_mesh2d_triangle;

struct MOTH_CORE moth_mesh2d_p final : public moth_p2d
{
    moth_mesh2d_triangle* pT{};
};  // struct moth_mesh2d_p

struct MOTH_CORE moth_mesh2d_e final
{
    moth_mesh2d_p* pP1{};
    moth_mesh2d_p* pP2{};
    moth_mesh2d_triangle* pT1{};
    moth_mesh2d_triangle* pT2{};
};  // struct moth_mesh2d_e

struct MOTH_CORE moth_mesh2d_triangle final
{
    moth_mesh2d_p* pP1{};
    moth_mesh2d_p* pP2{};
    moth_mesh2d_p* pP3{};
    moth_mesh2d_triangle* pT1{};
    moth_mesh2d_triangle* pT2{};
    moth_mesh2d_triangle* pT3{};
    bool bad{};
    bool visited{};

public:
    moth_triangle2d operator* () const
    {
        return {*pP1, *pP2, *pP3};
    }

public:
    void shift()
    {
        std::swap(pP1, pP2);
        std::swap(pP2, pP3);
        std::swap(pT1, pT2);
        std::swap(pT2, pT3);
    }
};  // struct moth_mesh2d_triangle

struct MOTH_CORE moth_mesh2d
{
    std::vector<moth_mesh2d_triangle*> pTriangles;
    std::vector<moth_mesh2d_e*> pEdges;
    std::vector<moth_mesh2d_p*> pPoints;

public:
    MOTH_HOST
    moth_mesh2d(moth_size_t capacity = 100000)
    {
        /* Create super triangle. */
        moth_mesh2d_triangle* pT = new moth_mesh2d_triangle{};
        pT->pP1 = new moth_mesh2d_p{{-4.0, -4.0}, pT};
        pT->pP2 = new moth_mesh2d_p{{+4.0, -4.0}, pT};
        pT->pP3 = new moth_mesh2d_p{{ 0.0, +4.0}, pT};

        /* Add super triangle. */
        pTriangles.reserve(capacity);
        pTriangles.push_back(pT);

        /* Add super triangle points. */
        pPoints.reserve(2 * capacity);
        pPoints.push_back(pT->pP1);
        pPoints.push_back(pT->pP2);
        pPoints.push_back(pT->pP3);

        /* Add super triangle edges. */
        pEdges.reserve(2 * capacity);
        pEdges.push_back(new moth_mesh2d_e{pT->pP2, pT->pP3, pT});
        pEdges.push_back(new moth_mesh2d_e{pT->pP3, pT->pP1, pT});
        pEdges.push_back(new moth_mesh2d_e{pT->pP1, pT->pP2, pT});
    }

public:
    MOTH_HOST
    void insert(const moth_p2d& p)
    {
        /* Insert the point. */
        pPoints.push_back(new moth_mesh2d_p{p});
        moth_mesh2d_p* pP = pPoints.back();

        /* Find bad triangles. */
        for (moth_mesh2d_triangle* pT : pTriangles) {
            pT->bad = moth_triangle2d::circle(**pT, *pP);
            pT->visited = false;
        }

        /* Find first border bad triangle and sort it. */
        moth_mesh2d_triangle* pT = nullptr;
        for (moth_size_t i = 0; i < pTriangles.size(); ++i) {
            pT = pTriangles[i];
            if (pT->bad) {
                if (pT->pT1 == nullptr || !pT->pT1->bad) {
                    break;
                } else if (pT->pT2 == nullptr || !pT->pT2->bad) {
                    pT->shift();
                    break;
                } else if (pT->pT3 == nullptr || !pT->pT3->bad) {
                    pT->shift();
                    pT->shift();
                    break;
                }
            }
        }
        if (pT == nullptr) {
            std::abort();
        }

        moth_mesh2d_p* pP_f = pT->pP2;
        moth_mesh2d_triangle* pT_f{};
        moth_mesh2d_triangle* pT_c{};
        moth_mesh2d_triangle* pT_p{};

        static int iii = 0;
        int kkk = 0;
        ++iii;

        do {
            if (pT->pT1 == nullptr || !pT->pT1->bad) {
                pT_p = pT_c;
                pT_c = new moth_mesh2d_triangle{pP, pT->pP2, pT->pP3};
                if (pT_f == nullptr) {
                    pT_f = pT_c;
                }

                /* Add new triangle and link it with outer
                 * neighbor and previous triangle. */
                pTriangles.push_back(pT_c);
                if (pT->pT1 != nullptr) {
                    if (pT->pT1->pT1 == pT) {
                        pT->pT1->pT1 = pT_c;
                    } else if (pT->pT1->pT2 == pT) {
                        pT->pT1->pT2 = pT_c;
                    } else if (pT->pT1->pT3 == pT) {
                        pT->pT1->pT3 = pT_c;
                    }
                    pT_c->pT1 = pT->pT1->pT1;
                }
                if (pT_p != nullptr) {
                    pT_p->pT2 = pT_c;
                    pT_c->pT3 = pT_p;
                }

                if (pT->pT2 != nullptr && pT->pT2->bad) {
                    /* CCW turn to the neighbor triangle. */
                    while (pT->pT2->pP1 != pT->pP1) {
                        pT->pT2->shift();
                    }
                    pT = pT->pT2;
                } else {
                    /* CCW turn inside current triangle. */
                    pT->shift();
                }
            } else {
                while (pT->pT1->pP1 != pT->pP3) {
                    pT->pT1->shift();
                }
                pT = pT->pT1;
            }

            ++kkk;
            if (iii == 3 && kkk == 3) break;
        } while (pT->pP2 != pP_f);

        /* Link first and last triangle. */
        pT_p = pT_c;
        pT_c = pT_f;
        pT_p->pT2 = pT_c;
        pT_c->pT3 = pT_p;

        /* Remove bad triangles. */
        pTriangles.erase(std::remove_if(pTriangles.begin(), pTriangles.end(), [](const moth_mesh2d_triangle* pT) {
            if (pT->bad) {
                delete pT;
                return true;
            } else {
                return false;
            }
        }), pTriangles.end());
    }

    void print()
    {
        std::string TrFilePath("res/tr-" + std::to_string(99999) + ".txt");
        std::ofstream TrFile(TrFilePath);

        for (moth_size_t i = 0; i < pTriangles.size(); ++i) {
            const moth_mesh2d_triangle* pT = pTriangles[i];


            TrFile << *pT->pP1 << std::endl
                   << *pT->pP2 << std::endl
                   << *pT->pP3 << std::endl
                   << *pT->pP1 << std::endl
                   << std::endl
                   << (*pT->pP1 + *pT->pP2 + *pT->pP3)/3 << std::endl
                   << (*pT->pP1 + *pT->pP2 + *pT->pP3)/3 << std::endl
                   << std::endl;
        }
    }
};  // struct moth_mesh2d



