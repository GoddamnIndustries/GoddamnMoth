#pragma once

#include "libGeometry2D/src/GeomTriangle.hh"

#include <vector>

struct MOTH_CORE moth_cell2d
{
    moth_size_t a{MOTH_NPOS}, b{MOTH_NPOS}, c{MOTH_NPOS};
    moth_size_t ab{MOTH_NPOS}, bc{MOTH_NPOS}, ca{MOTH_NPOS};
};  // moth_cell2d

MOTH_HOST MOTH_DEVICE MOTH_CORE
static std::vector<moth_cell2d> moth_triangulate(std::vector<moth_p2d> P)
{
    std::vector<moth_cell2d> Tr{};


    moth_tri2d Sx{{-100, -100}, {+100, -100}, {0, 200}};
    P.push_back(Sx.p1);
    P.push_back(Sx.p2);
    P.push_back(Sx.p3);

    moth_cell2d S{P.size() - 3, P.size() - 2, P.size() - 1};
    Tr.push_back(S);

    for (moth_size_t i = 0; i < P.size() - 3; ++i) {
        const moth_p2d& p = P[i];

        std::vector<moth_cell2d> BadTr{};
        for (const moth_cell2d& T : Tr) {
            moth_tri2d Tx{P[T.a], P[T.b], P[T.c]};
            if (moth_tri2d::circle(Tx, p)) {
                BadTr.push_back(T);
            }
        }
    }

    return Tr;
}
