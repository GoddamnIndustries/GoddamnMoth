
#include "libGeometry2D/src/GeomTetrahedron.hh"
#include "libGeometry2D/src/GeomDelaunay.hh"
#include "libGeometry2D/src/GeomTriangulation.hh"
#include "libGeometry2D/src/GeomMesh.hh"

moth_real_t f(const moth_p3d& x)
{
    return x.x*std::exp((x.x + x.y + x.z) / 0.4 * MOTH_PI_2);
}

moth_p3d df(const moth_p3d& x)
{
    double d = std::exp((x.x + x.y + x.z) / 0.4 * MOTH_PI) / 0.4 * MOTH_PI;
    double d1 = x.x * d;
    return {d1 + d, d1, d1};
}

moth_p3d moth_grad(const moth_tetrahedron& T, moth_real_t f1, moth_real_t f2,
                                              moth_real_t f3, moth_real_t f4)
{
    moth_tri3d t234{T.p2, T.p3, T.p4};
    moth_tri3d t143{T.p1, T.p4, T.p3};
    moth_tri3d t124{T.p1, T.p2, T.p4};
    moth_tri3d t132{T.p1, T.p3, T.p2};

    moth_p3d grad{};
    grad += (f2 + f3 + f4) / 3.0 * moth_tri3d::normal(t234) * moth_tri3d::area(t234);
    grad += (f1 + f4 + f3) / 3.0 * moth_tri3d::normal(t143) * moth_tri3d::area(t143);
    grad += (f1 + f2 + f4) / 3.0 * moth_tri3d::normal(t124) * moth_tri3d::area(t124);
    grad += (f1 + f3 + f2) / 3.0 * moth_tri3d::normal(t132) * moth_tri3d::area(t132);
    grad /= moth_tetrahedron::volume(T);

    return grad;
}

#include <ctime>
#include <random>

int main()
{
#if 1
    std::default_random_engine random_engine(228);
    std::uniform_real_distribution<moth_real_t> uniform_distribution(-1.0, 1.0);

    moth_mesh2d builder;

    auto c = clock();
    for (moth_size_t k = 0; k < 3; ++k) {
        moth_p2d p{uniform_distribution(random_engine),
                   uniform_distribution(random_engine)};
        builder.insert(p);
    }
    c = clock() - c;
    std::cerr << moth_real_t(c) / CLOCKS_PER_SEC << std::endl;
    builder.print();
#endif

#if 0
    moth_tri_grid2d_builder builder;
    std::vector<moth_p2d> P;
    for (moth_size_t i = 0; i < 40; ++i) {
        moth_real_t r = 5.0;
        moth_real_t phi = 2.0 * i / 40.0 * MOTH_PI;
        moth_real_t x = r * cos(phi);
        moth_real_t y = r * sin(phi) + x;
        builder.insert({x, y});
        P.push_back({x , y});
    }
    std::vector<moth_e2d> E;
    for (moth_size_t i = 0; i < P.size() - 1; ++i) {
        E.push_back({P[i], P[i + 1]});
    }
    E.push_back({P.back(), P.front()});
    builder.tri_edges = E;
    builder.refine();
    builder.print();
#endif

#if 0
    std::vector<moth_p2d> P;
    /*for (moth_size_t i = 1; i <= 5; ++i) {
        moth_real_t r = i / 5.0;
        for (moth_size_t j = 0; j < 20; ++j) {
            moth_real_t phi = 2.0 * j / 20.0 * MOTH_PI;
            moth_real_t x = r * cos(phi);
            moth_real_t y = r * sin(phi);
            P.push_back({x , y});
        }
    }*/
    for (moth_size_t i = 0; i < 10; ++i) {
        moth_real_t r = 5.0;
        moth_real_t phi = 2.0 * i / 10.0 * MOTH_PI;
        moth_real_t x = r * cos(phi);
        moth_real_t y = r * sin(phi);
        P.push_back({x , y});
    }
    //P.push_back({0,0});
    //moth_triangulate_bowyer_watson(P);

#if 0
    std::vector<moth_e2d> E;
    for (moth_size_t i = 0; i < P.size() - 1; ++i) {
        E.push_back({P[i], P[i + 1]});
    }
    E.push_back({P.back(), P.front()});

    moth_triangulate_ruppert(E);
#endif
#endif

    moth_tetrahedron T{
        {0.001, 0.0013, 0.001}, {-0.001, -0.001, 0.001}, {-0.001, 0.001, -0.0017}, {0.001, -0.0016, -0.001}
    };

    moth_real_t f1 = f(T.p1);
    moth_real_t f2 = f(T.p2);
    moth_real_t f3 = f(T.p3);
    moth_real_t f4 = f(T.p4);

    moth_p3d num_grad{moth_grad(T, f1, f2, f3, f4)};
    moth_p3d ext_grad{df(0.25*(T.p1 + T.p2 + T.p3 + T.p4))};

    return 0;
}
