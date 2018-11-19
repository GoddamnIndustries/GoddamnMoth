
#include "libGeometry2D/src/GeomTetrahedron.hh"
//#include "libGeometry2D/src/GeomMesh.hh"
#include "libGeometry2D/src/GeomMesh2.hh"
#include "libGeometry2D/src/GeomSort.hh"
#include "libGeometry2D/src/GeomPoly.hh"

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
    // Test Hilbert sort.
#if 0
    std::vector<moth_p2d> points;
    for (moth_real_t i = 0; i <= 31; ++i) {
        for (moth_real_t j = 0; j <= 31; ++j) {
            points.push_back({i, j});
        }
    }
    moth_sort(points.data(), points.data() + points.size());
#endif

    // Benchmark Hilbert sort.
#if 0
    std::default_random_engine random_engine;
    std::uniform_real_distribution<moth_real_t> uniform_distribution(-1.0, 1.0);

    std::vector<moth_p2d> points;
    points.reserve(100000000);

    for (moth_size_t m = 0; m < 10; ++m) {
        for (moth_size_t k = 0; k < 10000000; ++k) {
            moth_p2d p{uniform_distribution(random_engine),
                       uniform_distribution(random_engine)};
            points.push_back(p);
        }

        auto c = clock();
        moth_sort(points.data(), points.data() + points.size());
        c = clock() - c;
        std::cerr << m + 1 << " " << moth_real_t(c) / CLOCKS_PER_SEC << std::endl;
    }
#endif

    // Test incremental unconstrained triangulation.
#if 0
    std::default_random_engine random_engine;
    std::uniform_real_distribution<moth_real_t> uniform_distribution(-1.0, 1.0);

    moth_mesh2d builder;

    for (moth_size_t m = 0; m < 10; ++m) {
        auto c = clock();
        for (moth_size_t k = 0; k < 10000; ++k) {
            moth_p2d p{uniform_distribution(random_engine),
                       uniform_distribution(random_engine)};
            builder.insert(p);
        }

        c = clock() - c;
        std::cerr << m + 1 << " " << moth_real_t(c) / CLOCKS_PER_SEC << std::endl;
    }
    builder.print();
#endif

    // Test unconstrained triangulation.
#if 0
    std::default_random_engine random_engine;
    std::uniform_real_distribution<moth_real_t> uniform_distribution(-1.0, 1.0);

    moth_mesh2d builder;
    std::vector<moth_p2d> points;
    points.reserve(10000);

    for (moth_size_t m = 0; m < 1; ++m) {
        points.clear();
        for (moth_size_t k = 0; k < 10001; ++k) {
            moth_p2d p{uniform_distribution(random_engine),
                       uniform_distribution(random_engine)};
            points.push_back(p);
        }

        auto c = clock();
        builder.insert_unconstrained(points.data(), points.data() + points.size());
        c = clock() - c;
        std::cerr << m + 1 << " " << moth_real_t(c) / CLOCKS_PER_SEC << std::endl;
    }

#if 1
    builder.insert_unconstrained({-0.75, -0.75});
    builder.insert_unconstrained({+0.75, +0.75});
    builder.pConstraints.push_back({builder.pPoints.size() - 2, builder.pPoints.size() - 1});
    builder.apply_constrains_ignoring();
#endif
    builder.print();
#endif

    // Test constrained triangulation.
#if 1
    moth_poly2d rect{geom_poly2d_primitives::rect({0.0, 0.0}, {4.0, 4.0})};
    moth_poly2d crcl{geom_poly2d_primitives::circle({4.0/3.0, 2.0}, 4.0/5.0, 25)};
    //moth_poly2d crcl{geom_poly2d_primitives::star({4.0/3.0, 1.0}, 2.0/5.0, 1.5/5.0, 25)};

    moth_mesh2d builder;
    builder.insert_constrain(rect);
    builder.insert_constrain(crcl);

    //builder.peel();
    builder.refine();
    builder.print();

    //for (int m = 0; m < 100; ++m)
#endif

#if 0
    for (moth_size_t i = 0; i < 20; ++i) {
        std::cerr << i << std::endl;
        moth_real_t r = 1.0;
        moth_real_t phi = 2.0 * i / 20.0 * MOTH_PI;
        moth_real_t x = r * cos(phi);
        moth_real_t y = r * sin(phi);
        builder.insert({x, y});
    }
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
