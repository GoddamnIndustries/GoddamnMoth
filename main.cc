
#include "libGeometry2D/src/GeomEdgeList.hh"
#include <random>
#include <ctime>
#include <fstream>

struct StatsV
{
    float u1{}, v1{};
    float u2{}, v2{};
};

int main()
{
    geom_e2d_list* domain1 = geom_e2d_list_factory::new_rect_ccw({0.0, 0.0}, {10.0, 10.0});
    geom_e2d_list* domain2 = geom_e2d_list_factory::new_circle_ccw({5, 5}, 0.5);
    domain2 = geom_e2d_list_factory::copy_rev(domain2);

    std::vector<std::vector<StatsV>> estimator{150, std::vector<StatsV>{150}};

    std::default_random_engine gen_random;

    size_t m = 1;
    for (; m < 1000000; ++m) {

        auto cs = clock();

        std::vector<geom_p2d> particles;
        for (int i = 0; i < 70000; ++i) {
            while (true) {
                geom_p2d p;
                p.x = (geom_real_t) gen_random() / gen_random.max() * 10.0;
                p.y = (geom_real_t) gen_random() / gen_random.max() * 10.0;
                if (geom_e2d_list::contains(domain2, p)) continue;
                particles.push_back(p);
                break;
            }
        }

        geom_e2d_list* t = nullptr;
        for (int k = 0; k < particles.size(); ++k) {
            auto& p = particles[k];
            geom_real_t t0 = 0.0, t1 = 20.0;
            geom_size_t n = 1;
            {
                geom_p2d v;
                v.x = std::normal_distribution<double>(0.2, 0.09)(gen_random);
                v.y = std::normal_distribution<double>(0.0, 0.09)(gen_random);
                for (geom_size_t i = 0; i < n; ++i) {
                    geom_real_t dt = (t1 - t0) / n;
                    geom_p2d p1 = p + v * dt, p0 = p;
                    if (p.x > 9.5) {
                        p.x = 0.1;
                    }

                    geom_e2d e2;
                    if (geom_e2d_list::reflect(domain1, {p, p1}, e2)) {
                        p = e2.t;
                    } else {
                        if (geom_e2d_list::reflect(domain2, {p, p1}, e2)) {
                            p = e2.t;
                        } else {
                            p = p1;
                        }
                    }

                    v = (p - p0) / dt;

                    //v = geom_p2d::rotate(v, (geom_real_t) gen_random() / gen_random.max() * GEOM_PI * 2);
                    p.u = v.x;
                    p.v = v.y;
                    //geom_e2d_list::push(t, p);
                }
            }
            //geom_e2d_list::push(t, p);

            size_t ii = (size_t)round(p.x * (estimator.size() - 1) / 10);
            size_t jj = (size_t)round(p.y * (estimator[0].size() - 1) / 10);
            if (ii >= 0 && jj >= 0) {
                estimator[ii][jj].u1 += p.u;
                estimator[ii][jj].v1 += p.v;
                estimator[ii][jj].u2 += p.u * p.u;
                estimator[ii][jj].v2 += p.v * p.v;
            }
        }

        double max_dispersion = 0.0;
        for (auto ii = 1; ii < estimator.size(); ++ii) {
            for (auto jj = 1; jj < estimator[ii].size(); ++jj) {
                double expected_u = estimator[ii][jj].u1 / m;
                double expected_v = estimator[ii][jj].v1 / m;
                double dispersion_x = (estimator[ii][jj].u2 / m - expected_u * expected_u) / m;
                double dispersion_y = (estimator[ii][jj].v2 / m - expected_v * expected_v) / m;
                double dispersion = std::max(dispersion_x, dispersion_y);
                max_dispersion = std::max(max_dispersion, dispersion);
            }
        }

        cs = clock() - cs;
        std::cerr << "Iter took: " << cs / (float) CLOCKS_PER_SEC
                  << ", D = " << max_dispersion
                  << std::endl;

        if (max_dispersion < 0.0001) {
            break;
        }
    }

#if _WIN32
    std::ofstream output("H:\\GoddamnOpuwenijSolver\\data.dat");
#else
    std::ofstream output("data1.dat");
#endif
    for (auto ii = 0; ii < estimator.size(); ++ii) {
        for (auto jj = 0; jj < estimator[ii].size(); ++jj) {
            if (ii == 0 && jj == 0) continue;
            double expected_u = estimator[ii][jj].u1 / m;
            double expected_v = estimator[ii][jj].v1 / m;
            output << 10 * ii / float(estimator.size()) << " " << 10 * jj / float(estimator[ii].size())<< " "
            << expected_u << " " << expected_v << std::endl;
        }
    }

    //geom_e2d_list::plt(t);

    delete domain1;
    return 0;
}
