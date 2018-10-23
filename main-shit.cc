#include <cstdlib>
#include <ctime>

#include <algorithm>
#include <fstream>
#include <vector>
#include <random>


using phys_real_t = double;
using phys_size_t = size_t;

phys_real_t mass = 1.0;
phys_real_t bolts = 1.0;

struct phys_particle
{
    phys_real_t x{};
    phys_real_t v{};
    phys_real_t a{};
};  // struct phys_particle

struct phys_sample
{
};  // struct phys_sample1d

struct phys_grid_1d
{
    std::vector<phys_real_t> p_x;
    std::vector<phys_real_t> p_v;
    std::vector<phys_real_t> p_d;

public:
    void insert(phys_real_t x, phys_real_t v, phys_real_t rho)
    {
        p_x.push_back(x);
        p_v.push_back(v);
        p_d.push_back(rho);
    }
    void sort()
    {
        /* Generate the permutation. */
        std::vector<phys_size_t> p_x_perm(p_x.size());
        std::iota(p_x_perm.begin(), p_x_perm.end(), 0);
        std::sort(p_x_perm.begin(), p_x_perm.end(), [&](phys_size_t i, phys_size_t j){
            return p_x[i] < p_x[j];
        });

        std::vector<phys_real_t> p_x1(p_x);
        std::transform(p_x_perm.begin(), p_x_perm.end(), p_x.begin(), [&](phys_size_t i) {
            return p_x1[i];
        });

        std::vector<phys_real_t> p_v1(p_v);
        std::transform(p_x_perm.begin(), p_x_perm.end(), p_v.begin(), [&](phys_size_t i) {
            return p_v1[i];
        });

        std::vector<phys_real_t> p_d1(p_d);
        std::transform(p_x_perm.begin(), p_x_perm.end(), p_d.begin(), [&](phys_size_t i) {
            return p_d1[i];
        });
    }

public:
    phys_real_t V(phys_real_t x) const
    {
        for (phys_size_t i = 1; i < p_x.size(); ++i) {
            phys_real_t x0 = p_x[i - 1];
            phys_real_t x1 = p_x[i - 0];
            if (x0 < x && x <= x1) {
               phys_real_t delta = (x - x0) / (x1 - x0);
                phys_real_t v0 = p_v[i - 1];
                phys_real_t v1 = p_v[i - 0];
                return v0 * (1.0 - delta) + v1 * delta;
            }
        }
        return 0.0;
    }
    phys_real_t E(phys_real_t x) const
    {
        phys_real_t vel = V(x);
        return mass * vel * vel * 0.5;
    }
    phys_real_t T(phys_real_t x) const
    {
        phys_real_t enr = E(x);
        return enr / bolts / 3.0;
    }

    phys_real_t Rho(phys_real_t x) const
    {
        for (phys_size_t i = 1; i < p_x.size(); ++i) {
            phys_real_t x0 = p_x[i - 1];
            phys_real_t x1 = p_x[i - 0];
            if (x0 < x && x <= x1) {
                phys_real_t delta = (x - x0) / (x1 - x0);
                phys_real_t rho0 = p_d[i - 1];
                phys_real_t rho1 = p_d[i - 0];
                return rho0 * (1.0 - delta) + rho1 * delta;
            }
        }
        return 0.0;
    }
    phys_real_t P() const
    {
    }
};  // struct phys_grid_1d

int main()
{
    phys_real_t len = 10.0;
    phys_real_t x_0 = 5.0;
    phys_real_t rho_0 = 2.0, rho_1 = 0.5;
    phys_grid_1d domain;
    domain.insert(0.0, 0.0, rho_0);
    domain.insert(x_0, 0.0, rho_0);
    domain.insert(x_0, 0.0, rho_1);
    domain.insert(len, 0.0, rho_1);

    phys_real_t tau = 0.1;

    phys_grid_1d domain1;
    std::default_random_engine rand;
    for (phys_size_t n = 0; n < 10000; ++n) {
        std::piecewise_linear_distribution<phys_real_t> pld(domain.p_x.begin(), domain.p_x.end(), domain.p_d.begin());

        phys_particle p;
        p.x = pld(rand);
        p.v = domain.V(p.x);
        p.a = -(domain.Rho(p.x + 0.01) - domain.Rho(p.x - 0.01)) / 0.02 / domain.Rho(p.x);

        std::uniform_real_distribution<phys_real_t> ud(p.v, p.v * p.v / 3);
        p.v = ud(rand);
        p.v += tau * tau * p.a * 0.5;
        p.x += tau * p.v;
        domain1.insert(p.x, 0.0, p.v);

    }
    domain1.sort();

    std::ofstream file("test_dist.txt");
    for (phys_size_t n = 0; n < domain1.p_x.size(); ++n)
    {
        phys_size_t num = 1;
        phys_size_t m0 = (phys_size_t) std::max<int>(0, n - 10);
        for (phys_size_t m = n; m != m0; --m) {
            phys_real_t x0 = domain1.p_x[n];
            phys_real_t x1 = domain1.p_x[m];
            if ((x0 - x1) < 0.01) {
                ++num;
            } else {
                break;
            }
        }

        phys_size_t m1 = (phys_size_t) std::min<int>(domain1.p_x.size(), n + 10);
        for (phys_size_t m = n; m != m1; ++m) {
            phys_real_t x0 = domain1.p_x[n];
            phys_real_t x1 = domain1.p_x[m];
            if ((x1 - x0) < 0.01) {
                ++num;
            } else {
                break;
            }
        }

        domain1.p_d[n] = 1.0 * num * mass / 0.2;
        file << domain1.p_x[n] << " "
             << domain1.p_v[n] << " "
             << domain1.p_d[n] << " " << std::endl;
    }

    return 0;
}
