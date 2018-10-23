#include <vector>
#include <random>

#include <algorithm>
#include <iostream>
#include <fstream>

using moth_real_t = double;
using moth_size_t = size_t;

using moth_sample_t = std::vector<moth_real_t>;

struct moth_stat_sampler
{
private:
    std::vector<moth_real_t> m_dist_x;
    std::vector<moth_real_t> m_dist_y;
    mutable std::piecewise_linear_distribution<moth_real_t> m_dist;
    mutable std::default_random_engine m_rand;

public:
    moth_stat_sampler(const std::vector<moth_real_t>& dist_x, const std::vector<moth_real_t>& dist_y):
        m_dist_x(dist_x), m_dist_y(dist_y),
        m_dist(m_dist_x.begin(), m_dist_x.end(), m_dist_y.begin())
    {}

public:
    moth_real_t min_x() const
    {
        return m_dist_x.front();
    }
    moth_real_t max_x() const
    {
        return m_dist_x.back();
    }

    moth_sample_t gen_sample(moth_size_t n = 1000) const
    {
        moth_sample_t sample(n);
        for (moth_real_t& x : sample) {
            x = m_dist(m_rand);
        }
        std::sort(sample.begin(), sample.end());
        return sample;
    }
};

struct moth_stat_pdf_estimator
{
public:
    moth_stat_sampler m_sampler;
    std::vector<moth_real_t> m_pdf_x;
    std::vector<moth_real_t> m_pdf_y;

public:
    void estimate_histogram(moth_size_t n)
    {
        moth_real_t x_min{m_sampler.min_x()};
        moth_real_t x_max{m_sampler.max_x()};
        moth_real_t h = (x_max - x_min) / n;

        /* Generate the range. */
        m_pdf_x.clear();
        m_pdf_x.push_back(x_min);
        for (moth_size_t i = 0; i < n; ++i) {
            moth_real_t x = x_min + (i + 0.5) * h;
            m_pdf_x.push_back(x);
            m_pdf_x.push_back(x);
        }
        m_pdf_x.push_back(x_max);
        estimate();
    }

public:
    void estimate_adaptive()
    {
        moth_real_t x_min{m_sampler.min_x()};
        moth_real_t x_max{m_sampler.max_x()};
        moth_real_t x_mid{0.5 * (x_min + x_max)};
        m_pdf_x = {x_min, x_mid, x_mid, x_max};
        for (;;) {
            bool ff = false;
            estimate();
            moth_size_t n = m_pdf_x.size() / 2 - 1;
            for (moth_size_t i = 1; i <= n; ++i) {
                moth_real_t y0 = m_pdf_y[2 * i - 0];
                moth_real_t y1 = m_pdf_y[2 * i - 1];
                moth_real_t yy = fabs(y1 - y0) / fabs(y1 + y0);
                if (yy > 0.7) {
                    moth_real_t x0 = m_pdf_x[2 * i + 0];
                    moth_real_t x1 = m_pdf_x[2 * i + 1];
                    moth_real_t xx = 0.5 * (x0 + x1);
                    m_pdf_x.insert(m_pdf_x.begin() + 2 * i, xx);
                    m_pdf_x.insert(m_pdf_x.begin() + 2 * i, xx);
                    ff = true;
                }
            }
            if (!ff) {
                break;
            }
        }
    }

public:
    void estimate()
    {
        moth_size_t n = m_pdf_x.size() / 2 - 1;
        std::vector<moth_real_t> m_pdf_e(n + 1);
        std::vector<moth_real_t> m_pdf_d(n + 1);
        for (moth_size_t k = 1;; ++k) {
            moth_sample_t pdf_sample(m_sampler.gen_sample(10000));
            moth_real_t dspr_max = 0.0;
            for (moth_size_t i = 0; i <= n; ++i) {
                moth_real_t x0 = m_pdf_x[2 * i + 0];
                moth_real_t x1 = m_pdf_x[2 * i + 1];
                moth_real_t m = 0.0;
                for (moth_real_t x : pdf_sample) {
                    if (x0 <= x && x <= x1) {
                        m += 1.0;
                    }
                }

                moth_real_t y = m / pdf_sample.size() / (x1 - x0);
                m_pdf_e[i] += y;
                m_pdf_d[i] += y * y;

                moth_real_t mean = m_pdf_e[i] / k;
                moth_real_t dspr = mean * mean - m_pdf_d[i] / k;
                dspr_max = std::max(dspr_max, fabs(dspr) / k);
            }

            std::cerr << dspr_max << std::endl;
            if (0.0 < dspr_max && dspr_max < 1e-5) {
                for (moth_size_t i = 0; i <= n; ++i) {
                    moth_real_t mean = m_pdf_e[i] / k;
                    moth_real_t y = mean;
                    m_pdf_y.push_back(y);
                    m_pdf_y.push_back(y);
                }
                break;
            }
        }
    }
};

static void moth_print(const moth_sample_t& sample, const char* path)
{
    std::ofstream file(path);
    for (moth_real_t x : sample) {
        file << x << " " << x << std::endl;
    }
}

static void moth_print(const moth_stat_pdf_estimator& est, const char* path)
{
    std::ofstream file(path);
    for (moth_size_t i  = 0; i < est.m_pdf_x.size(); ++i) {
        moth_real_t x = est.m_pdf_x[i];
        moth_real_t y = est.m_pdf_y[i];
        file << x << " " << y << std::endl;
    }
}

int main()
{
    moth_stat_sampler density{{0.0, 0.7, 0.7, 1.0}, {0.1, 0.1, 0.9, 0.9}};
    moth_stat_pdf_estimator est{density};

    est.estimate_adaptive();
    moth_print(est, "sample.txt");

    return 0;
}
