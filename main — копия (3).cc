#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include <fstream>
#include <iostream>

static const double Gamma = 5.0 / 3.0;

std::vector<double> x_i;
std::vector<double> d_i;
std::vector<double> v_i;
std::vector<double> p_i;

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

/**
 * Abstract Euler equations solver.
 */
class euler_solver
{
public:
    virtual void integrate(double Tau) = 0;
};  // class euler_solver

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

/**
 * Exact solver by Toro.
 */
class euler_solver_toro : public euler_solver
{
public:
    void integrate(double Tau) override
    {
        rl = d_i.front();
        rr = d_i.back();
        ul = v_i.front();
        ur = v_i.back();
        pl = p_i.front();
        pr = p_i.back();

        g1 = (Gamma - 1.0) / (2.0 * Gamma);
        g2 = (Gamma + 1.0) / (2.0 * Gamma);
        g3 = 2.0 * Gamma / (Gamma - 1.0);
        g4 = 2.0/(Gamma - 1.0);
        g5 = 2.0/(Gamma + 1.0);
        g6 = (Gamma - 1.0) / (Gamma + 1.0);
        g7 = (Gamma - 1.0) / 2.0;

        // compute sound speeds
        cl = sqrt(Gamma * pl / rl);
        cr = sqrt(Gamma * pr / rr);

        // the pressure positivity condition is tested for
        if (g4*(cl+cr) <= (ur-ul)) {

            std::cerr << "the initial data is such that vacuum is generated"
                      << "\nstopping program" << std::endl;
            exit(1);
        }

        double pm, um;

        // exact solution for pressure and velocity in star region is found
        starpu(pm, um, 1.0);

        for(int i = 1; i < x_i.size(); ++i)
        {
            auto xpos = x_i[i];
            auto s = (xpos - 0.5) / Tau;

            double rs, us, ps;
            sample(pm, um, s, rs, us, ps);

            d_i[i] = rs;
            v_i[i] = us;
            p_i[i] = ps;
        }
    }

private:
    double g1, g2, g3, g4, g5, g6, g7;
    double rl, ul, pl, cl, rr, ur, pr, cr;

    void starpu(double& p, double& u, double pscale)
    {
        // purpose: to compute the solution for pressure and
        //          velocity in the Star Region

        const int nriter = 20;
        const double tolpre = 1.0e-6;
        double change, fl, fld, fr, frd, pold, pstart, udiff;

        // guessed value pstart is computed
        guessp(pstart);
        pold = pstart;
        udiff = ur - ul;

        std::cout << "----------------------------------------\n"
                  << "   Iteration number     Change\n"
                  << "----------------------------------------" << std::endl;

        int i = 1;
        for ( ; i <= nriter; i++) {
            prefun(fl, fld, pold, rl, pl, cl);
            prefun(fr, frd, pold, rr, pr, cr);
            p = pold - (fl + fr + udiff)/(fld + frd);
            change = 2.0*fabs((p - pold)/(p + pold));
            std::cout << '\t' << i <<  "\t\t" << change << std::endl;
            if (change <= tolpre)
                break;
            if (p < 0.0)
                p = tolpre;
            pold = p;
        }
        if (i > nriter) {
            std::cout << "divergence in Newton-Raphson iteration" << std::endl;
            return;
        }

        // compute velocity in star region
        u = 0.5*(ul + ur + fr - fl);
        std::cout << "----------------------------------------\n"
                  << "     Pressure           Velocity\n"
                  << "----------------------------------------\n"
                  << "     " << p/pscale << "\t\t" <<  u << '\n'
                  << "----------------------------------------" << std::endl;
    }

    void sample(double pm, double um, double s, double& d, double& u, double& p)
    {
        // purpose: to sample the solution throughout the wave
        //          pattern. Pressure pm and velocity um in the
        //          star region are known. Sampling is performed
        //          in terms of the 'speed' s = x/t. Sampled
        //          values are d, u, p
        double c, cml, cmr, pml, pmr, shl, shr, sl, sr, stl, str;

        if (s <= um) {
            // sampling point lies to the left of the contact discontinuity
            if (pm <= pl) {
                // left rarefaction
                shl = ul - cl;
                if (s <= shl) {
                    // sampled point is left data state
                    d = rl;
                    u = ul;
                    p = pl;
                } else {
                    cml = cl*pow(pm/pl, g1);
                    stl = um - cml;
                    if (s > stl) {
                        // sampled point is star left state
                        d = rl*pow(pm/pl, 1.0/Gamma);
                        u = um;
                        p = pm;
                    } else {
                        // sampled point is inside left fan
                        u = g5*(cl + g7*ul + s);
                        c = g5*(cl + g7*(ul - s));
                        d = rl*pow(c/cl, g4);
                        p = pl*pow(c/cl, g3);
                    }
                }
            } else {
                // left shock
                pml = pm/pl;
                sl = ul - cl*sqrt(g2*pml + g1);
                if (s <= sl) {
                    // sampled point is left data state
                    d = rl;
                    u = ul;
                    p = pl;
                } else {
                    // sampled point is star left state
                    d = rl*(pml + g6)/(pml*g6 + 1.0);
                    u = um;
                    p = pm;
                }
            }
        } else {
            // sampling point lies to the right of the contact discontinuity
            if (pm > pr) {
                // right shock
                pmr = pm/pr;
                sr  = ur + cr*sqrt(g2*pmr + g1);
                if (s >= sr) {
                    // sampled point is right data state
                    d = rr;
                    u = ur;
                    p = pr;
                } else {
                    // sampled point is star right state
                    d = rr*(pmr + g6)/(pmr*g6 + 1.0);
                    u = um;
                    p = pm;
                }
            } else {
                // right rarefaction
                shr = ur + cr;
                if (s >= shr) {
                    // sampled point is right data state
                    d = rr;
                    u = ur;
                    p = pr;
                } else {
                    cmr = cr*pow(pm/pr, g1);
                    str = um + cmr;
                    if (s <= str) {
                        // sampled point is star right state
                        d = rr*pow(pm/pr, 1.0/Gamma);
                        u = um;
                        p = pm;
                    } else {
                        // sampled point is inside left fan
                        u = g5*(-cr + g7*ur + s);
                        c = g5*(cr - g7*(ur - s));
                        d = rr*pow(c/cr, g4);
                        p = pr*pow(c/cr, g3);
                    }
                }
            }
        }
    }

    void guessp(double& pm)
    {
        // purpose: to provide a guessed value for pressure
        //          pm in the Star Region. The choice is made
        //          according to adaptive Riemann solver using
        //          the PVRS, TRRS and TSRS approximate
        //          Riemann solvers. See Sect. 9.5 of Chapt. 9 of Ref. 1

        double cup, gel, ger, pmax, pmin, ppv, pq, ptl, ptr,
                qmax, quser, um;

        quser = 2.0;

        // compute guess pressure from PVRS Riemann solver
        cup = 0.25*(rl + rr)*(cl + cr);
        ppv = 0.5*(pl + pr) + 0.5*(ul - ur)*cup;
        ppv = std::max(0.0, ppv);
        pmin = std::min(pl,  pr);
        pmax = std::max(pl,  pr);
        qmax = pmax/pmin;

        if (qmax <= quser && (pmin <= ppv && ppv <= pmax))
            pm = ppv;     // select PVRS Riemann solver
        else {
            if (ppv < pmin) {
                // select Two-Rarefaction Riemann solver
                pq = pow(pl/pr, g1);
                um = (pq*ul/cl + ur/cr + g4*(pq - 1.0))/(pq/cl + 1.0/cr);
                ptl = 1.0 + g7*(ul - um)/cl;
                ptr = 1.0 + g7*(um - ur)/cr;
                pm = 0.5*(pow(pl*ptl, g3) + pow(pr*ptr, g3));
            } else {
                // select Two-Shock Riemann solver with PVRS as estimate
                gel = sqrt((g5/rl)/(g6*pl + ppv));
                ger = sqrt((g5/rr)/(g6*pr + ppv));
                pm = (gel*pl + ger*pr - (ur - ul))/(gel + ger);
            }
        }
    }

    void prefun(double& f, double& fd, double& p, double& dk, double& pk, double& ck)
    {
        // purpose: to evaluate the pressure functions
        //          fl and fr in exact Riemann solver
        //          and their first derivatives

        double ak, bk, pratio, qrt;
        if (p <= pk) {
            // rarefaction wave
            pratio = p/pk;
            f = g4*ck*(pow(pratio, g1) - 1.0);
            fd = (1.0/(dk*ck))*pow(pratio, -g2);
        } else {
            //  shock wave
            ak = g5/dk;
            bk = g6*pk;
            qrt = sqrt(ak/(bk + p));
            f = (p - pk)*qrt;
            fd = (1.0 - 0.5*(p - pk)/(bk + p))*qrt;
        }
    }
};  // class euler_solver_toro

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

/**
 * Abstract solver in Lagrange variables.
 */
class euler_solver_invariants : public euler_solver
{
private:
    std::vector<std::pair<double, double>> r1_i;
    std::vector<std::pair<double, double>> r2_i;
    std::vector<std::pair<double, double>> r3_i;

private:
    /**
     * Interpolate Riemann invariants in given point.
     */
    virtual double interp_riemann_invariant(double x, size_t k) const
    {
        const decltype(r1_i)* rk_i;
        switch (k) {
            case 1:
                rk_i = &r1_i;
                break;
            case 2:
                rk_i = &r2_i;
                break;
            case 3:
                rk_i = &r3_i;
                break;
            default:
                abort();
        }

        for (std::size_t j = 1; j < r1_i.size(); ++j) {
            const double xrp = (*rk_i)[j - 0].first;
            const double xrm = (*rk_i)[j - 1].first;
            if (xrm <= x && x <= xrp) {
                const double rp = (*rk_i)[j - 0].second;
                const double rm = (*rk_i)[j - 1].second;

                /* Simple linear interpolation. */
                const double Delta = (x - xrm)/(xrp - xrm);
                return Delta*rp + (1.0 - Delta)*rm;
            }
        }
        abort();
    }

public:
    void integrate(double Tau) override
    {
        const std::size_t N = x_i.size() - 1;

        /* Calculate and advent the Riemann invariants. */
        r1_i.resize(N+1);
        r2_i.resize(N+1);
        r3_i.resize(N+1);
        for (std::size_t i = 0; i <= N; ++i) {
            const double x = x_i[i];
            const double d = d_i[i];
            const double v = v_i[i];
            const double p = p_i[i];

            /* Speed of sound and Entropy. */
            const double c = std::sqrt(Gamma * p / d);
            const double s = p / std::pow(d, Gamma);

            /* Invariants. */
            r1_i[i] = {x + Tau * (v), s};
            r2_i[i] = {x + Tau * (v + c), v + 2.0 * c / (Gamma - 1.0)};
            r3_i[i] = {x + Tau * (v - c), v - 2.0 * c / (Gamma - 1.0)};
        }

        /* Project the invariants on the Grid. */
        for (std::size_t i = 1; i < N; ++i) {
            const double x = x_i[i];

            /* Estimate the Riemann invariants on current grid. */
            double Rz = interp_riemann_invariant(x, 1);
            double Rp = interp_riemann_invariant(x, 2);
            double Rm = interp_riemann_invariant(x, 3);

            /* Recalculate the gas parameters. */
            const double v = 0.50 * (Rp + Rm);
            const double c = 0.25 * (Rp - Rm) * (Gamma - 1.0);

            const double s = Rz;
            const double d = std::pow(c * c / s / Gamma, 1.0 / (Gamma - 1.0));
            const double p = c * c * d / Gamma;

            v_i[i] = v;
            d_i[i] = d;
            p_i[i] = p;
        }
    }
};  // class euler_solver

int main()
{
    const std::size_t N = 200;
    const double x_0 = 0.0;
    const double x_l = 1.0;
    const double H = (x_l - x_0) / N;

    x_i.resize(N+1);
    d_i.resize(N+1);
    v_i.resize(N+1);
    p_i.resize(N+1);

    /* Initialize the initial uniform grid for
     * Riemann problem. */
    const double dm = 2.0, dp = 1.0;
    const double pm = 5.0, pp = 1.0;
    for (std::size_t i = 0; i <= N; ++i) {
        const double x = x_0 + i * H;
        x_i[i] = x;
        if (x <= 0.5 * (x_0 + x_l)) {
            d_i[i] = dm;
            p_i[i] = pm;
        } else {
            d_i[i] = dp;
            p_i[i] = pp;
        }
    }

#if 0
    euler_solver_toro solver_exact;
    solver_exact.integrate(0.001 * 50);
    std::string grid_file_path{"res/data-toro-" + std::to_string(50) + ".txt"};
    std::ofstream grid_file(grid_file_path);
    for (std::size_t i = 0; i <= N; ++i) {
        const double x = x_i[i];
        const double d = d_i[i];
        const double v = v_i[i];
        const double p = p_i[i];
        grid_file << x << "\t" << d << "\t" << v << "\t" << p << std::endl;
    }
#else
    /* Create the solver. */
    euler_solver_invariants solver;

    const double Tau = 0.001;
    for (std::size_t k = 0; k <= 50; ++k){

        /* Perform the time integrate. */
        solver.integrate(Tau);

        /* Output the result. */
        std::string grid_file_path{"res/data-" + std::to_string(k) + ".txt"};
        std::ofstream grid_file(grid_file_path);
        for (std::size_t i = 0; i <= N; ++i) {
            const double x = x_i[i];
            const double d = d_i[i];
            const double v = v_i[i];
            const double p = p_i[i];
            grid_file << x << "\t" << d << "\t" << v << "\t" << p << std::endl;
        }
    }
#endif

    return 0;
}
