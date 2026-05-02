#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

#include "RandomGenerator/LinearCongruential.h"
#include "RandomGenerator/Normal.h"
#include "RandomGenerator/AntitheticNormal.h"
#include "SDE/BSMilstein1D.h"
#include "SDE/BSMilsteinND.h"
#include "Payoffs/EuropeanCallPayoff.h"
#include "Payoffs/EuroCallBasketPayOff.h"
#include "Pricer/EuropeanMCPricer.h"
#include "Pricer/AntitheticMCPricer.h"

// nsteps=50: Milstein discretisation error is negligible, prices converge to
// true BS values. Both plain MC and antithetic use the same number of path
// evaluations (nbSim each), so the comparison is on equal footing.

static int failures = 0;

void check(bool condition, const std::string& name) {
    if (condition) std::cout << "[PASS] " << name << "\n";
    else { std::cout << "[FAIL] " << name << "\n"; ++failures; }
}

static const size_t SEED = 42, A = 16807, C = 0, M = 2147483647;
static const size_t NSIM   = 10000;
static const size_t NSTEPS = 50;
static const double SPOT = 100.0, STRIKE = 100.0, RATE = 0.05, T = 1.0;

static double popVar(double ci) {
    double se = ci / 1.96;
    return se * se * static_cast<double>(NSIM);
}

void compare(const std::string& label,
             double price_p, double ci_p,
             double price_a, double ci_a)
{
    double var_p = popVar(ci_p);
    double var_a = popVar(ci_a);
    double red   = (1.0 - var_a / var_p) * 100.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Plain MC    price=" << price_p
              << "  var=" << std::setprecision(1) << var_p
              << "  std=" << std::sqrt(var_p)
              << "  CI=+/-" << std::setprecision(4) << ci_p << "\n";
    std::cout << "  Antithetic  price=" << price_a
              << "  var=" << std::setprecision(1) << var_a
              << "  std=" << std::sqrt(var_a)
              << "  CI=+/-" << ci_a
              << "  red=" << std::setprecision(1) << red << "%\n";

    check(std::abs(price_p - price_a) < 0.5, label + ": antithetic price agrees with plain MC");
    check(ci_a < ci_p,                        label + ": antithetic reduces variance");
}

// ---------------------------------------------------------------------------

void test_1d() {
    const double vol = 0.2;
    std::cout << "\n[BSMilstein1D  K=" << STRIKE << "  sig=" << vol << "]\n";

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilstein1D proc_p(&norm_p, SPOT, RATE, vol);
    EuropeanCallPayoff payoff_p(STRIKE);
    EuropeanMCPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS);
    auto r_p = pricer_p.Price();

    LinearCongruential u_a(SEED, A, C, M);
    Normal norm_raw(0.0, 1.0, u_a);
    AntitheticNormal anti(norm_raw);
    BSMilstein1D proc_a(&anti, SPOT, RATE, vol);
    EuropeanCallPayoff payoff_a(STRIKE);
    AntitheticMCPricer pricer_a(&proc_a, &payoff_a, anti, RATE, T, NSIM, NSTEPS);
    auto r_a = pricer_a.Price();

    compare("1D sig=0.2", r_p.price, r_p.confidenceInterval,
                          r_a.price, r_a.confidenceInterval);
}

void test_1d_highvol() {
    const double vol = 0.4;
    std::cout << "\n[BSMilstein1D  K=" << STRIKE << "  sig=" << vol << "]\n";

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilstein1D proc_p(&norm_p, SPOT, RATE, vol);
    EuropeanCallPayoff payoff_p(STRIKE);
    EuropeanMCPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS);
    auto r_p = pricer_p.Price();

    LinearCongruential u_a(SEED, A, C, M);
    Normal norm_raw(0.0, 1.0, u_a);
    AntitheticNormal anti(norm_raw);
    BSMilstein1D proc_a(&anti, SPOT, RATE, vol);
    EuropeanCallPayoff payoff_a(STRIKE);
    AntitheticMCPricer pricer_a(&proc_a, &payoff_a, anti, RATE, T, NSIM, NSTEPS);
    auto r_a = pricer_a.Price();

    compare("1D sig=0.4", r_p.price, r_p.confidenceInterval,
                          r_a.price, r_a.confidenceInterval);
}

void test_2d_equal() {
    const double vol = 0.2;
    std::cout << "\n[BSMilsteinND N=2  K=" << STRIKE << "  w={0.5,0.5}  sig={" << vol << "," << vol << "}  rho=0.3]\n";

    std::vector<double> spots   = {SPOT, SPOT};
    std::vector<double> vols    = {vol, vol};
    std::vector<std::vector<double>> corr = {{1.0, 0.3}, {0.3, 1.0}};
    std::vector<double> weights = {0.5, 0.5};

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilsteinND proc_p(&norm_p, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_p(STRIKE, weights);
    EuropeanMCPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS);
    auto r_p = pricer_p.Price();

    LinearCongruential u_a(SEED, A, C, M);
    Normal norm_raw(0.0, 1.0, u_a);
    AntitheticNormal anti(norm_raw);
    BSMilsteinND proc_a(&anti, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_a(STRIKE, weights);
    AntitheticMCPricer pricer_a(&proc_a, &payoff_a, anti, RATE, T, NSIM, NSTEPS);
    auto r_a = pricer_a.Price();

    compare("2D w={0.5,0.5} sig={0.2,0.2}", r_p.price, r_p.confidenceInterval,
                                             r_a.price, r_a.confidenceInterval);
}

void test_2d_unequal_weights() {
    const double vol = 0.2;
    std::cout << "\n[BSMilsteinND N=2  K=" << STRIKE << "  w={0.7,0.3}  sig={" << vol << "," << vol << "}  rho=0.3]\n";

    std::vector<double> spots   = {SPOT, SPOT};
    std::vector<double> vols    = {vol, vol};
    std::vector<std::vector<double>> corr = {{1.0, 0.3}, {0.3, 1.0}};
    std::vector<double> weights = {0.7, 0.3};

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilsteinND proc_p(&norm_p, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_p(STRIKE, weights);
    EuropeanMCPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS);
    auto r_p = pricer_p.Price();

    LinearCongruential u_a(SEED, A, C, M);
    Normal norm_raw(0.0, 1.0, u_a);
    AntitheticNormal anti(norm_raw);
    BSMilsteinND proc_a(&anti, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_a(STRIKE, weights);
    AntitheticMCPricer pricer_a(&proc_a, &payoff_a, anti, RATE, T, NSIM, NSTEPS);
    auto r_a = pricer_a.Price();

    compare("2D w={0.7,0.3} sig={0.2,0.2}", r_p.price, r_p.confidenceInterval,
                                             r_a.price, r_a.confidenceInterval);
}

void test_2d_diffvols() {
    std::cout << "\n[BSMilsteinND N=2  K=" << STRIKE << "  w={0.5,0.5}  sig={0.1,0.3}  rho=0.3]\n";

    std::vector<double> spots   = {SPOT, SPOT};
    std::vector<double> vols    = {0.1, 0.3};
    std::vector<std::vector<double>> corr = {{1.0, 0.3}, {0.3, 1.0}};
    std::vector<double> weights = {0.5, 0.5};

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilsteinND proc_p(&norm_p, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_p(STRIKE, weights);
    EuropeanMCPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS);
    auto r_p = pricer_p.Price();

    LinearCongruential u_a(SEED, A, C, M);
    Normal norm_raw(0.0, 1.0, u_a);
    AntitheticNormal anti(norm_raw);
    BSMilsteinND proc_a(&anti, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_a(STRIKE, weights);
    AntitheticMCPricer pricer_a(&proc_a, &payoff_a, anti, RATE, T, NSIM, NSTEPS);
    auto r_a = pricer_a.Price();

    compare("2D w={0.5,0.5} sig={0.1,0.3}", r_p.price, r_p.confidenceInterval,
                                             r_a.price, r_a.confidenceInterval);
}

void test_3d_equal() {
    const double vol = 0.2;
    std::cout << "\n[BSMilsteinND N=3  K=" << STRIKE << "  w={1/3,1/3,1/3}  sig={" << vol << "," << vol << "," << vol << "}  rho=0.3]\n";

    std::vector<double> spots   = {SPOT, SPOT, SPOT};
    std::vector<double> vols    = {vol, vol, vol};
    std::vector<std::vector<double>> corr = {{1,0.3,0.3},{0.3,1,0.3},{0.3,0.3,1}};
    std::vector<double> weights = {1.0/3, 1.0/3, 1.0/3};

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilsteinND proc_p(&norm_p, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_p(STRIKE, weights);
    EuropeanMCPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS);
    auto r_p = pricer_p.Price();

    LinearCongruential u_a(SEED, A, C, M);
    Normal norm_raw(0.0, 1.0, u_a);
    AntitheticNormal anti(norm_raw);
    BSMilsteinND proc_a(&anti, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_a(STRIKE, weights);
    AntitheticMCPricer pricer_a(&proc_a, &payoff_a, anti, RATE, T, NSIM, NSTEPS);
    auto r_a = pricer_a.Price();

    compare("3D w={1/3,1/3,1/3} sig={0.2,0.2,0.2}", r_p.price, r_p.confidenceInterval,
                                                      r_a.price, r_a.confidenceInterval);
}

void test_3d_unequal_weights() {
    const double vol = 0.2;
    std::cout << "\n[BSMilsteinND N=3  K=" << STRIKE << "  w={0.6,0.3,0.1}  sig={" << vol << "," << vol << "," << vol << "}  rho=0.3]\n";

    std::vector<double> spots   = {SPOT, SPOT, SPOT};
    std::vector<double> vols    = {vol, vol, vol};
    std::vector<std::vector<double>> corr = {{1,0.3,0.3},{0.3,1,0.3},{0.3,0.3,1}};
    std::vector<double> weights = {0.6, 0.3, 0.1};

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilsteinND proc_p(&norm_p, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_p(STRIKE, weights);
    EuropeanMCPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS);
    auto r_p = pricer_p.Price();

    LinearCongruential u_a(SEED, A, C, M);
    Normal norm_raw(0.0, 1.0, u_a);
    AntitheticNormal anti(norm_raw);
    BSMilsteinND proc_a(&anti, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_a(STRIKE, weights);
    AntitheticMCPricer pricer_a(&proc_a, &payoff_a, anti, RATE, T, NSIM, NSTEPS);
    auto r_a = pricer_a.Price();

    compare("3D w={0.6,0.3,0.1} sig={0.2,0.2,0.2}", r_p.price, r_p.confidenceInterval,
                                                      r_a.price, r_a.confidenceInterval);
}

void test_3d_diffvols() {
    std::cout << "\n[BSMilsteinND N=3  K=" << STRIKE << "  w={1/3,1/3,1/3}  sig={0.1,0.2,0.35}  rho=0.3]\n";

    std::vector<double> spots   = {SPOT, SPOT, SPOT};
    std::vector<double> vols    = {0.1, 0.2, 0.35};
    std::vector<std::vector<double>> corr = {{1,0.3,0.3},{0.3,1,0.3},{0.3,0.3,1}};
    std::vector<double> weights = {1.0/3, 1.0/3, 1.0/3};

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilsteinND proc_p(&norm_p, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_p(STRIKE, weights);
    EuropeanMCPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS);
    auto r_p = pricer_p.Price();

    LinearCongruential u_a(SEED, A, C, M);
    Normal norm_raw(0.0, 1.0, u_a);
    AntitheticNormal anti(norm_raw);
    BSMilsteinND proc_a(&anti, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_a(STRIKE, weights);
    AntitheticMCPricer pricer_a(&proc_a, &payoff_a, anti, RATE, T, NSIM, NSTEPS);
    auto r_a = pricer_a.Price();

    compare("3D w={1/3,1/3,1/3} sig={0.1,0.2,0.35}", r_p.price, r_p.confidenceInterval,
                                                       r_a.price, r_a.confidenceInterval);
}

// ---------------------------------------------------------------------------

int main() {
    std::cout << "=== Antithetic Variables vs Plain MC ===\n";
    std::cout << "N=" << NSIM << " sims (pairs for antithetic), nsteps=" << NSTEPS
              << ", S0=K=" << SPOT << ", r=" << RATE << ", T=" << T << "\n";

    test_1d();
    test_1d_highvol();
    test_2d_equal();
    test_2d_unequal_weights();
    test_2d_diffvols();
    test_3d_equal();
    test_3d_unequal_weights();
    test_3d_diffvols();

    std::cout << "\n";
    if (failures == 0) std::cout << "All tests PASSED.\n";
    else               std::cout << failures << " test(s) FAILED.\n";
    return failures == 0 ? 0 : 1;
}
