#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

#include "RandomGenerator/LinearCongruential.h"
#include "RandomGenerator/Normal.h"
#include "SDE/BSMilstein1D.h"
#include "SDE/BSMilsteinND.h"
#include "Payoffs/EuropeanCallPayoff.h"
#include "Payoffs/EuroCallBasketPayOff.h"
#include "Pricer/EuropeanMCPricer.h"
#include "Pricer/VarRedMCPricer.h"
#include "VarRed/BasketGeomControlVariate.h"

// NOTE: QMC (Halton) vs plain MC is tested separately in test_qmc.cpp, where
// nsteps=1 is required to keep Halton bases small (large prime bases cause
// extreme Box-Muller values). Here we use nsteps=50 so that Milstein
// discretisation error is negligible and the plain MC price converges to the
// true BS price -- a prerequisite for a meaningful CV comparison.

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
             double price_v, double ci_v)
{
    double var_p = popVar(ci_p);
    double var_v = popVar(ci_v);
    double red   = (1.0 - var_v / var_p) * 100.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Plain MC  price=" << price_p
              << "  var=" << std::setprecision(1) << var_p
              << "  std=" << std::sqrt(var_p)
              << "  CI=+/-" << std::setprecision(4) << ci_p << "\n";
    std::cout << "  Static CV price=" << price_v
              << "  var=" << std::setprecision(1) << var_v
              << "  std=" << std::sqrt(var_v)
              << "  CI=+/-" << ci_v
              << "  red=" << std::setprecision(1) << red << "%\n";

    check(std::abs(price_p - price_v) < 0.5, label + ": CV price agrees with plain MC");
    check(ci_v < ci_p,                        label + ": static CV reduces variance");
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

    LinearCongruential u_v(SEED, A, C, M);
    Normal norm_v(0.0, 1.0, u_v);
    BSMilstein1D proc_v(&norm_v, SPOT, RATE, vol);
    EuropeanCallPayoff payoff_v(STRIKE);
    BasketGeomControlVariate cv({1.0}, {SPOT}, {{1.0}}, {vol}, T, STRIKE, RATE);
    VarRedMCPricer pricer_v(&proc_v, &payoff_v, &cv, RATE, T, NSIM, NSTEPS);
    auto r_v = pricer_v.Price();

    compare("1D sig=0.2", r_p.price, r_p.confidenceInterval,
                          r_v.price, r_v.confidenceInterval);
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

    LinearCongruential u_v(SEED, A, C, M);
    Normal norm_v(0.0, 1.0, u_v);
    BSMilsteinND proc_v(&norm_v, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_v(STRIKE, weights);
    BasketGeomControlVariate cv(weights, spots, corr, vols, T, STRIKE, RATE);
    VarRedMCPricer pricer_v(&proc_v, &payoff_v, &cv, RATE, T, NSIM, NSTEPS);
    auto r_v = pricer_v.Price();

    compare("2D w={0.5,0.5} sig={0.2,0.2}", r_p.price, r_p.confidenceInterval,
                                             r_v.price, r_v.confidenceInterval);
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

    LinearCongruential u_v(SEED, A, C, M);
    Normal norm_v(0.0, 1.0, u_v);
    BSMilsteinND proc_v(&norm_v, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_v(STRIKE, weights);
    BasketGeomControlVariate cv(weights, spots, corr, vols, T, STRIKE, RATE);
    VarRedMCPricer pricer_v(&proc_v, &payoff_v, &cv, RATE, T, NSIM, NSTEPS);
    auto r_v = pricer_v.Price();

    compare("2D w={0.7,0.3} sig={0.2,0.2}", r_p.price, r_p.confidenceInterval,
                                             r_v.price, r_v.confidenceInterval);
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

    LinearCongruential u_v(SEED, A, C, M);
    Normal norm_v(0.0, 1.0, u_v);
    BSMilsteinND proc_v(&norm_v, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_v(STRIKE, weights);
    BasketGeomControlVariate cv(weights, spots, corr, vols, T, STRIKE, RATE);
    VarRedMCPricer pricer_v(&proc_v, &payoff_v, &cv, RATE, T, NSIM, NSTEPS);
    auto r_v = pricer_v.Price();

    compare("2D w={0.5,0.5} sig={0.1,0.3}", r_p.price, r_p.confidenceInterval,
                                             r_v.price, r_v.confidenceInterval);
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

    LinearCongruential u_v(SEED, A, C, M);
    Normal norm_v(0.0, 1.0, u_v);
    BSMilsteinND proc_v(&norm_v, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_v(STRIKE, weights);
    BasketGeomControlVariate cv(weights, spots, corr, vols, T, STRIKE, RATE);
    VarRedMCPricer pricer_v(&proc_v, &payoff_v, &cv, RATE, T, NSIM, NSTEPS);
    auto r_v = pricer_v.Price();

    compare("3D w={1/3,1/3,1/3} sig={0.2,0.2,0.2}", r_p.price, r_p.confidenceInterval,
                                                      r_v.price, r_v.confidenceInterval);
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

    LinearCongruential u_v(SEED, A, C, M);
    Normal norm_v(0.0, 1.0, u_v);
    BSMilsteinND proc_v(&norm_v, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_v(STRIKE, weights);
    BasketGeomControlVariate cv(weights, spots, corr, vols, T, STRIKE, RATE);
    VarRedMCPricer pricer_v(&proc_v, &payoff_v, &cv, RATE, T, NSIM, NSTEPS);
    auto r_v = pricer_v.Price();

    compare("3D w={0.6,0.3,0.1} sig={0.2,0.2,0.2}", r_p.price, r_p.confidenceInterval,
                                                      r_v.price, r_v.confidenceInterval);
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

    LinearCongruential u_v(SEED, A, C, M);
    Normal norm_v(0.0, 1.0, u_v);
    BSMilsteinND proc_v(&norm_v, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_v(STRIKE, weights);
    BasketGeomControlVariate cv(weights, spots, corr, vols, T, STRIKE, RATE);
    VarRedMCPricer pricer_v(&proc_v, &payoff_v, &cv, RATE, T, NSIM, NSTEPS);
    auto r_v = pricer_v.Price();

    compare("3D w={1/3,1/3,1/3} sig={0.1,0.2,0.35}", r_p.price, r_p.confidenceInterval,
                                                       r_v.price, r_v.confidenceInterval);
}

// ---------------------------------------------------------------------------

int main() {
    std::cout << "=== Static Control Variate vs Plain MC ===\n";
    std::cout << "N=" << NSIM << " sims, nsteps=" << NSTEPS
              << ", S0=K=" << SPOT << ", r=" << RATE << ", T=" << T << "\n";
    std::cout << "(QMC comparison is in tests/run_tests_qmc — nsteps=1 required there)\n";

    test_1d();
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
