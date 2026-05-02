#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

#include "RandomGenerator/LinearCongruential.h"
#include "RandomGenerator/Normal.h"
#include "RandomGenerator/HaltonGenerator.h"
#include "RandomGenerator/AntitheticNormal.h"
#include "SDE/BSMilstein1D.h"
#include "SDE/BSMilsteinND.h"
#include "Payoffs/EuropeanCallPayoff.h"
#include "Payoffs/EuroCallBasketPayOff.h"
#include "Pricer/EuropeanMCPricer.h"
#include "Pricer/VarRedMCPricer.h"
#include "Pricer/AntitheticMCPricer.h"
#include "VarRed/BasketGeomControlVariate.h"

// All methods use nsteps=1 so that QMC (Halton) operates in the correct regime
// (small prime bases, no extreme Box-Muller values). With nsteps=1 the terminal
// law is exact for geometric Brownian motion, so there is no discretisation bias.
// This allows a direct, apples-to-apples variance comparison across all four methods.
//
// Base counts for Halton:  nsteps * dim * 2 (Box-Muller pairs)
//   1D -> 2 bases {2,3}    2D -> 4 bases {2,3,5,7}    3D -> 6 bases {2,3,5,7,11,13}

static int failures = 0;

void check(bool condition, const std::string& name) {
    if (condition) std::cout << "[PASS] " << name << "\n";
    else { std::cout << "[FAIL] " << name << "\n"; ++failures; }
}

static const size_t SEED = 42, A = 16807, C = 0, M = 2147483647;
static const size_t NSIM   = 10000;
static const size_t NSTEPS = 1;
static const double SPOT = 100.0, STRIKE = 100.0, RATE = 0.05, T = 1.0;

static double popVar(double ci) {
    double se = ci / 1.96;
    return se * se * static_cast<double>(NSIM);
}

struct Result { double price, ci; };

void printRow(const std::string& method, Result r, double var_plain) {
    double var = popVar(r.ci);
    double red = (1.0 - var / var_plain) * 100.0;
    std::cout << std::fixed << std::setprecision(4)
              << "  " << std::left << std::setw(14) << method
              << "  price=" << r.price
              << "  var=" << std::setprecision(1) << std::setw(7) << var
              << "  CI=+/-" << std::setprecision(4) << r.ci
              << "  red=" << std::setprecision(1) << red << "%\n";
}

// ---------------------------------------------------------------------------

void test_1d() {
    const double vol = 0.2;
    std::cout << "\n[1D  K=" << STRIKE << "  sig=" << vol << "  nsteps=" << NSTEPS << "]\n";

    static const std::vector<size_t> halton_bases = {2, 3};

    // Plain MC
    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilstein1D proc_p(&norm_p, SPOT, RATE, vol);
    EuropeanCallPayoff pay_p(STRIKE);
    EuropeanMCPricer pricer_p(&proc_p, &pay_p, RATE, T, NSIM, NSTEPS);
    Result r_p = {0, 0};
    { auto r = pricer_p.Price(); r_p = {r.price, r.confidenceInterval}; }

    // QMC (Halton)
    HaltonGenerator u_q(halton_bases);
    Normal norm_q(0.0, 1.0, u_q);
    BSMilstein1D proc_q(&norm_q, SPOT, RATE, vol);
    EuropeanCallPayoff pay_q(STRIKE);
    EuropeanMCPricer pricer_q(&proc_q, &pay_q, RATE, T, NSIM, NSTEPS);
    Result r_q = {0, 0};
    { auto r = pricer_q.Price(); r_q = {r.price, r.confidenceInterval}; }

    // Antithetic
    LinearCongruential u_a(SEED, A, C, M);
    Normal norm_raw_a(0.0, 1.0, u_a);
    AntitheticNormal anti_a(norm_raw_a);
    BSMilstein1D proc_a(&anti_a, SPOT, RATE, vol);
    EuropeanCallPayoff pay_a(STRIKE);
    AntitheticMCPricer pricer_a(&proc_a, &pay_a, anti_a, RATE, T, NSIM, NSTEPS);
    Result r_a = {0, 0};
    { auto r = pricer_a.Price(); r_a = {r.price, r.confidenceInterval}; }

    // Static CV (geometric basket, weight=1.0)
    LinearCongruential u_v(SEED, A, C, M);
    Normal norm_v(0.0, 1.0, u_v);
    BSMilstein1D proc_v(&norm_v, SPOT, RATE, vol);
    EuropeanCallPayoff pay_v(STRIKE);
    BasketGeomControlVariate cv({1.0}, {SPOT}, {{1.0}}, {vol}, T, STRIKE, RATE);
    VarRedMCPricer pricer_v(&proc_v, &pay_v, &cv, RATE, T, NSIM, NSTEPS);
    Result r_v = {0, 0};
    { auto r = pricer_v.Price(); r_v = {r.price, r.confidenceInterval}; }

    double var_p = popVar(r_p.ci);
    printRow("Plain MC",    r_p, var_p);
    printRow("QMC",         r_q, var_p);
    printRow("Antithetic",  r_a, var_p);
    printRow("Static CV",   r_v, var_p);

    check(r_q.ci < r_p.ci * 1.08, "1D: QMC CI not significantly worse than plain MC");
    check(r_a.ci < r_p.ci,        "1D: antithetic reduces variance");
    check(r_v.ci < r_p.ci,        "1D: static CV reduces variance");
    check(std::abs(r_q.price - r_p.price) < 1.0, "1D: QMC price agrees with plain MC");
    check(std::abs(r_a.price - r_p.price) < 0.5, "1D: antithetic price agrees with plain MC");
    check(std::abs(r_v.price - r_p.price) < 0.5, "1D: static CV price agrees with plain MC");
}

void test_2d_equal() {
    const double vol = 0.2;
    std::cout << "\n[2D  K=" << STRIKE << "  w={0.5,0.5}  sig={" << vol << "," << vol << "}  rho=0.3  nsteps=" << NSTEPS << "]\n";

    static const std::vector<size_t> halton_bases = {2, 3, 5, 7};
    std::vector<double> spots   = {SPOT, SPOT};
    std::vector<double> vols    = {vol, vol};
    std::vector<std::vector<double>> corr = {{1.0, 0.3}, {0.3, 1.0}};
    std::vector<double> weights = {0.5, 0.5};

    // Plain MC
    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilsteinND proc_p(&norm_p, spots, RATE, vols, corr);
    EuroCallBasketPayOff pay_p(STRIKE, weights);
    EuropeanMCPricer pricer_p(&proc_p, &pay_p, RATE, T, NSIM, NSTEPS);
    Result r_p = {0, 0};
    { auto r = pricer_p.Price(); r_p = {r.price, r.confidenceInterval}; }

    // QMC (Halton)
    HaltonGenerator u_q(halton_bases);
    Normal norm_q(0.0, 1.0, u_q);
    BSMilsteinND proc_q(&norm_q, spots, RATE, vols, corr);
    EuroCallBasketPayOff pay_q(STRIKE, weights);
    EuropeanMCPricer pricer_q(&proc_q, &pay_q, RATE, T, NSIM, NSTEPS);
    Result r_q = {0, 0};
    { auto r = pricer_q.Price(); r_q = {r.price, r.confidenceInterval}; }

    // Antithetic
    LinearCongruential u_a(SEED, A, C, M);
    Normal norm_raw_a(0.0, 1.0, u_a);
    AntitheticNormal anti_a(norm_raw_a);
    BSMilsteinND proc_a(&anti_a, spots, RATE, vols, corr);
    EuroCallBasketPayOff pay_a(STRIKE, weights);
    AntitheticMCPricer pricer_a(&proc_a, &pay_a, anti_a, RATE, T, NSIM, NSTEPS);
    Result r_a = {0, 0};
    { auto r = pricer_a.Price(); r_a = {r.price, r.confidenceInterval}; }

    // Static CV
    LinearCongruential u_v(SEED, A, C, M);
    Normal norm_v(0.0, 1.0, u_v);
    BSMilsteinND proc_v(&norm_v, spots, RATE, vols, corr);
    EuroCallBasketPayOff pay_v(STRIKE, weights);
    BasketGeomControlVariate cv(weights, spots, corr, vols, T, STRIKE, RATE);
    VarRedMCPricer pricer_v(&proc_v, &pay_v, &cv, RATE, T, NSIM, NSTEPS);
    Result r_v = {0, 0};
    { auto r = pricer_v.Price(); r_v = {r.price, r.confidenceInterval}; }

    double var_p = popVar(r_p.ci);
    printRow("Plain MC",    r_p, var_p);
    printRow("QMC",         r_q, var_p);
    printRow("Antithetic",  r_a, var_p);
    printRow("Static CV",   r_v, var_p);

    check(r_q.ci < r_p.ci * 1.08, "2D w={0.5,0.5}: QMC CI not significantly worse");
    check(r_a.ci < r_p.ci,        "2D w={0.5,0.5}: antithetic reduces variance");
    check(r_v.ci < r_p.ci,        "2D w={0.5,0.5}: static CV reduces variance");
    check(std::abs(r_q.price - r_p.price) < 1.0, "2D w={0.5,0.5}: QMC price agrees");
    check(std::abs(r_a.price - r_p.price) < 0.5, "2D w={0.5,0.5}: antithetic price agrees");
    check(std::abs(r_v.price - r_p.price) < 0.5, "2D w={0.5,0.5}: static CV price agrees");
}

void test_2d_diffvols() {
    std::cout << "\n[2D  K=" << STRIKE << "  w={0.5,0.5}  sig={0.1,0.3}  rho=0.3  nsteps=" << NSTEPS << "]\n";

    static const std::vector<size_t> halton_bases = {2, 3, 5, 7};
    std::vector<double> spots   = {SPOT, SPOT};
    std::vector<double> vols    = {0.1, 0.3};
    std::vector<std::vector<double>> corr = {{1.0, 0.3}, {0.3, 1.0}};
    std::vector<double> weights = {0.5, 0.5};

    // Plain MC
    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilsteinND proc_p(&norm_p, spots, RATE, vols, corr);
    EuroCallBasketPayOff pay_p(STRIKE, weights);
    EuropeanMCPricer pricer_p(&proc_p, &pay_p, RATE, T, NSIM, NSTEPS);
    Result r_p = {0, 0};
    { auto r = pricer_p.Price(); r_p = {r.price, r.confidenceInterval}; }

    // QMC (Halton)
    HaltonGenerator u_q(halton_bases);
    Normal norm_q(0.0, 1.0, u_q);
    BSMilsteinND proc_q(&norm_q, spots, RATE, vols, corr);
    EuroCallBasketPayOff pay_q(STRIKE, weights);
    EuropeanMCPricer pricer_q(&proc_q, &pay_q, RATE, T, NSIM, NSTEPS);
    Result r_q = {0, 0};
    { auto r = pricer_q.Price(); r_q = {r.price, r.confidenceInterval}; }

    // Antithetic
    LinearCongruential u_a(SEED, A, C, M);
    Normal norm_raw_a(0.0, 1.0, u_a);
    AntitheticNormal anti_a(norm_raw_a);
    BSMilsteinND proc_a(&anti_a, spots, RATE, vols, corr);
    EuroCallBasketPayOff pay_a(STRIKE, weights);
    AntitheticMCPricer pricer_a(&proc_a, &pay_a, anti_a, RATE, T, NSIM, NSTEPS);
    Result r_a = {0, 0};
    { auto r = pricer_a.Price(); r_a = {r.price, r.confidenceInterval}; }

    // Static CV
    LinearCongruential u_v(SEED, A, C, M);
    Normal norm_v(0.0, 1.0, u_v);
    BSMilsteinND proc_v(&norm_v, spots, RATE, vols, corr);
    EuroCallBasketPayOff pay_v(STRIKE, weights);
    BasketGeomControlVariate cv(weights, spots, corr, vols, T, STRIKE, RATE);
    VarRedMCPricer pricer_v(&proc_v, &pay_v, &cv, RATE, T, NSIM, NSTEPS);
    Result r_v = {0, 0};
    { auto r = pricer_v.Price(); r_v = {r.price, r.confidenceInterval}; }

    double var_p = popVar(r_p.ci);
    printRow("Plain MC",    r_p, var_p);
    printRow("QMC",         r_q, var_p);
    printRow("Antithetic",  r_a, var_p);
    printRow("Static CV",   r_v, var_p);

    check(r_q.ci < r_p.ci * 1.08, "2D sig={0.1,0.3}: QMC CI not significantly worse");
    check(r_a.ci < r_p.ci,        "2D sig={0.1,0.3}: antithetic reduces variance");
    check(r_v.ci < r_p.ci,        "2D sig={0.1,0.3}: static CV reduces variance");
    check(std::abs(r_q.price - r_p.price) < 1.0, "2D sig={0.1,0.3}: QMC price agrees");
    check(std::abs(r_a.price - r_p.price) < 0.5, "2D sig={0.1,0.3}: antithetic price agrees");
    check(std::abs(r_v.price - r_p.price) < 0.5, "2D sig={0.1,0.3}: static CV price agrees");
}

void test_3d_equal() {
    const double vol = 0.2;
    std::cout << "\n[3D  K=" << STRIKE << "  w={1/3,1/3,1/3}  sig={" << vol << "," << vol << "," << vol << "}  rho=0.3  nsteps=" << NSTEPS << "]\n";

    static const std::vector<size_t> halton_bases = {2, 3, 5, 7, 11, 13};
    std::vector<double> spots   = {SPOT, SPOT, SPOT};
    std::vector<double> vols    = {vol, vol, vol};
    std::vector<std::vector<double>> corr = {{1,0.3,0.3},{0.3,1,0.3},{0.3,0.3,1}};
    std::vector<double> weights = {1.0/3, 1.0/3, 1.0/3};

    // Plain MC
    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilsteinND proc_p(&norm_p, spots, RATE, vols, corr);
    EuroCallBasketPayOff pay_p(STRIKE, weights);
    EuropeanMCPricer pricer_p(&proc_p, &pay_p, RATE, T, NSIM, NSTEPS);
    Result r_p = {0, 0};
    { auto r = pricer_p.Price(); r_p = {r.price, r.confidenceInterval}; }

    // QMC (Halton)
    HaltonGenerator u_q(halton_bases);
    Normal norm_q(0.0, 1.0, u_q);
    BSMilsteinND proc_q(&norm_q, spots, RATE, vols, corr);
    EuroCallBasketPayOff pay_q(STRIKE, weights);
    EuropeanMCPricer pricer_q(&proc_q, &pay_q, RATE, T, NSIM, NSTEPS);
    Result r_q = {0, 0};
    { auto r = pricer_q.Price(); r_q = {r.price, r.confidenceInterval}; }

    // Antithetic
    LinearCongruential u_a(SEED, A, C, M);
    Normal norm_raw_a(0.0, 1.0, u_a);
    AntitheticNormal anti_a(norm_raw_a);
    BSMilsteinND proc_a(&anti_a, spots, RATE, vols, corr);
    EuroCallBasketPayOff pay_a(STRIKE, weights);
    AntitheticMCPricer pricer_a(&proc_a, &pay_a, anti_a, RATE, T, NSIM, NSTEPS);
    Result r_a = {0, 0};
    { auto r = pricer_a.Price(); r_a = {r.price, r.confidenceInterval}; }

    // Static CV
    LinearCongruential u_v(SEED, A, C, M);
    Normal norm_v(0.0, 1.0, u_v);
    BSMilsteinND proc_v(&norm_v, spots, RATE, vols, corr);
    EuroCallBasketPayOff pay_v(STRIKE, weights);
    BasketGeomControlVariate cv(weights, spots, corr, vols, T, STRIKE, RATE);
    VarRedMCPricer pricer_v(&proc_v, &pay_v, &cv, RATE, T, NSIM, NSTEPS);
    Result r_v = {0, 0};
    { auto r = pricer_v.Price(); r_v = {r.price, r.confidenceInterval}; }

    double var_p = popVar(r_p.ci);
    printRow("Plain MC",    r_p, var_p);
    printRow("QMC",         r_q, var_p);
    printRow("Antithetic",  r_a, var_p);
    printRow("Static CV",   r_v, var_p);

    check(r_q.ci < r_p.ci * 1.08, "3D equal: QMC CI not significantly worse");
    check(r_a.ci < r_p.ci,        "3D equal: antithetic reduces variance");
    check(r_v.ci < r_p.ci,        "3D equal: static CV reduces variance");
    check(std::abs(r_q.price - r_p.price) < 1.0, "3D equal: QMC price agrees");
    check(std::abs(r_a.price - r_p.price) < 0.5, "3D equal: antithetic price agrees");
    check(std::abs(r_v.price - r_p.price) < 0.5, "3D equal: static CV price agrees");
}

// ---------------------------------------------------------------------------

int main() {
    std::cout << "=== Full Variance Reduction Comparison: Plain MC / QMC / Antithetic / Static CV ===\n";
    std::cout << "N=" << NSIM << " sims, nsteps=" << NSTEPS
              << " (exact terminal law, QMC-safe prime bases)\n";
    std::cout << "S0=K=" << SPOT << ", r=" << RATE << ", T=" << T << "\n";

    test_1d();
    test_2d_equal();
    test_2d_diffvols();
    test_3d_equal();

    std::cout << "\n";
    if (failures == 0) std::cout << "All tests PASSED.\n";
    else               std::cout << failures << " test(s) FAILED.\n";
    return failures == 0 ? 0 : 1;
}
