#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

#include "RandomGenerator/LinearCongruential.h"
#include "RandomGenerator/Normal.h"
#include "RandomGenerator/HaltonGenerator.h"
#include "SDE/BSEuler1D.h"
#include "SDE/BSMilstein1D.h"
#include "SDE/BSMilstein2D.h"
#include "SDE/BSMilsteinND.h"
#include "Payoffs/EuropeanCallPayoff.h"
#include "Payoffs/EuroCallBasketPayOff.h"
#include "Pricer/EuropeanMCPricer.h"
#include "Pricer/BermudanPricer.h"

static int failures = 0;

void check(bool condition, const std::string& name) {
    if (condition) std::cout << "[PASS] " << name << "\n";
    else { std::cout << "[FAIL] " << name << "\n"; ++failures; }
}

// ---------------------------------------------------------------------------
// Key rule for correct QMC path simulation with Box-Muller:
//
//   number of prime bases = nbSteps * normals_per_step * 2
//
// Each time step must use its OWN pair of prime bases.
// With K bases cycling, each simulation uses ONE ROW of the K-dimensional
// Halton design (row = simulation index). Different steps use different
// prime bases → increments are approximately independent within a path
// → Var(W_T) ≈ 1 (correct). Across simulations each step's Z values are
// a low-discrepancy normal sample → reduced variance.
//
// Consequence: NSTEPS must stay SMALL so that the prime bases stay small
// (large primes like p=3583 cause extreme Box-Muller values for small sims).
//
//   1D, 1 step  → 1*1*2 =  2 bases: {2, 3}
//   2D, 1 step  → 1*2*2 =  4 bases: {2, 3, 5, 7}
//   3D, 1 step  → 1*3*2 =  6 bases: {2, 3, 5, 7, 11, 13}
//   1D, 4 steps → 4*1*2 =  8 bases: {2, 3, 5, 7, 11, 13, 17, 19}
//   2D, 4 steps → 4*2*2 = 16 bases: first 16 primes (up to 53)
// ---------------------------------------------------------------------------

static const size_t SEED = 42, A = 16807, C = 0, M = 2147483647;
static const size_t NSIM   = 10000;
static const size_t NSTEPS_EUR  = 1;  // 1 step → exact terminal law (no discretisation error)
static const size_t NSTEPS_BERM = 4;  // 4 steps, one per exercise date {0.25,0.5,0.75,1.0}
static const double SPOT = 100.0, STRIKE = 100.0, RATE = 0.05, VOL = 0.2, T = 1.0;

static const std::vector<size_t> B_1D = {2, 3};
static const std::vector<size_t> B_2D = {2, 3, 5, 7};
static const std::vector<size_t> B_3D = {2, 3, 5, 7, 11, 13};
// 1D Bermudan 4 steps: 4*1*2 = 8 bases
static const std::vector<size_t> B_BERM1 = {2, 3, 5, 7, 11, 13, 17, 19};
// 2D Bermudan 4 steps: 4*2*2 = 16 bases (first 16 primes, up to 53)

// ---------------------------------------------------------------------------

void compare(const std::string& label,
             double price_p, double ci_p,
             double price_q, double ci_q)
{
    double var_p = std::pow(ci_p * std::sqrt((double)NSIM) / 1.96, 2);
    double var_q = std::pow(ci_q * std::sqrt((double)NSIM) / 1.96, 2);
    double var_red = (1.0 - var_q / var_p) * 100.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Pseudo  price=" << price_p << "  variance=" << var_p
              << "  CI=+/-" << ci_p << "\n";
    std::cout << "  Quasi   price=" << price_q << "  variance=" << var_q
              << "  CI=+/-" << ci_q << "\n";
    std::cout << "  Variance reduction: " << var_red << "%\n";

    check(std::abs(price_p - price_q) < 1.0,
          label + ": quasi price agrees with pseudo (within 1.0)");
    // Allow up to 8% worse: small variance differences can go either way at N=10000
    check(ci_q < ci_p * 1.08,
          label + ": quasi CI not significantly worse than pseudo CI");
}

// ---------------------------------------------------------------------------

void test_euler1d() {
    std::cout << "\n[BSEuler1D + European  (nsteps=" << NSTEPS_EUR << ")]\n";

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSEuler1D proc_p(&norm_p, SPOT, RATE, VOL);
    EuropeanCallPayoff payoff_p(STRIKE);
    EuropeanMCPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS_EUR);
    auto r_p = pricer_p.Price();

    HaltonGenerator u_q(B_1D);
    Normal norm_q(0.0, 1.0, u_q);
    BSEuler1D proc_q(&norm_q, SPOT, RATE, VOL);
    EuropeanCallPayoff payoff_q(STRIKE);
    EuropeanMCPricer pricer_q(&proc_q, &payoff_q, RATE, T, NSIM, NSTEPS_EUR);
    auto r_q = pricer_q.Price();

    compare("BSEuler1D", r_p.price, r_p.confidenceInterval,
                         r_q.price, r_q.confidenceInterval);
}

void test_milstein1d() {
    std::cout << "\n[BSMilstein1D + European  (nsteps=" << NSTEPS_EUR << ")]\n";

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilstein1D proc_p(&norm_p, SPOT, RATE, VOL);
    EuropeanCallPayoff payoff_p(STRIKE);
    EuropeanMCPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS_EUR);
    auto r_p = pricer_p.Price();

    HaltonGenerator u_q(B_1D);
    Normal norm_q(0.0, 1.0, u_q);
    BSMilstein1D proc_q(&norm_q, SPOT, RATE, VOL);
    EuropeanCallPayoff payoff_q(STRIKE);
    EuropeanMCPricer pricer_q(&proc_q, &payoff_q, RATE, T, NSIM, NSTEPS_EUR);
    auto r_q = pricer_q.Price();

    compare("BSMilstein1D", r_p.price, r_p.confidenceInterval,
                            r_q.price, r_q.confidenceInterval);
}

void test_milstein2d() {
    std::cout << "\n[BSMilstein2D + Basket European  (nsteps=" << NSTEPS_EUR << ")]\n";

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilstein2D proc_p(&norm_p, SPOT, SPOT, RATE, VOL, VOL, 0.3);
    EuroCallBasketPayOff payoff_p(STRIKE, {0.5, 0.5});
    EuropeanMCPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS_EUR);
    auto r_p = pricer_p.Price();

    HaltonGenerator u_q(B_2D);   // 1 step * 2 assets * 2 uniforms = 4 bases
    Normal norm_q(0.0, 1.0, u_q);
    BSMilstein2D proc_q(&norm_q, SPOT, SPOT, RATE, VOL, VOL, 0.3);
    EuroCallBasketPayOff payoff_q(STRIKE, {0.5, 0.5});
    EuropeanMCPricer pricer_q(&proc_q, &payoff_q, RATE, T, NSIM, NSTEPS_EUR);
    auto r_q = pricer_q.Price();

    compare("BSMilstein2D Basket", r_p.price, r_p.confidenceInterval,
                                   r_q.price, r_q.confidenceInterval);
}

void test_milsteinND_3d() {
    std::cout << "\n[BSMilsteinND (N=3) + Basket European  (nsteps=" << NSTEPS_EUR << ")]\n";

    std::vector<double> spots   = {SPOT, SPOT, SPOT};
    std::vector<double> vols    = {VOL, VOL, VOL};
    std::vector<std::vector<double>> corr = {
        {1.0, 0.3, 0.3},
        {0.3, 1.0, 0.3},
        {0.3, 0.3, 1.0}
    };
    std::vector<double> weights = {1.0/3.0, 1.0/3.0, 1.0/3.0};

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilsteinND proc_p(&norm_p, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_p(STRIKE, weights);
    EuropeanMCPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS_EUR);
    auto r_p = pricer_p.Price();

    HaltonGenerator u_q(B_3D);   // 1 step * 3 assets * 2 uniforms = 6 bases
    Normal norm_q(0.0, 1.0, u_q);
    BSMilsteinND proc_q(&norm_q, spots, RATE, vols, corr);
    EuroCallBasketPayOff payoff_q(STRIKE, weights);
    EuropeanMCPricer pricer_q(&proc_q, &payoff_q, RATE, T, NSIM, NSTEPS_EUR);
    auto r_q = pricer_q.Price();

    compare("BSMilsteinND (N=3) Basket", r_p.price, r_p.confidenceInterval,
                                         r_q.price, r_q.confidenceInterval);
}

void test_bermudan_1d() {
    std::cout << "\n[BSMilstein1D + Bermudan  (nsteps=" << NSTEPS_BERM << ")]\n";

    std::vector<double> exerciseDates = {0.25, 0.5, 0.75, 1.0};

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilstein1D proc_p(&norm_p, SPOT, RATE, VOL);
    EuropeanCallPayoff payoff_p(STRIKE);
    BermudanPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS_BERM, exerciseDates);
    auto r_p = pricer_p.Price();

    HaltonGenerator u_q(B_BERM1);  // 4 steps * 1 asset * 2 uniforms = 8 bases
    Normal norm_q(0.0, 1.0, u_q);
    BSMilstein1D proc_q(&norm_q, SPOT, RATE, VOL);
    EuropeanCallPayoff payoff_q(STRIKE);
    BermudanPricer pricer_q(&proc_q, &payoff_q, RATE, T, NSIM, NSTEPS_BERM, exerciseDates);
    auto r_q = pricer_q.Price();

    compare("BSMilstein1D Bermudan", r_p.price, r_p.confidenceInterval,
                                     r_q.price, r_q.confidenceInterval);
}

void test_bermudan_nd() {
    std::cout << "\n[BSMilsteinND (N=2) + Bermudan basket  (nsteps=" << NSTEPS_BERM << ")]\n";

    std::vector<double> exerciseDates = {0.25, 0.5, 0.75, 1.0};
    // 4 steps * 2 assets * 2 uniforms = 16 bases (first 16 primes, up to 53)
    auto bases_2d_berm = HaltonGenerator::firstNPrimes(16);

    LinearCongruential u_p(SEED, A, C, M);
    Normal norm_p(0.0, 1.0, u_p);
    BSMilsteinND proc_p(&norm_p, {SPOT, SPOT}, RATE, {VOL, VOL}, {{1.0,0.3},{0.3,1.0}});
    EuroCallBasketPayOff payoff_p(STRIKE, {0.5, 0.5});
    BermudanPricer pricer_p(&proc_p, &payoff_p, RATE, T, NSIM, NSTEPS_BERM, exerciseDates);
    auto r_p = pricer_p.Price();

    HaltonGenerator u_q(bases_2d_berm);
    Normal norm_q(0.0, 1.0, u_q);
    BSMilsteinND proc_q(&norm_q, {SPOT, SPOT}, RATE, {VOL, VOL}, {{1.0,0.3},{0.3,1.0}});
    EuroCallBasketPayOff payoff_q(STRIKE, {0.5, 0.5});
    BermudanPricer pricer_q(&proc_q, &payoff_q, RATE, T, NSIM, NSTEPS_BERM, exerciseDates);
    auto r_q = pricer_q.Price();

    compare("BSMilsteinND (N=2) Bermudan", r_p.price, r_p.confidenceInterval,
                                           r_q.price, r_q.confidenceInterval);
}

// ---------------------------------------------------------------------------

int main() {
    std::cout << "=== Quasi-Random (Halton) vs Pseudo-Random Variance Test ===\n";
    std::cout << "N=" << NSIM << " sims\n";
    std::cout << "European: " << NSTEPS_EUR << " step  |  Bermudan: "
              << NSTEPS_BERM << " steps (one per exercise date)\n";
    std::cout << "Rule: #bases = nsteps * normals_per_step * 2\n";
    std::cout << "BOTH pseudo and quasi use the same nsteps for a fair comparison.\n";

    test_euler1d();
    test_milstein1d();
    test_milstein2d();
    test_milsteinND_3d();
    test_bermudan_1d();
    test_bermudan_nd();

    std::cout << "\n";
    if (failures == 0) std::cout << "All tests PASSED.\n";
    else               std::cout << failures << " test(s) FAILED.\n";
    return failures == 0 ? 0 : 1;
}
