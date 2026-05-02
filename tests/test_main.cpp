#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

// RandomGenerator
#include "RandomGenerator/LinearCongruential.h"
#include "RandomGenerator/Normal.h"

// SDE
#include "SDE/BSMilstein1D.h"
#include "SDE/BSMilstein2D.h"
#include "SDE/SinglePath.h"
#include "SDE/RandomProcess.h"

// Payoffs
#include "Payoffs/EuropeanCallPayoff.h"
#include "Payoffs/EuroCallBasketPayOff.h"

// Pricer
#include "Pricer/EuropeanMCPricer.h"
#include "Pricer/BermudanPricer.h"

// ---------------------------------------------------------------------------
// Test harness
// ---------------------------------------------------------------------------
static int failures = 0;

void check(bool condition, const std::string& name) {
    if (condition) {
        std::cout << "[PASS] " << name << "\n";
    } else {
        std::cout << "[FAIL] " << name << "\n";
        ++failures;
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

void test_linear_congruential() {
    const size_t seed = 42, a = 16807, c = 0, m = 2147483647;
    LinearCongruential gen(seed, a, c, m);

    // First value: new_seed = (a*seed + c) % m, then normalized
    size_t expected_seed = (a * seed + c) % m;
    double expected_val = static_cast<double>(expected_seed) / m;
    double val = gen.Generate();
    check(std::abs(val - expected_val) < 1e-15, "LinearCongruential: first value matches formula");

    // Output in [0, 1]
    bool all_in_range = true;
    for (int i = 0; i < 1000; ++i) {
        double v = gen.Generate();
        if (v < 0.0 || v > 1.0) { all_in_range = false; break; }
    }
    check(all_in_range, "LinearCongruential: output in [0,1]");

    // Same seed => same sequence
    LinearCongruential gen2(seed, a, c, m);
    LinearCongruential gen3(seed, a, c, m);
    bool same_seq = true;
    for (int i = 0; i < 100; ++i) {
        if (gen2.Generate() != gen3.Generate()) { same_seq = false; break; }
    }
    check(same_seq, "LinearCongruential: same seed gives same sequence");
}

void test_normal() {
    const size_t seed = 123;
    const size_t a = 16807, c = 0, m = 2147483647;
    LinearCongruential unifGen(seed, a, c, m);
    Normal normGen(0.0, 1.0, unifGen);

    const int N = 10000;
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < N; ++i) {
        double v = normGen.Generate();
        sum += v;
        sum2 += v * v;
    }
    double mean = sum / N;
    double variance = sum2 / N - mean * mean;

    check(std::abs(mean) < 0.05, "Normal: mean close to 0");
    check(std::abs(variance - 1.0) < 0.1, "Normal: variance close to 1");
}

void test_bsmilstein1d() {
    const size_t seed = 42, a = 16807, c = 0, m = 2147483647;
    LinearCongruential unifGen(seed, a, c, m);
    Normal normGen(0.0, 1.0, unifGen);
    BSMilstein1D process(&normGen, 100.0, 0.05, 0.2);

    const size_t nbSteps = 50;
    process.Simulate(0.0, 1.0, nbSteps);

    const std::vector<double>& path = process.GetPath(0)->GetAllValues();
    check(path.size() == nbSteps + 1, "BSMilstein1D: path length = nbSteps+1");

    bool all_positive = true;
    for (double v : path) { if (v <= 0.0) { all_positive = false; break; } }
    check(all_positive, "BSMilstein1D: all path values > 0");

    // Determinism: same seed => same path
    LinearCongruential unifGen2(seed, a, c, m);
    Normal normGen2(0.0, 1.0, unifGen2);
    BSMilstein1D process2(&normGen2, 100.0, 0.05, 0.2);
    process2.Simulate(0.0, 1.0, nbSteps);
    const std::vector<double>& path2 = process2.GetPath(0)->GetAllValues();

    bool same_path = (path.size() == path2.size());
    if (same_path) {
        for (size_t i = 0; i < path.size(); ++i) {
            if (std::abs(path[i] - path2[i]) > 1e-12) { same_path = false; break; }
        }
    }
    check(same_path, "BSMilstein1D: same seed gives same path");
}

void test_bsmilstein2d() {
    const size_t seed = 77, a = 16807, c = 0, m = 2147483647;
    LinearCongruential unifGen(seed, a, c, m);
    Normal normGen(0.0, 1.0, unifGen);
    BSMilstein2D process(&normGen, 100.0, 100.0, 0.05, 0.2, 0.25, 0.3);

    const size_t nbSteps = 30;
    process.Simulate(0.0, 1.0, nbSteps);

    check(process.GetPath(0)->GetAllValues().size() == nbSteps + 1, "BSMilstein2D: path[0] length correct");
    check(process.GetPath(1)->GetAllValues().size() == nbSteps + 1, "BSMilstein2D: path[1] length correct");

    bool all_positive = true;
    for (double v : process.GetPath(0)->GetAllValues()) { if (v <= 0.0) { all_positive = false; break; } }
    for (double v : process.GetPath(1)->GetAllValues()) { if (v <= 0.0) { all_positive = false; break; } }
    check(all_positive, "BSMilstein2D: all values > 0");
}

void test_european_call_payoff() {
    // Build a trivial single-step path to test payoff logic
    const double strike = 100.0;
    EuropeanCallPayoff payoff(strike);

    // ITM: S_T = 120 => payoff = 20
    {
        SinglePath path(0.0, 1.0, 1);
        path.InsertValue(120.0);
        std::vector<SinglePath*> paths = {&path};
        double p = payoff(paths);
        check(std::abs(p - 20.0) < 1e-12, "EuropeanCallPayoff: ITM payoff = 20");
    }
    // OTM: S_T = 80 => payoff = 0
    {
        SinglePath path(0.0, 1.0, 1);
        path.InsertValue(80.0);
        std::vector<SinglePath*> paths = {&path};
        double p = payoff(paths);
        check(std::abs(p) < 1e-12, "EuropeanCallPayoff: OTM payoff = 0");
    }
    // ATM: S_T = 100 => payoff = 0
    {
        SinglePath path(0.0, 1.0, 1);
        path.InsertValue(100.0);
        std::vector<SinglePath*> paths = {&path};
        double p = payoff(paths);
        check(std::abs(p) < 1e-12, "EuropeanCallPayoff: ATM payoff = 0");
    }
}

void test_european_mc_pricer() {
    // ATM call, BS analytical price ~ 10.45
    const size_t seed = 12345, a = 16807, c = 0, m = 2147483647;
    LinearCongruential unifGen(seed, a, c, m);
    Normal normGen(0.0, 1.0, unifGen);
    BSMilstein1D process(&normGen, 100.0, 0.05, 0.2);
    EuropeanCallPayoff payoff(100.0);
    EuropeanMCPricer pricer(&process, &payoff, 0.05, 1.0, 10000, 252);

    double price = pricer.Price();
    check(price > 9.0 && price < 12.0, "EuropeanMCPricer: ATM call price in [9.0, 12.0]");
}

void test_bermudan_pricer() {
    const size_t seed = 99, a = 16807, c = 0, m = 2147483647;
    LinearCongruential unifGen(seed, a, c, m);
    Normal normGen(0.0, 1.0, unifGen);
    BSMilstein1D process(&normGen, 100.0, 0.05, 0.2);
    EuropeanCallPayoff payoff(100.0);
    std::vector<double> exerciseDates = {0.25, 0.5, 0.75, 1.0};
    BermudanPricer pricer(&process, &payoff, 0.05, 1.0, 5000, 252, exerciseDates);

    double price = pricer.Price();
    check(price > 0.0 && price < 25.0, "BermudanPricer: price > 0 and in plausible range");
}

void test_basket_payoff() {
    // With weight=1.0 on a single path, basket behaves like vanilla call
    const double strike = 100.0;
    EuroCallBasketPayOff basketPayoff(strike, {1.0});
    EuropeanCallPayoff vanillaPayoff(strike);

    const size_t seed = 55, a = 16807, c = 0, m = 2147483647;
    LinearCongruential unifGen1(seed, a, c, m);
    Normal normGen1(0.0, 1.0, unifGen1);
    BSMilstein1D process1(&normGen1, 100.0, 0.05, 0.2);
    process1.Simulate(0.0, 1.0, 50);
    std::vector<SinglePath*> paths1 = {process1.GetPath(0)};

    LinearCongruential unifGen2(seed, a, c, m);
    Normal normGen2(0.0, 1.0, unifGen2);
    BSMilstein1D process2(&normGen2, 100.0, 0.05, 0.2);
    process2.Simulate(0.0, 1.0, 50);
    std::vector<SinglePath*> paths2 = {process2.GetPath(0)};

    double basketVal = basketPayoff(paths1);
    double vanillaVal = vanillaPayoff(paths2);
    check(std::abs(basketVal - vanillaVal) < 1e-10,
          "EuroCallBasketPayOff: weight=1.0 matches vanilla call with same seed");
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main() {
    std::cout << "=== Pricer-Lib Test Suite ===\n\n";

    test_linear_congruential();
    test_normal();
    test_bsmilstein1d();
    test_bsmilstein2d();
    test_european_call_payoff();
    test_european_mc_pricer();
    test_bermudan_pricer();
    test_basket_payoff();

    std::cout << "\n";
    if (failures == 0) {
        std::cout << "All tests PASSED.\n";
        return 0;
    } else {
        std::cout << failures << " test(s) FAILED.\n";
        return 1;
    }
}
