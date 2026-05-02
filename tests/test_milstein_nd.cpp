#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "RandomGenerator/LinearCongruential.h"
#include "RandomGenerator/Normal.h"
#include "SDE/BSMilstein1D.h"
#include "SDE/BSMilstein2D.h"
#include "SDE/BSMilsteinND.h"
#include "Payoffs/EuroCallBasketPayOff.h"
#include "Pricer/EuropeanMCPricer.h"
#include "Pricer/BermudanPricer.h"

static int failures = 0;

void check(bool condition, const std::string& name) {
    if (condition)
        std::cout << "[PASS] " << name << "\n";
    else {
        std::cout << "[FAIL] " << name << "\n";
        ++failures;
    }
}

// ---------------------------------------------------------------------------
// helpers to build fresh generators with the same seed
// ---------------------------------------------------------------------------
static const size_t SEED = 42;
static const size_t A    = 16807;
static const size_t C    = 0;
static const size_t M    = 2147483647;

// ---------------------------------------------------------------------------

void test_n1_path_properties() {
    LinearCongruential u(SEED, A, C, M);
    Normal norm(0.0, 1.0, u);

    std::vector<double> spots = {100.0};
    std::vector<double> vols  = {0.2};
    std::vector<std::vector<double>> corr = {{1.0}};

    BSMilsteinND proc(&norm, spots, 0.05, vols, corr);
    proc.Simulate(0.0, 1.0, 50);

    const auto& path = proc.GetPath(0)->GetAllValues();
    check(path.size() == 51, "N=1: path length = 51");

    bool all_pos = true;
    for (double v : path) if (v <= 0.0) { all_pos = false; break; }
    check(all_pos, "N=1: all values > 0");
}

void test_n1_matches_milstein1d() {
    // BSMilstein1D and BSMilsteinND(N=1) must produce the exact same path
    // because both delegate to BrownianD1 / BrownianND(1) which draw one
    // Generate() call per step.

    LinearCongruential u1(SEED, A, C, M);
    Normal norm1(0.0, 1.0, u1);
    BSMilstein1D ref(&norm1, 100.0, 0.05, 0.2);
    ref.Simulate(0.0, 1.0, 50);
    const auto& ref_path = ref.GetPath(0)->GetAllValues();

    LinearCongruential u2(SEED, A, C, M);
    Normal norm2(0.0, 1.0, u2);
    BSMilsteinND nd(&norm2, {100.0}, 0.05, {0.2}, {{1.0}});
    nd.Simulate(0.0, 1.0, 50);
    const auto& nd_path = nd.GetPath(0)->GetAllValues();

    bool same = (ref_path.size() == nd_path.size());
    for (size_t i = 0; same && i < ref_path.size(); ++i)
        if (std::abs(ref_path[i] - nd_path[i]) > 1e-12) same = false;
    check(same, "N=1: identical path to BSMilstein1D with same seed");
}

void test_n2_path_properties() {
    LinearCongruential u(SEED, A, C, M);
    Normal norm(0.0, 1.0, u);

    std::vector<double> spots = {100.0, 80.0};
    std::vector<double> vols  = {0.2, 0.25};
    std::vector<std::vector<double>> corr = {{1.0, 0.3}, {0.3, 1.0}};

    BSMilsteinND proc(&norm, spots, 0.05, vols, corr);
    proc.Simulate(0.0, 1.0, 50);

    check(proc.GetPath(0)->GetAllValues().size() == 51, "N=2: path[0] length = 51");
    check(proc.GetPath(1)->GetAllValues().size() == 51, "N=2: path[1] length = 51");

    bool all_pos = true;
    for (double v : proc.GetPath(0)->GetAllValues()) if (v <= 0.0) { all_pos = false; break; }
    for (double v : proc.GetPath(1)->GetAllValues()) if (v <= 0.0) { all_pos = false; break; }
    check(all_pos, "N=2: all values > 0");
}

void test_n2_matches_milstein2d() {
    // BSMilstein2D and BSMilsteinND(N=2) both internally create BrownianND(2)
    // with the same correlation matrix, so they must produce identical paths.
    const double rho = 0.3;

    LinearCongruential u1(SEED, A, C, M);
    Normal norm1(0.0, 1.0, u1);
    BSMilstein2D ref(&norm1, 100.0, 80.0, 0.05, 0.2, 0.25, rho);
    ref.Simulate(0.0, 1.0, 50);

    LinearCongruential u2(SEED, A, C, M);
    Normal norm2(0.0, 1.0, u2);
    BSMilsteinND nd(&norm2, {100.0, 80.0}, 0.05, {0.2, 0.25},
                    {{1.0, rho}, {rho, 1.0}});
    nd.Simulate(0.0, 1.0, 50);

    bool same0 = true, same1 = true;
    const auto& r0 = ref.GetPath(0)->GetAllValues();
    const auto& r1 = ref.GetPath(1)->GetAllValues();
    const auto& n0 = nd.GetPath(0)->GetAllValues();
    const auto& n1 = nd.GetPath(1)->GetAllValues();

    for (size_t i = 0; i < r0.size(); ++i)
        if (std::abs(r0[i] - n0[i]) > 1e-12) { same0 = false; break; }
    for (size_t i = 0; i < r1.size(); ++i)
        if (std::abs(r1[i] - n1[i]) > 1e-12) { same1 = false; break; }

    check(same0, "N=2: path[0] identical to BSMilstein2D with same seed");
    check(same1, "N=2: path[1] identical to BSMilstein2D with same seed");
}

void test_n3_path_properties() {
    LinearCongruential u(SEED, A, C, M);
    Normal norm(0.0, 1.0, u);

    std::vector<double> spots = {100.0, 80.0, 120.0};
    std::vector<double> vols  = {0.2, 0.25, 0.15};
    std::vector<std::vector<double>> corr = {
        {1.0,  0.3,  0.1},
        {0.3,  1.0,  0.2},
        {0.1,  0.2,  1.0}
    };

    BSMilsteinND proc(&norm, spots, 0.05, vols, corr);
    proc.Simulate(0.0, 1.0, 50);

    check(proc.GetPath(0)->GetAllValues().size() == 51, "N=3: path[0] length = 51");
    check(proc.GetPath(1)->GetAllValues().size() == 51, "N=3: path[1] length = 51");
    check(proc.GetPath(2)->GetAllValues().size() == 51, "N=3: path[2] length = 51");

    bool all_pos = true;
    for (int d = 0; d < 3; ++d)
        for (double v : proc.GetPath(d)->GetAllValues())
            if (v <= 0.0) { all_pos = false; break; }
    check(all_pos, "N=3: all values > 0");
}

void test_n3_basket_pricing() {
    // Equal-weight basket call on 3 assets, ATM
    LinearCongruential u(SEED, A, C, M);
    Normal norm(0.0, 1.0, u);

    std::vector<double> spots = {100.0, 100.0, 100.0};
    std::vector<double> vols  = {0.2, 0.2, 0.2};
    std::vector<std::vector<double>> corr = {
        {1.0, 0.3, 0.3},
        {0.3, 1.0, 0.3},
        {0.3, 0.3, 1.0}
    };

    BSMilsteinND proc(&norm, spots, 0.05, vols, corr);
    EuroCallBasketPayOff payoff(100.0, {1.0/3.0, 1.0/3.0, 1.0/3.0});
    EuropeanMCPricer pricer(&proc, &payoff, 0.05, 1.0, 10000, 252);

    double price = pricer.Price().price;
    // ATM call, correlation reduces vol slightly vs independent, price should
    // still be in a plausible range around the 1D BS price (~10.45)
    check(price > 5.0 && price < 15.0, "N=3 basket: ATM price in [5.0, 15.0]");
}

void test_bermudan_nd_2assets() {
    // 2-asset basket Bermudan call, equal weights, ATM, quarterly exercise dates
    // Tests that BermudanPricer works end-to-end with a multi-dimensional process:
    // - the D-dimensional polynomial basis is exercised (D=2 here)
    // - price should be >= European basket price (more rights = more value)
    // - price should be in a plausible range

    const double rate = 0.05, maturity = 1.0, strike = 100.0;
    std::vector<double> exerciseDates = {0.25, 0.5, 0.75, 1.0};

    // European basket price for reference
    LinearCongruential u1(SEED, A, C, M);
    Normal norm1(0.0, 1.0, u1);
    BSMilsteinND proc1(&norm1, {100.0, 100.0}, rate, {0.2, 0.2}, {{1.0, 0.3},{0.3, 1.0}});
    EuroCallBasketPayOff payoff1(strike, {0.5, 0.5});
    EuropeanMCPricer european(&proc1, &payoff1, rate, maturity, 5000, 252);
    double europeanPrice = european.Price().price;

    // Bermudan basket price
    LinearCongruential u2(SEED, A, C, M);
    Normal norm2(0.0, 1.0, u2);
    BSMilsteinND proc2(&norm2, {100.0, 100.0}, rate, {0.2, 0.2}, {{1.0, 0.3},{0.3, 1.0}});
    EuroCallBasketPayOff payoff2(strike, {0.5, 0.5});
    BermudanPricer bermudan(&proc2, &payoff2, rate, maturity, 5000, 252, exerciseDates);
    double bermudanPrice = bermudan.Price().price;

    check(bermudanPrice > 0.0,               "Bermudan N=2: price > 0");
    check(bermudanPrice < 25.0,              "Bermudan N=2: price in plausible range");
    check(bermudanPrice >= europeanPrice - 1.0,
          "Bermudan N=2: price >= European basket price (within MC noise)");
}

void test_bermudan_nd_3assets() {
    // 3-asset basket Bermudan — mainly checks it runs without crashing
    // and produces a sensible number with the 3D polynomial basis
    LinearCongruential u(SEED, A, C, M);
    Normal norm(0.0, 1.0, u);

    BSMilsteinND proc(&norm,
        {100.0, 100.0, 100.0}, 0.05, {0.2, 0.2, 0.2},
        {{1.0, 0.3, 0.3},{0.3, 1.0, 0.3},{0.3, 0.3, 1.0}});
    EuroCallBasketPayOff payoff(100.0, {1.0/3.0, 1.0/3.0, 1.0/3.0});
    BermudanPricer bermudan(&proc, &payoff, 0.05, 1.0, 5000, 252,
                            {0.25, 0.5, 0.75, 1.0});
    double price = bermudan.Price().price;

    check(price > 0.0 && price < 25.0, "Bermudan N=3: price in plausible range");
}

// ---------------------------------------------------------------------------

int main() {
    std::cout << "=== BSMilsteinND Test Suite ===\n\n";

    test_n1_path_properties();
    test_n1_matches_milstein1d();
    test_n2_path_properties();
    test_n2_matches_milstein2d();
    test_n3_path_properties();
    test_n3_basket_pricing();
    test_bermudan_nd_2assets();
    test_bermudan_nd_3assets();

    std::cout << "\n";
    if (failures == 0) {
        std::cout << "All tests PASSED.\n";
        return 0;
    } else {
        std::cout << failures << " test(s) FAILED.\n";
        return 1;
    }
}
