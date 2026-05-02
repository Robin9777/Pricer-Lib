#!/usr/bin/env bash
set -e

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

OUT_TESTS="$ROOT/tests/run_tests"
OUT_ND="$ROOT/tests/run_tests_nd"
OUT_QMC="$ROOT/tests/run_tests_qmc"
OUT_VR="$ROOT/tests/run_tests_varred_compare"
OUT_AT="$ROOT/tests/run_tests_antithetic"
OUT_FC="$ROOT/tests/run_tests_full_comparison"

if [ "${1:-}" = "clean" ]; then
    rm -f "$OUT_TESTS" "$OUT_ND" "$OUT_QMC" "$OUT_VR" "$OUT_AT" "$OUT_FC"
    echo "Cleaned."
    exit 0
fi

# ---------------------------------------------------------------------------
# Library sources (shared by both executables)
# ---------------------------------------------------------------------------
LIB_SOURCES=(
    # RandomGenerator
    RandomGenerator/RandomGenerator.cpp
    RandomGenerator/UniformGenerator.cpp
    RandomGenerator/PseudoGenerator.cpp
    RandomGenerator/LinearCongruential.cpp
    RandomGenerator/EcuyerCombined.cpp
    RandomGenerator/ContinuousGenerator.cpp
    RandomGenerator/Normal.cpp
    RandomGenerator/Exponential.cpp
    RandomGenerator/DiscreteGenerator.cpp
    RandomGenerator/Bernoulli.cpp
    RandomGenerator/Binomial.cpp
    RandomGenerator/Poisson.cpp
    RandomGenerator/FiniteSet.cpp
    RandomGenerator/HeadTail.cpp
    RandomGenerator/HaltonGenerator.cpp
    RandomGenerator/AntitheticNormal.cpp

    # SDE
    SDE/SinglePath.cpp
    SDE/RandomProcess.cpp
    SDE/BlackScholes1D.cpp
    SDE/BlackScholes2D.cpp
    SDE/BSEuler1D.cpp
    SDE/BSMilstein1D.cpp
    SDE/BSMilstein2D.cpp
    SDE/BrownianD1.cpp
    SDE/BrownianND.cpp
    SDE/Heston.cpp
    SDE/BlackScholesND.cpp
    SDE/BSMilsteinND.cpp

    # Payoffs
    Payoffs/PayOff.cpp
    Payoffs/EuropeanCallPayoff.cpp
    Payoffs/EuroCallBasketPayOff.cpp

    # Pricer
    Pricer/MCPricer.cpp
    Pricer/EuropeanMCPricer.cpp
    Pricer/BermudanPricer.cpp
    Pricer/VarRedMCPricer.cpp
    Pricer/AntitheticMCPricer.cpp

    # VarRed
    VarRed/BSClosedForm.cpp
    VarRed/BasketGeomControlVariate.cpp
    VarRed/ControlVariate.cpp
)

FLAGS="-std=c++17 -O2 -I."

# ---------------------------------------------------------------------------
# Build 1: main test suite -> tests/run_tests
# ---------------------------------------------------------------------------
echo "Compiling tests/run_tests..."
clang++ $FLAGS "${LIB_SOURCES[@]}" tests/test_main.cpp -o "$OUT_TESTS"

# ---------------------------------------------------------------------------
# Build 2: BSMilsteinND tests -> tests/run_tests_nd
# ---------------------------------------------------------------------------
echo "Compiling tests/run_tests_nd..."
clang++ $FLAGS "${LIB_SOURCES[@]}" tests/test_milstein_nd.cpp -o "$OUT_ND"

# ---------------------------------------------------------------------------
# Run both
# ---------------------------------------------------------------------------
echo ""
echo "--- tests/run_tests ---"
"$OUT_TESTS"

echo ""
echo "--- tests/run_tests_nd ---"
"$OUT_ND"

# ---------------------------------------------------------------------------
# Build 3: QMC variance reduction tests -> tests/run_tests_qmc
# ---------------------------------------------------------------------------
echo "Compiling tests/run_tests_qmc..."
clang++ $FLAGS "${LIB_SOURCES[@]}" tests/test_qmc.cpp -o "$OUT_QMC"

echo ""
echo "--- tests/run_tests_qmc ---"
"$OUT_QMC"

# ---------------------------------------------------------------------------
# Build 4: Variance reduction comparison -> tests/run_tests_varred_compare
# ---------------------------------------------------------------------------
echo "Compiling tests/run_tests_varred_compare..."
clang++ $FLAGS "${LIB_SOURCES[@]}" tests/test_varred_compare.cpp -o "$OUT_VR"

echo ""
echo "--- tests/run_tests_varred_compare ---"
"$OUT_VR"

# ---------------------------------------------------------------------------
# Build 5: Antithetic variables -> tests/run_tests_antithetic
# ---------------------------------------------------------------------------
echo "Compiling tests/run_tests_antithetic..."
clang++ $FLAGS "${LIB_SOURCES[@]}" tests/test_antithetic.cpp -o "$OUT_AT"

echo ""
echo "--- tests/run_tests_antithetic ---"
"$OUT_AT"

# ---------------------------------------------------------------------------
# Build 6: Full variance reduction comparison -> tests/run_tests_full_comparison
# ---------------------------------------------------------------------------
echo "Compiling tests/run_tests_full_comparison..."
clang++ $FLAGS "${LIB_SOURCES[@]}" tests/test_full_comparison.cpp -o "$OUT_FC"

echo ""
echo "--- tests/run_tests_full_comparison ---"
"$OUT_FC"
