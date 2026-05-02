# Pricer-Lib — Implementation Notes

## What was built

This library is a C++ Monte Carlo option pricer. It is structured in four layers:

```
RandomGenerator  →  SDE (stochastic processes)  →  Payoffs  →  Pricer
```

---

## 1. Mac Build Support

The project was originally Windows/MSVC only. The following changes make it compile on Mac with clang++:

**`RandomGenerator/framework.h`** — guarded `<windows.h>` with `#ifdef _WIN32`.

**`RandomGenerator/RandomGenerator.h`** — added `#include <cstddef>` so that `size_t` is available on Mac (Windows provided it transitively via `<windows.h>`).

**`build.sh`** — shell script to compile and run all three test executables:
```bash
./build.sh        # build + run all tests
./build.sh clean  # remove compiled binaries
```

Manual compilation (full command):
```bash
clang++ -std=c++17 -O2 -I. \
  RandomGenerator/*.cpp SDE/*.cpp Payoffs/*.cpp Pricer/*.cpp \
  tests/test_main.cpp -o tests/run_tests
```

---

## 2. Architecture Overview

### RandomGenerator layer
Abstract base `RandomGenerator` with a single `Generate() → double` method.

- `UniformGenerator` → `LinearCongruential`, `EcuyerCombined`, **`HaltonGenerator`**
- `ContinuousGenerator` → `Normal` (Box-Muller / CLT / rejection), `Exponential`

### SDE layer
`RandomProcess` holds a `RandomGenerator*` and simulates paths via `Simulate(t0, T, nbSteps)`.

- `BrownianD1` — 1D Brownian motion
- `BrownianND` — N-dimensional correlated Brownian (Cholesky decomposition)
- `BSEuler1D` — Euler-Maruyama scheme on 1 asset
- `BSMilstein1D` — Milstein scheme on 1 asset (higher-order correction: `+0.5σ²S(dW²-dt)`)
- `BSMilstein2D` — Milstein on 2 correlated assets
- **`BSMilsteinND`** — Milstein on N correlated assets (new)

### Payoffs layer
`PayOff` takes a `vector<SinglePath*>` and returns a `double`.

- `EuropeanCallPayoff` — `max(S_T - K, 0)` on path[0]
- `EuroCallBasketPayOff` — `max(Σ w_i * S_T^{(i)} - K, 0)` weighted average over N assets

### Pricer layer
`MCPricer` runs simulations and returns `PriceResult{double price, double confidenceInterval}`.

- `EuropeanMCPricer` — standard MC, discounted average payoff
- `BermudanPricer` — Longstaff-Schwartz: polynomial regression for early exercise value; uses Eigen

---

## 3. BSMilsteinND — N-Asset Milstein

**Motivation:** `BSMilstein2D` was hardcoded for exactly 2 assets. `BSMilsteinND` generalises it to any N.

**Implementation** (`SDE/BSMilsteinND.h/.cpp`):
- Stores `vector<double> Spots, Vols` and an `N×N` correlation matrix
- `Simulate()` creates a `BrownianND(N, correlationMatrix)` internally
- Applies the Milstein correction per asset per step:
  ```
  S_i += r*S_i*dt + σ_i*S_i*dW_i + 0.5*σ_i²*S_i*(dW_i² - dt)
  ```
- Verified identical to `BSMilstein1D` for N=1 and `BSMilstein2D` for N=2

**Usage:**
```cpp
BSMilsteinND proc(&normGen,
    {100.0, 100.0, 100.0},   // spots
    0.05,                     // rate
    {0.2, 0.2, 0.2},          // vols
    {{1.0,0.3,0.3},{0.3,1.0,0.3},{0.3,0.3,1.0}}); // correlation
```

**Basket weights** are passed to `EuroCallBasketPayOff`:
```cpp
EuroCallBasketPayOff payoff(100.0, {1.0/3.0, 1.0/3.0, 1.0/3.0});  // equal weight
EuroCallBasketPayOff payoff(100.0, {0.7, 0.3});                     // custom weight
```

---

## 4. Bermudan Pricer (Longstaff-Schwartz)

**`BermudanPricer`** prices options with early exercise rights at specified dates.

**Algorithm:**
1. Forward pass: simulate all paths
2. Backward pass (from last exercise date to first):
   - At each date, compute intrinsic value (immediate exercise)
   - Fit a polynomial regression (degree 3) of continuation value on current asset values
   - Exercise if intrinsic > estimated continuation value

**Multi-asset support:** The polynomial basis includes cross terms for all dimensions.

**Usage:**
```cpp
BermudanPricer pricer(&proc, &payoff, 0.05, 1.0, 10000, 252,
                      {0.25, 0.5, 0.75, 1.0});  // exercise dates
```

---

## 5. Quasi-Random Numbers (Halton Sequences)

### Theory

The **Van der Corput sequence** in base b maps integer n to a number in (0,1) by reflecting the base-b digits of n across the decimal point:

```
n = 1 → Φ_2(1) = 0.1 (binary) = 0.5
n = 2 → Φ_2(2) = 0.01         = 0.25
n = 3 → Φ_2(3) = 0.11         = 0.75
n = 4 → Φ_2(4) = 0.001        = 0.125
```

The **Halton sequence** extends this to d dimensions using d different prime bases:
- Dimension 1: Van der Corput base 2
- Dimension 2: Van der Corput base 3
- Dimension k: Van der Corput base p_k (k-th prime)

These sequences have **low discrepancy**: they fill [0,1]^d more uniformly than pseudo-random numbers, which tend to cluster and leave gaps.

### Implementation (`RandomGenerator/HaltonGenerator`)

```cpp
// Single base (Van der Corput)
HaltonGenerator gen(2);        // base-2 sequence: 0.5, 0.25, 0.75, 0.125, ...

// Multiple bases (cycling)
HaltonGenerator gen({2, 3});   // cycles: Φ_2(1), Φ_3(1), Φ_2(2), Φ_3(2), ...

// Helper to generate prime bases
auto bases = HaltonGenerator::firstNPrimes(504);  // first 504 primes
```

`HaltonGenerator` inherits from `UniformGenerator`, so it can replace `LinearCongruential` anywhere a uniform source is needed — including as the uniform input for `Normal` (Box-Muller).

---

## 6. The Critical QMC Implementation Mistake (and fix)

### What went wrong first

The first attempt used `{2, 3}` Halton for all configurations, including 252-step simulations:

```
Result: pseudo price = 10.55, quasi price = 2.91  ← WRONG
        CI reduction: 90%                          ← MISLEADING
```

This looked like excellent variance reduction, but the price was completely wrong.

**Root cause:** With only 2 bases for 252 time steps, each simulation uses 252 *consecutive rows* of the same 2D Halton sequence. These consecutive rows are anti-correlated by design (that is what low-discrepancy means in 2D). The sum of 252 anti-correlated increments nearly cancels out:

```
Var(W_T) with 2 bases, 252 steps = 0.005   (should be 1.0)
```

`W_T ≈ 0` for every simulation → all paths cluster near S_0 × e^{drift} ≈ 103 → payoff ≈ 3. The estimator was consistently wrong with very low variance — a biased estimator with fake precision.

### The fix

**Rule:** the number of prime bases must equal `nsteps × normals_per_step × 2`.

Each time step must have its own pair of prime bases. Then simulation i uses row i of the full K-dimensional Halton design, and different steps use different primes → increments within a path are approximately independent.

```
# uniforms per simulation:
1D, 1 step:   1 step × 1 normal × 2 uniforms = 2  → {2, 3}
2D, 1 step:   1 step × 2 normals × 2 uniforms = 4  → {2, 3, 5, 7}
3D, 1 step:   1 step × 3 normals × 2 uniforms = 6  → {2, 3, 5, 7, 11, 13}
1D, 4 steps:  4 steps × 1 normal × 2 uniforms = 8  → {2, 3, 5, 7, 11, 13, 17, 19}
2D, 4 steps:  4 steps × 2 normals × 2 uniforms = 16 → firstNPrimes(16)
```

**Diagnostic confirmation:**

| Configuration | Var(W_T) | Expected |
|---|---|---|
| 2 bases, 252 steps | 0.005 | 1.0 — **broken** |
| 2 bases, 1 step | 1.000 | 1.0 — correct |
| 8 bases, 4 steps | 1.003 | 1.0 — correct |

### Why NSTEPS must stay small

Using many steps requires many prime bases. Large primes (e.g., p = 3583, the 504th prime) cause extreme Box-Muller values for the first simulation: `Φ_{3583}(1) = 1/3583 ≈ 0.00028` → `sqrt(-2 ln(0.00028)) ≈ 4.0`. Attempting 504 bases gave prices of ~4000 due to these outliers dominating the average with N=10000 simulations.

The practical limit with N=10000: use primes up to ~53 (16th prime) → at most 16 base cycling → at most 4 steps per 1D simulation, or 2 steps per 2D simulation.

---

## 7. QMC Test Results

The `tests/run_tests_qmc` executable compares pseudo-random vs Halton quasi-random. Both use the same number of steps for a fair comparison.

**Parameters:** N = 10,000 simulations, S₀ = K = 100, r = 5%, σ = 20%, T = 1 year.

### European options (1 step)

| Configuration | Pseudo price | Pseudo variance | Quasi price | Quasi variance | Reduction |
|---|---|---|---|---|---|
| BSEuler1D | 10.24 | 163.4 | 10.20 | 160.9 | **1.5%** |
| BSMilstein1D | 10.10 | 199.8 | 10.05 | 195.6 | **2.1%** |
| BSMilstein2D basket | 8.77 | 130.6 | 8.67 | 129.2 | **1.1%** |
| BSMilsteinND (N=3) basket | 8.19 | 108.7 | 8.14 | 107.2 | **1.3%** |

### Bermudan options (4 steps, exercise at 0.25, 0.5, 0.75, 1.0)

| Configuration | Pseudo price | Pseudo variance | Quasi price | Quasi variance | Reduction |
|---|---|---|---|---|---|
| BSMilstein1D | 10.51 | 213.5 | 10.35 | 209.7 | **1.8%** |
| BSMilsteinND (N=2) basket | 8.90 | 135.4 | 8.95 | 143.1 | -5.6% (noise) |

### Why the reduction is modest (1–5%)

The theoretical advantage of QMC is `O(log^d(N) / N)` error vs `O(1/sqrt(N))` for MC. For N=10,000 and d=2 this gives roughly a 13% error reduction, corresponding to ~25% variance reduction.

In practice, for a European call payoff, the reduction is much smaller for two reasons:

1. **Kink in the payoff.** The function `max(S_T - K, 0)` has a discontinuous first derivative at S_T = K. The Koksma-Hlawka bound (which guarantees QMC improvement) requires the integrand to have bounded variation. A kinked function has infinite total variation in the classical sense, eliminating the theoretical guarantee.

2. **Box-Muller transformation.** Box-Muller maps 2 uniforms to 1 normal via a non-linear, non-monotone transformation `sqrt(-2 ln U₁) × cos(2π U₂)`. This transformation partially destroys the low-discrepancy structure of the Halton sequence in [0,1]², so the normals produced don't fully inherit the equidistribution property. Using the inverse CDF (probit function) instead of Box-Muller would preserve equidistribution, but gives similar empirical results for this payoff.

**When QMC works better:** smooth integrals (e.g., digital options, path-average Asian options), low effective dimension problems, or when combined with Brownian Bridge construction which concentrates QMC advantage on the most important dimensions.

---

## 8. Antithetic Variables

**`AntitheticMCPricer`** reduces variance by pairing each simulation with its antithetic counterpart.

**Algorithm:**
For each of N pairs:
1. Simulate a path using Z₁, Z₂, ..., Z_k (the normal draws) → payoff Y
2. Simulate the antithetic path using −Z₁, −Z₂, ..., −Z_k → payoff Y_anti
3. Estimator for this pair: Ỹ = (Y + Y_anti) / 2

**Why it works:**
```
Var(Ỹ) = ¼ · (Var(Y) + Var(Y_anti) + 2·Cov(Y, Y_anti))
       = ½ · Var(Y) · (1 + ρ)
```
where ρ = Corr(Y(Z), Y(−Z)). For a monotone increasing payoff (European call), Y(Z) and Y(−Z) are negatively correlated (ρ < 0), so Var(Ỹ) < ½·Var(Y). In practice ρ ≈ −0.5 to −0.6 for ATM calls → 65–80% variance reduction.

**Implementation:**

`AntitheticNormal` wraps a `Normal` generator and acts as the `RandomGenerator*` for the process:
- **Normal pass** (`antithetic=false`): calls the inner `Normal`, stores each Z in a buffer, returns Z
- **Antithetic pass** (`antithetic=true`): replays the buffer returning −Z values

`AntitheticMCPricer` controls the mode between passes:
```cpp
// For each of nbSim pairs:
antiNorm.ResetBuffer();   antiNorm.SetAntithetic(false);
proc.Simulate(...);       // draws Z values → stored in buffer
double y = discount * payoff(paths);

antiNorm.ResetIndex();    antiNorm.SetAntithetic(true);
proc.Simulate(...);       // replays -Z values from buffer
double y_anti = discount * payoff(paths);

avg = 0.5 * (y + y_anti);
```

**Observed variance reduction (N=10,000, nsteps=50):**

| Configuration | Plain MC var | Antithetic var | Reduction |
|---|---|---|---|
| 1D call σ=0.2 | 225 | 56 | **75%** |
| 1D call σ=0.4 | 1031 | 343 | **67%** |
| 2D basket equal w | 146 | 32 | **78%** |
| 2D basket w={0.7,0.3} | 157 | 36 | **77%** |
| 3D basket equal w | 120 | 24 | **80%** |
| 3D basket w={0.6,0.3,0.1} | 137 | 30 | **79%** |

Higher correlation between assets → more effective antithetic (reduction improves for baskets vs single asset).

**Comparison of all three variance reduction techniques:**

| Technique | Typical reduction | Works for | Requirement |
|---|---|---|---|
| QMC (Halton) | 1–3% | Any payoff | nsteps must be small |
| Antithetic variables | 65–80% | Monotone payoffs | Payoff monotone in Z |
| Static control variate | 95–99% | Basket calls | Needs analytical E[CV] |

---

## 9. Test Executables

| Executable | Source | What it tests |
|---|---|---|
| `tests/run_tests` | `tests/test_main.cpp` | LinearCongruential, Normal, BSMilstein1D/2D, payoffs, EuropeanMCPricer, BermudanPricer |
| `tests/run_tests_nd` | `tests/test_milstein_nd.cpp` | BSMilsteinND: N=1,2,3 path properties, identity with BSMilstein1D/2D, basket pricing, Bermudan N=2,3 |
| `tests/run_tests_qmc` | `tests/test_qmc.cpp` | Quasi vs pseudo variance comparison (nsteps=1) |
| `tests/run_tests_varred_compare` | `tests/test_varred_compare.cpp` | Static CV vs plain MC (nsteps=50, 1D/2D/3D, various weights and vols) |
| `tests/run_tests_antithetic` | `tests/test_antithetic.cpp` | Antithetic variables vs plain MC (nsteps=50, 1D/2D/3D, various configs) |

Run all:
```bash
./build.sh
```

Run one:
```bash
./tests/run_tests_antithetic
```

---

## 10. Files Added / Modified

| File | Change |
|---|---|
| `RandomGenerator/framework.h` | Guarded `<windows.h>` with `#ifdef _WIN32` |
| `RandomGenerator/RandomGenerator.h` | Added `#include <cstddef>` |
| `RandomGenerator/HaltonGenerator.h` | New class (+ `firstNPrimes` static helper) |
| `RandomGenerator/HaltonGenerator.cpp` | Radical inverse implementation, cycling over bases |
| `SDE/BlackScholesND.h/.cpp` | Abstract base for N-asset Black-Scholes |
| `SDE/BSMilsteinND.h/.cpp` | N-asset Milstein (reuses BrownianND + Cholesky) |
| `Payoffs/EuroCallBasketPayOff.h/.cpp` | N-weight basket call payoff |
| `tests/test_main.cpp` | Core test suite (17 tests) |
| `tests/test_milstein_nd.cpp` | BSMilsteinND test suite (17 tests) |
| `tests/test_qmc.cpp` | QMC variance comparison (12 tests) |
| `RandomGenerator/AntitheticNormal.h/.cpp` | Antithetic normal generator (wraps Normal, buffers/negates Z) |
| `Pricer/AntitheticMCPricer.h/.cpp` | Antithetic variables pricer |
| `VarRed/BasketGeomControlVariate.cpp` | Fixed: sqrt(varGeom) and drift adjustment in AnalyticalExpectation |
| `VarRed/BSClosedForm.cpp` | Fixed: added `<algorithm>` for std::max on Mac |
| `Pricer/VarRedMCPricer.cpp` | Fixed: discount cvValue to match scale of cvExpectation |
| `tests/test_varred_compare.cpp` | Static CV vs plain MC (7 basket configs, nsteps=50) |
| `tests/test_antithetic.cpp` | Antithetic vs plain MC (8 configs, nsteps=50) |
| `build.sh` | Mac build script |
| `README.md` | Updated with Mac build instructions |
