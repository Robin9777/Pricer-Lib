# Pricer-Lib — Numerical Finance C++ Library

Modular C++ library for Monte Carlo and PDE option pricing.

---

## Table of Contents

1. [Why Abstract Classes?](#why-abstract-classes)
<<<<<<< HEAD
2. [Project Structure](#project-structure)
3. [Full Workflow](#full-workflow)
4. [Layer 1 — RandomGenerator](#layer-1--randomgenerator)
5. [Layer 2 — SDE](#layer-2--sde)
6. [Layer 3 — Payoffs](#layer-3--payoffs)
7. [Layer 4 — Pricer](#layer-4--pricer)
8. [PDE Module](#pde)
9. [VolSurf](#volsurf)
10. [VarRed](#varred)
11. [Eigen](#eigen)
12. [Unit Tests](#unit-tests)
13. [Compilation & Build](#compilation--build)
14. [Full Code Examples](#full-usage-examples)
=======
2. [Architecture](#architecture)
3. [Layer 1 — RandomGenerator](#layer-1--randomgenerator)
4. [Layer 2 — SDE](#layer-2--sde)
5. [Layer 3 — Payoffs](#layer-3--payoffs)
6. [Layer 4 — Pricer](#layer-4--pricer)
7. [PDE Module](#pde-module)
8. [VarRed](#varred)
9. [Unit Tests](#unit-tests)
10. [Build](#build)
11. [Demo](#demo)
12. [Code Examples](#code-examples)
>>>>>>> 663106143a3d19ab6d4b0eed0aea2db5bd2303fc

---

## Why Abstract Classes?

The pricer should not care what process or payoff you plug in. Abstract classes enforce that.

```cpp
// BAD — tightly coupled, must rewrite for every new model
class EuropeanMCPricer { BSMilstein1D* process; EuropeanCallPayoff* payoff; };

// GOOD — decoupled, works with anything
class MCPricer { RandomProcess* Process; PayOff* Payoff; };
```

At runtime C++ dispatches the right `Simulate()` or `operator()` through the vtable. Swap `BSMilstein1D` for `Heston` — the pricer code is unchanged.

| Abstract class | Pure virtual | Contract |
|---|---|---|
| `RandomGenerator` | `Generate()` | produces a `double` |
| `RandomProcess` | `Simulate(t0, T, n)` | fills its paths |
| `PayOff` | `operator()(paths)` | evaluates the payoff |
| `MCPricer` | `Price()` | returns the price |
| `R2R1Function` | `operator()(x, t)` | PDE coefficient |
| `R1R1Function` | `operator()(x)` | boundary/terminal condition |

Intermediate abstract layers like `UniformGenerator` exist so `Normal` can require a uniform source specifically — not just any `RandomGenerator`.

---

## Architecture

```
RandomGenerator  →  RandomProcess  →  SinglePath[]  →  PayOff  →  MCPricer
                                                                   (or PDEGrid)
```

```
Pricer-Lib/
├── RandomGenerator/    Layer 1 — RNG
├── SDE/                Layer 2 — Stochastic processes
├── Payoffs/            Layer 3 — Payoff functions
├── Pricer/             Layer 4 — Monte Carlo pricers
├── VarRed/             Variance reduction (control variates)
├── PDE/                Finite-difference PDE solver
├── VolSurf/            Local volatility (Dupire)
├── Eigen/              Bundled header-only linear algebra
├── tests/              Mac test runner
└── demo/               Interactive pricer + Jupyter analysis
```

---

## Layer 1 — RandomGenerator

```
RandomGenerator              (abstract — Generate())
├── UniformGenerator         (abstract)
│   └── PseudoGenerator      (adds seed)
│       ├── LinearCongruential   xn+1 = (a*xn + c) mod m
│       ├── EcuyerCombined       combines two LCGs, period ~2.3e18
│       └── HaltonGenerator      quasi-random low-discrepancy sequence
├── ContinuousGenerator      (abstract)
│   ├── Normal               Box-Muller on top of a UniformGenerator
│   ├── Exponential
│   └── AntitheticNormal     wraps Normal, alternates Z and -Z
└── DiscreteGenerator        (abstract)
    ├── Bernoulli, Binomial, Poisson, FiniteSet, HeadTail
```

`Normal` constructor: `(mu, sigma, UniformGenerator&)`. It needs a uniform source specifically because Box-Muller calls `Ugen.Generate()` twice per sample.

---

## Layer 2 — SDE

```
RandomProcess                (abstract — Simulate(t0, T, n), GetPath(d))
├── BrownianD1 / BrownianND
├── BlackScholes1D           (abstract — stores Spot/Rate/Vol)
│   ├── BSEuler1D            Euler-Maruyama
│   └── BSMilstein1D         Milstein (higher order, adds 0.5σ²S(dW²-dt))
├── BlackScholes2D / BSMilstein2D   2 correlated assets via Cholesky
├── BlackScholesND / BSMilsteinND   N assets
└── Heston                   spot + CIR variance process, correlation rho
```

`SinglePath` stores a time series as `vector<double>`. After `Simulate()`, use `GetPath(d)->GetAllValues()` to read asset `d`.

**Euler vs Milstein**: Euler is O(dt) error, Milstein is O(dt²). For a European call with 1 step, both are exact under GBM — the step size doesn't matter.

---

## Layer 3 — Payoffs

```
PayOff                       (abstract — operator()(vector<SinglePath*>))
├── EuropeanCallPayoff        max(S_T - K, 0) on path[0]
└── EuroCallBasketPayOff      max(Σ wᵢ·Sᵢ_T - K, 0) across all paths
```

---

## Layer 4 — Pricer

```
MCPricer                     (abstract — Price() → PriceResult{price, ci})
├── EuropeanMCPricer          plain MC average of discounted payoffs
├── BermudanPricer            Longstaff-Schwartz (needs exerciseDates vector)
├── VarRedMCPricer            control variate MC (needs a ControlVariate*)
└── AntitheticMCPricer        antithetic variables (needs AntitheticNormal&)
```

**Longstaff-Schwartz**: simulates all paths forward, then works backwards from the last exercise date. At each date, fits a degree-3 polynomial regression (via Eigen) of continuation values against spot states. If immediate exercise > fitted continuation → exercise.

**VarRed**: `VarRedMCPricer` takes a `ControlVariate*` that provides both a simulated sample and its known analytical price. The estimator corrects each payoff by `β*(cv_sample - cv_price)` where `β` is the optimal coefficient.

---

## PDE Module

Solves `dV/dt + a·d²V/dS² + b·dV/dS - r·V = 0` on a (spot × time) grid.

```
PDEGrid2D                    (abstract — FillNodes(), GetValue(t, S))
├── PDEGridExplicit           stable only with fine grid (dt < dS²)
└── PDEGridImplicit           unconditionally stable, tridiagonal solve per step

R2R1Function  (abstract)  →  BSVariance, BSTrend, BSActualization, LVVariance
R1R1Function  (abstract)  →  CallTerminalCondition, PutTerminalCondition, boundaries
```

`GetValue(t, S)` interpolates the grid at any (time, spot). The PDE and MC approaches converge to the same price for European calls.

---

## VarRed

```
ControlVariate               (abstract — simulate() + analytical_price())
├── BasketGeomControlVariate  geometric basket closed-form (exact for log-normal)
└── BSClosedForm              Black-Scholes formula
```

For a 1D European call, `BasketGeomControlVariate` with weight=1 IS the Black-Scholes price — the residual variance collapses to near-zero.

---

## Unit Tests

### Mac — `./build.sh`

Compiles and runs 6 test suites:

| Binary | What it covers |
|---|---|
| `tests/run_tests` | LCG formula, Normal moments, BSMilstein1D/2D paths, payoffs, MC prices, Bermudan |
| `tests/run_tests_nd` | BSMilsteinND for 1/2/3 assets, basket pricing, Bermudan ND |
| `tests/run_tests_qmc` | Halton vs pseudo-random agreement |
| `tests/run_tests_varred_compare` | VarRedMCPricer variance < plain |
| `tests/run_tests_antithetic` | AntitheticMCPricer variance < plain |
| `tests/run_tests_full_comparison` | All 4 methods head-to-head, 24 configs |

All print `[PASS]` / `[FAIL]` and exit 1 if anything fails.

### Windows

Open `Numerical Finance VS.sln` in Visual Studio → `F7`, or:

```
msbuild "Numerical Finance VS.sln" /p:Configuration=Debug /p:Platform=x64
```

Tests via `Test > Test Explorer`.

---

## Compilation & Build

<<<<<<< HEAD
The project uses **CMake** as its primary, cross-platform build system. This allows the same commands to work on Windows, Linux, and macOS.

### 1. Prerequisites
- **CMake** (3.10 or higher)
- A C++ compiler supporting **C++17** (GCC, Clang, or MSVC)

### 2. Universal Build Commands
From the root of the repository, run:

```bash
# 1. Configure the project (creates the build directory)
cmake -B build

# 2. Compile everything (Demo, Studies, Tests)
cmake --build build
```

The executables will be generated in the `build/` directory (e.g., `build/Debug/run_demo.exe` on Windows or `build/run_demo` on Unix).

### 3. Alternative Build Methods (Legacy)

#### Windows (Visual Studio)
Open `Numerical Finance VS.sln` in Visual Studio.
- Press **F7** to build the entire solution.
- Right-click a project (e.g., `Numerical Finance VS` or `UnitTestPricer`) and select **Set as Startup Project**.
- Press **F5** to run.

#### Windows (Command Line / MinGW)
Use the provided batch script:
```powershell
.\build.bat
```

#### macOS / Linux
Use the provided shell script:
```bash
chmod +x build.sh
./build.sh
```

=======
### Mac

```bash
xcode-select --install   # once, installs clang++
chmod +x build.sh        # once
./build.sh               # compiles everything + runs all tests
./build.sh clean         # remove compiled binaries
```

### Windows — full solution

```
msbuild "Numerical Finance VS.sln" /p:Configuration=Release /p:Platform=x64
```

### Windows — compile a single .cpp

Open a **Developer Command Prompt for VS** (search it in Start menu — it sets up `cl.exe` and the include paths):

```bat
rem MSVC
cl /std:c++17 /O2 /EHsc /I. MyFile.cpp /Fe:MyFile.exe

rem or with MinGW/MSYS2 (g++ must be on PATH)
g++ -std=c++17 -O2 -I. MyFile.cpp -o MyFile.exe
```

The `-I.` / `/I.` flag is important — it lets `#include "Eigen/Dense"` and `#include "RandomGenerator/Normal.h"` resolve from the repo root.

### Adding a new source file

Add the `.cpp` to the relevant section in `build.sh`:

```bash
Payoffs/PutPayoff.cpp   # just add a line
```

Do NOT include `pch.h` (MSVC only). Include paths are relative to the repo root.

>>>>>>> 663106143a3d19ab6d4b0eed0aea2db5bd2303fc
---

## Demo

Two tools in `demo/`. Both are built by `./build.sh`.

---

### Interactive pricer — `run_demo`

```bash
./demo/run_demo
```

Asks you for everything step by step — option type, number of assets, spots, vols, strike, rate, maturity, paths, variance reduction. All inputs are validated.

**European or Bermudan?**
European = exercised at maturity. Bermudan = early exercise at dates you specify (Longstaff-Schwartz).

**1 asset or multiple?**
1 asset = vanilla call. More = basket call with weights and a correlation rho.

**Variance reduction (European only):**

| Option | Method | Typical reduction |
|---|---|---|
| 0 | None — plain MC | baseline |
| 1 | QMC (Halton) — forces nbSteps=1 | ~1% |
| 2 | Antithetic — pairs Z and -Z | ~75% |
| 3 | Static CV — geometric basket closed form | ~99% for 1D |
| 4 | QMC + Static CV | best of both |

For Bermudan: only `none` and `qmc` available.

When you pick a VR method the output always shows plain MC alongside so you can see the variance reduction %.

---

### Config file mode

```bash
./demo/run_demo demo/example_config.json
```

Prices all options in the JSON file in one shot. `example_config.json` has 5 examples (1D plain, antithetic, 2D basket CV, Bermudan, 3D basket).

```json
[
  {
    "type": "european",
    "spots": [100.0],
    "vols": [0.2],
    "weights": [1.0],
    "strike": 100.0,
    "rate": 0.05,
    "maturity": 1.0,
    "nbSim": 10000,
    "nbSteps": 1,
    "varReduction": "none",
    "seed": 42
  }
]
```

| Field | Notes |
|---|---|
| `type` | `"european"` / `"bermudan"` |
| `spots`, `vols`, `weights` | arrays, one per asset. `weights` defaults to equal if omitted. |
| `rho` | correlation, single value applied to all pairs |
| `nbSteps` | for European, 1 is exact. For Bermudan use 50+ |
| `varReduction` | `"none"` / `"qmc"` / `"antithetic"` / `"staticcv"` / `"qmc_cv"` |
| `exerciseDates` | required for Bermudan, array of dates ≤ maturity |
| `seed` | optional, default 42 |

---

### Jupyter analysis — `run_study` + `analysis.ipynb`

For graphs: variance by method, CI width vs N, convergence, etc.

```
study_config.json  →  run_study  →  results/*.csv  →  analysis.ipynb  →  plots
```

**Step 1 — edit `demo/study_config.json`**

Set your option in the `"option"` block. The file has 5 example option blocks at the bottom (European, Bermudan, 2D basket, 3D basket, high-vol OTM) — copy the one you want.

`nbSim_values` drives the convergence sweep — each value is one data point on the graphs. `nbSim` inside `option` is ignored.

**Step 2 — run the study**

```bash
./demo/run_study demo/study_config.json
```

Writes to `results/`:
- `convergence.csv` — price, CI, variance for all 5 methods × all N values
- `terminal_plain.csv` — terminal spot prices and payoffs
- `antithetic_pairs.csv` — paired (Z, -Z) paths

**Step 3 — open the notebook**

```bash
jupyter notebook demo/analysis.ipynb
```

In the first cell, update to match your config:

```python
S0 = 100.0;  K = 100.0;  r = 0.05;  sigma = 0.2;  T = 1.0
option_label = 'European call, S=100, K=100, σ=20%, T=1y'
```

Hit **Run All**. Six plots:

| # | Plot |
|---|---|
| 1 | $S_T$ distribution vs theoretical log-normal |
| 2 | Antithetic pairing — $S_T$ correlation + convergence comparison |
| 3 | Price convergence vs N for all methods |
| 4 | CI width vs N log-log — $1/\sqrt{N}$ reference slope |
| 5 | Variance by method — bar chart, log scale, % reduction labeled |
| 6 | Paths needed to hit a target CI |

Plots also saved as PNGs in `results/`.

> **Why Static CV has no line in graph 4**: for a 1D European call, the geometric basket CV is the Black-Scholes price itself — residual variance is ~0. `log(0)` is undefined so there's nothing to plot. This is correct. Use a 2D/3D basket config to see a real line.

---

## Code Examples

<<<<<<< HEAD
### Adding a new source file

If you use **CMake**, adding a new file is **automatic**. Just create your `.cpp` file in one of the core folders (e.g. `Payoffs/`, `Pricer/`), and CMake will detect it automatically the next time you run:
```bash
cmake --build build
```

If you use the legacy scripts (`build.sh` or `build.bat`), you must manually add the path of the new file to the `SOURCES` or `LIB_SOURCES` list inside the script.

**Rules for new `.cpp` files:**

- Do NOT include `pch.h` — that's MSVC only. Just include what you need directly.
- Include paths are relative to the repo root (e.g. `#include "RandomGenerator/Normal.h"`).
- Do not add `dllmain.cpp`, `pch.cpp`, or `SDE/SDE.cpp` — those are MSVC stubs with no content.

---

### Adding a new test

Open `tests/test_main.cpp` and:

1. Write a function `void test_something() { ... }` using the `check(condition, "name")` helper:
=======
### European call (Monte Carlo)
>>>>>>> 663106143a3d19ab6d4b0eed0aea2db5bd2303fc

```cpp
LinearCongruential unifGen(1234, 16807, 0, 2147483647);
Normal normGen(0.0, 1.0, unifGen);
BSMilstein1D process(&normGen, 100.0, 0.05, 0.2);
EuropeanCallPayoff payoff(100.0);
EuropeanMCPricer pricer(&process, &payoff, 0.05, 1.0, 10000, 1);
double price = pricer.Price().price;  // ~10.45
```

### Bermudan call (quarterly exercise)

```cpp
std::vector<double> dates = {0.25, 0.5, 0.75, 1.0};
BermudanPricer pricer(&process, &payoff, 0.05, 1.0, 10000, 50, dates);
double price = pricer.Price().price;  // >= European price
```

### Basket call — 2 correlated assets

```cpp
BSMilsteinND proc(&normGen, {100.0, 100.0}, 0.05, {0.2, 0.25},
                  {{1.0, 0.3}, {0.3, 1.0}});
EuroCallBasketPayOff payoff(100.0, {0.5, 0.5});
EuropeanMCPricer pricer(&proc, &payoff, 0.05, 1.0, 10000, 1);
```

### With control variate

```cpp
LinearCongruential ugen(42, 16807, 0, 2147483647);
Normal normGen(0.0, 1.0, ugen);
BSMilstein1D proc(&normGen, 100.0, 0.05, 0.2);
EuropeanCallPayoff payoff(100.0);
BasketGeomControlVariate cv({1.0}, {100.0}, {{1.0}}, {0.2}, 1.0, 100.0, 0.05);
VarRedMCPricer pricer(&proc, &payoff, &cv, 0.05, 1.0, 10000, 1);
```

### European call via PDE

```cpp
BSVariance a(0.2);  BSTrend b(0.05);  BSActualization r(0.05);  NullFunction f;
CallTerminalCondition terminal(100.0);
CallTopBoundary top(100.0, 0.05, 0.2);  CallBottomBoundary bottom;
PDEGridImplicit grid(1.0, 0.0, 300.0, 100, 200, &a, &b, &r, &f, &top, &bottom, &terminal);
grid.FillNodes();
double price = grid.GetValue(0.0, 100.0);
```

---

## Data Flow

```
RandomGenerator  →  RandomProcess  →  SinglePath[]  →  PayOff  →  MCPricer → price
    (Normal)         (BSMilstein)      (time series)    (call)    (average)

Or without MC:
R2R1Function + R1R1Function  →  PDEGrid2D.FillNodes()  →  GetValue(0, spot)
```
