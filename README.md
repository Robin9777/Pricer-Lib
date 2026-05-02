# Pricer-Lib — Numerical Finance C++ Library

A modular C++ library for Monte Carlo and PDE-based option pricing. Built as a Visual Studio (MSVC) solution targeting Windows x64.

---

## Table of Contents

1. [Why Abstract Classes?](#why-abstract-classes)
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
13. [Build](#build)
14. [Full Code Examples](#full-usage-examples)

---

## Why Abstract Classes?

This is the central design question of the entire library. The short answer: **the pricer should not care what process or payoff you give it**. Abstract classes make this possible.

### The problem without abstract classes

Imagine you write `EuropeanMCPricer` to work specifically with `BSMilstein1D`. Your pricer would look like:

```cpp
// BAD: tightly coupled
class EuropeanMCPricer {
    BSMilstein1D* process;   // what if we want Heston? or Brownian?
    EuropeanCallPayoff* payoff; // what if we want a basket? or a put?
    ...
};
```

Every time you want a new model or a new payoff, you would have to rewrite the pricer. The code is not reusable.

### The solution: abstract classes + polymorphism

Instead, the pricer holds **abstract pointers**:

```cpp
// GOOD: decoupled via abstract classes
class MCPricer {
    RandomProcess* Process;  // could be BSMilstein1D, Heston, BrownianND...
    PayOff* Payoff;          // could be EuropeanCall, Basket, Put...
};
```

At runtime, C++ calls the correct `Simulate()` or `operator()` through the **vtable** (virtual dispatch). The pricer calls `Process->Simulate()` without knowing whether it's Black-Scholes or Heston — it just knows that *something* will simulate a path.

### What each abstract class enforces

| Abstract Class | Pure Virtual Method | What it guarantees |
|---|---|---|
| `RandomGenerator` | `Generate()` | Any generator can produce a `double` on demand |
| `RandomProcess` | `Simulate(t0, T, n)` | Any process can fill its paths given a time range and step count |
| `PayOff` | `operator()(paths)` | Any payoff can compute a value from a set of paths |
| `MCPricer` | `Price()` | Any pricer exposes a single `Price()` entry point |
| `R2R1Function` | `operator()(x, t)` | Any PDE coefficient can be evaluated at any (spot, time) point |
| `R1R1Function` | `operator()(x)` | Any boundary/terminal condition can be evaluated at any spot |

### Concrete example: swapping a model

```cpp
// Price a European call under Black-Scholes
BSMilstein1D bsProcess(&normGen, 100.0, 0.05, 0.2);
EuropeanCallPayoff callPayoff(100.0);
EuropeanMCPricer pricer(&bsProcess, &callPayoff, 0.05, 1.0, 10000, 252);
double bsPrice = pricer.Price();

// Now price the same option under Heston — pricer code unchanged
Heston hestonProcess(&normGen, 100.0, 0.04, 0.05, 0.04, 2.0, 0.3, -0.7);
EuropeanMCPricer pricer2(&hestonProcess, &callPayoff, 0.05, 1.0, 10000, 252);
double hestonPrice = pricer2.Price();
```

The pricer is identical in both cases. Only the concrete object passed as `RandomProcess*` changes.

### Intermediate abstract layers (UniformGenerator, ContinuousGenerator, etc.)

These exist for a second reason: **some classes need a more specific type than the root abstract**.

For example, `Normal` needs a `UniformGenerator&` — not just any `RandomGenerator` — because the Box-Muller algorithm internally needs uniform samples. So `UniformGenerator` acts as a contract: *"I am a generator that produces uniform values"*. This prevents you from accidentally passing a `Normal` generator to another `Normal` generator.

---

## Project Structure

```
Pricer-Lib/
├── RandomGenerator/        # Layer 1 — Random number generation
├── SDE/                    # Layer 2 — Stochastic processes
├── Payoffs/                # Layer 3 — Payoff functions
├── Pricer/                 # Layer 4 — Monte Carlo pricers
├── PDE/                    # PDE solver (Black-Scholes finite differences)
├── VolSurf/                # Volatility surface (Dupire local vol)
├── VarRed/                 # Variance reduction (placeholder)
├── Eigen/                  # Bundled Eigen header-only library
├── UnitTestRandomNumber/   # Tests for RandomGenerator
├── UnitTestSDE/            # Tests for SDE
├── UnitTestPayoffs/        # Tests for Payoffs
├── UnitTestPricer/         # Tests for Pricer
├── UnitTestPDE/            # Tests for PDE
├── Numerical Finance VS/   # Executable entry point (main.cpp)
└── Numerical Finance VS.sln
```

---

## Full Workflow

### Monte Carlo path (the main approach)

The full pipeline from random numbers to an option price:

```
Step 1: You create a random generator
        LinearCongruential -> gives uniform numbers in [0,1]
        Normal (wraps the uniform generator) -> gives Gaussian samples

Step 2: You create a stochastic process
        BSMilstein1D (wraps the Normal generator)
        When Simulate() is called, it draws Gaussian samples at each time step
        and builds a discretized price path S(t0), S(t1), ..., S(T)
        stored inside a SinglePath object

Step 3: You create a payoff
        EuropeanCallPayoff (knows the strike K)
        When called with a SinglePath, it reads the final value S(T)
        and returns max(S(T) - K, 0)

Step 4: You create a pricer
        EuropeanMCPricer (wraps the process and the payoff)
        When Price() is called:
          - it repeats steps 2-3 nbSim times (e.g. 10,000 simulations)
          - each time it gets a different path because the generator
            advances its internal state (it's pseudo-random, deterministic
            but the sequence is long enough to behave randomly)
          - it discounts each payoff: exp(-r * T) * payoff
          - it averages everything: the result is the Monte Carlo estimate
            of the option price, which converges to the true price
            as nbSim -> infinity (law of large numbers)
```

### Why Monte Carlo works financially

By the **risk-neutral pricing theorem**, the fair price of an option is:

```
Price = E^Q [ exp(-r*T) * Payoff(S_T) ]
```

Where the expectation is taken under the risk-neutral measure Q (drift = r, not mu). Monte Carlo estimates this expectation by averaging many samples. The more simulations, the more accurate. The standard error decreases as `1/sqrt(nbSim)`.

### PDE path (alternative approach)

Instead of simulating random paths, you can solve the Black-Scholes PDE directly on a grid:

```
Step 1: Define the PDE coefficients as R2R1Function objects
        BSVariance  -> a(S,t) = 0.5 * sigma^2 * S^2  (diffusion term)
        BSTrend     -> b(S,t) = r * S                  (drift term)
        BSActualization -> r(S,t) = r                  (discount term)

Step 2: Define terminal and boundary conditions as R1R1Function objects
        CallTerminalCondition -> V(S, T) = max(S - K, 0)
        Boundary conditions at S_min and S_max

Step 3: Create the grid and fill it
        PDEGridImplicit grid(T, Smin, Smax, Nt, Ns, ...)
        grid.FillNodes()  <- solves the PDE backwards in time

Step 4: Read the price
        double price = grid.GetValue(0.0, spot)
        <- interpolates the grid at time 0 and the desired spot
```

The PDE and MC approaches give the same result for simple options (European calls/puts). PDE is faster and exact for 1D problems; MC is more flexible and scales to high-dimensional problems (baskets, Heston, etc.).

---

## Layer 1 — `RandomGenerator/`

Abstract base class and full hierarchy for random number generation.

### Class Hierarchy

```
RandomGenerator              (abstract)
│   Generate() = 0           <- must be implemented by every concrete class
│   Mean(n)                  <- calls Generate() n times, returns average
│   Variance(n)              <- calls Generate() n times, returns variance
│
├── UniformGenerator         (abstract — marks generators that produce uniform [0,1] values)
│   └── PseudoGenerator      (adds a seed member)
│       ├── LinearCongruential   — xn+1 = (a*xn + c) mod m, normalized to [0,1]
│       └── EcuyerCombined       — combines two LinearCongruential for a much longer period
│
├── ContinuousGenerator      (abstract — marks generators for continuous distributions)
│   ├── Normal               — N(mu, sigma^2), built on top of a UniformGenerator
│   │                          three internal algorithms: BoxMuller, CLT, RejectionSampling
│   └── Exponential          — exponential distribution
│
└── DiscreteGenerator        (abstract — marks generators for discrete distributions)
    ├── Bernoulli
    ├── Binomial
    ├── Poisson
    ├── FiniteSet
    └── HeadTail
```

### Why this layered hierarchy?

- `RandomGenerator` defines the contract for everything: *"produce a number"*.
- `UniformGenerator` is a sub-contract: *"produce a uniform number in [0,1]"*. This is needed because `Normal` requires uniform inputs — it needs to know it has a uniform source, not just any generator.
- `ContinuousGenerator` and `DiscreteGenerator` are semantic tags — they group generators by the type of distribution they produce, making the code self-documenting and allowing future type-based constraints.

### Key Classes

| Class | File | Description |
|---|---|---|
| `RandomGenerator` | `RandomGenerator.h` | Pure abstract base: `Generate()`, `Mean(n)`, `Variance(n)` |
| `LinearCongruential` | `LinearCongruential.h` | Formula: `xn+1 = (a*xn + c) mod m`. Constructor: `(seed, multiplier, increment, modulus)`. Common params: a=16807, c=0, m=2147483647 (Park-Miller) |
| `EcuyerCombined` | `EcuyerCombined.h` | Takes two `LinearCongruential&`. Subtracts their outputs and normalizes, achieving a period ~2.3 * 10^18 |
| `Normal` | `Normal.h` | Takes `(mu, sigma, UniformGenerator&)`. Default `Generate()` = Box-Muller. `Generate(NormalAlgo::RejectionSampling)` lets you pick |

### How `Normal::BoxMuller` works

Box-Muller transforms two independent uniforms U1, U2 into two independent Gaussians:

```
Z1 = sqrt(-2 * ln(U1)) * cos(2*pi*U2)
Z2 = sqrt(-2 * ln(U1)) * sin(2*pi*U2)
```

This is why `Normal` needs a `UniformGenerator` reference — it calls `Ugen.Generate()` twice per sample.

---

## Layer 2 — `SDE/`

Stochastic processes that simulate price paths using a `RandomGenerator`.

### Class Hierarchy

```
RandomProcess                (abstract)
│   Generator*               <- pointer to any RandomGenerator
│   vector<SinglePath*>      <- one path per dimension
│   Dimension                <- number of assets/processes
│   Simulate(t0, T, n) = 0   <- must fill the paths
│   GetPath(d)               <- returns path d after simulation
│
├── BrownianD1               — 1D standard Brownian motion W(t)
├── BrownianND               — N-dimensional independent Brownian motions
├── BlackScholes1D           (abstract base, stores Spot/Rate/Vol)
│   ├── BSEuler1D            — Euler-Maruyama: S(t+dt) = S(t) + r*S*dt + sigma*S*dW
│   └── BSMilstein1D         — Milstein: adds correction term 0.5*sigma^2*S*(dW^2 - dt)
├── BlackScholes2D           (abstract base, stores 2 spots, correlation rho)
│   └── BSMilstein2D         — Milstein on 2 correlated assets via Cholesky decomposition
└── Heston                   — spot + variance (CIR process), with correlation rho
```

### Why separate `BlackScholes1D` from `BSMilstein1D`?

`BlackScholes1D` is an abstract intermediate that stores the financial parameters (Spot, Rate, Vol) shared by all 1D Black-Scholes schemes. `BSEuler1D` and `BSMilstein1D` only differ in how they discretize the SDE — they inherit the parameters and implement their own `Simulate()`. This avoids duplicating the constructor and parameter storage.

### Discretization schemes explained

**Euler-Maruyama** (first order):
```
S(t + dt) = S(t) * [1 + r*dt + sigma*sqrt(dt)*Z]
            where Z ~ N(0,1)
```

**Milstein** (second order, higher accuracy):
```
S(t + dt) = S(t) * [1 + r*dt + sigma*sqrt(dt)*Z + 0.5*sigma^2*(Z^2*dt - dt)]
```
The extra term corrects for the curvature of the SDE and reduces the discretization error from O(dt) to O(dt^2).

**Heston model** (stochastic volatility):
```
dS = mu * S * dt + sqrt(V) * S * dW1
dV = kappa * (theta - V) * dt + sigma_v * sqrt(V) * dW2
corr(dW1, dW2) = rho
```
The variance V follows a mean-reverting CIR process. Heston requires 2 correlated Brownian motions — the correlation is applied via a Cholesky-like decomposition: `dW2 = rho*dW1 + sqrt(1-rho^2)*dW_independent`.

### Key Classes

| Class | File | Description |
|---|---|---|
| `RandomProcess` | `RandomProcess.h` | Abstract base. Holds `Generator*`, `vector<SinglePath*>`, `Dimension` |
| `SinglePath` | `SinglePath.h/.cpp` | Stores one time series as a `vector<double>`. `GetAllValues()` returns all steps. `InsertValue(v)` appends a step |
| `BSMilstein1D` | `BSMilstein1D.h` | Constructor: `(Generator*, spot, rate, vol)`. After `Simulate()`, `GetPath(0)` holds S(0)..S(T) |
| `BSMilstein2D` | `BSMilstein2D.h` | Constructor: `(Generator*, spot1, spot2, rate, vol1, vol2, rho)`. After `Simulate()`, `GetPath(0)` = asset 1, `GetPath(1)` = asset 2 |
| `Heston` | `Heston.h` | Constructor: `(Generator*, spot, initVariance, mu, theta, kappa, sigma, rho)`. `GetPath(0)` = spot, `GetPath(1)` = variance |

---

## Layer 3 — `Payoffs/`

Payoff functions that take a set of simulated paths and return the cash flow of the option.

### Class Hierarchy

```
PayOff                         (abstract)
│   operator()(paths) = 0      <- takes vector<SinglePath*>, returns double
│
├── EuropeanCallPayoff         — max(S_T - K, 0)  on path[0]
└── EuroCallBasketPayOff       — max(sum_i(w_i * S_T^i) - K, 0)  across all paths
```

### Why abstract PayOff?

The pricer calls `(*Payoff)(paths)` through a pointer. It does not know — and should not know — whether the payoff is a call, a put, a basket, a digital, or anything else. Adding a new exotic payoff is as simple as creating a new class that inherits `PayOff` and implementing `operator()`. The pricer, process, and generator do not change at all.

### Key Classes

| Class | File | Description |
|---|---|---|
| `PayOff` | `PayOff.h` | Abstract base. Pure virtual `operator()(const vector<SinglePath*>&) const` |
| `EuropeanCallPayoff` | `EuropeanCallPayoff.cpp` | Constructor: `(strike)`. Reads `Paths[0]->GetAllValues().back()` = final spot price, returns `max(S_T - K, 0)` |
| `EuroCallBasketPayOff` | `EuroCallBasketPayOff.cpp` | Constructor: `(strike, weights)`. Computes `sum(w_i * Paths[i]->last())` across all assets, returns `max(basket - K, 0)` |

---

## Layer 4 — `Pricer/`

Monte Carlo pricers that tie together a process, a payoff, and simulation parameters.

### Class Hierarchy

```
MCPricer                     (abstract)
│   Process*                 <- any RandomProcess
│   Payoff*                  <- any PayOff
│   rate, maturity           <- financial parameters
│   nbSim, nbSteps           <- simulation parameters
│   Price() = 0              <- must return the option price
│
├── EuropeanMCPricer         — plain Monte Carlo averaging
└── BermudanPricer           — Longstaff-Schwartz with early exercise
```

### Key Classes

| Class | File | Description |
|---|---|---|
| `MCPricer` | `MCPricer.cpp` | Abstract base. Constructor stores all parameters |
| `EuropeanMCPricer` | `EuropeanMCPricer.cpp` | Runs `nbSim` simulations, averages `exp(-r*T) * payoff` |
| `BermudanPricer` | `BermudanPricer.cpp` | Longstaff-Schwartz. Constructor adds `vector<double> exerciseDates`. Uses Eigen for regression. `PolynomialDegree = 3` |

### `EuropeanMCPricer::Price()` — step by step

```
for sim = 1 to nbSim:
    Process->Simulate(0, T, nbSteps)      // generate a new path
    paths = [Process->GetPath(0), ...]    // collect all dimensions
    payoff = (*Payoff)(paths)             // evaluate the payoff at maturity
    discounted = exp(-r * T) * payoff     // bring value back to today
    accumulate sum and sum-of-squares

return sum / nbSim                        // Monte Carlo estimate
```

The variance of the estimate is `Var[payoff] / nbSim`. To halve the error, you need 4x more simulations.

### `BermudanPricer::Price()` — Longstaff-Schwartz algorithm

A Bermudan option can be exercised at specific dates (not just at maturity). The challenge is deciding at each date: *should I exercise now or wait?*. Longstaff-Schwartz solves this by backward induction with regression.

```
Phase 1 — Forward simulation
  For each simulation, run the full path and record:
    - the spot state at each exercise date
    - the payoff at each exercise date

Phase 2 — Backward induction (from last date to first)
  Start: V[sim] = payoff at last exercise date for each sim

  For each exercise date k going backwards:
    1. Find in-the-money paths (payoff > 0 at date k)
    2. For these paths, regress V[sim] * exp(-r*dt) against the spot states
       using a polynomial basis of degree 3 (via Eigen least-squares)
       -> this gives an estimate of E[continuation value | S_k]
    3. For each in-the-money path:
       if immediate exercise >= estimated continuation: exercise (V[sim] = payoff_k)
       else: continue  (V[sim] = discounted future value)

Phase 3 — Price
  return average of V[sim] * exp(-r * first_exercise_date)
```

The regression step is why Eigen is needed: it solves the least-squares system `Phi * beta = Y` where `Phi` is the matrix of polynomial basis functions evaluated at each path's spot value.

---

## `PDE/`

Finite-difference solver for the Black-Scholes PDE. An alternative to Monte Carlo that works directly on a grid of (spot, time) values.

### PDE solved

```
dV/dt + a(S,t) * d²V/dS² + b(S,t) * dV/dS - r(S,t) * V + f(S,t) = 0
```

For standard Black-Scholes: `a = 0.5*sigma^2*S^2`, `b = r*S`, `r(S,t) = r`, `f = 0`.

### Class Hierarchy

```
PDEGrid2D                    (abstract)
│   Nodes[height][width]     <- the 2D grid of option values V(S, t)
│   h0 = space step (dS)
│   h1 = time step (dt)
│   a, b, r, f               <- R2R1Function* (PDE coefficients)
│   Top/Bottom/Right         <- R1R1Function* (boundary/terminal conditions)
│   FillNodes() = 0          <- must fill the grid
│   GetValue(t, S)           <- interpolates the grid
│
├── PDEGridExplicit          — computes each node from already-known future nodes
└── PDEGridImplicit          — solves a tridiagonal linear system at each time step

R2R1Function                 (abstract: operator()(spot, time) -> double)
├── BSVariance               — 0.5 * sigma^2 * S^2
├── BSTrend                  — r * S
├── BSActualization          — r
├── LVVariance               — local vol variance from DupireSurface
└── NullFunction             — 0.0

R1R1Function                 (abstract: operator()(spot) -> double)
├── CallTerminalCondition    — max(S - K, 0) at t = T
├── PutTerminalCondition     — max(K - S, 0) at t = T
├── VanillaTerminalCondition — generic vanilla terminal payoff
├── CallTopBoundary          — V(Smax, t) = Smax - K*exp(-r*(T-t))
├── CallBottomBoundary       — V(0, t) = 0
├── PutTopBoundary           — V(Smax, t) = 0
└── PutBottomBoundary        — V(0, t) = K*exp(-r*(T-t))
```

### Why R2R1Function and R1R1Function as abstract classes?

The grid solver is written generically — it calls `(*a)(S, t)` at each node without knowing whether that's Black-Scholes variance, local vol, or anything else. This lets you plug any model coefficients into the same solver. Swapping from constant vol to local vol only means passing a different `R2R1Function*` — the grid solving code is untouched.

### Explicit vs Implicit scheme

**Explicit**: each node `V(S, t)` is computed directly from `V(S, t+dt)`, `V(S+dS, t+dt)`, `V(S-dS, t+dt)` — values already known. Simple but **conditionally stable**: requires small time steps (`dt < dS^2`).

**Implicit**: the unknowns at time `t` are coupled — solving requires inverting a tridiagonal system at each time step. More expensive per step but **unconditionally stable**: works with larger time steps.

### Key Classes

| Class | File | Description |
|---|---|---|
| `PDEGrid2D` | `PDEGrid2D.h` | Abstract base. Stores the grid and all function pointers. `GetValue(t, S)` interpolates |
| `PDEGridExplicit` | `PDEGridExplicit.h` | Explicit scheme. Fast setup, requires fine grid for stability |
| `PDEGridImplicit` | `PDEGridImplicit.h` | Implicit scheme. Tridiagonal solve per time step. Stable for coarse grids |
| `R2R1Function` | `R2R1Function.h` | Abstract base for `(spot, time) -> double` coefficient functions |
| `R1R1Function` | `R1R1Function.h` | Abstract base for `(spot) -> double` boundary/terminal functions |

---

## `VolSurf/`

Local volatility surface support for Dupire's model.

| Class | File | Description |
|---|---|---|
| `DupireSurface` | `DupireSurface.h` | Holds a 2D grid of implied/local vols indexed by `(spot, time)`. `GetVariance(spot, time)` interpolates the surface |

`DupireSurface` is used by `LVVariance` (in `PDE/`) as the coefficient function `a(S,t)` in the PDE. This allows the PDE solver to price options under a local volatility model (where vol is a function of both spot and time) rather than constant vol.

---

## `VarRed/`

Placeholder project for variance reduction techniques. Currently empty. Intended for:
- **Antithetic variates**: for each path with draw `Z`, also compute the path with `-Z`. Averaging the two reduces variance.
- **Control variates**: use a correlated quantity with a known analytical price to reduce MC error.
- **Importance sampling**: shift the sampling distribution to focus on paths that contribute more to the payoff.

---

## `Eigen/`

Bundled header-only linear algebra library (version 3.x). Required by `BermudanPricer` for the least-squares polynomial regression step in Longstaff-Schwartz.

The key operation is:
```cpp
Eigen::VectorXd beta = Phi.colPivHouseholderQr().solve(Y);
// Phi: matrix of basis functions evaluated at each path state
// Y:   vector of discounted continuation values
// beta: regression coefficients
```

Include with:
```cpp
#include "Eigen/Dense"
```

---

## Unit Tests

### On Mac — `tests/test_main.cpp`

The file `tests/test_main.cpp` is a standalone test runner that works on Mac with `clang++`. It does not use any MSVC/Windows test framework. You run it via `build.sh` (see [Build](#build) below).

**What it tests:**

| Test | What it checks |
|---|---|
| `LinearCongruential: first value matches formula` | `Generate()` returns `(a*seed + c) % m` normalized — verifies the LCG formula is implemented correctly |
| `LinearCongruential: output in [0,1]` | 1000 samples all fall in `[0.0, 1.0]` |
| `LinearCongruential: same seed gives same sequence` | Two generators with the same seed produce identical output — verifies determinism |
| `Normal: mean close to 0` | Mean of 10,000 samples is within `[-0.05, 0.05]` |
| `Normal: variance close to 1` | Variance of 10,000 samples is within `[0.9, 1.1]` |
| `BSMilstein1D: path length = nbSteps+1` | After `Simulate(0, T, N)`, the path vector has `N+1` values (initial spot + N steps) |
| `BSMilstein1D: all path values > 0` | Stock price never goes negative (geometric SDE property) |
| `BSMilstein1D: same seed gives same path` | Two processes with the same seed produce identical paths |
| `BSMilstein2D: path[0] and path[1] length correct` | Both asset paths have the right number of steps |
| `BSMilstein2D: all values > 0` | Both assets stay positive |
| `EuropeanCallPayoff: ITM payoff = 20` | `S_T = 120, K = 100` → payoff = `20.0` |
| `EuropeanCallPayoff: OTM payoff = 0` | `S_T = 80, K = 100` → payoff = `0.0` |
| `EuropeanCallPayoff: ATM payoff = 0` | `S_T = 100, K = 100` → payoff = `0.0` (not in the money) |
| `EuropeanMCPricer: ATM call price in [9.0, 12.0]` | 10,000 sims, S=K=100, r=5%, vol=20%, T=1y → MC price should be near the Black-Scholes analytical value of ~10.45 |
| `BermudanPricer: price > 0 and in plausible range` | Bermudan call with quarterly exercise dates prices between 0 and 25 |
| `EuroCallBasketPayOff: weight=1.0 matches vanilla call` | A basket payoff with a single asset and weight=1.0 should be identical to a vanilla European call on the same path |

**What is NOT tested (Mac):**

- The PDE module (`PDE/`) — finite-difference solver. It compiles and is tested separately on Windows via `UnitTestPDE/`. It has no cross-platform build yet.
- `Heston` process pricing — the process compiles and the path simulation works, but there is no dedicated pricer test for it in `test_main.cpp`.

**Output format:**

Every test prints `[PASS] <name>` or `[FAIL] <name>`. If any test fails, the program exits with code 1 (useful for CI). Example:

```
=== Pricer-Lib Test Suite ===

[PASS] LinearCongruential: first value matches formula
[PASS] LinearCongruential: output in [0,1]
[PASS] LinearCongruential: same seed gives same sequence
[PASS] Normal: mean close to 0
[PASS] Normal: variance close to 1
[PASS] BSMilstein1D: path length = nbSteps+1
[PASS] BSMilstein1D: all path values > 0
[PASS] BSMilstein1D: same seed gives same path
[PASS] BSMilstein2D: path[0] length correct
[PASS] BSMilstein2D: path[1] length correct
[PASS] BSMilstein2D: all values > 0
[PASS] EuropeanCallPayoff: ITM payoff = 20
[PASS] EuropeanCallPayoff: OTM payoff = 0
[PASS] EuropeanCallPayoff: ATM payoff = 0
[PASS] EuropeanMCPricer: ATM call price in [9.0, 12.0]
[PASS] BermudanPricer: price > 0 and in plausible range
[PASS] EuroCallBasketPayOff: weight=1.0 matches vanilla call with same seed

All tests PASSED.
```

### On Windows — MSVC test projects

Each module has a corresponding test project using **Microsoft CppUnitTestFramework**.

| Project | Test files | What is tested |
|---|---|---|
| `UnitTestRandomNumber/` | `UnitTestLinearCongruential.cpp`, `UnitTestEcuyerCombined.cpp`, `UnitTestNormal.cpp`, `UnitTestExponential.cpp` | Statistical properties (mean, variance), sequence correctness |
| `UnitTestSDE/` | `TestBrownian1D.cpp`, `TestBrownianND.cpp`, `TestBSEuler1D.cpp`, `TestBSMilstein.cpp`, `TestHeston.cpp` | Path length, positivity, moment properties |
| `UnitTestPayoffs/` | `UnitTestPayoffs.cpp` | Payoff values for known inputs (at-the-money, in-the-money, out-of-the-money) |
| `UnitTestPricer/` | `UnitTestPricer.cpp` | MC prices compared to Black-Scholes analytical formula |
| `UnitTestPDE/` | `UnitTestPDE.cpp` | PDE prices compared to analytical formula, convergence with grid refinement |

Run via Visual Studio Test Explorer: `Test > Test Explorer`.

---

## Build

### Files added for Mac support

```
build.sh            — shell script: compiles everything and runs the tests
tests/
└── test_main.cpp   — standalone test runner (no MSVC dependencies)
```

Two small header changes were also made to fix Windows-only includes:
- `RandomGenerator/framework.h` — `#include <windows.h>` is now guarded with `#ifdef _WIN32`
- `RandomGenerator/RandomGenerator.h` — added `#include <cstddef>` so `size_t` is always available

These changes do not affect the Windows build.

### Mac

**Prerequisites**: Xcode command line tools (provides `clang++`):

```bash
xcode-select --install
```

**Option 1 — script (compile + run in one shot):**

```bash
chmod +x build.sh   # only needed once
./build.sh
```

**Option 2 — by hand:**

Compile (run from the repo root):

```bash
clang++ -std=c++17 -O2 -I. \
  RandomGenerator/RandomGenerator.cpp \
  RandomGenerator/UniformGenerator.cpp \
  RandomGenerator/PseudoGenerator.cpp \
  RandomGenerator/LinearCongruential.cpp \
  RandomGenerator/EcuyerCombined.cpp \
  RandomGenerator/ContinuousGenerator.cpp \
  RandomGenerator/Normal.cpp \
  RandomGenerator/Exponential.cpp \
  RandomGenerator/DiscreteGenerator.cpp \
  RandomGenerator/Bernoulli.cpp \
  RandomGenerator/Binomial.cpp \
  RandomGenerator/Poisson.cpp \
  RandomGenerator/FiniteSet.cpp \
  RandomGenerator/HeadTail.cpp \
  SDE/SinglePath.cpp \
  SDE/RandomProcess.cpp \
  SDE/BlackScholes1D.cpp \
  SDE/BlackScholes2D.cpp \
  SDE/BSEuler1D.cpp \
  SDE/BSMilstein1D.cpp \
  SDE/BSMilstein2D.cpp \
  SDE/BrownianD1.cpp \
  SDE/BrownianND.cpp \
  SDE/Heston.cpp \
  Payoffs/PayOff.cpp \
  Payoffs/EuropeanCallPayoff.cpp \
  Payoffs/EuroCallBasketPayOff.cpp \
  Pricer/MCPricer.cpp \
  Pricer/EuropeanMCPricer.cpp \
  Pricer/BermudanPricer.cpp \
  tests/test_main.cpp \
  -o tests/run_tests
```

Then run it:

```bash
./tests/run_tests
```

**Clean:**

```bash
./build.sh clean
# or by hand:
rm tests/run_tests
```

---

### Adding a new source file

If you create a new `.cpp` file (e.g. `Payoffs/PutPayoff.cpp`), open `build.sh` and add it to the `SOURCES` array in the relevant section:

```bash
# Payoffs
Payoffs/PayOff.cpp
Payoffs/EuropeanCallPayoff.cpp
Payoffs/EuroCallBasketPayOff.cpp
Payoffs/PutPayoff.cpp        # <-- add it here
```

That's it. Run `./build.sh` again and it will be compiled in.

**Rules for new `.cpp` files:**

- Do NOT include `pch.h` — that's MSVC only. Just include what you need directly.
- Include paths are relative to the repo root (e.g. `#include "RandomGenerator/Normal.h"`).
- Do not add `dllmain.cpp`, `pch.cpp`, or `SDE/SDE.cpp` — those are MSVC stubs with no content.

---

### Adding a new test

Open `tests/test_main.cpp` and:

1. Write a function `void test_something() { ... }` using the `check(condition, "name")` helper:

```cpp
void test_put_payoff() {
    PutPayoff payoff(100.0);
    SinglePath path(0.0, 1.0, 1);
    path.InsertValue(80.0);
    std::vector<SinglePath*> paths = {&path};
    check(std::abs(payoff(paths) - 20.0) < 1e-12, "PutPayoff: ITM = 20");
}
```

2. Call it from `main()`:

```cpp
int main() {
    // existing tests...
    test_put_payoff();   // <-- add your call here
    ...
}
```

`check(condition, name)` prints `[PASS]` or `[FAIL]` and counts failures. The program exits with code 1 if anything fails.

### Windows

Open `Numerical Finance VS.sln` in Visual Studio and press `F7`, or via MSBuild:

```bash
msbuild "Numerical Finance VS.sln" /p:Configuration=Debug /p:Platform=x64
```

Run tests via Visual Studio Test Explorer: `Test > Test Explorer`.

---

## Full Usage Examples

### European call under Black-Scholes (Monte Carlo)

```cpp
// Step 1: random generator
LinearCongruential unifGen(1234, 16807, 0, 2147483647);
Normal normGen(0.0, 1.0, unifGen);

// Step 2: stochastic process
BSMilstein1D process(&normGen, /*spot=*/100.0, /*rate=*/0.05, /*vol=*/0.2);

// Step 3: payoff
EuropeanCallPayoff payoff(/*strike=*/100.0);

// Step 4: pricer
EuropeanMCPricer pricer(&process, &payoff, /*rate=*/0.05, /*maturity=*/1.0,
                        /*nbSim=*/10000, /*nbSteps=*/252);
double price = pricer.Price();
// Expected: ~10.45 (Black-Scholes analytical value)
```

### Bermudan option (early exercise at quarterly dates)

```cpp
std::vector<double> exerciseDates = {0.25, 0.5, 0.75, 1.0};
BermudanPricer bermPricer(&process, &payoff, 0.05, 1.0, 10000, 252, exerciseDates);
double bermPrice = bermPricer.Price();
// Bermudan price >= European price (more rights = more value)
```

### Basket option on 2 correlated assets

```cpp
// 2 assets: S1=100, S2=100, rate=5%, vol1=20%, vol2=25%, correlation=0.3
BSMilstein2D process2D(&normGen, 100.0, 100.0, 0.05, 0.2, 0.25, /*rho=*/0.3);

// Equal-weight basket call, strike 100
EuroCallBasketPayOff basketPayoff(100.0, {0.5, 0.5});

EuropeanMCPricer basketPricer(&process2D, &basketPayoff, 0.05, 1.0, 10000, 252);
double basketPrice = basketPricer.Price();
```

### European call under Heston (stochastic volatility)

```cpp
// spot=100, V0=0.04 (vol=20%), mu=5%, theta=0.04, kappa=2.0, sigma_v=0.3, rho=-0.7
Heston hestonProcess(&normGen, 100.0, 0.04, 0.05, 0.04, 2.0, 0.3, -0.7);

EuropeanMCPricer hestonPricer(&hestonProcess, &payoff, 0.05, 1.0, 10000, 252);
double hestonPrice = hestonPricer.Price();
```

### European call via PDE (finite differences)

```cpp
BSVariance    a(/*sigma=*/0.2);
BSTrend       b(/*rate=*/0.05);
BSActualization r(0.05);
NullFunction  f;
CallTerminalCondition terminal(/*strike=*/100.0);
CallTopBoundary  top(100.0, 0.05, 0.2);
CallBottomBoundary bottom;

PDEGridImplicit grid(/*T=*/1.0, /*Smin=*/0.0, /*Smax=*/300.0,
                     /*Nt=*/100, /*Ns=*/200,
                     &a, &b, &r, &f, &top, &bottom, &terminal);
grid.FillNodes();
double price = grid.GetValue(/*time=*/0.0, /*spot=*/100.0);
```

---

## Data Flow Diagram

```
                     +------------------+
                     |  RandomGenerator  |  Generate() -> double
                     |  (LinearCongr,    |
                     |   Normal, ...)    |
                     +--------+---------+
                              |
                              | feeds random numbers
                              v
                     +------------------+
                     |  RandomProcess    |  Simulate(t0, T, n)
                     |  (BSMilstein1D,   |  -> fills SinglePath[]
                     |   Heston, ...)    |
                     +--------+---------+
                              |
                              | produces paths S(t0)...S(T)
                              v
                     +------------------+
                     |  SinglePath[]     |  GetAllValues() -> vector<double>
                     +--------+---------+
                              |
                              | paths passed to payoff
                              v
                     +------------------+
                     |  PayOff           |  operator()(paths) -> double
                     |  (EuropeanCall,   |
                     |   Basket, ...)    |
                     +--------+---------+
                              |
                              | raw payoff value
                              v
                     +------------------+
                     |  MCPricer         |  Price() -> double
                     |  (European,       |  discounts and averages
                     |   Bermudan)       |  over nbSim simulations
                     +------------------+


  Alternatively, without Monte Carlo:

  R2R1Function (a,b,r,f)  +  R1R1Function (terminal, boundaries)
                              |
                              v
                     +------------------+
                     |  PDEGrid2D        |  FillNodes() solves the PDE
                     |  (Explicit or     |  backwards from T to 0
                     |   Implicit)       |
                     +--------+---------+
                              |
                              v
                     GetValue(0, spot) -> option price
```
