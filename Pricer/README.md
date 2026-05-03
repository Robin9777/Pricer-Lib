# Pricer

Monte Carlo pricing layer. All pricers derive from `MCPricer` and return a `PriceResult`.

## Class hierarchy

```
MCPricer                 (abstract — Price() → PriceResult)
├── EuropeanMCPricer
├── BermudanPricer
├── VarRedMCPricer       (static control variate)
└── AntitheticMCPricer
```

---

## `PriceResult`

```cpp
struct PriceResult {
    double price;
    double confidenceInterval;   // 95% CI half-width: 1.96 * std(mean)
};
```

---

## Classes

### `MCPricer`
Abstract base. Holds `RandomProcess*`, `PayOff*`, rate, maturity, nbSim, nbSteps.

---

### `EuropeanMCPricer`
Standard Monte Carlo: average discounted payoff over nbSim paths.

```
price = e^{-rT} · (1/N) · Σ payoff(path_i)
```

```cpp
EuropeanMCPricer pricer(&process, &payoff, rate, maturity, nbSim, nbSteps);
auto r = pricer.Price();
```

---

### `BermudanPricer`
Longstaff-Schwartz algorithm for options with early exercise rights at specified dates.

**Algorithm:**
1. Simulate all paths forward
2. Backward pass from last to first exercise date:
   - Compute intrinsic value (immediate exercise payoff)
   - Fit polynomial regression (degree 3) of discounted continuation value on current spot values
   - Exercise if intrinsic > fitted continuation value
3. Average discounted payoffs across paths

Uses **Eigen** for least-squares regression. Supports multi-asset (polynomial basis includes cross terms).

```cpp
BermudanPricer pricer(&process, &payoff, rate, maturity, nbSim, nbSteps,
                      {0.25, 0.5, 0.75, 1.0});   // exercise dates
auto r = pricer.Price();
```

---

### `VarRedMCPricer`
Monte Carlo with a static control variate.

```
adjusted_payoff = disc · arith_payoff − disc · CV_simulated + CV_analytical
```

The control variate must provide `SimulatedValue(paths)` and `AnalyticalExpectation()` (see VarRed module).

```cpp
BasketGeomControlVariate cv(weights, spots, corr, vols, T, K, r);
VarRedMCPricer pricer(&process, &payoff, &cv, rate, maturity, nbSim, nbSteps);
auto r = pricer.Price();
```

Typical variance reduction: **95–99%** for basket calls.

---

### `AntitheticMCPricer`
Monte Carlo with antithetic variables. For each of nbSim pairs:
1. Simulate with Z draws → payoff Y
2. Simulate with −Z draws → payoff Y_anti
3. Estimator = 0.5 · (Y + Y_anti)

Requires an `AntitheticNormal` generator (which must also be the process's generator).

```cpp
LinearCongruential ugen(seed, 16807, 0, 2147483647);
Normal normRaw(0.0, 1.0, ugen);
AntitheticNormal anti(normRaw);
BSMilstein1D process(&anti, spot, rate, vol);
AntitheticMCPricer pricer(&process, &payoff, anti, rate, maturity, nbSim, nbSteps);
auto r = pricer.Price();
```

Typical variance reduction: **65–80%** for European calls.
