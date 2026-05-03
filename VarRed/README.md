# VarRed

Variance reduction helpers for the control variate method. Used by `VarRedMCPricer`.

## Class hierarchy

```
ControlVariate                    (abstract)
└── BasketGeomControlVariate

BSClosedForm                      (standalone — Black-Scholes analytical pricer)
```

---

## Classes

### `ControlVariate`
Abstract base for any control variate. Two methods:

- `SimulatedValue(paths)` — the CV value computed from the simulated path (same randomness as the main payoff)
- `AnalyticalExpectation()` — the known closed-form expectation of the CV

`VarRedMCPricer` computes: `adjusted = disc·payoff − disc·CV_simulated + CV_analytical`

---

### `BasketGeomControlVariate`
Control variate based on the **geometric basket**: `G_T = Π S_T^(i)^{wᵢ}` (weighted geometric mean).

- `SimulatedValue(paths)` — `max(G_T − K, 0)` using last path values
- `AnalyticalExpectation()` — exact Black-Scholes price of the geometric basket call

The geometric basket follows a log-normal law with:
- Geometric vol: `σ_G² = Σᵢ Σⱼ wᵢ wⱼ σᵢ σⱼ ρᵢⱼ`
- Drift adjustment: `G₀_eff = G₀ · exp(0.5 · (σ_G² − Σᵢ wᵢ σᵢ²) · T)`

```cpp
BasketGeomControlVariate cv(
    {0.5, 0.5},                   // weights
    {100.0, 100.0},               // spots
    {{1,0.3},{0.3,1}},            // correlation matrix
    {0.2, 0.2},                   // vols
    1.0,                          // maturity
    100.0,                        // strike
    0.05                          // rate
);
double simVal  = cv.SimulatedValue(paths);
double anaVal  = cv.AnalyticalExpectation();
```

**Why it works:** arithmetic ≥ geometric (AM-GM), so G_T is a near-perfect proxy for the arithmetic basket A_T. Their difference has near-zero variance → almost all Monte Carlo variance is explained by the CV.

---

### `BSClosedForm`
Standalone Black-Scholes analytical pricer. Used internally by `BasketGeomControlVariate`.

```cpp
BSClosedForm bs;
double call = bs.CallPrice(spot, strike, rate, vol, maturity);
double put  = bs.PutPrice(spot, strike, rate, vol, maturity);
```
