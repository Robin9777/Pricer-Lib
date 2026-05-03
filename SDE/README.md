# SDE

Stochastic process simulation layer. All concrete processes derive from `RandomProcess`.

## Class hierarchy

```
RandomProcess            (abstract — Simulate(t0, T, nbSteps))
├── Brownian1D
├── BrownianND
├── BlackScholes1D       (abstract — adds Spot, Rate, Vol)
│   ├── BSEuler1D
│   └── BSMilstein1D
├── BlackScholes2D       (abstract — 2 correlated assets)
│   └── BSMilstein2D
├── BlackScholesND       (abstract — N correlated assets)
│   └── BSMilsteinND
└── Heston
```

---

## Core classes

### `SinglePath`
Stores a simulated path as a time series of doubles.

- `InsertValue(val)` — append a value
- `GetState(time)` — retrieve value at a given time
- `GetAllValues()` — full path as `vector<double>`

---

### `RandomProcess`
Abstract base for all processes. Holds a `RandomGenerator*`, a dimension count, and a `vector<SinglePath*>`.

- `Simulate(t0, T, nbSteps)` — fills the paths (pure virtual)
- `GetPath(d)` — returns `SinglePath*` for dimension d
- `GetDimension()` — number of assets/dimensions

---

### `Brownian1D` / `BrownianND`
Pure Brownian motion with no drift.

`BrownianND` uses Cholesky decomposition to produce correlated increments:

```cpp
BrownianND proc(&gen, N, &correlationMatrix);
proc.Simulate(0.0, 1.0, 52);
```

---

### `BSEuler1D`
Euler-Maruyama scheme on one Black-Scholes asset:

```
S_{t+dt} = S_t + r·S_t·dt + σ·S_t·dW
```

```cpp
BSEuler1D proc(&normGen, spot, rate, vol);
proc.Simulate(0.0, 1.0, 50);
```

---

### `BSMilstein1D`
Milstein scheme on one asset (higher-order correction):

```
S_{t+dt} = S_t + r·S_t·dt + σ·S_t·dW + 0.5·σ²·S_t·(dW² − dt)
```

```cpp
BSMilstein1D proc(&normGen, spot, rate, vol);
```

---

### `BSMilstein2D`
Milstein scheme on 2 correlated assets. Uses a single correlation coefficient `ρ`.

```cpp
BSMilstein2D proc(&normGen, spot1, spot2, rate, vol1, vol2, rho);
```

---

### `BSMilsteinND`
Milstein scheme on N correlated assets. Generalises `BSMilstein2D` to any N using `BrownianND` (Cholesky) internally.

```cpp
BSMilsteinND proc(&normGen,
    {100.0, 100.0, 100.0},             // spots
    0.05,                               // rate
    {0.2, 0.2, 0.2},                   // vols
    {{1,0.3,0.3},{0.3,1,0.3},{0.3,0.3,1}});  // correlation
proc.Simulate(0.0, 1.0, 50);
// paths: proc.GetPath(0), GetPath(1), GetPath(2)
```

---

### `Heston`
Stochastic volatility model. Simulates 2 correlated processes: asset price S and variance V.

```
dS = μ·S·dt + √V·S·dW₁
dV = κ(θ − V)dt + σ·√V·dW₂     corr(dW₁, dW₂) = ρ
```

```cpp
Heston proc(&normGen, spot, initVariance, mu, theta, kappa, sigma, rho);
```

Parameters:
- `initVariance` — initial variance V₀
- `mu` — drift of the asset
- `theta` — long-run variance (mean reversion target)
- `kappa` — mean reversion speed
- `sigma` — vol of vol
- `rho` — correlation between asset and variance processes
