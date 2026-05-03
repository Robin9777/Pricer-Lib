# Payoffs

Payoff functions applied to simulated paths. All classes derive from the abstract base `PayOff`.

## Class hierarchy

```
PayOff                    (abstract — operator()(paths) → double)
├── EuropeanCallPayoff
└── EuroCallBasketPayOff
```

---

## Classes

### `PayOff`
Abstract base. One method:

```cpp
virtual double operator()(const std::vector<SinglePath*>& paths) const = 0;
```

Takes the full set of simulated paths (one per asset) and returns the undiscounted payoff.

---

### `EuropeanCallPayoff`
Vanilla European call on a single asset.

```
payoff = max(S_T − K, 0)
```

Uses the last value of `paths[0]`.

```cpp
EuropeanCallPayoff payoff(strike);
double p = payoff(paths);
```

---

### `EuroCallBasketPayOff`
Weighted basket European call on N assets.

```
payoff = max(Σ wᵢ · S_T^(i) − K, 0)
```

Uses the last value of each path. Weights must sum to 1.

```cpp
EuroCallBasketPayOff payoff(strike, {0.5, 0.5});       // 2-asset equal weight
EuroCallBasketPayOff payoff(strike, {0.7, 0.2, 0.1});  // 3-asset custom weight
double p = payoff(paths);
```
