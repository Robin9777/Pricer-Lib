# RandomGenerator

Random number generation layer. All classes derive from the abstract base `RandomGenerator`.

## Class hierarchy

```
RandomGenerator          (abstract — Generate() → double)
├── UniformGenerator     (abstract — uniform [0,1] sources)
│   ├── PseudoGenerator  (abstract — adds seed management)
│   │   ├── LinearCongruential
│   │   └── EcuyerCombined
│   └── HaltonGenerator
├── ContinuousGenerator  (abstract)
│   ├── Normal
│   └── Exponential
└── DiscreteGenerator    (abstract)
    ├── Bernoulli
    ├── Binomial
    ├── Poisson
    ├── FiniteSet
    └── HeadTail
```

Also: `AntitheticNormal` — wraps a `Normal`, used for antithetic variance reduction.

---

## Classes

### `RandomGenerator`
Abstract base. Every generator implements one method:
- `Generate() → double` — draw one sample
- `Mean(nbSim)`, `Variance(nbSim)` — empirical stats over nbSim draws

---

### `LinearCongruential`
Standard LCG: `x_{n+1} = (a·x_n + c) mod m`, output normalized to [0,1].

```cpp
LinearCongruential gen(seed, multiplier, increment, modulus);
// Park-Miller: LinearCongruential gen(42, 16807, 0, 2147483647);
double u = gen.Generate();
```

---

### `EcuyerCombined`
L'Ecuyer combined generator — combines two `LinearCongruential` instances for a longer period.

```cpp
LinearCongruential g1(42, 40014, 0, 2147483563);
LinearCongruential g2(42, 40692, 0, 2147483399);
EcuyerCombined gen(g1, g2);
```

---

### `HaltonGenerator`
Quasi-random (low-discrepancy) generator based on Van der Corput sequences in prime bases. Replaces `LinearCongruential` as the uniform source for reduced variance in Monte Carlo.

```cpp
HaltonGenerator gen({2, 3});           // 2D Halton, cycles bases
HaltonGenerator gen(2);               // single base
auto bases = HaltonGenerator::firstNPrimes(8);  // helper
```

**Important:** the number of bases must equal `nsteps × assets × 2` (Box-Muller pairs). Keep nsteps small (1–4) to avoid large prime bases producing extreme values.

---

### `Normal`
Transforms a `UniformGenerator` into N(μ, σ²) samples. Three algorithms selectable at call time:

```cpp
Normal gen(mu, sigma, uniformGen);
double z = gen.Generate();                        // default: Box-Muller
double z = gen.Generate(NormalAlgo::BoxMuller);
double z = gen.Generate(NormalAlgo::CentralLimitTheorem);
double z = gen.Generate(NormalAlgo::RejectionSampling);
```

---

### `Exponential`
Transforms a `UniformGenerator` into Exp(λ) samples.

```cpp
Exponential gen(lambda, uniformGen);
double x = gen.Generate();
double x = gen.Generate(ExpoAlgo::InverseTransform);
double x = gen.Generate(ExpoAlgo::RejectionSampling);
```

---

### `Bernoulli`, `Binomial`, `Poisson`
Discrete distributions. Each wraps a `UniformGenerator`.

```cpp
Bernoulli  bern(p, uniformGen);
Binomial   bino(n, bern);
Poisson    pois(lambda, uniformGen);
```

---

### `AntitheticNormal`
Wraps a `Normal` generator. On the normal pass, records each draw Z into a buffer. On the antithetic pass, replays −Z. Used by `AntitheticMCPricer`.

```cpp
Normal inner(0.0, 1.0, uniformGen);
AntitheticNormal anti(inner);

anti.ResetBuffer();       anti.SetAntithetic(false);
// ... simulate normal path ...
anti.ResetIndex();        anti.SetAntithetic(true);
// ... simulate antithetic path (−Z values) ...
```
