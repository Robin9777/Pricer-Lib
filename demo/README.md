# demo/

Two tools here: an interactive pricer (`run_demo`) and a study runner that feeds a Jupyter notebook (`run_study` + `analysis.ipynb`).

Build both with `./build.sh` from the project root.

---

## 1. Interactive pricer — `run_demo`

Just run it and answer the questions.

```bash
./demo/run_demo
```

It asks you: option type, number of assets, spots, vols, strike, rate, maturity, number of paths, variance reduction method. All inputs are validated so you can't crash it with bad values.

At the end it prints the price with a 95% confidence interval. If you picked a variance reduction method it also shows the plain MC baseline so you can compare.

**European or Bermudan?**

European = standard, exercised at maturity only.
Bermudan = can be exercised early at specific dates (Longstaff-Schwartz under the hood).

**1 asset or multiple?**

1 asset = vanilla call. More than 1 = basket call (weighted sum of spots). You set the weights and a correlation rho between assets.

**Variance reduction methods (European only):**

- `0. None` — plain Monte Carlo
- `1. QMC` — Halton low-discrepancy sequence instead of pseudo-random. Forces nbSteps=1.
- `2. Antithetic` — for each path with draw Z, also computes -Z. Reduces variance ~75%.
- `3. Static CV` — uses the geometric basket closed-form price as a control variate. Reduces variance ~99% for 1D, less for baskets.
- `4. QMC + Static CV` — combines both. Usually the best.

For Bermudan only `none` and `qmc` are available (the regression makes antithetic/CV harder to combine cleanly).

---

## 2. Config file mode — `run_demo <file.json>`

If you don't want to go through the menu every time, write a JSON file and pass it:

```bash
./demo/run_demo demo/example_config.json
```

It prices all the options in the file one by one. Useful when you want to compare multiple configs in one shot.

`example_config.json` already has 5 examples: 1D plain, 1D antithetic, 2D basket with CV, Bermudan, 3D basket.

**JSON format:**

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

It's an array so you can put as many options as you want in one file.

| field | what it does |
|---|---|
| `type` | `"european"` or `"bermudan"` |
| `spots` | array of initial spot prices, one per asset |
| `vols` | array of volatilities, same length as spots |
| `weights` | basket weights, must sum to 1. If omitted, equal weights. |
| `rho` | correlation between assets (single value applied to all pairs) |
| `strike` | strike price |
| `rate` | risk-free rate, e.g. `0.05` |
| `maturity` | in years |
| `nbSim` | number of Monte Carlo paths |
| `nbSteps` | time steps. For European, 1 is exact (no discretization error under GBM). For Bermudan, use more. |
| `varReduction` | `"none"`, `"qmc"`, `"antithetic"`, `"staticcv"`, `"qmc_cv"` |
| `exerciseDates` | array of exercise dates, required for Bermudan |
| `seed` | random seed for reproducibility |

---

## 3. Jupyter analysis — `run_study` + `analysis.ipynb`

This is for when you want graphs: variance by method, CI width vs N, convergence plots, etc.

**Step 1 — edit `study_config.json`**

Set your option parameters in the `"option"` block. The file has examples at the bottom for European, Bermudan, 2D basket, 3D basket, high-vol OTM — just copy the one you want into the active `"option"` block.

`nbSim_values` is the list of path counts for the convergence sweep. Each value becomes one data point on the graphs. Don't put too many or it'll take forever.

`nSamples` is how many paths to use for the terminal distribution and antithetic plots (separate from the sweep).

**Step 2 — run the study**

```bash
./demo/run_study demo/study_config.json
```

This creates `results/` with three CSV files:
- `convergence.csv` — price, CI, variance for every method × every N
- `terminal_plain.csv` — terminal spot prices and payoffs (for the distribution plot)
- `antithetic_pairs.csv` — paired paths (Z and -Z) for the antithetic plot

Takes maybe 30 seconds depending on your `nbSim_values`.

**Step 3 — open the notebook**

```bash
jupyter notebook demo/analysis.ipynb
```

In the first cell, update the parameters to match your config:

```python
S0    = 100.0
K     = 100.0
r     = 0.05
sigma = 0.2
T     = 1.0
option_label = 'European call, S=100, K=100, σ=20%, T=1y'
```

`option_label` is just what shows up in the plot titles.

Then hit **Run All**. You get 6 plots:

1. Distribution of $S_T$ vs theoretical log-normal — checks the simulation is correct
2. Antithetic pairing — left: scatter showing negative correlation between paired paths. Right: running price estimate, antithetic (blue) oscillates less than plain MC (grey)
3. Price convergence vs N — all methods converging to the analytical Black-Scholes price
4. CI width vs N on log-log scale — the 1/sqrt(N) reference slope, each method sits lower = better
5. Variance by method — bar chart on log scale, % reduction vs plain MC labeled
6. Paths needed to hit a target CI — shows the practical speedup of each method

All plots are also saved as PNGs in `results/`.

**Note on Static CV in the CI plot**

For a 1D European call, the geometric basket control variate is essentially a perfect hedge — the variance collapses to near zero. So Static CV and QMC+CV don't appear as lines on the CI plot (log of zero is undefined). They show in the legend as `CI ≈ 0` and a note explains it. This is correct behavior, not a bug.
