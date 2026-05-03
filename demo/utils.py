import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import norm

plt.style.use('dark_background')

METHOD_LABELS = {
    'plain':      'Plain MC',
    'antithetic': 'Antithetic',
    'staticcv':   'Static CV',
    'qmc':        'QMC (Halton)',
    'qmc_cv':     'QMC + Static CV',
}

METHOD_COLORS = {
    'plain':      '#555555',
    'antithetic': '#2196F3',
    'staticcv':   '#4CAF50',
    'qmc':        '#FF9800',
    'qmc_cv':     '#9C27B0',
}


def load_data(results_dir='results'):
    conv = pd.read_csv(os.path.join(results_dir, 'convergence.csv'))
    term = pd.read_csv(os.path.join(results_dir, 'terminal_plain.csv'))
    anti = pd.read_csv(os.path.join(results_dir, 'antithetic_pairs.csv'))
    return conv, term, anti


def bs_price(S0, K, r, sigma, T):
    d1 = (np.log(S0 / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    return np.exp(-r * T) * (S0 * np.exp(r * T) * norm.cdf(d1) - K * norm.cdf(d2))


def _methods_in(conv):
    order = list(METHOD_LABELS.keys())
    present = conv['method'].unique()
    return [m for m in order if m in present]


def plot_distribution(term, S0, K, r, sigma, T, save_dir='results', option_label=None):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    if option_label:
        fig.suptitle(option_label, fontsize=11, alpha=0.7)

    ax = axes[0]
    ax.hist(term['sT'], bins=80, density=True, alpha=0.6, color='steelblue', label='Simulated')
    mu_log  = np.log(S0) + (r - 0.5 * sigma**2) * T
    sig_log = sigma * np.sqrt(T)
    x = np.linspace(term['sT'].quantile(0.001), term['sT'].quantile(0.999), 400)
    ax.plot(x, stats.lognorm.pdf(x, s=sig_log, scale=np.exp(mu_log)),
            color='tomato', lw=2, label='Theoretical log-normal')
    ax.axvline(K, color='white', ls='--', lw=1.5, label=f'Strike K={K}')
    ax.set_xlabel('$S_T$')
    ax.set_ylabel('Density')
    ax.set_title('Distribution of $S_T$ (plain MC)')
    ax.legend()

    ax = axes[1]
    ax.hist(term['payoff'], bins=80, density=True, alpha=0.6, color='darkorange')
    otm_pct = (term['payoff'] == 0).mean() * 100
    ax.text(0.55, 0.80,
            f'{otm_pct:.1f}% of paths OTM\n(discounted payoff = 0)',
            transform=ax.transAxes, fontsize=10,
            bbox=dict(boxstyle='round', facecolor='#222222', alpha=0.8))
    ax.set_xlabel('Discounted payoff $e^{-rT}\\max(S_T - K, 0)$')
    ax.set_ylabel('Density')
    ax.set_title('Discounted payoff distribution')

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'fig_distribution.png'), dpi=150, bbox_inches='tight')
    plt.show()


def plot_antithetic(anti, K, save_dir='results', option_label=None, analytical_price=None):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    if option_label:
        fig.suptitle(option_label, fontsize=11, alpha=0.7)

    # Left: S_T scatter — negative correlation between paired terminal prices
    ax = axes[0]
    both_otm = (anti['payoff_z'] == 0) & (anti['payoff_negz'] == 0)
    one_itm  = (anti['payoff_z'] == 0) ^ (anti['payoff_negz'] == 0)
    both_itm = (anti['payoff_z'] >  0) & (anti['payoff_negz'] >  0)
    ax.scatter(anti['sT_z'][both_otm], anti['sT_negz'][both_otm],
               s=2, alpha=0.5, color='#4FC3F7', label='Both OTM')
    ax.scatter(anti['sT_z'][one_itm],  anti['sT_negz'][one_itm],
               s=4, alpha=0.4, color='orange',  label='One ITM')
    ax.scatter(anti['sT_z'][both_itm], anti['sT_negz'][both_itm],
               s=8, alpha=0.8, color='tomato',  label='Both ITM')
    ax.axhline(K, color='white', ls=':', lw=1, alpha=0.4)
    ax.axvline(K, color='white', ls=':', lw=1, alpha=0.4)
    corr_st = np.corrcoef(anti['sT_z'], anti['sT_negz'])[0, 1]
    ax.set_xlabel('$S_T(Z)$')
    ax.set_ylabel('$S_T(-Z)$')
    ax.set_title(f'$S_T$ pairs — corr = {corr_st:.3f}')
    ax.legend(markerscale=3, fontsize=9)

    # Right: running price estimate — antithetic converges faster
    ax = axes[1]
    plain_pay = anti['payoff_z'].values
    anti_avg  = (anti['payoff_z'].values + anti['payoff_negz'].values) / 2.0
    ns        = np.arange(1, len(plain_pay) + 1)

    plain_running = np.cumsum(plain_pay) / ns
    anti_running  = np.cumsum(anti_avg)  / ns

    ax.plot(ns, plain_running, color='#777777', lw=1,   alpha=0.9, label='Plain MC')
    ax.plot(ns, anti_running,  color='#2196F3', lw=1.2, alpha=0.95, label='Antithetic')

    ref = analytical_price if analytical_price is not None else anti_running[-1]
    ax.axhline(ref, color='#FFEB3B', ls='--', lw=1.5, label=f'True price ({ref:.4f})')

    var_plain = plain_pay.var()
    var_anti  = anti_avg.var()
    reduction = (1 - var_anti / var_plain) * 100 if var_plain > 0 else 0
    ax.text(0.97, 0.05,
            f'Var reduction: {reduction:.1f}%\nPlain std: {plain_pay.std():.3f}\nAnti std:  {anti_avg.std():.3f}',
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='#1a1a1a', alpha=0.85))

    ax.set_xlabel('Number of pairs')
    ax.set_ylabel('Running price estimate')
    ax.set_title('Convergence — antithetic oscillates less')
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'fig_antithetic.png'), dpi=150, bbox_inches='tight')
    plt.show()


def plot_convergence(conv, analytical_price=None, save_dir='results', option_label=None):
    fig, ax = plt.subplots(figsize=(10, 5))
    if option_label:
        fig.suptitle(option_label, fontsize=11, alpha=0.7)

    for method in _methods_in(conv):
        sub = conv[conv['method'] == method].sort_values('nbSim')
        ax.plot(sub['nbSim'], sub['price'], 'o-',
                label=METHOD_LABELS[method], color=METHOD_COLORS[method], lw=1.5, ms=4)

    if analytical_price is not None:
        ax.axhline(analytical_price, color='white', ls='--', lw=1.5,
                   label=f'Analytical ({analytical_price:.4f})')

    ax.set_xscale('log')
    ax.set_xlabel('Number of paths N')
    ax.set_ylabel('Price estimate')
    ax.set_title('Price convergence vs N')
    ax.legend()
    ax.grid(True, alpha=0.2, which='major')

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'fig_convergence.png'), dpi=150, bbox_inches='tight')
    plt.show()


def plot_ci_vs_n(conv, save_dir='results', option_label=None):
    fig, ax = plt.subplots(figsize=(10, 5))
    if option_label:
        fig.suptitle(option_label, fontsize=11, alpha=0.7)

    zero_var_labels = []

    for method in _methods_in(conv):
        sub   = conv[conv['method'] == method].sort_values('nbSim')
        valid = sub[sub['ci'] > 1e-12]

        if len(valid) == 0:
            # CI is numerically zero for all N — variance ≈ 0
            zero_var_labels.append(METHOD_LABELS[method])
            # Add a phantom line to the legend with a ≈0 annotation
            ax.plot([], [], color=METHOD_COLORS[method], lw=2,
                    label=f'{METHOD_LABELS[method]} (CI ≈ 0)')
            continue

        ax.loglog(valid['nbSim'], valid['ci'], 'o-',
                  label=METHOD_LABELS[method], color=METHOD_COLORS[method], lw=1.5, ms=4)

    plain = conv[conv['method'] == 'plain'].sort_values('nbSim')
    c0    = plain['ci'].iloc[-1] * np.sqrt(plain['nbSim'].iloc[-1])
    ns    = np.array([conv['nbSim'].min(), conv['nbSim'].max()])
    ax.loglog(ns, c0 / np.sqrt(ns), color='#FFEB3B', ls='--', lw=1.5,
              alpha=0.7, label='$1/\\sqrt{N}$ ref slope')

    if zero_var_labels:
        note = ', '.join(zero_var_labels)
        ax.text(0.98, 0.97,
                f'{note}:\nvariance ≈ 0, CI not plottable\n(control variate nearly perfect)',
                transform=ax.transAxes, fontsize=8.5, color='#cccccc',
                ha='right', va='top',
                bbox=dict(boxstyle='round', facecolor='#1a1a1a', alpha=0.85))

    ax.set_xlabel('Number of paths N')
    ax.set_ylabel('95% CI half-width')
    ax.set_title('CI width vs N — log-log')
    ax.legend(loc='lower left')
    ax.grid(True, alpha=0.2, which='major')

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'fig_ci_vs_n.png'), dpi=150, bbox_inches='tight')
    plt.show()


def plot_variance_bars(conv, save_dir='results', option_label=None):
    n_ref     = conv['nbSim'].max()
    at_n      = conv[conv['nbSim'] == n_ref].copy()
    plain_var = at_n[at_n['method'] == 'plain']['variance'].values[0]
    at_n['reduction'] = (1 - at_n['variance'] / plain_var) * 100
    at_n = at_n[at_n['method'].isin(METHOD_LABELS)].copy()

    # Clamp zero/NaN variance to a small positive floor so log scale works
    floor = plain_var * 1e-8
    at_n['var_plot'] = at_n['variance'].fillna(0).clip(lower=floor)

    fig, ax = plt.subplots(figsize=(9, 5))
    if option_label:
        fig.suptitle(option_label, fontsize=11, alpha=0.7)
    bars = ax.bar(
        [METHOD_LABELS[m] for m in at_n['method']],
        at_n['var_plot'],
        color=[METHOD_COLORS[m] for m in at_n['method']],
        alpha=0.85, edgecolor='white', linewidth=0.6
    )
    for bar, (_, row) in zip(bars, at_n.iterrows()):
        h = bar.get_height()
        if not np.isfinite(h) or h <= 0:
            continue
        ax.text(bar.get_x() + bar.get_width() / 2,
                h * 2,
                f"{row['reduction']:.1f}%",
                ha='center', va='bottom', fontsize=10, fontweight='bold', color='white')

    ax.set_yscale('log')
    ax.set_ylabel('Variance of the estimator (log scale)')
    ax.set_title(f'Variance by method — N = {n_ref:,}')
    ax.grid(True, alpha=0.2, axis='y', which='major')

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'fig_variance_bars.png'), dpi=150, bbox_inches='tight')
    plt.show()


def plot_required_n(conv, target_ci=0.10, save_dir='results', option_label=None):
    large   = conv[conv['nbSim'] >= 5000]
    var_est = large.groupby('method')['variance'].mean()

    rows = []
    for method in _methods_in(conv):
        if method not in var_est:
            continue
        var = var_est[method]
        if not np.isfinite(var) or var <= 0:
            # Variance ≈ 0: n_needed = 1 (can't get lower)
            rows.append({'method': method, 'label': METHOD_LABELS[method], 'var': var if np.isfinite(var) else 0.0, 'n_needed': 1})
            continue
        n = int(np.ceil((1.96 * np.sqrt(var) / target_ci) ** 2))
        rows.append({'method': method, 'label': METHOD_LABELS[method], 'var': var, 'n_needed': n})

    df       = pd.DataFrame(rows)
    plain_n  = df[df['method'] == 'plain']['n_needed'].values[0]
    df['speedup'] = plain_n / df['n_needed']

    print(f'Target: 95% CI ≤ {target_ci}\n')
    print(df[['label', 'var', 'n_needed', 'speedup']]
          .rename(columns={'label': 'Method', 'var': 'Variance',
                           'n_needed': 'Paths needed', 'speedup': 'Speedup vs plain'})
          .to_string(index=False, float_format=lambda x: f'{x:.1f}'))

    fig, ax = plt.subplots(figsize=(9, 5))
    if option_label:
        fig.suptitle(option_label, fontsize=11, alpha=0.7)
    bars = ax.bar(
        df['label'], df['n_needed'],
        color=[METHOD_COLORS[m] for m in df['method']],
        alpha=0.85, edgecolor='white', linewidth=0.6
    )
    for bar, (_, row) in zip(bars, df.iterrows()):
        h = bar.get_height()
        if not np.isfinite(h) or h <= 0:
            continue
        speedup_str = f"{row['speedup']:.1f}x" if np.isfinite(row['speedup']) else '∞x'
        ax.text(bar.get_x() + bar.get_width() / 2,
                h * 1.5,
                f"{row['n_needed']:,}\n({speedup_str})",
                ha='center', va='bottom', fontsize=9, color='white')

    ax.set_yscale('log')
    ax.set_ylabel('Paths required (log scale)')
    ax.set_title(f'Paths needed — 95% CI ≤ {target_ci}')
    ax.grid(True, alpha=0.2, axis='y', which='major')

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'fig_required_n.png'), dpi=150, bbox_inches='tight')
    plt.show()

    return df
