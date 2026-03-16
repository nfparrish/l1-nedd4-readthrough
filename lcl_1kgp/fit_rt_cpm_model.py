#!/usr/bin/env python3
"""
Histogram + 3-genotype Gaussian mixture fit for GSE226189 RT_CPM values.

Model: log(RT_CPM) ~ mixture of 3 log-normals
  Mixing weights: Hardy-Weinberg equilibrium using q from 1000 Genomes
  (q ≈ 0.302, from 1232 0/0 : 1028 0/1 : 242 1/1 in PRJNA851328)

Tests two models vs. a single-Gaussian null:
  Free 3G:    3 free log-means + 1 shared log-sigma  (k=4)
  Dosage 3G:  1 base log-mean, μ_01 = μ_00+log(2), μ_11 = μ_00+log(4)  (k=2)

Likelihood ratio test (LRT): dosage constraint vs. free means → χ²(2)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm, chi2
from scipy.optimize import minimize
from pathlib import Path

RNG = np.random.default_rng(42)
SCRIPT_DIR = Path(__file__).parent

# ── Load data ──────────────────────────────────────────────────────────────
import csv
rows = list(csv.DictReader(open(SCRIPT_DIR / "GSE226189_rt_cpm.tsv"), delimiter='\t'))
srr_ids = [r['sample'] for r in rows]
cpm = np.array([float(r['RT_CPM']) for r in rows])
log_cpm = np.log(cpm)
n = len(cpm)

# ── HWE genotype frequencies from 1000 Genomes ────────────────────────────
# PRJNA851328 selection metadata: 0/0=1232, 0/1=1028, 1/1=242 (N=2502)
n_00, n_01, n_11 = 1232, 1028, 242
N_1kg = n_00 + n_01 + n_11
q = (n_01 + 2 * n_11) / (2 * N_1kg)   # insertion allele frequency
p = 1 - q
w = np.array([p**2, 2*p*q, q**2])     # [0/0, 0/1, 1/1] HWE weights

print(f"L1 insertion allele freq (1KGP):  q = {q:.4f}")
print(f"HWE genotype weights:  0/0={w[0]:.4f}  0/1={w[1]:.4f}  1/1={w[2]:.4f}")
print(f"Expected counts in n={n}:  0/0≈{w[0]*n:.1f}  0/1≈{w[1]*n:.1f}  1/1≈{w[2]*n:.1f}")

# ── model helpers ──────────────────────────────────────────────────────────
def mixture_pdf(y, mus, s, weights):
    pdf = sum(weights[i] * norm.pdf(y, mus[i], s) for i in range(3))
    return np.clip(pdf, 1e-300, None)

def nll_free(params, y, weights):
    """3-component, 3 free means, 1 shared sigma; weights fixed."""
    mu0, mu1, mu2, log_s = params
    s = np.exp(log_s)
    pdf = mixture_pdf(y, [mu0, mu1, mu2], s, weights)
    return -np.sum(np.log(pdf))

def nll_dosage(params, y, weights):
    """Dosage constraint: μ_01 = μ_00+log2, μ_11 = μ_00+log4."""
    mu0, log_s = params
    s = np.exp(log_s)
    mus = [mu0, mu0 + np.log(2), mu0 + np.log(4)]
    pdf = mixture_pdf(y, mus, s, weights)
    return -np.sum(np.log(pdf))

def nll_null(params, y):
    """Single Gaussian."""
    mu, log_s = params
    return -np.sum(norm.logpdf(y, mu, np.exp(log_s)))

# ── Fit all three models ───────────────────────────────────────────────────
# Initialization: spread initial means across log-CPM range
lq = np.percentile(log_cpm, [20, 55, 92])
log_s_init = np.log(np.std(log_cpm) * 0.45)

res_free = minimize(
    nll_free, [lq[0], lq[1], lq[2], log_s_init],
    args=(log_cpm, w), method='Nelder-Mead',
    options={'maxiter': 20000, 'xatol': 1e-7, 'fatol': 1e-7}
)
res_dosage = minimize(
    nll_dosage, [lq[0], log_s_init],
    args=(log_cpm, w), method='Nelder-Mead',
    options={'maxiter': 10000, 'xatol': 1e-7, 'fatol': 1e-7}
)
res_null = minimize(
    nll_null, [np.mean(log_cpm), np.log(np.std(log_cpm))],
    args=(log_cpm,), method='Nelder-Mead'
)

# ── Extract + sort parameters ──────────────────────────────────────────────
mu_free_raw = np.array(res_free.x[:3])
order = np.argsort(mu_free_raw)
mu_free = mu_free_raw[order]
s_free = np.exp(res_free.x[3])

mu0_d = res_dosage.x[0]
s_d = np.exp(res_dosage.x[1])
mu_dosage = np.array([mu0_d, mu0_d + np.log(2), mu0_d + np.log(4)])

mu_null_fit, s_null_fit = res_null.x[0], np.exp(res_null.x[1])

# ── AIC / BIC ──────────────────────────────────────────────────────────────
ll_free   = -res_free.fun
ll_dosage = -res_dosage.fun
ll_null   = -res_null.fun

k_free, k_dosage, k_null = 4, 2, 2

def aic(ll, k): return 2 * k - 2 * ll
def bic(ll, k): return k * np.log(n) - 2 * ll

models = {
    'Null (1G)':   dict(ll=ll_null,   k=k_null,   aic=aic(ll_null, k_null),     bic=bic(ll_null, k_null)),
    'Dosage 3G':   dict(ll=ll_dosage, k=k_dosage, aic=aic(ll_dosage, k_dosage), bic=bic(ll_dosage, k_dosage)),
    'Free 3G':     dict(ll=ll_free,   k=k_free,   aic=aic(ll_free, k_free),     bic=bic(ll_free, k_free)),
}

print(f"\n{'Model':<14} {'k':>3} {'LogLik':>10} {'AIC':>10} {'BIC':>10}")
print("─" * 52)
for name, m in models.items():
    print(f"{name:<14} {m['k']:>3} {m['ll']:>10.2f} {m['aic']:>10.2f} {m['bic']:>10.2f}")

# ── Likelihood ratio test: dosage vs free ────────────────────────────────-
lrt_stat = 2 * (ll_free - ll_dosage)
lrt_df = k_free - k_dosage   # 2
lrt_p = chi2.sf(lrt_stat, lrt_df)
print(f"\nLRT (Free 3G vs Dosage 3G):  χ²({lrt_df}) = {lrt_stat:.2f},  p = {lrt_p:.4f}")
if lrt_p > 0.05:
    print("  → Dosage constraint NOT rejected (data consistent with linear dosage model)")
else:
    print("  → Dosage constraint REJECTED (expression does not scale linearly with copy number)")

# ── Component summary ──────────────────────────────────────────────────────
print(f"\nFitted component means (Free 3G model):")
for i, (gt, mu) in enumerate(zip(['0/0', '0/1', '1/1'], mu_free)):
    print(f"  {gt}:  log_μ={mu:.3f}  →  CPM={np.exp(mu):.3f}  (σ={s_free:.3f} on log scale)")
print(f"  Fold-change 0/1 vs 0/0: {np.exp(mu_free[1]-mu_free[0]):.2f}×  (linear dosage predicts 2.00×)")
print(f"  Fold-change 1/1 vs 0/0: {np.exp(mu_free[2]-mu_free[0]):.2f}×  (linear dosage predicts 4.00×)")

# ── Posterior MAP genotype assignments ────────────────────────────────────
def posteriors(y, mus, s, weights):
    liks = np.array([w_i * norm.pdf(y, mu, s) for w_i, mu in zip(weights, mus)])
    return (liks / liks.sum(axis=0)).T   # shape (n, 3)

posts_free = posteriors(log_cpm, mu_free, s_free, w)
map_gt = np.argmax(posts_free, axis=1)
gt_labels = ['0/0', '0/1', '1/1']
n_assigned = [np.sum(map_gt == i) for i in range(3)]

print(f"\nMAP genotype assignments (Free 3G posterior):")
for i, (gt, ni) in enumerate(zip(gt_labels, n_assigned)):
    print(f"  {gt}: {ni:2d} samples  (HWE expected {w[i]*n:.1f})")

# Identify top outlier
top5_idx = np.argsort(cpm)[::-1][:5]
print(f"\nTop-5 samples by RT_CPM:")
print(f"  {'SRR':<15} {'RT_CPM':>8}  {'MAP_gt':<6}  {'P(1/1)':>8}")
for idx in top5_idx:
    print(f"  {srr_ids[idx]:<15} {cpm[idx]:>8.3f}  {gt_labels[map_gt[idx]]:<6}  {posts_free[idx,2]:>8.4f}")

# ── Plot ───────────────────────────────────────────────────────────────────
gt_colors = ['#2196F3', '#FF9800', '#E53935']
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# ── Panel 1: Raw CPM histogram ─────────────────────────────────────────────
ax = axes[0]
ax.hist(cpm, bins=np.linspace(0, cpm.max() * 1.02, 35), density=True,
        color='steelblue', alpha=0.75, edgecolor='white')
ax.set_xlabel('RT_CPM', fontsize=11)
ax.set_ylabel('Density', fontsize=11)
ax.set_title('Observed RT_CPM distribution\nGSE226189 fibroblasts (n=82)', fontsize=10)
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

# ── Panel 2: log-CPM histogram + mixture overlay ──────────────────────────
ax = axes[1]
ax.hist(log_cpm, bins=22, density=True,
        color='steelblue', alpha=0.75, edgecolor='white', label='Observed')

y_grid = np.linspace(log_cpm.min() - 0.3, log_cpm.max() + 0.3, 400)
total_pdf = np.zeros_like(y_grid)
for i, (mu, col, lab) in enumerate(zip(mu_free, gt_colors, gt_labels)):
    comp = w[i] * norm.pdf(y_grid, mu, s_free)
    ax.plot(y_grid, comp, color=col, lw=2, linestyle='--',
            label=f'{lab}  μ={np.exp(mu):.2f} CPM  (n={n_assigned[i]})')
    total_pdf += comp
ax.plot(y_grid, total_pdf, 'k-', lw=2.5, label='Mixture total')

ax.set_xlabel('log(RT_CPM)', fontsize=11)
ax.set_ylabel('Density', fontsize=11)
ax.set_title(f'log-CPM + 3-Genotype GMM\n(HWE weights, q={q:.3f}, free μ)', fontsize=10)
ax.legend(fontsize=8)
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

# Annotate LRT result
pval_str = f"LRT p={'<0.001' if lrt_p < 0.001 else f'{lrt_p:.3f}'}"
ax.text(0.03, 0.97, pval_str, transform=ax.transAxes,
        fontsize=8, va='top', style='italic',
        color='darkgreen' if lrt_p > 0.05 else 'firebrick')

# ── Panel 3: AIC/BIC bar chart ────────────────────────────────────────────
ax = axes[2]
model_names = list(models.keys())
aics = [m['aic'] for m in models.values()]
bics = [m['bic'] for m in models.values()]

x_pos = np.arange(len(model_names))
dx = 0.32
best_aic, best_bic = min(aics), min(bics)

b1 = ax.bar(x_pos - dx/2, aics, dx, label='AIC', color='#1976D2', alpha=0.82)
b2 = ax.bar(x_pos + dx/2, bics, dx, label='BIC', color='#F57C00', alpha=0.82)

for bar, val in zip(b1, aics):
    delta = val - best_aic
    txt = 'best' if delta < 0.1 else f'Δ{delta:.0f}'
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
            txt, ha='center', va='bottom', fontsize=8, color='#1565C0')

for bar, val in zip(b2, bics):
    delta = val - best_bic
    txt = 'best' if delta < 0.1 else f'Δ{delta:.0f}'
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
            txt, ha='center', va='bottom', fontsize=8, color='#E65100')

ax.set_xticks(x_pos)
ax.set_xticklabels(model_names, fontsize=9)
ax.set_ylabel('Information Criterion (lower = better)', fontsize=10)
ax.set_title('Model Fit Comparison\n(null vs dosage vs free 3-genotype)', fontsize=10)
ax.legend()
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

plt.suptitle('L1-NEDD4 Readthrough — GSE226189 Skin Fibroblasts\n'
             '3-Genotype Dosage Model (HWE weights from 1000 Genomes)',
             fontweight='bold', fontsize=12)
plt.tight_layout()

out_png = SCRIPT_DIR / 'GSE226189_rt_cpm_histogram.png'
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"\nFigure saved: {out_png}")
