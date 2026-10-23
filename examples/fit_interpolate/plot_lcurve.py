"""L-curve alpha selection for mesh fitting.

Demonstrates the L-curve criterion used by ANUGA's Fit class to automatically
choose the smoothing parameter alpha when fitting scattered point data onto an
unstructured triangular mesh.

The L-curve is a log-log plot of the residual norm (data misfit) against the
smoothing norm (roughness of the fitted surface) as alpha varies.  The optimal
alpha sits at the corner of the L — balancing fit quality against smoothness.

The example uses two overdetermined setups (more data than mesh nodes) with
different noise levels so that the L-curve corner is clearly visible.

Run::

    python plot_lcurve.py

Requires matplotlib and scipy.
"""

import numpy as num
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

import anuga.fit_interpolate.fitsmooth_ext as fitsmooth
from anuga.fit_interpolate.fit import Fit

# ---------------------------------------------------------------------------
# Build a coarse 7x7 triangular mesh on [0,1]x[0,1]  (49 nodes, 72 triangles)
# ---------------------------------------------------------------------------
nx, ny = 7, 7
xs = num.linspace(0, 1, nx)
ys = num.linspace(0, 1, ny)
XX, YY = num.meshgrid(xs, ys)
vertices = num.column_stack([XX.ravel(), YY.ravel()])

triangles = []
for j in range(ny - 1):
    for i in range(nx - 1):
        bl = j * nx + i
        br = bl + 1
        tl = bl + nx
        tr = tl + 1
        triangles.append([bl, br, tr])
        triangles.append([bl, tr, tl])
triangles = num.array(triangles, dtype=num.int64)
n_nodes = nx * ny                      # 49

# ---------------------------------------------------------------------------
# Scattered data: true surface + two noise levels
# Both are overdetermined (n_pts >> n_nodes) so the residual is genuinely > 0
# ---------------------------------------------------------------------------
n_pts = 500
rng = num.random.default_rng(42)
xy = rng.uniform(0.05, 0.95, (n_pts, 2))
z_true = num.sin(3 * xy[:, 0]) * num.cos(3 * xy[:, 1])

noise_low  = 0.1
noise_high = 1.0
z_low  = z_true + rng.standard_normal(n_pts) * noise_low
z_high = z_true + rng.standard_normal(n_pts) * noise_high


# ---------------------------------------------------------------------------
# Helper: extract AtA and D from a Fit object as scipy sparse matrices
# ---------------------------------------------------------------------------
def extract_scipy_matrices(fit_obj):
    m = fit_obj.mesh.number_of_nodes

    def _to_sp(cap):
        raw = fitsmooth.dok_to_csr(cap)
        data_   = num.array(raw[0], dtype=float)
        colind_ = num.array(raw[1], dtype=num.int64)
        rp      = num.array(raw[2], dtype=num.int64)
        n_ent = len(data_)
        if len(rp) < m + 1:
            rp = num.concatenate(
                [rp, num.full(m + 1 - len(rp), n_ent, dtype=num.int64)])
        rp = num.clip(rp, 0, n_ent)
        for k in range(len(rp) - 2, -1, -1):
            if rp[k] > rp[k + 1]:
                rp[k] = rp[k + 1]
        return sp.csr_matrix((data_, colind_, rp), shape=(m, m))

    return _to_sp(fit_obj.AtA), _to_sp(fit_obj.D)


# ---------------------------------------------------------------------------
# L-curve computation over a custom alpha grid
# ---------------------------------------------------------------------------
# For overdetermined systems the L-corner typically sits at alpha ~ O(1).
# A log-spaced grid from 1e-3 to 10 captures the full curve.
alphas = num.logspace(-4, 4, 25)


def compute_lcurve(z):
    f_obj = Fit(vertices, triangles, alpha='auto')
    f_obj._build_matrix_AtA_Atz(xy, z)
    AtA_sp, D_sp = extract_scipy_matrices(f_obj)

    Atz   = num.asarray(f_obj.Atz)
    z_sq  = float(f_obj.z_sq)

    log_rss = num.zeros(len(alphas))
    log_s   = num.zeros(len(alphas))

    for i, alpha in enumerate(alphas):
        B  = AtA_sp + float(alpha) * D_sp
        f  = spla.spsolve(B, Atz)
        # Numerically stable RSS via the normal equations:
        # (AtA + alpha*D) f = Atz  =>  f^T AtA f = f·Atz - alpha*f_D_f
        # RSS = z_sq - 2*f·Atz + f^T AtA f = z_sq - f·Atz - alpha*f_D_f
        Atz_f = float(Atz @ f)
        f_D_f = float(f @ (D_sp @ f))
        rss   = max(z_sq - Atz_f - float(alpha) * f_D_f, 1e-300)
        log_rss[i] = num.log10(rss)
        log_s[i]   = num.log10(max(f_D_f, 1e-300))

    # L-curve corner: maximum curvature
    x, y  = log_rss, log_s
    dx    = num.gradient(x)
    ddx   = num.gradient(dx)
    dy    = num.gradient(y)
    ddy   = num.gradient(dy)
    denom = num.where((dx**2 + dy**2)**1.5 < 1e-300, 1e-300,
                      (dx**2 + dy**2)**1.5)
    kappa = (ddx * dy - dx * ddy) / denom
    best  = int(num.argmax(kappa))
    opt   = float(alphas[best]) if kappa[best] > 0 else None

    # Also run ANUGA's built-in auto selection
    anuga_opt = f_obj.select_alpha()

    return log_rss, log_s, kappa, opt, anuga_opt


log_rss_low,  log_s_low,  kappa_low,  opt_low,  anuga_low  = compute_lcurve(z_low)
log_rss_high, log_s_high, kappa_high, opt_high, anuga_high = compute_lcurve(z_high)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(13, 9))
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.38)

log_a = num.log10(alphas)


def plot_lcurve_panel(ax, log_rss, log_s, kappa, opt, anuga_opt, title, cmap_name):
    colors = cm.get_cmap(cmap_name)(num.linspace(0.2, 0.9, len(alphas)))
    for i in range(len(alphas) - 1):
        ax.plot(log_rss[i:i+2], log_s[i:i+2], '-', color=colors[i], lw=1.5)
    sc = ax.scatter(log_rss, log_s, c=log_a, cmap=cmap_name,
                    vmin=log_a[0], vmax=log_a[-1], s=40, zorder=5)
    cbar = plt.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label(r'$\log_{10}(\alpha)$', fontsize=9)

    if opt is not None:
        idx = int(num.argmin(num.abs(alphas - opt)))
        ax.plot(log_rss[idx], log_s[idx], '*', color='red', ms=16, zorder=6,
                label=f'L-corner α≈{opt:.2g}')
        ax.legend(fontsize=9)

    ax.set_xlabel(r'$\log_{10}$ residual norm  $\log\,\|z-Af\|^2$', fontsize=10)
    ax.set_ylabel(r'$\log_{10}$ smoothing norm  $\log\,f^\top D f$', fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.grid(True, alpha=0.3)


def plot_kappa_panel(ax, kappa, opt, anuga_opt, title, color):
    ax.plot(log_a, kappa, 'o-', color=color, lw=1.8, ms=5)
    ax.axhline(0, color='k', lw=0.7, ls=':')
    if opt is not None:
        ax.axvline(num.log10(opt), color='red', ls='--', lw=1.5,
                   label=f'L-corner α≈{opt:.2g}')
    if anuga_opt is not None:
        ax.axvline(num.log10(anuga_opt), color='green', ls=':', lw=1.5,
                   label=f'ANUGA auto α={anuga_opt:.0e}')
    ax.legend(fontsize=8)
    ax.set_xlabel(r'$\log_{10}(\alpha)$', fontsize=10)
    ax.set_ylabel(r'Curvature $\kappa$', fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.grid(True, alpha=0.3)


plot_lcurve_panel(
    fig.add_subplot(gs[0, 0]),
    log_rss_low, log_s_low, kappa_low, opt_low, anuga_low,
    f'L-curve — low noise (σ={noise_low})\n{n_pts} pts, {n_nodes} nodes',
    'Blues')

plot_lcurve_panel(
    fig.add_subplot(gs[0, 1]),
    log_rss_high, log_s_high, kappa_high, opt_high, anuga_high,
    f'L-curve — high noise (σ={noise_high})\n{n_pts} pts, {n_nodes} nodes',
    'Oranges')

plot_kappa_panel(
    fig.add_subplot(gs[1, 0]),
    kappa_low, opt_low, anuga_low,
    'Curvature κ vs α — low noise', 'steelblue')

plot_kappa_panel(
    fig.add_subplot(gs[1, 1]),
    kappa_high, opt_high, anuga_high,
    'Curvature κ vs α — high noise', 'darkorange')

fig.suptitle(
    'L-curve criterion for smoothing parameter selection\n'
    r'Fitting $z = \sin(3x)\cos(3y)$ + noise onto a $7\times7$ mesh'
    f' ({n_nodes} nodes, {n_pts} data points)',
    fontsize=12, y=1.01)

plt.savefig('lcurve_example.png', dpi=150, bbox_inches='tight')
print(f'Low-noise  L-corner α ≈ {opt_low}   ANUGA auto: {anuga_low:.2e}')
print(f'High-noise L-corner α ≈ {opt_high}   ANUGA auto: {anuga_high:.2e}')
print('Plot saved to lcurve_example.png')
plt.show()
