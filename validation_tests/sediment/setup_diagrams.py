"""Schematic diagrams of the flumes the experimental sediment cases reproduce.

Each function draws one apparatus to scale in the vertical or the plan, with
the dimensions and flow conditions the case uses, and writes `setup.png` in
the case directory. Called by each case's `plot_results.py`, so the figure is
regenerated with the rest of the report and never committed.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle

WATER = '#cfe3f2'
SAND = '#e8d8a8'
RIGID = '#9a9a9a'
ARROW = dict(arrowstyle='->', lw=1.2, color='0.25')


def _finish(fig, ax, path, xlabel='distance (m)', ylabel='height (m)', equal=False):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if equal:
        ax.set_aspect('equal')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _dim(ax, x0, x1, y, label, colour='0.25', fs=8):
    """A dimension line with a centred label."""
    ax.annotate('', xy=(x0, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle='<->', lw=0.9, color=colour))
    ax.text(0.5 * (x0 + x1), y, label, ha='center', va='bottom', fontsize=fs, color=colour)


def pickup_flume(path='setup.png'):
    """Van Rijn's pick-up flume: clear water onto a sand bed (side view)."""
    L, h, d = 12.0, 0.25, 0.23e-3
    fig, ax = plt.subplots(figsize=(9, 3.0))
    ax.add_patch(Rectangle((0, 0), L, h, fc=WATER, ec='none'))
    ax.add_patch(Rectangle((0, -0.05), 2.0, 0.05, fc=RIGID, ec='k', lw=0.6))
    ax.add_patch(Rectangle((2.0, -0.05), L - 2.0, 0.05, fc=SAND, ec='k', lw=0.6))
    ax.plot([0, L], [h, h], color='#2b6ca3', lw=1.4)
    for y in (0.05, 0.11, 0.17):
        ax.annotate('', xy=(1.75, y), xytext=(0.25, y), arrowprops=ARROW)
    ax.text(0.25, h - 0.028, r'clear water in at $\bar u = 0.67$ m/s', fontsize=9, color='#2b6ca3')
    ax.text(1.0, -0.025, 'rigid', ha='center', va='center', fontsize=8)
    ax.text(7.0, -0.025, 'sand bed, $d_{50}$ = 230 μm', ha='center', va='center', fontsize=8)
    for x, lab in ((2.0 + 4 * 0.25, '4d'), (2.0 + 10 * 0.25, '10d'), (2.0 + 20 * 0.25, '20d'), (2.0 + 40 * 0.25, '40d')):
        ax.plot([x, x], [0, h], color='0.5', ls=':', lw=0.8)
        ax.text(x, h + 0.012, lab, ha='center', fontsize=8, color='0.35')
    _dim(ax, 0, 2.0, -0.085, 'rigid approach 2 m')
    _dim(ax, 2.0, L, -0.085, 'erodible bed')
    ax.annotate('', xy=(11.6, 0.0), xytext=(11.6, h), arrowprops=dict(arrowstyle='<->', lw=0.9, color='#2b6ca3'))
    ax.text(11.5, 0.5 * h, 'h = 0.25 m', fontsize=8, color='#2b6ca3', va='center', ha='right')
    ax.set_xlim(-0.3, L + 0.3)
    ax.set_ylim(-0.12, h + 0.05)
    ax.set_title("Van Rijn's pick-up flume: sediment entrained into clear water", fontsize=10)
    _finish(fig, ax, path)


def trench(path='setup.png', slope=10.0, depth=0.15, test=3):
    """Van Rijn's migrating trench (side view, vertical exaggerated)."""
    h0 = 0.39
    x_lip = 1.5 if test == 3 else (1.0 if test == 1 else 1.2)
    run = depth * slope
    if test == 1:
        x0, x1 = 1.0, 5.5
        depth, run = 0.175, 0.175 * 3
    elif test == 2:
        x0, x1 = 1.2, 3.8
        depth, run = 0.170, 0.170 * 7
    else:
        x0, x1 = 1.5, 4.5
        depth, run = 0.150, 0.150 * 10
    L = 11.0
    bed = [(0, 0), (x0 - depth * (run / depth) / (run / depth), 0)]
    prof = [(0, 0), (x0, -depth), (x1, -depth), (x1 + run, 0), (L, 0)]
    fig, ax = plt.subplots(figsize=(9, 3.2))
    ax.add_patch(Rectangle((0, -0.30), L, h0 + 0.30, fc=WATER, ec='none'))
    ax.add_patch(Polygon(prof + [(L, -0.30), (0, -0.30)], closed=True, fc=SAND, ec='k', lw=0.8))
    ax.plot([0, L], [h0, h0], color='#2b6ca3', lw=1.4)
    for y in (0.08, 0.20, 0.32):
        ax.annotate('', xy=(1.05, y), xytext=(0.15, y), arrowprops=ARROW)
    ax.text(0.15, h0 + 0.016, r'$h_0$ = 0.39 m at $\bar u_0$ = 0.51 m/s', fontsize=9, color='#2b6ca3')
    ax.text(0.5 * (x0 + x1), -depth - 0.045,
            'trench, 1:%d side slopes' % (1, 3, 7, 10)[test], ha='center', fontsize=9)
    ax.text(1.6, 0.30, 'suspended 0.03 kg/s/m\nbedload 0.01 kg/s/m', fontsize=8, color='0.3')
    ax.text(8.0, -0.075, '$d_{50}$ = 160 μm, $w_s$ = 0.013 m/s', ha='center', fontsize=8)
    _dim(ax, x0, x1, -depth - 0.10, 'floor %.1f m' % (x1 - x0))
    ax.set_xlim(-0.3, L + 0.3)
    ax.set_ylim(-0.30, h0 + 0.09)
    ax.set_title("Van Rijn's trench, test %d: suspended sand settling into a dredged trench" % test,
                 fontsize=10)
    _finish(fig, ax, path)


def settling_flume(path='setup.png'):
    """Wang and Ribberink's perforated-bed flume (side view)."""
    h = 0.2155
    fig, ax = plt.subplots(figsize=(9, 3.0))
    ax.add_patch(Rectangle((-10, 0), 30, h, fc=WATER, ec='none'))
    ax.add_patch(Rectangle((-10, -0.045), 10, 0.045, fc=RIGID, ec='k', lw=0.6))
    ax.add_patch(Rectangle((0, -0.045), 16, 0.045, fc='none', ec='k', lw=0.8, hatch='||'))
    ax.add_patch(Rectangle((16, -0.045), 4, 0.045, fc=RIGID, ec='k', lw=0.6))
    ax.add_patch(Rectangle((0, -0.115), 16, 0.07, fc='0.93', ec='k', lw=0.6))
    ax.text(8, -0.082, 'compartmented chamber: sand that reaches the bed is lost',
            ha='center', va='center', fontsize=8)
    ax.plot([-10, 20], [h, h], color='#2b6ca3', lw=1.4)
    for y in (0.05, 0.11, 0.17):
        ax.annotate('', xy=(-7.6, y), xytext=(-9.6, y), arrowprops=ARROW)
    ax.text(-6.9, 0.105, r'$\bar u = 0.558$ m/s', fontsize=9, va='center', color='#2b6ca3')
    ax.annotate('', xy=(-9.0, h), xytext=(-9.0, h + 0.05),
                arrowprops=dict(arrowstyle='->', lw=1.4, color='#8a6d1f'))
    ax.text(-8.8, h + 0.036, 'sand fed 70.8 kg/h', fontsize=8, color='#8a6d1f')
    ax.text(-5.0, -0.022, 'rigid, roughened (10 m)', ha='center', va='center', fontsize=8)
    ax.text(18.0, -0.022, 'out', ha='center', va='center', fontsize=8)
    _dim(ax, 0, 16, -0.142, 'perforated plate, 16 m')
    for x in (0, 1, 3, 6, 11, 15.5):
        ax.plot([x, x], [0, h], color='0.5', ls=':', lw=0.8)
    ax.text(7.0, h + 0.012, 'siphon profiles, 8 points each', fontsize=8, color='0.35')
    ax.annotate('', xy=(-9.95, 0.0), xytext=(-9.95, h), arrowprops=dict(arrowstyle='<->', lw=0.9, color='#2b6ca3'))
    ax.text(-10.3, 0.5 * h, 'h = 0.216 m', fontsize=8, color='#2b6ca3', rotation=90,
            va='center', ha='center')
    ax.text(16.6, 0.12, '$d_{50}$ = 100 μm\n$w_s$ = 0.7 cm/s', fontsize=8)
    ax.set_xlim(-10.9, 20.5)
    ax.set_ylim(-0.17, h + 0.075)
    ax.set_title('Wang and Ribberink: a suspension settling out through a porous bed', fontsize=10)
    _finish(fig, ax, path, xlabel='x from the start of the perforated bed (m)')


def pier_dambreak(path='setup.png', case='P1'):
    """The Zaragoza dam-break flume: side view above, plan of the pier below."""
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(9, 5.0),
                                 gridspec_kw=dict(height_ratios=[1.0, 1.15]))
    # --- side view ---
    ax.add_patch(Rectangle((-1.57, 0), 1.57, 0.08, fc=WATER, ec='none'))
    ax.plot([-1.57, 0], [0.08, 0.08], color='#2b6ca3', lw=1.4)
    ax.plot([0, 0], [0, 0.12], color='k', lw=2.5)
    ax.text(0.02, 0.125, 'pneumatic gate', fontsize=8)
    ax.add_patch(Rectangle((-1.57, -0.02), 1.57, 0.02, fc=RIGID, ec='k', lw=0.6))
    ax.add_patch(Rectangle((0, -0.02), 1.0, 0.02, fc=RIGID, ec='k', lw=0.6))
    ax.add_patch(Rectangle((1.0, -0.05), 1.5, 0.05, fc=SAND, ec='k', lw=0.8))
    ax.add_patch(Rectangle((2.5, -0.02), 0.5, 0.02, fc=RIGID, ec='k', lw=0.6))
    ax.plot([3.0, 3.5], [-0.02, -0.04], color=RIGID, lw=3)
    ax.text(3.25, -0.055, '4 % to traps', fontsize=8, ha='center')
    ax.plot([1.6, 1.6], [0.0, 0.055], color='k', lw=2)
    ax.text(1.6, 0.062, 'pier', ha='center', fontsize=8)
    ax.text(-0.785, 0.038, 'reservoir, 8 cm deep', ha='center', fontsize=8, color='#2b6ca3')
    _dim(ax, -1.57, 0.0, -0.085, '157 cm')
    _dim(ax, 0.0, 1.0, -0.085, '100 cm')
    _dim(ax, 1.0, 2.5, -0.085, 'erodible sand, 150 cm, 5 cm deep')
    ax.set_xlim(-1.66, 3.6)
    ax.set_ylim(-0.108, 0.16)
    ax.set_ylabel('height (m)')
    ax.set_title('The Zaragoza dam-break flume, case %s' % case, fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    # --- plan ---
    W = 0.24
    bx.add_patch(Rectangle((1.1, 0), 0.95, W, fc=SAND, ec='none'))
    bx.plot([1.1, 2.05], [0, 0], color='k', lw=1.2)
    bx.plot([1.1, 2.05], [W, W], color='k', lw=1.2)
    if case == 'P1':
        bx.add_patch(plt.Circle((1.6, 0.5 * W), 0.015, fc='0.4', ec='k'))
        bx.text(1.6, 0.5 * W + 0.028, 'cylinder, 3 cm', ha='center', fontsize=8)
    else:
        bx.add_patch(Rectangle((1.6 - 0.0175, 0.5 * W - 0.015), 0.035, 0.03, fc='0.4', ec='k'))
        bx.add_patch(plt.Circle((1.6 - 0.0175, 0.5 * W), 0.015, fc='0.4', ec='k'))
        bx.add_patch(plt.Circle((1.6 + 0.0175, 0.5 * W), 0.015, fc='0.4', ec='k'))
        bx.text(1.6, 0.5 * W + 0.028, 'rounded pier, 3 x 6.5 cm', ha='center', fontsize=8)
    for y in (0.05, 0.12, 0.19):
        bx.annotate('', xy=(1.32, y), xytext=(1.14, y), arrowprops=ARROW)
    bx.add_patch(Rectangle((1.12, 0.005), 0.93, W - 0.01, fill=False, ec='#2b6ca3', ls='--', lw=1.0))
    bx.text(1.95, W - 0.03, 'Kinect window', fontsize=8, color='#2b6ca3', ha='right')
    _dim(bx, 1.1, 1.6, -0.028, '1.6 m from the gate')
    bx.annotate('', xy=(2.03, 0), xytext=(2.03, W), arrowprops=dict(arrowstyle='<->', lw=0.9, color='0.25'))
    bx.text(2.045, 0.5 * W, '24 cm', fontsize=8, va='center')
    bx.set_xlim(1.05, 2.12)
    bx.set_ylim(-0.045, W + 0.05)
    bx.set_aspect('equal')
    bx.set_xlabel('x from the gate (m)')
    bx.set_ylabel('y (m)')
    bx.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def mesh_figure(sww, path='mesh.png', zoom=None, title=None, zoom_title='detail',
                counts=None):
    """Draw the mesh an sww file carries, optionally with a zoom panel.

    `zoom` is (x0, x1, y0, y1) in the sww's coordinates; `counts` is a list of
    (label, n_cells, side_mm) printed under the whole-domain panel.
    """
    import numpy as np
    from netCDF4 import Dataset
    with Dataset(sww) as f:
        x = f.variables['x'][:].astype(float)
        y = f.variables['y'][:].astype(float)
        vols = f.variables['volumes'][:]
    n = len(vols)
    if zoom is None:
        fig, axes = plt.subplots(1, 1, figsize=(11, 3.4))
        axes = [axes]
    else:
        fig, axes = plt.subplots(2, 1, figsize=(11, 6.6),
                                 gridspec_kw=dict(height_ratios=[1.0, 1.25]))
    ax = axes[0]
    ax.triplot(x, y, vols, lw=0.15, color='0.35')
    ax.set_aspect('equal')
    ax.set_title(title or ('mesh, %d triangles' % n), fontsize=10)
    ax.set_ylabel('y (m)')
    if counts:
        ax.text(0.01, -0.34, '   '.join('%s: %d cells, ~%d mm' % c for c in counts),
                transform=ax.transAxes, fontsize=8, color='0.3')
    if zoom is not None:
        bx = axes[1]
        bx.triplot(x, y, vols, lw=0.4, color='0.35')
        bx.set_xlim(zoom[0], zoom[1])
        bx.set_ylim(zoom[2], zoom[3])
        bx.set_aspect('equal')
        bx.set_title(zoom_title, fontsize=10)
        bx.set_ylabel('y (m)')
        ax.add_patch(Rectangle((zoom[0], zoom[2]), zoom[1] - zoom[0], zoom[3] - zoom[2],
                               fill=False, ec='#c1462f', lw=1.0))
        axes[-1].set_xlabel('x from the gate (m)')
    else:
        ax.set_xlabel('x (m)')
    for a in axes:
        a.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
