"""Figures for the EA Test 6A report: set-up, levels and speeds, the limiter
study (run_limiter_study.py) and checks.

Exits non-zero if a check fails, so the validation report marks the case.
"""

import os
import sys

import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib.path import Path  # noqa: E402
import netCDF4  # noqa: E402
import numpy as np  # noqa: E402
import anuga  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import ea_common as ea  # noqa: E402
import checks  # noqa: E402
import ea_dambreak  # noqa: E402
from run_limiter_study import BETAS, FINAL_TIME  # noqa: E402

BETA_COLOUR = {1.0: 'C0', 0.5: 'C1', 0.0: 'C2'}
SNAPSHOTS = (6, 12, 18, 24, 30)                  # seconds
MEASURED_JUMP_G2 = (13.0, 18.0)                  # the report's Figure 4.25


def beta_label(beta):
    return {0.0: r'$\beta = 0$', 0.5: r'$\beta = 0.5$ (DE0)',
            1.0: r'$\beta = 1$ (DE1, DE_ader2)'}[beta]


def plot_limiter_gauges(filename):
    """Level at G1-G6 over the first 40 s for each beta."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 6.5), sharex=True)
    for beta in BETAS:
        g = ea.load_gauges(ea_dambreak.run_name('6A', beta)[0] + '_gauges.csv')
        for ax, name in zip(axes.flat, checks.NAMES):
            ax.plot(g['time'], g['stage_' + name], color=BETA_COLOUR[beta],
                    lw=1.2, label=beta_label(beta))
    for ax, name in zip(axes.flat, checks.NAMES):
        ax.set_title(name, fontsize=10)
        ax.grid(alpha=0.3)
        ax.set_xlim(0, FINAL_TIME)
    axes[0][1].axvspan(*MEASURED_JUMP_G2, color='0.85', zorder=0,
                       label='jump measured')
    for ax in axes[:, 0]:
        ax.set_ylabel('Level (m)')
    for ax in axes[-1]:
        ax.set_xlabel('Time (s)')
    axes[0][1].legend(fontsize=8, loc='lower right')
    fig.suptitle('EA Test 6A: water level by slope limiter, first %g s' % FINAL_TIME,
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(filename, dpi=120)
    plt.close(fig)


def plot_limiter_field(filename):
    """Depth, flow direction and the Froude = 1 line around the building."""
    geo = json.load(open('test6a_geometry.json'))
    building = Path(geo['building'])
    xlim, ylim = (-0.8, 7.0), (-1.88, 1.78)
    xg, yg = np.meshgrid(np.arange(-0.6, 7.0, 0.3), np.arange(-1.65, 1.7, 0.3))
    fig, axes = plt.subplots(len(SNAPSHOTS), len(BETAS), figsize=(13, 9.4),
                             sharex=True, sharey=True)
    for c, beta in enumerate(sorted(BETAS)):
        s = netCDF4.Dataset(ea_dambreak.run_name('6A', beta)[0] + '.sww')
        x = s['x'][:] + s.xllcorner
        y = s['y'][:] + s.yllcorner
        v = s['volumes'][:]
        mesh = mtri.Triangulation(x, y, v)
        # The centroids, for contours and arrows, without the building.
        cen = mtri.Triangulation(x[v].mean(1), y[v].mean(1))
        tc = np.column_stack((cen.x[cen.triangles].mean(1), cen.y[cen.triangles].mean(1)))
        cen.set_mask(building.contains_points(tc))
        z = s['elevation_c'][:]
        t = s['time'][:]
        for r, tt in enumerate(SNAPSHOTS):
            k = np.argmin(abs(t - tt))
            h = np.maximum(s['stage_c'][k] - z, 0.0)
            hs = np.maximum(h, 1.0e-3)
            u, w = s['xmomentum_c'][k] / hs, s['ymomentum_c'][k] / hs
            fr = np.where(h > 0.005, np.hypot(u, w) / np.sqrt(9.8 * hs), 0.0)
            ax = axes[r, c]
            im = ax.tripcolor(mesh, facecolors=h, vmin=0.0, vmax=0.2, cmap='Blues')
            ax.tricontour(cen, fr, levels=[1.0], colors='#c0392b', linewidths=0.9)
            ug = mtri.LinearTriInterpolator(cen, u)(xg, yg)
            wg = mtri.LinearTriInterpolator(cen, w)(xg, yg)
            ax.quiver(xg, yg, ug, wg, scale=40, width=0.0025, color='0.15')
            ax.fill(*np.array(geo['building']).T, color='0.35')
            for name in ('G1', 'G2', 'G3', 'G4', 'G5'):
                px, py = GAUGES[name]
                ax.plot(px, py, 'k+', ms=6)
                ax.text(px + 0.08, py + 0.08, name, fontsize=7)
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_aspect('equal')
            if r == 0:
                ax.set_title(beta_label(beta), fontsize=11)
            if c == 0:
                ax.set_ylabel('t = %d s' % tt, fontsize=11)
        s.close()
    fig.subplots_adjust(left=0.06, right=0.88, top=0.96, bottom=0.04,
                        wspace=0.05, hspace=0.12)
    cax = fig.add_axes([0.9, 0.35, 0.015, 0.4])
    fig.colorbar(im, cax=cax, label='Depth (m)')
    fig.text(0.895, 0.3, 'Red line:\nFroude number 1\n\nArrows:\nvelocity',
             fontsize=9, va='top')
    fig.savefig(filename, dpi=110)
    plt.close(fig)

args = anuga.get_args()

ea_dambreak.plot_setup('6A', 'setup.png')

g = ea.load_gauges('test6a_gauges.csv')
ea.plot_gauges(g, checks.NAMES, 'stage', 'stage_plot.png', 'EA Test 6A: water level, first 60 s',
               'Level (m)', ncol=3, tmax=60.0,
               label='ANUGA (%s)' % args.alg)
ea.plot_gauges(g, checks.NAMES, 'speed', 'speed_plot.png', 'EA Test 6A: speed, first 60 s',
               'Speed (m/s)', ncol=3, tmax=60.0,
               label='ANUGA (%s)' % args.alg)

GAUGES = ea.read_points(ea.data_path('Test6', 'Test6Aoutput.csv'))
plot_limiter_gauges('limiter_stage.png')
plot_limiter_field('limiter_field.png')

c = checks.evaluate(args.alg)
c.figure('checks.png', 'EA Test 6A checks (%s)' % args.alg)
c.report('    ')
sys.exit(1 if c.failed() else 0)
