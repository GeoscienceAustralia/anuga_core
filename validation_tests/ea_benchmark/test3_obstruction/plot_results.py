"""Figures for the EA Test 3 report: set-up, levels and speeds, checks.

Exits non-zero if a check fails, so the validation report marks the case.
"""

import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import anuga  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import ea_common as ea  # noqa: E402
import checks  # noqa: E402

args = anuga.get_args()

# Set-up: the long profile through both points, and the inflow.
x_dem, y_dem, z_dem = ea.read_asc(ea.data_path('Test3', 'test3DEM.asc'))
t_bc, q_bc = ea.read_inflow(ea.data_path('Test3', 'Test3BC.csv'))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3.4), gridspec_kw={'width_ratios': [2, 1]})
inside = (x_dem >= 0.0) & (x_dem <= 300.0)
ax1.plot(x_dem[inside], z_dem[np.argmin(np.abs(y_dem - 50.0))][inside], 'k-')
for n, x in (('1', 150.0), ('2', 250.0)):
    ax1.axvline(x, color='C0', ls=':')
    ax1.text(x + 3, 10.4, 'point %s' % n, color='C0')
ax1.annotate('inflow', xy=(2.5, 10.5), xytext=(30, 10.62), arrowprops={'arrowstyle': '->'})
ax1.set_xlabel('x (m)')
ax1.set_ylabel('Bed level (m)')
ax1.set_title('Bed along the channel centreline (100 m wide)', fontsize=10)
ax1.grid(alpha=0.3)
ax2.plot(t_bc, q_bc, 'k-')
ax2.set_xlim(0, 60)
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Inflow (m$^3$/s)')
ax2.set_title('Inflow over the full width (zero after 40 s)', fontsize=10)
ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig('setup.png', dpi=120)
plt.close(fig)

# Levels, with the report's summary of the full SWE packages.
g = ea.load_gauges('test3_gauges.csv')
fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
for ax, n in zip(axes, ('1', '2')):
    ax.plot(g['time'], g['stage_' + n], 'C0-', label='ANUGA (%s)' % args.alg)
    ax.set_title('Point %s' % n, fontsize=10)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Level (m)')
    ax.grid(alpha=0.3)
axes[0].fill_between([110, 900], 9.98, 9.99, color='0.7', alpha=0.6,
                     label='report: ~9.98-9.99 m from ~120 s')
axes[1].fill_between([150, 900], 9.80, 9.81, color='0.7', alpha=0.6,
                     label='report: rises ~5-6 cm after ~150 s')
for ax in axes:
    ax.legend(fontsize=8, loc='lower right')
fig.suptitle('EA Test 3: water level', fontsize=11)
fig.tight_layout()
fig.savefig('stage_plot.png', dpi=120)
plt.close(fig)

ea.plot_gauges(g, ['1', '2'], 'speed', 'speed_plot.png', 'EA Test 3: speed', 'Speed (m/s)',
               label='ANUGA (%s)' % args.alg)

c = checks.evaluate(args.alg)
c.figure('checks.png', 'EA Test 3 checks (%s)' % args.alg)
c.report('    ')
sys.exit(1 if c.failed() else 0)
