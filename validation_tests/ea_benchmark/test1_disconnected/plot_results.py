"""Figures for the EA Test 1 report: set-up, levels, checks.

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

# Set-up: the ground profile with the output points, and the imposed level.
x_dem, y_dem, z_dem = ea.read_asc(ea.data_path('Test1', 'test1DEM.asc'))
t_bc, h_bc = ea.read_inflow(ea.data_path('Test1', 'Test1BC.csv'))
points = ea.read_points(ea.data_path('Test1', 'Test1output.csv'))
row = z_dem[np.argmin(np.abs(y_dem - 50.0))]
keep = (x_dem >= 0.0) & (x_dem <= 700.0)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3.8))
ax1.fill_between(x_dem[keep], 9.5, row[keep], color='0.8')
ax1.plot(x_dem[keep], row[keep], 'k-', lw=1)
ax1.axhline(10.35, color='C0', ls='--', lw=1, label='peak level 10.35 m')
ax1.axhline(10.25, color='C0', ls=':', lw=1, label='sill 10.25 m')
for n, (x, y) in points.items():
    ax1.plot(x, np.interp(x, x_dem, row), 'o', color='C3', ms=5)
    ax1.text(x, np.interp(x, x_dem, row) - 0.07, 'Point ' + n, ha='center', fontsize=8)
ax1.set_xlim(0, 700)
ax1.set_ylim(9.55, 10.45)
ax1.set_xlabel('x (m)')
ax1.set_ylabel('Level (m)')
ax1.set_title('Ground profile (uniform across the 100 m width)', fontsize=10)
ax1.legend(fontsize=8, loc='lower left')
ax2.plot(t_bc / 3600.0, h_bc, 'k-')
ax2.set_xlim(0, 20)
ax2.set_xlabel('Time (h)')
ax2.set_ylabel('Level (m)')
ax2.set_title('Water level imposed along x = 0', fontsize=10)
ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig('setup.png', dpi=120)
plt.close(fig)

g = ea.load_gauges('test1_gauges.csv')
ea.plot_gauges(g, checks.NAMES, 'stage', 'stage_plot.png', 'EA Test 1: water level',
               'Level (m)', time_scale=3600.0, time_label='Time (h)',
               label='ANUGA (%s)' % args.alg)

c = checks.evaluate(args.alg)
c.figure('checks.png', 'EA Test 1 checks (%s)' % args.alg)
c.report('    ')
sys.exit(1 if c.failed() else 0)
