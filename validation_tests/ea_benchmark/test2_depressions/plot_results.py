"""Figures for the EA Test 2 report: set-up, levels, checks.

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

# Set-up: the DEM with the inflow line and the output points, and the
# inflow hydrograph.
x_dem, y_dem, z_dem = ea.read_asc(ea.data_path('Test2', 'test2DEM.asc'))
t_bc, q_bc = ea.read_inflow(ea.data_path('Test2', 'Test2_BC.csv'))
points = ea.read_points(ea.data_path('Test2', 'Test2output.csv'))
keep_x = (x_dem >= 0.0) & (x_dem <= 2000.0)
keep_y = (y_dem >= 0.0) & (y_dem <= 2000.0)
z = z_dem[np.ix_(keep_y, keep_x)]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={'width_ratios': [1.2, 1]})
im = ax1.imshow(z[::4, ::4], origin='lower', cmap='terrain', extent=[0, 2000, 0, 2000])
ax1.contour(x_dem[keep_x][::4], y_dem[keep_y][::4], z[::4, ::4],
            levels=np.arange(7.4, 10.3, 0.1), colors='k', linewidths=0.3, alpha=0.5)
ax1.plot([0, 0], [1900, 2000], 'r-', lw=4, label='inflow', clip_on=False)
for n, (x, y) in points.items():
    ax1.plot(x, y, '+', color='k', ms=7)
    ax1.text(x + 30, y + 30, n, fontsize=9)
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.set_title('Ground level, inflow and the 16 depression centres', fontsize=10)
ax1.legend(fontsize=8, loc='lower left')
fig.colorbar(im, ax=ax1, shrink=0.8, label='Ground level (m)')
ax2.plot(t_bc / 60.0, q_bc, 'k-')
ax2.set_xlim(0, 120)
ax2.set_xlabel('Time (min)')
ax2.set_ylabel('Inflow (m$^3$/s)')
ax2.set_title('Inflow (97,200 m$^3$ in all)', fontsize=10)
ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig('setup.png', dpi=120)
plt.close(fig)

# One panel per depression, laid out as on the map (north at the top).
MAP_ORDER = [str(4 * col + row) for row in (4, 3, 2, 1) for col in range(4)]
g = ea.load_gauges('test2_gauges.csv')
ea.plot_gauges(g, MAP_ORDER, 'stage', 'stage_plot.png',
               'EA Test 2: water level at the depression centres (laid out as on the map)',
               'Level (m)', time_scale=3600.0, time_label='Time (h)', ncol=4,
               label='ANUGA (%s)' % args.alg)

c = checks.evaluate(args.alg)
c.figure('checks.png', 'EA Test 2 checks (%s)' % args.alg)
c.report('    ')
sys.exit(1 if c.failed() else 0)
