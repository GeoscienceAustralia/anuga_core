"""Figures for the EA Test 5 report: set-up, levels and speeds, checks.

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

# Set-up: the valley with the modelled area, the inflow line and the
# output points, and the inflow hydrograph.
x_dem, y_dem, z_dem = ea.read_asc(ea.data_path('Test5', 'Test5DEM.asc'))
t_bc, q_bc = ea.read_inflow(ea.data_path('Test5', 'Test5BC.csv'))
points = ea.read_points(ea.data_path('Test5', 'Test5output.csv'))
outline = np.array(ea.read_polygon(ea.data_path('Test5', 'Test5ActiveArea_region.shp')))
km = 1000.0
x0, y0 = x_dem[0], y_dem[0]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.2), gridspec_kw={'width_ratios': [1.6, 1]})
z = np.ma.masked_less(z_dem[::4, ::4], -1000.0)
im = ax1.imshow(z, origin='lower', cmap='terrain', vmin=140, vmax=260,
                extent=[0, (x_dem[-1] - x0) / km, 0, (y_dem[-1] - y0) / km])
ax1.plot(np.r_[outline[:, 0], outline[0, 0]] / km - x0 / km,
         np.r_[outline[:, 1], outline[0, 1]] / km - y0 / km, 'k-', lw=1)
ax1.plot([(232595 - x0) / km, (232785 - x0) / km], [(830480 - y0) / km, (830290 - y0) / km],
         'r-', lw=3, label='inflow')
for n, (x, y) in points.items():
    ax1.plot((x - x0) / km, (y - y0) / km, 'o', color='k', ms=4)
    ax1.text((x - x0) / km + 0.2, (y - y0) / km - 0.5, n, fontsize=9)
ax1.set_xlabel('x (km)')
ax1.set_ylabel('y (km)')
ax1.set_title('Valley, modelled area and output points', fontsize=10)
ax1.legend(fontsize=8, loc='lower right')
fig.colorbar(im, ax=ax1, shrink=0.7, label='Ground level (m)')
ax2.plot(t_bc / 60.0, q_bc, 'k-')
ax2.set_xlim(0, 120)
ax2.set_xlabel('Time (min)')
ax2.set_ylabel('Inflow (m$^3$/s)')
ax2.set_title('Inflow from the dam failure (zero after 100 min)', fontsize=10)
ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig('setup.png', dpi=120)
plt.close(fig)

g = ea.load_gauges('test5_gauges.csv')
ea.plot_gauges(g, checks.NAMES, 'stage', 'stage_plot.png', 'EA Test 5: water level',
               'Level (m)', time_scale=3600.0, time_label='Time (h)', ncol=4,
               label='ANUGA (%s)' % args.alg)
ea.plot_gauges(g, checks.NAMES, 'speed', 'speed_plot.png', 'EA Test 5: speed',
               'Speed (m/s)', time_scale=3600.0, time_label='Time (h)', ncol=4,
               label='ANUGA (%s)' % args.alg)

c = checks.evaluate(args.alg)
c.figure('checks.png', 'EA Test 5 checks (%s)' % args.alg)
c.report('    ')
sys.exit(1 if c.failed() else 0)
