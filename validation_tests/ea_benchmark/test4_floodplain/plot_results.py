"""Figures for the EA Test 4 report: set-up, levels and speeds, checks.

Exits non-zero if a check fails, so the validation report marks the case.
"""

import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import anuga  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import ea_common as ea  # noqa: E402
import checks  # noqa: E402

args = anuga.get_args()

# Set-up: the floodplain with the breach and the output points, and the inflow.
t_bc, q_bc = ea.read_inflow(ea.data_path('Test4', 'Test4BC.csv'))
points = ea.read_points(ea.data_path('Test4', 'Test4output.csv'))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.6), gridspec_kw={'width_ratios': [1, 1.4]})
ax1.plot([0, 1000, 1000, 0, 0], [0, 0, 2000, 2000, 0], 'k-')
ax1.plot([0, 0], [990, 1010], color='C3', lw=5, label='20 m breach')
for n, (x, y) in points.items():
    ax1.plot(x, y, 'o', color='C0', ms=4)
    up = n in '135' or n == '6'
    ax1.text(x, y + (40 if up else -40), n, color='C0', ha='center',
             va='bottom' if up else 'top', fontsize=9)
ax1.set_aspect('equal')
ax1.set_xlim(-60, 1060)
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.set_title('Flat floodplain, bed 0 m, n = 0.05', fontsize=10)
ax1.legend(fontsize=8, loc='upper right')
ax2.plot(t_bc / 3600.0, q_bc, 'k-')
ax2.set_xlim(0, 5)
ax2.set_xlabel('Time (h)')
ax2.set_ylabel('Inflow (m$^3$/s)')
ax2.set_title('Inflow through the breach', fontsize=10)
ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig('setup.png', dpi=120)
plt.close(fig)

g = ea.load_gauges('test4_gauges.csv')
ea.plot_gauges(g, checks.NAMES, 'depth', 'stage_plot.png', 'EA Test 4: depth (bed at 0 m)',
               'Depth (m)', time_scale=3600.0, time_label='Time (h)', ncol=3,
               label='ANUGA (%s)' % args.alg)
ea.plot_gauges(g, checks.NAMES, 'speed', 'speed_plot.png', 'EA Test 4: speed',
               'Speed (m/s)', time_scale=3600.0, time_label='Time (h)', ncol=3,
               label='ANUGA (%s)' % args.alg)

c = checks.evaluate(args.alg)
c.figure('checks.png', 'EA Test 4 checks (%s)' % args.alg)
c.report('    ')
sys.exit(1 if c.failed() else 0)
