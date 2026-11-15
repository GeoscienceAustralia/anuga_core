"""Figures for the EA Test 6B report: set-up, levels and speeds, checks.

Exits non-zero if a check fails, so the validation report marks the case.
"""

import os
import sys

import matplotlib
matplotlib.use('Agg')
import anuga  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import ea_common as ea  # noqa: E402
import checks  # noqa: E402
import ea_dambreak  # noqa: E402

args = anuga.get_args()

ea_dambreak.plot_setup('6B', 'setup.png')

g = ea.load_gauges('test6b_gauges.csv')
ea.plot_gauges(g, checks.NAMES, 'stage', 'stage_plot.png', 'EA Test 6B: water level, first 10 min',
               'Level (m)', ncol=3, tmax=600.0,
               label='ANUGA (%s)' % args.alg)
ea.plot_gauges(g, checks.NAMES, 'speed', 'speed_plot.png', 'EA Test 6B: speed, first 10 min',
               'Speed (m/s)', ncol=3, tmax=600.0,
               label='ANUGA (%s)' % args.alg)

c = checks.evaluate(args.alg)
c.figure('checks.png', 'EA Test 6B checks (%s)' % args.alg)
c.report('    ')
sys.exit(1 if c.failed() else 0)
