import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import setup_diagrams          # noqa: E402  the shared flume schematics

setup_diagrams.trench('setup.png', test=3)

"""
    Plot the trench bed profiles against van Rijn's measurements
"""
import json
import os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
import numpy
import flume_tools as F

init_all = numpy.loadtxt('van_rijn_1986_fig16_initial.csv', delimiter=',', comments='#')
meas_all = numpy.loadtxt('van_rijn_1986_fig16_measured.csv', delimiter=',', comments='#')

for test in (1, 2, 3):
    if not os.path.exists('trench_test%d.sww' % test):
        continue
    with open('trench_test%d_parameters.json' % test) as f:
        prm = json.load(f)
    t, c, z, w, uh, cc = F.read_centroid_series('trench_test%d.sww' % test)
    t = t * prm.get('morfac', 1.0)                   # morphological time
    xf = cc[:, 0] - prm['X_OFF']
    plane = prm['S'] * (prm['L'] - cc[:, 0])
    bins = numpy.arange(-0.05, 11.05, prm['dx'])
    init = init_all[init_all[:, 0] == test][:, 1:]
    meas = meas_all[meas_all[:, 0] == test][:, 1:]

    pyplot.clf()
    pyplot.plot(init[:, 0], -init[:, 1], 'k-', lw=1, label='initial bed')
    for k in (len(t) // 3, 2 * len(t) // 3, len(t) - 1):
        xm, zm = F.along_flume(xf, plane - z[k], bins)
        pyplot.plot(xm, -zm, '-', label='ANUGA, t = %.1f h' % (t[k] / 3600.0))
    pyplot.plot(meas[:, 0], -meas[:, 1], 'kx', label='measured, 15 h (van Rijn 1986b Fig. 16)')
    pyplot.xlim(0, 11)
    pyplot.title('Trench test %d: bed level below the flume bed' % test)
    pyplot.xlabel('x (m)')
    pyplot.ylabel('bed level (m)')
    pyplot.legend(loc='lower right', fontsize='small')
    pyplot.savefig('bed_profile_test%d.png' % test)

    pyplot.clf()
    xm, qm = F.along_flume(xf, prm['rho_s'] * c[-1] * uh[-1], bins)
    pyplot.plot(xm, qm / prm['qs_in'], 'b-')
    pyplot.xlim(0, 11)
    pyplot.title('Trench test %d: suspended load along the flume at t = %.1f h' % (test, t[-1] / 3600.0))
    pyplot.xlabel('x (m)')
    pyplot.ylabel('$q_s / q_{s,in}$')
    pyplot.savefig('suspended_load_test%d.png' % test)
