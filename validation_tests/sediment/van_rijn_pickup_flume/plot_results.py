import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import setup_diagrams          # noqa: E402  the shared flume schematics

setup_diagrams.pickup_flume('setup.png')

"""
    Plot the pick-up flume outputs against van Rijn's measurements
"""
import json
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
import numpy
import flume_tools as F

with open('pickup_flume_parameters.json') as f:
    prm = json.load(f)

t, c, z, w, uh, cc = F.read_centroid_series('pickup_flume.sww')
x = cc[:, 0]
h0, rho_s, qs_inf = prm['h0'], prm['rho_s'], prm['qs_inf']
xd_meas, ratio_meas = F.read_measured('van_rijn_1986_fig3.csv')

bins = numpy.arange(0.0, prm['L'] + prm['dx'], prm['dx'])
xm, qsm = F.along_flume(x, rho_s * c[-1] * uh[-1], bins)

pyplot.clf()
pyplot.plot(xm / h0, qsm / qs_inf, 'b-', label='ANUGA, t = %g s' % t[-1])
pyplot.plot(xd_meas, ratio_meas, 'ko', label='measured, van Rijn (1986c) Fig. 3')
pyplot.axhline(1.0, color='k', ls='--', lw=0.8, label='$q_{s,\\infty}$ = %.3f kg/s/m' % qs_inf)
pyplot.xlim(0, 45)
pyplot.title('Depth-integrated suspended load along the flume')
pyplot.xlabel('x / d')
pyplot.ylabel('$q_s / q_{s,\\infty}$')
pyplot.legend(loc='best')
pyplot.savefig('suspended_load_plot.png')

pyplot.clf()
xm, cm = F.along_flume(x, c[-1], bins)
pyplot.plot(xm / h0, cm * rho_s * 1000.0, 'b-')
pyplot.xlim(0, 45)
pyplot.title('Depth-averaged concentration along the flume')
pyplot.xlabel('x / d')
pyplot.ylabel('concentration (mg/L)')
pyplot.savefig('concentration_plot.png')

pyplot.clf()
xm, dm = F.along_flume(x, w[-1] - z[-1], bins)
pyplot.plot(xm / h0, (dm - h0) * 1.0e6, 'b-')
pyplot.axhline(0.0, color='r', lw=0.8)
pyplot.xlim(0, 45)
pyplot.title('Depth error relative to the normal depth d = %g m at the end of the run' % h0)
pyplot.xlabel('x / d')
pyplot.ylabel('depth error (micrometres)')
pyplot.savefig('depth_plot.png')
