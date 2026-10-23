"""
    Quick plot of the bed-hump outputs against the reference solution
"""
import json
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
import analytical_bed_hump as analytic
import numpy

with open('bed_hump_parameters.json') as f:
    prm = json.load(f)
hump = prm['hump']
w0, q, A_g, m, lam = prm['w0'], prm['q'], prm['A_g'], prm['m'], prm['porosity']

t, z, w, uh, cc = analytic.read_centroid_series('bed_hump.sww')
x = cc[:, 0]
order = numpy.argsort(x)
xx = numpy.linspace(200.0, 700.0, 2000)

pyplot.clf()
for k in (0, len(t) // 2, len(t) - 1):
    line, = pyplot.plot(xx, analytic.bed(xx, t[k], w0, q, A_g, m, lam, **hump), '-', lw=1.0)
    pyplot.plot(x[order], z[k][order], '.', ms=3, color=line.get_color(),
                label='t = %.0f s' % t[k])
pyplot.xlim(200.0, 700.0)
pyplot.title('Bed elevation: numerical (dots) and characteristic solution (lines)')
pyplot.xlabel('x (m)')
pyplot.ylabel('z (m)')
pyplot.legend(loc='best')
pyplot.savefig('bed_plot.png')

pyplot.clf()
z_ref = analytic.bed(x, t[-1], w0, q, A_g, m, lam, **hump)
pyplot.plot(x[order], (z[-1] - z_ref)[order], 'b.', ms=3)
pyplot.xlim(200.0, 700.0)
pyplot.title('Bed error at t = %.0f s' % t[-1])
pyplot.xlabel('x (m)')
pyplot.ylabel('z - z_ref (m)')
pyplot.savefig('bed_error_plot.png')

pyplot.clf()
pyplot.plot(x[order], w[-1][order], 'b.', ms=3, label='stage')
pyplot.axhline(w0, color='r', lw=0.8, label='boundary stage')
pyplot.title('Free surface at t = %.0f s' % t[-1])
pyplot.xlabel('x (m)')
pyplot.ylabel('w (m)')
pyplot.legend(loc='best')
pyplot.savefig('stage_plot.png')
