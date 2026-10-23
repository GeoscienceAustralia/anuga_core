"""
    Quick plot of the settling-basin outputs against the reference solution
"""
import json
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
import analytical_settling_basin as analytic
import numpy

with open('settling_basin_parameters.json') as f:
    prm = json.load(f)

t, c, z, w, uh, cc = analytic.read_centroid_series('settling_basin.sww')
x = cc[:, 0]
order = numpy.argsort(x)
q, v_s, d_star = prm['q'], prm['v_s'], prm['d_star']
xx = numpy.linspace(0.0, prm['L'], 400)

pyplot.clf()
pyplot.semilogy(x[order], c[-1][order], 'b.', label='numerical, t = %g s' % t[-1])
pyplot.semilogy(xx, analytic.concentration(xx, prm['c0'], q, v_s, d_star), 'r-',
                label='reference $c_0 e^{-x/L_s}$, $L_s$ = %.1f m'
                % analytic.settling_length(q, v_s, d_star))
pyplot.title('Steady concentration along the channel')
pyplot.xlabel('x (m)')
pyplot.ylabel('Concentration (volumetric)')
pyplot.legend(loc='best')
pyplot.savefig('concentration_plot.png')

pyplot.clf()
rate = (z[-1] - z[-3]) / (t[-1] - t[-3])
pyplot.semilogy(x[order], rate[order], 'b.', label='numerical, last %g s' % (t[-1] - t[-3]))
pyplot.semilogy(xx, analytic.bed_rise_rate(xx, prm['c0'], q, v_s, d_star, prm['porosity']),
                'r-', label=r'reference $d^* v_s c(x) / (1-\lambda)$')
pyplot.title('Bed-rise rate along the channel')
pyplot.xlabel('x (m)')
pyplot.ylabel('dz/dt (m/s)')
pyplot.legend(loc='best')
pyplot.savefig('bed_rise_plot.png')

pyplot.clf()
for k in (0, len(order) // 4, len(order) // 2, 3 * len(order) // 4):
    i = order[k]
    pyplot.plot(t, c[:, i] / prm['c0'], label='x = %.0f m' % x[i])
pyplot.title('Approach to the steady state')
pyplot.xlabel('Time (s)')
pyplot.ylabel('$c / c_0$')
pyplot.legend(loc='best')
pyplot.savefig('approach_plot.png')
