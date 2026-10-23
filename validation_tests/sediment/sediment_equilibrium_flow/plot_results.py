"""
    Quick plot of the equilibrium-flow outputs against the reference solution
"""
import json
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
import analytical_equilibrium_flow as analytic
import numpy

with open('equilibrium_flow_parameters.json') as f:
    prm = json.load(f)

t, c, z, w, uh, cc = analytic.read_centroid_series('equilibrium_flow.sww')
x = cc[:, 0]
order = numpy.argsort(x)
q, v_s, d_star, c_eq = prm['q'], prm['v_s'], prm['d_star'], prm['c_eq']
xx = numpy.linspace(0.0, prm['L'], 400)

pyplot.clf()
pyplot.plot(x[order], c[-1][order], 'b.', label='numerical, t = %g s' % t[-1])
pyplot.plot(xx, analytic.concentration(xx, c_eq, q, v_s, d_star), 'r-',
            label='reference $c_{eq}(1 - e^{-x/L_s})$, $L_s$ = %.1f m'
            % analytic.settling_length(q, v_s, d_star))
pyplot.axhline(c_eq, color='k', ls='--', lw=0.8, label='$c_{eq}$ = %.4f' % c_eq)
pyplot.title('Steady concentration along the channel')
pyplot.xlabel('x (m)')
pyplot.ylabel('Concentration (volumetric)')
pyplot.legend(loc='best')
pyplot.savefig('concentration_plot.png')

pyplot.clf()
pyplot.plot(x[order], (w[-1] - z[-1])[order], 'b.', label='numerical depth')
pyplot.axhline(prm['h0'], color='r', label='normal depth')
pyplot.title('Depth along the channel at the end of the run')
pyplot.xlabel('x (m)')
pyplot.ylabel('Depth (m)')
pyplot.legend(loc='best')
pyplot.savefig('depth_plot.png')

pyplot.clf()
for k in (0, len(order) // 4, len(order) // 2, 3 * len(order) // 4):
    i = order[k]
    pyplot.plot(t, c[:, i] / c_eq, label='x = %.0f m' % x[i])
pyplot.title('Approach to the steady state')
pyplot.xlabel('Time (s)')
pyplot.ylabel('$c / c_{eq}$')
pyplot.legend(loc='best')
pyplot.savefig('approach_plot.png')
