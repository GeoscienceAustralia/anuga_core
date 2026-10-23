"""
    Quick plot of the settling-tank outputs against the reference solution
"""
import json
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
import analytical_sediment_settling as analytic
import numpy

with open('settling_parameters.json') as f:
    prm = json.load(f)

t, c_all, z_all, w_all = analytic.read_centroid_series('sediment_settling.sww')
c_num = c_all.mean(axis=1)
dz_num = (z_all - z_all[0]).mean(axis=1)

tt = numpy.linspace(0.0, t.max(), 400)
c_ref, z_ref, h_ref = analytic.settling(tt, prm['h0'], prm['c0'], prm['v_s'],
                                        prm['d_star'], prm['porosity'])
c_exp, z_exp = analytic.settling_no_feedback(tt, prm['h0'], prm['c0'], prm['v_s'],
                                             prm['d_star'], prm['porosity'])

pyplot.clf()
pyplot.plot(t, c_num, 'b.', label='numerical')
pyplot.plot(tt, c_ref, 'r-', label='reference (with bed feedback)')
pyplot.plot(tt, c_exp, 'k--', lw=0.8, label='exponential (fixed depth)')
pyplot.title('Depth-averaged concentration in the tank')
pyplot.xlabel('Time (s)')
pyplot.ylabel('Concentration (volumetric)')
pyplot.legend(loc='best')
pyplot.savefig('concentration_plot.png')

pyplot.clf()
pyplot.plot(t, dz_num, 'b.', label='numerical')
pyplot.plot(tt, z_ref - z_ref[0], 'r-', label='reference (with bed feedback)')
pyplot.plot(tt, z_exp - z_exp[0], 'k--', lw=0.8, label='exponential (fixed depth)')
pyplot.title('Bed rise from deposition')
pyplot.xlabel('Time (s)')
pyplot.ylabel('Bed rise (m)')
pyplot.legend(loc='best')
pyplot.savefig('bed_rise_plot.png')
