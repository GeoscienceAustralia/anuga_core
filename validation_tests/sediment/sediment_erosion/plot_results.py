"""
    Quick plot of the erosion-tank outputs against the reference solution
"""
import json
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
import analytical_sediment_erosion as analytic
import numpy

with open('erosion_parameters.json') as f:
    prm = json.load(f)

t, c_all, z_all, w_all, uh, vh = analytic.read_centroid_series('sediment_erosion.sww')
_, c_ctl, _, _, _, _ = analytic.read_centroid_series('sediment_erosion.sww', name='control')
h_cell = numpy.array(prm['h_cell'])
x_cell = numpy.array(prm['x_cell'])
law = (prm['v_s'], prm['R'], prm['diameter'], prm['tau_c_star'],
       prm['gamma0'], prm['d_star'])

tt = numpy.linspace(0.0, t.max(), 400)
c_ref, z_ref, h_ref = analytic.erosion_with_feedback(
    tt, h_cell, prm['S_bed'], *law, porosity=prm['porosity'], z0=z_all[0])
c_fixed = analytic.erosion(tt, h_cell, prm['S_bed'], *law)

pyplot.clf()
pyplot.plot(t, c_all.mean(axis=1), 'b.', label='numerical, sand')
pyplot.plot(tt, c_ref.mean(axis=1), 'r-', label='reference (bed evolving)')
pyplot.plot(tt, c_fixed.mean(axis=1), 'k--', lw=0.8, label='fixed-bed exponential')
pyplot.plot(t, c_ctl.mean(axis=1), 'g.', label='numerical, control fraction (below threshold)')
pyplot.title('Depth-averaged concentration, tank mean')
pyplot.xlabel('Time (s)')
pyplot.ylabel('Concentration (volumetric)')
pyplot.legend(loc='best')
pyplot.savefig('concentration_plot.png')

pyplot.clf()
pyplot.plot(t, (z_all - z_all[0]).mean(axis=1), 'b.', label='numerical')
pyplot.plot(tt, (z_ref - z_ref[0]).mean(axis=1), 'r-', label='reference')
pyplot.title('Bed lowering from entrainment, tank mean')
pyplot.xlabel('Time (s)')
pyplot.ylabel('Bed change (m)')
pyplot.legend(loc='best')
pyplot.savefig('bed_lowering_plot.png')

pyplot.clf()
order = numpy.argsort(x_cell)
c_ref_end = analytic.erosion_with_feedback(t[-1:], h_cell, prm['S_bed'], *law,
                                           porosity=prm['porosity'], z0=z_all[0])[0][0]
pyplot.plot(x_cell[order], c_all[-1, order], 'b.', label='numerical, t = %g s' % t[-1])
pyplot.plot(x_cell[order], c_ref_end[order], 'r-', lw=0.8, label='reference')
pyplot.title('Concentration in every cell along the tank at the end of the run')
pyplot.xlabel('x (m)')
pyplot.ylabel('Concentration (volumetric)')
pyplot.legend(loc='best')
pyplot.savefig('profile_plot.png')
