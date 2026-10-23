"""Example 1 -- sediment-laden flow down a short channel.

Ported from anugaSed's `run_simple_sed_transport.py` (Perignon 2016, MIT).
Same geometry, slope, boundaries and duration; our API and, deliberately, the
non-cohesive erosion route -- see README.md.

Demonstrates
  * add_sediment_fraction()          registers the grain size and its operator
  * set_tracer_boundary()     prescribed inflow concentration, per boundary tag
  * entrainment [E-1], deposition [D-1], Exner bed change [G-4]
"""
import argparse

import numpy as np
import anuga

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--bed', default='noncohesive',
                choices=['noncohesive', 'cohesive'],
                help="bed material (spec 4.1.1). 'cohesive' reproduces "
                     "anugaSed's own erosion route [E-3] with their "
                     "tau_crit = 0.088 Pa; 'noncohesive' uses the Shields "
                     "route [E-1] this work targets.")
args = ap.parse_args()

# --- their setup ------------------------------------------------------------
length = 5.0
width = 5.0
dx = dy = 2.0
INFLOW_CONC = 0.01          # their domain.set_quantity('concentration', 0.01)
D50 = 65.0e-6               # their grain_size, Griffin et al. 2014
FINALTIME = 30.0


def topography(x, y):
    return -x / 100.0


points, vertices, boundary = anuga.rectangular_cross(
    int(length / dx), int(width / dy), len1=length, len2=width)
domain = anuga.Domain(points, vertices, boundary)
domain.set_name('sed_channel')
domain.set_flow_algorithm('DE0')
domain.set_quantity('elevation', topography)
domain.set_quantity('stage', expression='elevation')      # dry start
domain.set_quantity('friction', 0.03)

max_elev = domain.quantities['elevation'].vertex_values.max()
min_elev = domain.quantities['elevation'].vertex_values.min()
domain.set_boundary({
    'left':  anuga.Dirichlet_boundary([max_elev + 0.5, 0.0, 0.0]),   # inflow
    'right': anuga.Dirichlet_boundary([min_elev - 1.0, 0.0, 0.0]),   # outflow
    'top':   anuga.Reflective_boundary(domain),
    'bottom': anuga.Reflective_boundary(domain),
})

# --- our API ----------------------------------------------------------------
domain.sediment_porosity = 0.3
if args.bed == 'cohesive':
    # anugaSed's own configuration: Hanson & Simon [E-3] with tau_c = 0.088 Pa.
    domain.set_bed_material('cohesive', tau_crit=0.088)
domain.add_sediment_fraction(name='sand', diameter=D50, rho_s=2650.0,
                      initial_concentration=INFLOW_CONC)
# The domain starts DRY, so the initial concentration washes out immediately;
# what actually matters is the concentration of the water flowing in.
domain.set_tracer_boundary('sand', 'left', INFLOW_CONC)   # 'left' is the inflow

areas = domain.areas
z0 = domain.quantities['elevation'].centroid_values.copy()
lam = domain.sediment_porosity

print('bed material: %s   %s'
      % (args.bed,
         ('[E-3] Hanson & Simon, tau_c = %.3f Pa, K_e = %.4e m3/N/s'
          % (domain.sediment_tau_crit, domain.sediment_K_e))
         if args.bed == 'cohesive' else '[E-1] Shields, tau_c* = 0.04'))
print('v_s = %.3e m/s   D50 = %.0f um   inflow c = %.3f'
      % (domain.sediment_settling_velocity[0], D50 * 1e6, INFLOW_CONC))
print()
print('   time    water vol    suspended     bed change     max |dz|')
for t in domain.evolve(yieldstep=5.0, finaltime=FINALTIME):
    z = domain.quantities['elevation'].centroid_values
    h = np.maximum(domain.quantities['stage'].centroid_values - z, 0.0)
    m = float((domain.tracer_conserved_values[0] * areas).sum())
    bed = float(((1.0 - lam) * (z - z0) * areas).sum())
    print('  %5.1f s   %8.4f m3  %.6e   %+.4e   %.3e'
          % (t, float((h * areas).sum()), m, bed, np.abs(z - z0).max()))

z = domain.quantities['elevation'].centroid_values
dz = z - z0
print()
eroded = dz[dz < 0]
deposited = dz[dz > 0]
print('Bed change over %.0f s:  min %+.3e m   max %+.3e m' % (FINALTIME, dz.min(), dz.max()))
print('  cells eroding: %d   cells accreting: %d   cells unchanged: %d'
      % (eroded.size, deposited.size, dz.size - eroded.size - deposited.size))
print()
print('The whole reach erodes here: the imposed 0.5 m head over a 5 m channel')
print('drives a fast, under-loaded flow, so entrainment [E-1] dominates')
print('deposition everywhere. Bed volume is NOT conserved in this example --')
print('sediment enters with the inflow and leaves through the outflow.')
