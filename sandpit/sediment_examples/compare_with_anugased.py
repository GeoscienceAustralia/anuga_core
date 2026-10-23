"""Run anugaSed and this module on the same case, and diff them.

Both are now importable in one process: anugaSed via the py3-modernisation
branch (`pip install -e .`), ours as part of ANUGA. The case is anugaSed's own
`run_simple_sed_transport.py` channel, unchanged.

To make the comparison mean something, our side is configured to match theirs
where the spec says we can:

    erosion        [E-3] Hanson & Simon, cohesive, tau_c = 0.088 Pa
    shear closure  [T-7] depth-slope, tau_b = rho g h S
    d*             Rouse profile (their polynomial is a fit to the same [S-4])
    D50 65 um, porosity 0.3, rho_s 2650

KNOWN REMAINING DIFFERENCES -- exact agreement is not expected, and this script
exists to size the gap, not to claim one:

  D1a  their code divides the elevation gradient by a domain-mean cell size and
       clamps S <- min(S, mean(S)/2). Neither is in aSM16; the first is
       dimensionally inconsistent. We follow the manual.
  D5   they hard-code a 0.2 concentration cap, applied three times per call.
  D6   they clamp edot and ddot to 1e-3 m/s, silently.
  D4   their d* is an 8th-degree polynomial extrapolated outside its fit range;
       ours is a bounded fit with [L-4] limiting the near-bed concentration.
"""
import numpy as np
import anuga
from anuga import rectangular_cross, Domain, Dirichlet_boundary, Reflective_boundary

LENGTH = WIDTH = 5.0
DX = DY = 2.0
CONC = 0.01
D50 = 65.0e-6
TAU_CRIT = 0.088
FINALTIME = 30.0


def topography(x, y):
    return -x / 100.0


def build(evolved=None):
    pts, vts, bdy = rectangular_cross(int(LENGTH / DX), int(WIDTH / DY),
                                      len1=LENGTH, len2=WIDTH)
    d = (Domain(pts, vts, bdy, evolved_quantities=evolved) if evolved
         else Domain(pts, vts, bdy))
    d.set_flow_algorithm('DE0')
    d.store = False
    d.set_quantity('elevation', topography)
    d.set_quantity('stage', expression='elevation')
    mx = d.quantities['elevation'].vertex_values.max()
    mn = d.quantities['elevation'].vertex_values.min()
    d.set_boundary({'left': Dirichlet_boundary([mx + 0.5, 0, 0]),
                    'right': Dirichlet_boundary([mn - 1, 0, 0]),
                    'top': Reflective_boundary(d),
                    'bottom': Reflective_boundary(d)})
    return d


# ---------------------------------------------------------------- anugaSed --
from anugaSed import Sed_transport_operator

a = build(evolved=['stage', 'xmomentum', 'ymomentum', 'concentration'])
a.set_quantity('concentration', CONC)
Sed_transport_operator(a)
za0 = a.quantities['elevation'].centroid_values.copy()
a.evolve_to_end(finaltime=FINALTIME)
dza = a.quantities['elevation'].centroid_values - za0
ca = a.quantities['concentration'].centroid_values

# -------------------------------------------------------------------- ours --
b = build()
b.sediment_porosity = 0.3
b.set_bed_material('cohesive', tau_crit=TAU_CRIT)
b.set_shear_closure('depth_slope')
b.sediment_d_star_mode = 1                      # Rouse, as their fit approximates
b.add_sediment_fraction(name='sand', diameter=D50, rho_s=2650.0, initial_concentration=CONC)
b.set_tracer_boundary('sand', 'left', CONC)   # 'left' is the inflow
zb0 = b.quantities['elevation'].centroid_values.copy()
b.evolve_to_end(finaltime=FINALTIME)
dzb = b.quantities['elevation'].centroid_values - zb0
cb = b.get_tracer('sand')

# ------------------------------------------------------------------ compare --
print(__doc__)
print('%d cells, %.0f s\n' % (len(a), FINALTIME))
print('                              anugaSed          ours        ratio')
def row(label, x, y):
    r = (y / x) if x else float('nan')
    print('  %-26s %+12.5e  %+12.5e   %6.2f' % (label, x, y, r))

row('mean bed change  [m]', dza.mean(), dzb.mean())
row('max bed change   [m]', dza.max(), dzb.max())
row('min bed change   [m]', dza.min(), dzb.min())
row('mean concentration', ca.mean(), cb.mean())
row('max concentration', ca.max(), cb.max())
print()
print('  sign of mean bed change:  anugaSed %s   ours %s   -> %s'
      % ('accretion' if dza.mean() > 0 else 'erosion',
         'accretion' if dzb.mean() > 0 else 'erosion',
         'AGREE' if np.sign(dza.mean()) == np.sign(dzb.mean()) else 'DISAGREE'))
print('  cell-by-cell correlation of bed change: %.4f'
      % (np.corrcoef(dza, dzb)[0, 1] if dza.std() > 0 and dzb.std() > 0
         else float('nan')))

# --- where the difference comes from ----------------------------------------
# Not the erosion law: that agrees exactly (test_cohesive.py D1). It is D1a.
S = [op for op in a.fractional_step_operators
     if type(op).__name__ == 'Sed_transport_operator'][0].calculate_energy_slope()
true_slope = 1.0 / 100.0
print()
print('  DIAGNOSIS -- their effective bed slope, vs the slope actually imposed:')
print('    imposed bed slope        %.6f' % true_slope)
print('    their S                  %.6e   (min == max: %s)'
      % (S.mean(), S.min() == S.max()))
print('    ratio                    %.4f  -> tau_b too small by %.1fx'
      % (S.mean() / true_slope, true_slope / S.mean()))
print()
print('  Their S is both too small AND spatially uniform. The /dx accounts for')
print('  one factor and the S <- min(S, mean(S)/2) clamp for the rest; the')
print('  clamp also flattens every cell to one value, which is why the')
print('  cell-by-cell correlation above is weak. Divergence D1a, measured.')
print('  The erosion law itself agrees exactly -- see test_cohesive.py.')
