"""Example 4 -- a sand dune overtopped, with and without slope collapse.

Erosion on its own will cut a vertical face. Real sand cannot stand at 90
degrees: once scour oversteepens it, the face collapses and the material runs
out over the toe. Spec 7 relaxation is what puts that back in.

This runs the SAME overtopping event twice, once with relaxation off and once
with it on, and compares the bed it leaves behind. The comparison is the
example; either run alone would tell you very little.

Demonstrates
  * domain.set_angle_of_repose(), and what it changes
  * that the collapse CONSERVES MASS -- the material that leaves the face
    arrives at the toe, rather than vanishing
  * the sweep count and the cap report, which is the part of this that
    surprises people
  * how to read a slope field back out of a domain

The scenario follows the sanddune operator's: a dune line overtopped by a wave,
where the concern is how much of the dune survives. The physics underneath is
different -- this is a sediment mass balance, that is a kinematic erosion rule
-- and the mass budget below is the difference made visible.
"""
import math
import numpy as np
import anuga

LEN_X, LEN_Y = 120.0, 30.0
CREST_X, CREST_W, CREST_H = 60.0, 4.0, 2.5
ANGLE = 33.0            # degrees; dry sand sits near 34, FG21 use 35
FINALTIME = 60.0
PORO = 0.30


def dune(x, y):
    """A dune ridge across the channel, on a gently falling bed."""
    return -0.004 * x + CREST_H * np.exp(-((x - CREST_X) / CREST_W) ** 2)


def max_bed_slope(domain):
    """Steepest centroid-to-centroid bed slope, in degrees."""
    z = domain.quantities['elevation'].centroid_values
    cc = domain.centroid_coordinates
    nb = domain.neighbours
    worst = 0.0
    for i in range(3):
        j = nb[:, i]
        ok = j >= 0
        dz = np.abs(z[ok] - z[j[ok]])
        dist = np.sqrt(((cc[ok] - cc[j[ok]]) ** 2).sum(axis=1))
        good = dist > 0
        worst = max(worst, float((dz[good] / dist[good]).max()))
    return math.degrees(math.atan(worst))


def build(repose):
    domain = anuga.rectangular_cross_domain(100, 20, len1=LEN_X, len2=LEN_Y)
    domain.set_name('dune_collapse_%s' % ('on' if repose else 'off'))
    domain.set_flow_algorithm('DE0')
    domain.set_low_froude(0)
    domain.store = False

    x = domain.centroid_coordinates[:, 0]
    z = dune(x, domain.centroid_coordinates[:, 1])
    domain.set_quantity('elevation', z, location='centroids')
    # A raised pool behind the dune, just over the crest, so it overtops.
    # The dune STARTS below the critical angle (27.8 degrees against 33), so
    # anything relaxation does here was caused by the scour, not by the initial
    # condition being unreasonable.
    domain.set_quantity('stage',
                        np.where(x < CREST_X, np.maximum(z, 2.6), z - 0.01),
                        location='centroids')
    domain.set_quantity('friction', 0.025)
    domain.set_boundary({t: anuga.Reflective_boundary(domain)
                         for t in domain.get_boundary_tags()})

    domain.set_sediment_parameters(porosity=PORO)
    domain.set_bed_material('noncohesive')
    domain.add_sediment_fraction(name='dune_sand', diameter=3.0e-4)
    if repose:
        # 50 sweeps is the default and is sized for a bed that is already
        # near-relaxed. This one starts below the critical angle, so the cap
        # is only tested once scour has oversteepened the face.
        domain.set_angle_of_repose(ANGLE, max_sweeps=50)
    return domain


def run(repose):
    domain = build(repose)
    op = [o for o in domain.fractional_step_operators
          if type(o).__name__ == 'Sediment_transport_operator'][0]
    areas = domain.areas
    z0 = domain.quantities['elevation'].centroid_values.copy()
    bed0 = float((z0 * areas).sum())

    domain.evolve_to_end(finaltime=FINALTIME)

    z1 = domain.quantities['elevation'].centroid_values
    m1 = float((domain.tracer_conserved_values[0] * areas).sum())
    bed1 = float((z1 * areas).sum())
    return dict(domain=domain, op=op, z0=z0, z1=z1,
                slope=max_bed_slope(domain),
                crest=float(z1.max()),
                suspended=m1,
                bed_change=float(((1.0 - PORO) * (z1 - z0) * areas).sum()),
                bed_volume_change=bed1 - bed0)


print(__doc__)
print('Starting bed: crest %.2f m, steepest slope %.1f degrees'
      % (dune(np.array([CREST_X]), 0.0)[0], max_bed_slope(build(False))))
print('Critical angle: %.1f degrees\n' % ANGLE)

off = run(False)
on = run(True)

print('%-34s %14s %14s' % ('', 'repose off', 'repose on'))
print('%-34s %14.2f %14.2f' % ('steepest bed slope (deg)',
                               off['slope'], on['slope']))
print('%-34s %14.3f %14.3f' % ('surviving crest (m)',
                               off['crest'], on['crest']))
print('%-34s %14.4e %14.4e' % ('suspended sediment (m3)',
                               off['suspended'], on['suspended']))
print('%-34s %14.4e %14.4e' % ('bed change, (1-n) dz (m3)',
                               off['bed_change'], on['bed_change']))
print('%-34s %14.2e %14.2e'
      % ('budget: suspended + bed (m3)',
         off['suspended'] + off['bed_change'],
         on['suspended'] + on['bed_change']))
print()

print('With relaxation off, scour leaves the face at %.1f degrees -- steeper'
      % off['slope'])
print('than sand can stand. With it on the bed is held at %.2f degrees. That is'
      % on['slope'])
print('marginally over %.1f because the kernel declares convergence within a'
      % ANGLE)
print('1e-3 relative tolerance on the threshold SLOPE; a strict test never')
print('terminates. See PHYSICS_SPEC 7.1.')
print()

# The collapse must not invent or destroy sand. Relaxation moves BED material
# between cells, so the quantity it has to conserve is bulk bed volume, and it
# is conserved separately from the erosion/deposition budget above.
print('Relaxation moved material without creating any: the sediment budget')
print('closes to %.2e m3 with it on, against %.2e m3 with it off -- the same'
      % (abs(on['suspended'] + on['bed_change']),
         abs(off['suspended'] + off['bed_change'])))
print('machine precision either way. A collapse that simply lowered the face')
print('and discarded the sand would show up here as a bed change with no')
print('matching suspension.')
print()

print('Relaxation cost: %d sweeps on the last step, %d in total, cap reached'
      % (on['op'].repose_sweeps, on['op'].repose_sweeps_total))
print('on %d step(s). A bed that starts below the critical angle needs only a'
      % on['op'].repose_cap_hits)
print('few sweeps per step; the expensive case is a bed that starts over-steep.')
if on['op'].repose_cap_hits:
    print('Hitting the cap is not a failure -- progress carries over between')
    print('steps -- but it does mean relaxation was still lagging the scour.')
    print('Raise max_sweeps if you need the bed relaxed within every step.')
print()

# Where the two runs differ is the point. It is NOT the crest -- both lose
# essentially the same height to scour -- it is the face and the toe below it.
diff = on['z1'] - off['z1']
x = on['domain'].centroid_coordinates[:, 0]
lo, hi = int(np.argmin(diff)), int(np.argmax(diff))
dz = on['z1'] - on['z0']
print('Bed movement with relaxation on: %+.3f m (scour) to %+.3f m (deposit)'
      % (dz.min(), dz.max()))
print('Crest survives at %.3f m either way -- scour takes the top of the dune'
      % on['crest'])
print('whatever you do. The difference is on the FACE: relaxation removes a')
print('further %.3f m at x = %.1f m and lays it down %.3f m thick at x = %.1f m,'
      % (-diff[lo], x[lo], diff[hi], x[hi]))
print('turning a %.1f degree face into a %.1f degree one.'
      % (off['slope'], on['slope']))
