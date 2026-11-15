"""The checks on EA Test 2, shared by validate_test2.py and plot_results.py.

No exact solution exists, and the report publishes the other packages'
results only as figures (Neelz and Pender 2013, Figures 4.5-4.6). The
case is checked against

* the inflow volume (closed domain, so all of it must still be there);
* the final levels. Depressions 1-4 and 6-8 overflow and drain back to
  the lowest sill out of each (Section 4.3.3), whose level is a property
  of the DEM: the final level must lie within 2 cm above it (still
  draining). The full SWE packages differed by a few cm at most, so the
  final levels at points 1-4, 7, 8, 11 and 12 must also lie within 0.03 m
  of the range of TUFLOW's five current solvers (read from the charts in
  BMT 2017; the chart for point 6 did not read reliably). Point 5 is still
  filling at 48 h and point 10 holds little water, so neither is checked;
* the transient peak above the final level near the inflow: ~0.2 m at
  point 4, ~0.1 m at points 3 and 8, ~0.05 m at points 2 and 7;
* the flow reaching points 10-12 by its momentum (the 2- and 0-term
  packages left them dry) and points 13-16 staying dry;
* a stored baseline per flow algorithm (baseline_test2_<alg>.csv). Set
  EA_REGENERATE_BASELINE=1 to write or rewrite it.
"""

import ea_common as ea

NAMES = [str(n) for n in range(1, 17)]
# The lowest sill out of each depression that overflows [m]: the highest
# ground on the lowest path from its floor to lower ground, found by a
# priority flood over the 2 m DEM from the depression's lowest cell.
SILL = {'1': 8.579, '2': 8.868, '3': 9.201, '4': 9.535, '6': 8.618, '7': 8.951,
        '8': 9.285}
SILL_TOL = 0.02
# Final levels of TUFLOW's five current solvers [m] (BMT 2017, Table 7)
TUFLOW_FINAL = {'1': (8.588, 8.596), '2': (8.872, 8.881), '3': (9.205, 9.214),
                '4': (9.539, 9.545), '7': (8.955, 8.962),
                '8': (9.289, 9.295), '11': (8.530, 8.546), '12': (8.886, 8.894)}
FINAL_TOL = 0.03
# Transient peak above the final level (report, Section 4.3.3) [m]
PEAK_RISE = {'4': (0.15, 0.25), '3': (0.05, 0.15), '8': (0.05, 0.15),
             '2': (0.02, 0.08), '7': (0.02, 0.08)}
WET = 0.01                 # depth counted as wet [m]


def evaluate(alg):
    g = ea.load_gauges('test2_gauges.csv')
    vol = ea.read_volume('test2_volume.txt')
    vol_err = abs(vol['final'] - vol['inflow']) / vol['inflow']

    c = ea.Checks()
    c.add('Volume kept (closed domain)', '%.3f%% error' % (100 * vol_err), '< 0.1%',
          vol_err < 0.001)
    above = {n: g['stage_' + n][-1] - level for n, level in SILL.items()}
    c.add('Final levels against sill (1-4, 6-8)',
          '%.3f-%.3f m above' % (min(above.values()), max(above.values())),
          '0-%g m above (drains to the sill)' % SILL_TOL,
          all(-0.002 <= a <= SILL_TOL for a in above.values()))
    off = {}
    for n, (lo, hi) in TUFLOW_FINAL.items():
        s = g['stage_' + n][-1]
        off[n] = max(lo - s, s - hi, 0.0)
    worst = max(off, key=off.get)
    c.add('Final levels against TUFLOW (%s)' % ', '.join(TUFLOW_FINAL),
          '%.3f m outside (point %s)' % (off[worst], worst),
          'within %g m of TUFLOW' % FINAL_TOL, off[worst] <= FINAL_TOL)
    for n, (lo, hi) in PEAK_RISE.items():
        rise = g['stage_' + n].max() - g['stage_' + n][-1]
        c.add('Transient peak at point %s' % n, '%.2f m above final' % rise,
              '%g-%g m (report ~%g)' % (lo, hi, round((lo + hi) / 2, 2)), lo <= rise <= hi)
    wet_east = {n: g['depth_' + n].max() for n in ('10', '11', '12')}
    c.add('Flow reaches points 10-12', ', '.join('%.2f' % v for v in wet_east.values()) + ' m deep',
          '> %g m (carried east)' % WET, min(wet_east.values()) > WET)
    dry = max(g['depth_' + n].max() for n in ('13', '14', '15', '16'))
    c.add('Points 13-16 stay dry', 'max depth %.3f m' % dry, '< %g m (report)' % WET, dry < WET)
    ea.baseline_check(c, 'test2_gauges.csv', 'test2', alg, ['stage_' + n for n in NAMES],
                      {'stage': 0.01})
    return c
