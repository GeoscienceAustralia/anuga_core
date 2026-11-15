"""The checks on EA Test 5, shared by validate_test5.py and plot_results.py.

No exact solution exists, and the report publishes the other packages'
results only as figures (Neelz and Pender 2013, Figures 4.16-4.23). The
case is checked against

* the inflow volume (closed domain, so all of it must still be there);
* the geometry: the flood reaches the points in order of distance down
  the valley (1, 6, 2, 3, 7, 4, 5);
* the report's summary of the full SWE packages (Section 4.6.4): a travel
  time of almost 3 h to point 5, the packages within ~30 min of each other
  (in 2009 ANUGA was ~1.5 h late), and peak speeds of 1.5-2.2 m/s at
  point 4 and 1.2-1.6 m/s at point 7;
* the final level in the pond at point 5, where the packages agreed within
  0.4 m (ANUGA was ~0.5 m low in 2009): TUFLOW's five solvers give
  148.6-148.9 m (BMT 2017);
* a stored baseline per flow algorithm (baseline_test5_<alg>.csv). Set
  EA_REGENERATE_BASELINE=1 to write or rewrite it.
"""

import numpy as np

import ea_common as ea

NAMES = ['1', '2', '3', '4', '5', '6', '7']
DOWN_VALLEY = ['1', '6', '2', '3', '7', '4', '5']
WET = 0.05                 # depth marking the flood's arrival [m]


def arrivals(g):
    """Flood arrival (depth > WET) at each point [h]."""
    t = g['time']
    out = {}
    for n in NAMES:
        wet = g['depth_' + n] > WET
        out[n] = t[np.argmax(wet)] / 3600.0 if wet.any() else np.inf
    return out


def evaluate(alg):
    g = ea.load_gauges('test5_gauges.csv')
    vol = ea.read_volume('test5_volume.txt')
    vol_err = abs(vol['final'] - vol['inflow']) / vol['inflow']
    a = arrivals(g)
    order = [a[n] for n in DOWN_VALLEY]
    final_5 = g['stage_5'][-1]

    c = ea.Checks()
    c.add('Volume kept (closed domain)', '%.3f%% error' % (100 * vol_err), '< 1%',
          vol_err < 0.01)
    c.add('Flood reaches 1, 6, 2, 3, 7, 4, 5 in order', ', '.join('%.1f' % x for x in order) + ' h',
          'increasing', all(x < y for x, y in zip(order, order[1:])))
    c.add('Flood reaches point 5', '%.2f h' % a['5'], '2.5-3.5 h (report: almost 3 h)',
          2.5 < a['5'] < 3.5)
    c.add('Peak speed at point 4', '%.2f m/s' % g['speed_4'].max(), '1.5-2.2 m/s (report)',
          1.5 <= g['speed_4'].max() <= 2.2)
    c.add('Peak speed at point 7', '%.2f m/s' % g['speed_7'].max(), '1.2-1.6 m/s (report)',
          1.2 <= g['speed_7'].max() <= 1.6)
    c.add('Final level at point 5', '%.2f m' % final_5, '148.4-149.1 m (TUFLOW 148.6-148.9)',
          148.4 < final_5 < 149.1)
    ea.baseline_check(c, 'test5_gauges.csv', 'test5', alg, ['stage_' + n for n in NAMES],
                      {'stage': 0.02})
    return c
