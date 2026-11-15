"""The checks on EA Test 4, shared by validate_test4.py and plot_results.py.

No exact solution exists, and the report publishes the other packages'
results only as figures (Neelz and Pender 2013, Figures 4.10-4.11). The
case is checked against

* the inflow volume (closed domain, so all of it must still be there);
* the report's summary (Section 4.5.4): the full SWE packages agree within
  a few per cent, with a travel time of ~1 h to points 5 and 6;
* the geometry: the front reaches points 1-5 in order of distance from
  the breach, and point 6 (424 m away) shortly after point 5 (400 m);
* a stored baseline per flow algorithm (baseline_test4_<alg>.csv). Set
  EA_REGENERATE_BASELINE=1 to write or rewrite it.
"""

import numpy as np

import ea_common as ea

NAMES = ['1', '2', '3', '4', '5', '6']
WET = 0.01                 # depth marking the front's arrival [m]


def arrivals(g):
    """Front arrival (depth > WET) at each point [min]."""
    t = g['time']
    out = {}
    for n in NAMES:
        wet = g['depth_' + n] > WET
        out[n] = t[np.argmax(wet)] / 60.0 if wet.any() else np.inf
    return out


def evaluate(alg):
    g = ea.load_gauges('test4_gauges.csv')
    vol = ea.read_volume('test4_volume.txt')
    vol_err = abs(vol['final'] - vol['inflow']) / vol['inflow']
    arrival = arrivals(g)
    order = [arrival[n] for n in NAMES[:5]]

    c = ea.Checks()
    c.add('Volume kept (closed domain)', '%.3f%% error' % (100 * vol_err), '< 1%',
          vol_err < 0.01)
    c.add('Front reaches points 1-5 in order', ', '.join('%.0f' % a for a in order) + ' min',
          'increasing', all(a < b for a, b in zip(order, order[1:])))
    c.add('Front at point 5', '%.1f min' % arrival['5'], '45-80 min (report ~1 h)',
          45.0 < arrival['5'] < 80.0)
    c.add('Front at point 6', '%.1f min' % arrival['6'], 'after point 5, within x1.25',
          arrival['5'] < arrival['6'] < 1.25 * arrival['5'])
    ea.baseline_check(c, 'test4_gauges.csv', 'test4', alg, ['stage_' + n for n in NAMES],
                      {'stage': 0.01})
    return c
