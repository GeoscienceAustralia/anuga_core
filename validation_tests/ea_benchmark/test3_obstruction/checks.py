"""The checks on EA Test 3, shared by validate_test3.py and plot_results.py.

No exact solution exists, and the report publishes the other packages'
results only as figures (Neelz and Pender 2013, Figure 4.8). The case is
checked against

* the inflow volume (closed domain, so all of it must still be there);
* the report's summary of the full SWE packages (Section 4.4.3): point 1
  rises from 9.75 m to ~9.98-9.99 m by t ~ 120 s, and after t ~ 150 s the
  level at point 2 rises by ~5-6 cm as water crosses the crest;
* a stored baseline per flow algorithm (baseline_test3_<alg>.csv), so that
  a change that moves the answer is caught. Set EA_REGENERATE_BASELINE=1
  to write or rewrite the baseline for the algorithm being run.
"""

import numpy as np

import ea_common as ea

BED = 9.75                 # bed level at both points


def evaluate(alg):
    g = ea.load_gauges('test3_gauges.csv')
    t = g['time']
    vol = ea.read_volume('test3_volume.txt')

    vol_err = abs(vol['final'] - vol['inflow']) / vol['inflow']
    level_1_120 = np.interp(120.0, t, g['stage_1'])
    level_1_end = g['stage_1'][-1]
    rise_2 = g['stage_2'][-1] - BED
    wet = g['depth_2'] > 0.01
    wet_2 = t[np.argmax(wet)] if wet.any() else np.inf

    c = ea.Checks()
    c.add('Volume kept (closed domain)', '%.2f%% error' % (100 * vol_err), '< 1%',
          vol_err < 0.01)
    c.add('Point 1 level at 120 s', '%.3f m' % level_1_120, '9.94-10.01 m (report ~9.98-9.99)',
          9.94 < level_1_120 < 10.01)
    c.add('Point 1 level at 900 s', '%.3f m' % level_1_end, '9.95-10.0 m',
          9.95 < level_1_end < 10.0)
    c.add('Point 2 wets (depth > 1 cm)', '%.0f s' % wet_2, '100-250 s (report ~150 s)',
          100.0 < wet_2 < 250.0)
    c.add('Point 2 final rise', '%.3f m' % rise_2, '0.03-0.09 m (report ~0.05-0.06)',
          0.03 < rise_2 < 0.09)
    ea.baseline_check(c, 'test3_gauges.csv', 'test3', alg, ['stage_1', 'stage_2'],
                      {'stage': 0.01})
    return c
