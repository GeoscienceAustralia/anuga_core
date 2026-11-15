"""The checks on EA Test 1, shared by validate_test1.py and plot_results.py.

The expected outcome is stated in the report (Neelz and Pender 2013,
Section 4.2): a peak level of ~10.35 m and a final level of ~10.25 m (the
sill) at both points, which every package reproduced; differences of up
to +0.01 m in the final level were put down to the wet/dry threshold.
TUFLOW's five current solvers end at 10.251-10.258 m (BMT 2017). The case
is checked against

* the peak and final levels at points 1 and 2;
* the water surface in the depression being flat: the report found the
  two points within a few mm of each other after t = 2 h;
* the sheet flow over point 1 (on the far slope of the sill, above the
  depression's floor) starting at about 1 h, before the depression fills
  to point 2;
* a stored baseline per flow algorithm (baseline_test1_<alg>.csv). Set
  EA_REGENERATE_BASELINE=1 to write or rewrite it.
"""

import numpy as np

import ea_common as ea

NAMES = ['1', '2']
WET = 0.001                # depth marking the water's arrival [m]


def evaluate(alg):
    g = ea.load_gauges('test1_gauges.csv')
    t = g['time'] / 3600.0
    c = ea.Checks()
    for n in NAMES:
        peak = g['stage_' + n].max()
        c.add('Peak level at point %s' % n, '%.3f m' % peak, '10.34-10.36 m (report ~10.35)',
              10.34 <= peak <= 10.36)
    for n in NAMES:
        final = g['stage_' + n][-1]
        c.add('Final level at point %s' % n, '%.3f m' % final,
              '10.25-10.26 m (TUFLOW 10.251-10.258)', 10.25 <= final <= 10.26)
    diff = np.abs(g['stage_1'] - g['stage_2'])[t > 2.0].max()
    c.add('Flat pond: |level 1 - level 2|, t > 2 h', '%.4f m' % diff, '< 0.01 m (report: a few mm)',
          diff < 0.01)
    arrive = {n: t[np.argmax(g['depth_' + n] > WET)] for n in NAMES}
    c.add('Water reaches point 1 before point 2', '%.2f h, %.2f h' % (arrive['1'], arrive['2']),
          'point 1 at 0.8-1.3 h (report ~1 h)', 0.8 <= arrive['1'] <= 1.3 and arrive['1'] < arrive['2'])
    ea.baseline_check(c, 'test1_gauges.csv', 'test1', alg, ['stage_' + n for n in NAMES],
                      {'stage': 0.005})
    return c
