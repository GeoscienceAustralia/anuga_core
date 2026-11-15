"""The checks on EA Test 6B, shared by validate_test6b.py and plot_results.py.

No exact solution exists, and the report publishes the other packages'
results only as figures (Neelz and Pender 2013, Figures 4.27-4.28). The
case is checked against

* volume: the channel is closed, so the initial volume must be kept;
* the geometry: the bore reaches G1/G2 (x = 53 m) before G3/G4 (80 m)
  and G5 (104 m), and the reservoir gauge G6 drains;
* the report's 2013 concern about ANUGA's high peak velocities (Tests 5
  and 6B): the largest speed anywhere in water deeper than 5 cm must stay
  below the dam-break bound 2 sqrt(g h0) = 17.7 m/s;
* the end state: by 30 min the levels at the gauges agree;
* a stored baseline per flow algorithm (baseline_test6b_<alg>.csv). Set
  EA_REGENERATE_BASELINE=1 to write or rewrite it.
"""

import numpy as np

import ea_common as ea

NAMES = ['G1', 'G2', 'G3', 'G4', 'G5', 'G6']
H_DOWN = 0.4
SPEED_BOUND = 2.0 * np.sqrt(9.8 * 8.0)


def arrivals(g):
    """Bore arrival (level 5 cm above the initial 0.4 m) at G1-G5 [s]."""
    t = g['time']
    out = {}
    for n in NAMES[:5]:
        risen = g['stage_' + n] > H_DOWN + 0.05
        out[n] = t[np.argmax(risen)] if risen.any() else np.inf
    return out


def evaluate(alg):
    g = ea.load_gauges('test6b_gauges.csv')
    vol = ea.read_volume('test6b_volume.txt')
    vol_err = abs(vol['final'] - vol['initial']) / vol['initial']
    a = arrivals(g)
    final = np.array([g['stage_' + n][-1] for n in NAMES])

    c = ea.Checks()
    c.add('Volume kept (closed channel)', '%.1e relative' % vol_err, '< 1e-3', vol_err < 1.0e-3)
    c.add('Bore order G1/G2, G3/G4, G5', '%.0f/%.0f, %.0f/%.0f, %.0f s'
          % (a['G1'], a['G2'], a['G3'], a['G4'], a['G5']), 'increasing',
          max(a['G1'], a['G2']) < min(a['G3'], a['G4']) < a['G5'])
    c.add('Reservoir (G6) drains', '%.2f m at 30 min' % g['stage_G6'][-1], '< 4 m',
          g['stage_G6'][-1] < 0.5 * g['stage_G6'][0])
    c.add('Largest speed (depth > 5 cm)', '%.2f m/s' % vol['max_speed'],
          '< 2 sqrt(g h0) = %.1f m/s' % SPEED_BOUND, vol['max_speed'] < SPEED_BOUND)
    c.add('Gauge levels agree at 30 min', 'spread %.3f m' % (final.max() - final.min()),
          '< 0.05 m', final.max() - final.min() < 0.05)
    ea.baseline_check(c, 'test6b_gauges.csv', 'test6b', alg, ['stage_' + n for n in NAMES],
                      {'stage': 0.02})
    return c
