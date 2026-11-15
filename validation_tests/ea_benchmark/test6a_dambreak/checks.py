"""The checks on EA Test 6A, shared by validate_test6a.py and plot_results.py.

Test 6A is the laboratory flume, for which the report shows the UCL
measurements against the packages (Neelz and Pender 2013, Figures
4.25-4.26) but publishes them only as figures. The case is checked against

* volume: the flume is closed, so the initial volume must be kept;
* the geometry: the bore reaches G1/G2 (x = 2.65 m) before G3/G4 (4 m)
  and G5 (5.2 m);
* the report's summary (Section 4.7.1): every package matched the
  measured level at G6, upstream of the gate, which reflects the periodic
  flow through the gate with a period of ~9 s;
* a bound on the largest speed in water deeper than 2.5 mm: the dam-break
  bound 2 sqrt(g h0) = 4.0 m/s;
* a stored baseline per flow algorithm (baseline_test6a_<alg>.csv). Set
  EA_REGENERATE_BASELINE=1 to write or rewrite it.
"""

import numpy as np

import ea_common as ea

NAMES = ['G1', 'G2', 'G3', 'G4', 'G5', 'G6']
H_DOWN = 0.02
SPEED_BOUND = 2.0 * np.sqrt(9.8 * 0.4)


def arrivals(g):
    """Bore arrival (level 2.5 mm above the initial 0.02 m) at G1-G5 [s]."""
    t = g['time']
    out = {}
    for n in NAMES[:5]:
        risen = g['stage_' + n] > H_DOWN + 0.0025
        out[n] = t[np.argmax(risen)] if risen.any() else np.inf
    return out


def gate_period(g):
    """Dominant period of the flow past G6 between 2 and 40 s [s]."""
    t = g['time']
    m = (t > 2.0) & (t <= 40.0)
    y = g['speed_G6'][m]
    y = y - np.polyval(np.polyfit(t[m], y, 2), t[m])
    dt = t[1] - t[0]
    f = np.fft.rfftfreq(len(y), dt)
    power = np.abs(np.fft.rfft(y * np.hanning(len(y)))) ** 2
    band = (f > 1.0 / 20.0) & (f < 1.0 / 4.0)
    return 1.0 / f[band][np.argmax(power[band])]


def evaluate(alg):
    g = ea.load_gauges('test6a_gauges.csv')
    vol = ea.read_volume('test6a_volume.txt')
    vol_err = abs(vol['final'] - vol['initial']) / vol['initial']
    a = arrivals(g)
    period = gate_period(g)

    c = ea.Checks()
    c.add('Volume kept (closed flume)', '%.1e relative' % vol_err, '< 1e-3', vol_err < 1.0e-3)
    c.add('Bore order G1/G2, G3/G4, G5', '%.1f/%.1f, %.1f/%.1f, %.1f s'
          % (a['G1'], a['G2'], a['G3'], a['G4'], a['G5']), 'increasing',
          max(a['G1'], a['G2']) < min(a['G3'], a['G4']) < a['G5'])
    c.add('Period of the gate flow (G6)', '%.1f s' % period, '7-11 s (report ~9 s)',
          7.0 < period < 11.0)
    c.add('Largest speed (depth > 2.5 mm)', '%.2f m/s' % vol['max_speed'],
          '< 2 sqrt(g h0) = %.1f m/s' % SPEED_BOUND, vol['max_speed'] < SPEED_BOUND)
    ea.baseline_check(c, 'test6a_gauges.csv', 'test6a', alg, ['stage_' + n for n in NAMES],
                      {'stage': 0.005})
    return c
