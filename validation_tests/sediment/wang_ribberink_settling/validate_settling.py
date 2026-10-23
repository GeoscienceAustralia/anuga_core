"""Score the settling-flume runs against Table 1 of Wang and Ribberink (1986)."""
import os
import unittest

import numpy as np

import anuga
from anuga import get_args

args = get_args()
indent = anuga.indent

RUNS = (1, 2) if args.long else (1,)
# the case validates the two-layer closure with the velocity split, the one
# that reproduces the measured decay; the band is the RMS of
# ln(model/measured) over the stations from 1 m on. The instantaneous
# exchange gives 0.73 and 0.58 on the two runs (see results.tex).
ADAPT = 'two_layer'
LOG_RMS = 0.15


def measured(run):
    d = np.loadtxt('wang_ribberink_1986_table1.csv', delimiter=',', comments='#')
    m = d[d[:, 0] == run]
    return m[:, 1], m[:, -1]


def compare(run, adaptation=ADAPT):
    r = np.load('settling_run%d_%s.npz' % (run, adaptation))
    x, c = r['x'], r['c'] * 2650.0 * 1.0e3            # volume -> ppm (mg/l)
    xm, cm = measured(run)
    bins = np.arange(-0.05, 20.05, 0.1)
    idx = np.digitize(x, bins) - 1
    xc = np.array([x[idx == k].mean() for k in range(len(bins) - 1) if (idx == k).any()])
    cc = np.array([c[idx == k].mean() for k in range(len(bins) - 1) if (idx == k).any()])
    model = np.interp(xm, xc, cc)
    sel = xm >= 1.0
    lr = np.log(model[sel] / cm[sel])
    k_meas = -np.polyfit(xm[sel], np.log(cm[sel]), 1)[0]
    k_model = -np.polyfit(xm[sel], np.log(model[sel]), 1)[0]
    return dict(x=xm, meas=cm, model=model, log_rms=float(np.sqrt(np.mean(lr ** 2))),
                bias=float(lr.mean()), k_meas=float(k_meas), k_model=float(k_model))


def simulate(run, adaptation=ADAPT):
    """Run numerical_settling.py for one run, so the score is never read from
    a stale (or, in a clean checkout, missing) result file."""
    out = 'settling_run%d_%s.npz' % (run, adaptation)
    if os.path.exists(out):
        os.remove(out)
    env = {'WANG_RUN': str(run), 'WANG_ADAPTATION': adaptation}
    saved = {k: os.environ.get(k) for k in env}
    os.environ.update(env)
    try:
        if args.verbose:
            print(indent + 'Running simulation script for run %d (%s)' % (run, adaptation))
        res = anuga.run_anuga_script('numerical_settling.py', args=args)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    assert res == 0, 'numerical_settling.py failed for run %d' % run


class Test_results(unittest.TestCase):

    def test_settling(self):
        for run in RUNS:
            simulate(run)
            r = compare(run)
            if args.verbose:
                print(indent + 'run %d: decay rate %.4f /m (measured %.4f); ln-ratio RMS %.3f, bias %+.3f over x >= 1 m'
                      % (run, r['k_model'], r['k_meas'], r['log_rms'], r['bias']))
            assert r['log_rms'] < LOG_RMS, 'run %d: ln-ratio RMS %.3f' % (run, r['log_rms'])


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
