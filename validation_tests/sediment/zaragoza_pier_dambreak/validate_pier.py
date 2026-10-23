"""Score a pier dam-break run against the Kinect bed and water-surface data."""
import json
import os
import sys
import unittest

import numpy as np
from scipy.interpolate import LinearNDInterpolator

import anuga
from anuga import get_args

args = get_args()
indent = anuga.indent

CASES = ('P1', 'P2') if args.long else ('P1',)
# pass band (provisional): bed change after the first dam-break, RMS over the
# Kinect window. The Kinect field itself has an RMS of 3.6 (P1) and 4.0 mm
# (P2) about zero, so 6 mm is "no worse than the scatter". The scour depth
# at the pier is reported, not tested: ANUGA gets a third of it (see
# results.tex), which is the open item of the case.
RMS_MM = 6.0


def load_bed(case):
    d = np.loadtxt('kinect_%s_bed.csv' % case, delimiter=',', comments='#')
    return d[:, 0], d[:, 1], d[:, 2:]


def load_wse(case):
    d = np.loadtxt('kinect_%s_wse.csv' % case, delimiter=',', comments='#')
    return d[:, 0], d[:, 1], d[:, 2]


def model_dz(case, event=1):
    """ANUGA bed change after `event` dam-breaks, interpolated from the
    centroids (mm), as a callable on (x, y)."""
    r = np.load('pier_%s_bed.npz' % case)
    dz = 1000.0 * (r['beds'][event - 1] - r['z0'])
    return LinearNDInterpolator(np.column_stack([r['x'], r['y']]), dz), r


def compare(case, event=1):
    x, y, meas = load_bed(case)
    f, r = model_dz(case, event)
    mod = f(x, y)
    ok = np.isfinite(mod)
    m = meas[:, event - 1][ok]
    a = mod[ok]
    err = a - m
    # the deepest scour within 15 cm of the pier, and the mean over a ring
    # 2-4 cm from its centre (the scour hole) and over the wake
    near = (np.abs(x[ok] - 1.6) < 0.15)
    rad = np.hypot(x[ok] - 1.6, y[ok] - 0.12)
    ring = (rad > 0.02) & (rad < 0.04)
    wake = (x[ok] > 1.63) & (x[ok] < 1.75) & (np.abs(y[ok] - 0.12) < 0.03)
    return dict(n=int(ok.sum()), rms=float(np.sqrt(np.mean(err ** 2))),
                bias=float(err.mean()),
                ring_meas=float(m[ring].mean()), ring_model=float(a[ring].mean()),
                wake_meas=float(m[wake].mean()), wake_model=float(a[wake].mean()),
                scour_meas=float(-m[near].min()), scour_model=float(-a[near].min()),
                fill_meas=float(m[near].max()), fill_model=float(a[near].max()),
                eroded_meas=float(np.sum(np.minimum(m, 0.0)) * 25e-6 / 1000.0),   # m^3 over 5 mm blocks
                eroded_model=float(np.sum(np.minimum(a, 0.0)) * 25e-6 / 1000.0))


def simulate(case):
    """Run numerical_pier.py for one pier, so the score is never read from a
    stale (or, in a clean checkout, missing) result file. One dam-break is
    enough: the score is the bed after the first."""
    out = 'pier_%s_bed.npz' % case
    if os.path.exists(out):
        os.remove(out)
    env = {'ZARAGOZA_CASE': case, 'ZARAGOZA_EVENTS': '1', 'ZARAGOZA_TAG': ''}
    saved = {k: os.environ.get(k) for k in env}
    os.environ.update(env)
    try:
        if args.verbose:
            print(indent + 'Running simulation script for %s' % case)
        res = anuga.run_anuga_script('numerical_pier.py', args=args)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    assert res == 0, 'numerical_pier.py failed for %s' % case


class Test_results(unittest.TestCase):

    def test_pier_scour(self):
        for case in CASES:
            simulate(case)
            r = compare(case)
            if args.verbose:
                print(indent + '%s after dam-break 1: bed change RMS %.1f mm, bias %.1f mm over %d points; '
                      'scour ring 2-4 cm %.1f mm (measured %.1f), deepest %.1f (%.1f); wake %+.1f (%+.1f)'
                      % (case, r['rms'], r['bias'], r['n'], r['ring_model'], r['ring_meas'],
                         r['scour_model'], r['scour_meas'], r['wake_model'], r['wake_meas']))
            assert r['rms'] < RMS_MM, '%s: RMS %.1f mm' % (case, r['rms'])


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
