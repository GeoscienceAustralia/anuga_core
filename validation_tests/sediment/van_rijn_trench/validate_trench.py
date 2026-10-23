"""Automatic verification of ANUGA flows.
See functions exercised by this wrapper for more details
"""
import json
import os
import sys
import unittest

import numpy
import anuga

args = anuga.get_args()
indent = anuga.indent
verbose = args.verbose

# The gentle 1:10 trench, the one a depth-averaged model is expected to get,
# runs routinely (about 5 min with the morphological factor); the 1:3 and 1:7
# trenches join it under -l/--long.
TESTS = (1, 2, 3) if args.long else (3,)


def compare(test):
    """RMS and max bed-level error at the measured points, m, plus the
    position and depth of the deepest point of the fill."""
    import flume_tools as F
    with open('trench_test%d_parameters.json' % test) as f:
        prm = json.load(f)
    t, c, z, w, uh, cc = F.read_centroid_series('trench_test%d.sww' % test)
    t = t * prm.get('morfac', 1.0)                   # morphological time
    xf = cc[:, 0] - prm['X_OFF']                     # figure coordinate
    S, L = prm['S'], prm['L']
    below = S * (L - cc[:, 0]) - z[-1]               # depth below the plane
    meas = numpy.loadtxt('van_rijn_1986_fig16_measured.csv', delimiter=',', comments='#')
    meas = meas[meas[:, 0] == test][:, 1:]
    bins = numpy.arange(-0.05, 11.05, prm['dx'])
    xm, zm = F.along_flume(xf, below, bins)
    z_num = numpy.interp(meas[:, 0], xm, zm)
    err = z_num - meas[:, 1]
    k = numpy.argmax(zm[(xm > 4.0) & (xm < 11.0)])
    xsel = xm[(xm > 4.0) & (xm < 11.0)]
    zsel = zm[(xm > 4.0) & (xm < 11.0)]
    km = numpy.argmax(meas[:, 1])
    return dict(t=t[-1], rms=float(numpy.sqrt(numpy.mean(err ** 2))),
                emax=float(numpy.abs(err).max()),
                x_deep=float(xsel[k]), z_deep=float(zsel[k]),
                x_deep_meas=float(meas[km, 0]), z_deep_meas=float(meas[km, 1]),
                scale=prm['scale'], calibration=prm['calibration'])


class Test_results(unittest.TestCase):

    def setUp(self):
        for file in os.listdir('.'):
            if file.endswith('.stdout') or \
                    file.endswith('.sww') or \
                    file.endswith('.msh') or \
                    file.endswith('.png'):
                os.remove(file)

    def tearDown(self):
        pass

    def test_trench(self):
        results = {}
        for test in TESTS:
            if verbose:
                print()
                print(indent + 'Running simulation script, test %d' % test)
            os.environ['VAN_RIJN_TRENCH_TEST'] = str(test)
            res = anuga.run_anuga_script('numerical_trench.py', args=args)
            assert res == 0
            results[test] = compare(test)

        if verbose:
            print(indent + 'Testing accuracy')
        print()
        for test, r in results.items():
            print(indent + 'test %d (t = %.0f h): entrainment scale %.3f; bed after 15 h vs measured:'
                  % (test, r['t'] / 3600.0, r['scale']))
            print(indent + '   RMS error %.4f m, max error %.4f m' % (r['rms'], r['emax']))
            print(indent + '   deepest point of the remaining trench: ANUGA x = %.2f m, %.3f m; '
                  'measured x = %.2f m, %.3f m'
                  % (r['x_deep'], r['z_deep'], r['x_deep_meas'], r['z_deep_meas']))

        # The measured fill is 0.08-0.10 m deep over 7 m of flume. Van Rijn's
        # own model matched it to about 0.01 m; a depth-averaged model with
        # an instantaneous near-bed adjustment is held to 0.03 m RMS and to
        # placing the deepest remaining point within 1.5 m of the measured one.
        for test, r in results.items():
            assert r['rms'] < 0.03, 'test %d: bed RMS error %.3f m' % (test, r['rms'])
            assert abs(r['x_deep'] - r['x_deep_meas']) < 1.5, \
                'test %d: deepest point at %.2f m, measured %.2f m' % (test, r['x_deep'], r['x_deep_meas'])


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
