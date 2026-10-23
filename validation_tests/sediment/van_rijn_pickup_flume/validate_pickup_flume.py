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

    def test_pickup_flume(self):
        if verbose:
            print()
            print(indent + 'Running simulation script')

        s = 'numerical_pickup_flume.py'
        res = anuga.run_anuga_script(s, args=args)

        # Test that script runs ok
        assert res == 0

        if verbose:
            print(indent + 'Testing accuracy')

        import flume_tools as F

        with open('pickup_flume_parameters.json') as f:
            prm = json.load(f)

        t, c, z, w, uh, cc = F.read_centroid_series('pickup_flume.sww')
        x = cc[:, 0]
        h0, q, rho_s, qs_inf = prm['h0'], prm['q'], prm['rho_s'], prm['qs_inf']

        # The flow must be the normal flow throughout, and the bed fixed
        edepth = numpy.abs(w[-1] - z[-1] - h0).max() / h0
        eq = numpy.abs(uh[-1] - q).max() / q
        ebed = numpy.abs(z[-1] - z[0]).max()
        steady = numpy.abs(c[-1] - c[-3]).max() / max(c[-1].max(), 1e-30)

        # Depth-integrated suspended load rho_s c uh (kg/s/m) at the stations
        xd_meas, ratio_meas = F.read_measured('van_rijn_1986_fig3.csv')
        qs = rho_s * c[-1] * uh[-1]
        xs = xd_meas * h0
        qs_num = numpy.array([qs[numpy.abs(x - xi) < prm['dx']].mean() for xi in xs])
        ratio_num = qs_num / qs_inf
        # discrepancy ratio, van Rijn's own measure of a transport prediction
        r = ratio_num / ratio_meas

        print()
        print(indent + 'Rouse number %.2f, v_s %.4f m/s, u* %.4f m/s'
              % (prm['rouse'], prm['v_s'], prm['u_star']))
        print(indent + 'Max relative depth / discharge deviation:   %.3e / %.3e' % (edepth, eq))
        print(indent + 'Max bed movement (m):                       %.3e' % ebed)
        print(indent + 'Change over the last two outputs / max c:   %.3e' % steady)
        print(indent + 'Suspended load q_s / q_s_inf (q_s_inf = %.3f kg/s/m):' % qs_inf)
        for xd, m, n_, ri in zip(xd_meas, ratio_meas, ratio_num, r):
            print(indent + '   x/d = %4.0f   measured %.2f   ANUGA %.3f   ratio %.2f'
                  % (xd, m, n_, ri))

        assert edepth < 1.0e-3, 'the depth is not the normal depth: %.2e' % edepth
        assert eq < 1.0e-3, 'the discharge is not the normal discharge: %.2e' % eq
        assert ebed == 0.0, 'the fixed bed moved: %.2e' % ebed
        assert steady < 1.0e-3, 'not steady: %.2e' % steady
        # The load must grow along the flume, and the prediction is scored
        # the way van Rijn (1986c) scores a transport formula, by its
        # discrepancy ratio: over the four stations the geometric mean must
        # lie in his class 0.5-2, and no station outside 0.33-3. Measured
        # with de Leeuw entrainment: 0.74, 0.66, 0.57, 0.46 (mean 0.60);
        # van Rijn's own SUTRENCH model gives 0.68 at x = 40 d.
        gm = float(numpy.exp(numpy.mean(numpy.log(r))))
        print(indent + 'Geometric-mean discrepancy ratio:           %.2f' % gm)
        assert numpy.all(numpy.diff(qs_num) >= -1e-6 * qs_inf), \
            'suspended load not increasing along the flume'
        assert 0.5 < gm < 2.0, 'suspended load off by more than a factor of two: %.2f' % gm
        assert numpy.all((r > 1.0 / 3.0) & (r < 3.0)), \
            'a station is off by more than a factor of three: %s' % r


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
