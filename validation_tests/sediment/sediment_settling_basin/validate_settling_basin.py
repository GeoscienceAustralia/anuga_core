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

    def test_settling_basin(self):
        if verbose:
            print()
            print(indent + 'Running simulation script')

        s = 'numerical_settling_basin.py'
        res = anuga.run_anuga_script(s, args=args)

        # Test that script runs ok
        assert res == 0

        if verbose:
            print(indent + 'Testing accuracy')

        import analytical_settling_basin as analytic

        with open('settling_basin_parameters.json') as f:
            prm = json.load(f)

        t, c, z, w, uh, cc = analytic.read_centroid_series('settling_basin.sww')
        x = cc[:, 0]
        q, v_s, d_star = prm['q'], prm['v_s'], prm['d_star']

        # The flow must stay the uniform steady state. The bed rises by the
        # deposit (a few mm here) and the depth falls with it while the stage
        # and the discharge per unit width are held by the boundaries; the
        # velocity rises by the same fraction, and the dynamic head that
        # costs (u du / g ~ 5e-5 m) is the only stage change to expect.
        estage = numpy.abs(w - prm['h0']).max() / prm['h0']
        eq = numpy.abs(uh - q).max() / q

        # Steady concentration at the end, and steady means steady
        c_ref = analytic.concentration(x, prm['c0'], q, v_s, d_star)
        ec = numpy.sum(numpy.abs(c[-1] - c_ref)) / numpy.sum(c_ref)
        steady = numpy.abs(c[-1] - c[-3]).max() / prm['c0']
        # Cross-channel uniformity: cells at (nearly) the same x agree
        order = numpy.argsort(x)
        spread = 0.0
        for k in range(0, len(order), 4):
            grp = order[k:k + 4]
            spread = max(spread, numpy.ptp(c[-1][grp]) / c_ref[grp].mean())

        # Bed-rise rate over the last two outputs against d* v_s c / (1 - lambda)
        rate_num = (z[-1] - z[-3]) / (t[-1] - t[-3])
        rate_ref = analytic.bed_rise_rate(x, prm['c0'], q, v_s, d_star, prm['porosity'])
        ez = numpy.sum(numpy.abs(rate_num - rate_ref)) / numpy.sum(rate_ref)
        assert (z[-1] > z[0]).all(), 'the bed did not rise everywhere'

        # Sediment budget: in - out = water column + bed
        budget = prm['boundary_flux_in'] - prm['water_column_change'] - prm['bed_sediment_volume']
        ebudget = abs(budget) / prm['boundary_flux_in']

        print()
        print(indent + 'Settling length L_s = %.1f m, channel %.0f m' % (
            analytic.settling_length(q, v_s, d_star), prm['L']))
        print(indent + 'Max relative stage / discharge deviation:   %.3e / %.3e' % (estage, eq))
        print(indent + 'Max bed rise (m):                           %.3e' % (z[-1] - z[0]).max())
        print(indent + 'Relative L1 error in steady concentration:  %.3e' % ec)
        print(indent + 'Change over the last two outputs / c0:      %.3e' % steady)
        print(indent + 'Max cross-channel spread / local c:         %.3e' % spread)
        print(indent + 'Relative L1 error in bed-rise rate:         %.3e' % ez)
        print(indent + 'Sediment budget closure (relative):         %.3e' % ebudget)
        print(indent + 'Concentration at the outflow: numerical %.3e, reference %.3e'
              % (c[-1][order[-4:]].mean(), c_ref[order[-4:]].mean()))

        # The tracer advection is first-order upwind, so the steady profile
        # carries an O(dx / L_s) error; measured 5e-4 at dx = 2.5 m and
        # L_s = 62 m. The bed-rise rate inherits it. The concentration keeps
        # drifting at the 1e-5 level as the bed rises, and the cross-channel
        # spread reaches 1e-3 only in the last row of cells at the outflow
        # boundary, where c is 1% of c0.
        assert estage < 1.0e-3, 'the stage moved: %.2e' % estage
        assert eq < 1.0e-3, 'the discharge is not the boundary value: %.2e' % eq
        assert ec < 0.01, 'steady concentration off by %.2e' % ec
        assert steady < 1.0e-4, 'not steady: %.2e' % steady
        assert spread < 5.0e-3, 'concentration varies across the channel: %.2e' % spread
        assert ez < 0.01, 'bed-rise rate off by %.2e' % ez
        assert ebudget < 1.0e-3, 'sediment budget does not close: %.2e' % ebudget


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
