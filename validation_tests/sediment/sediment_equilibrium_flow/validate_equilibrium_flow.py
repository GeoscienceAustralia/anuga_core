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

    def test_equilibrium_flow(self):
        if verbose:
            print()
            print(indent + 'Running simulation script')

        s = 'numerical_equilibrium_flow.py'
        res = anuga.run_anuga_script(s, args=args)

        # Test that script runs ok
        assert res == 0

        if verbose:
            print(indent + 'Testing accuracy')

        import analytical_equilibrium_flow as analytic

        with open('equilibrium_flow_parameters.json') as f:
            prm = json.load(f)

        t, c, z, w, uh, cc = analytic.read_centroid_series('equilibrium_flow.sww')
        x = cc[:, 0]
        q, v_s, d_star, c_eq = prm['q'], prm['v_s'], prm['d_star'], prm['c_eq']
        h0 = prm['h0']

        # The flow must be the normal flow throughout, and the bed fixed
        edepth = numpy.abs(w[-1] - z[-1] - h0).max() / h0
        eq = numpy.abs(uh[-1] - q).max() / q
        ebed = numpy.abs(z[-1] - z[0]).max()

        # Steady concentration at the end, and steady means steady
        c_ref = analytic.concentration(x, c_eq, q, v_s, d_star)
        ec = numpy.sum(numpy.abs(c[-1] - c_ref)) / numpy.sum(c_ref)
        steady = numpy.abs(c[-1] - c[-3]).max() / c_eq
        # Cross-channel uniformity: cells at (nearly) the same x agree
        order = numpy.argsort(x)
        spread = 0.0
        for k in range(0, len(order), 4):
            grp = order[k:k + 4]
            spread = max(spread, numpy.ptp(c[-1][grp]) / c_ref[grp].mean())
        # The outflow end has all but reached equilibrium: L / L_s = 5, so
        # c_ref there is 0.6% below c_eq, and the numerical value must agree
        # with c_ref, not with c_eq.
        c_out = c[-1][order[-4:]].mean()
        c_out_ref = c_ref[order[-4:]].mean()
        eeq = abs(c_out - c_out_ref) / c_out_ref

        print()
        print(indent + 'Shields stress tau* = %.3g, c_eq = %.4g, L_s = %.1f m, channel %.0f m' % (
            analytic.shields_stress(h0, prm['S'], prm['R'], prm['diameter']), c_eq,
            analytic.settling_length(q, v_s, d_star), prm['L']))
        print(indent + 'Max relative depth / discharge deviation:   %.3e / %.3e' % (edepth, eq))
        print(indent + 'Max bed movement (m):                       %.3e' % ebed)
        print(indent + 'Relative L1 error in steady concentration:  %.3e' % ec)
        print(indent + 'Change over the last two outputs / c_eq:    %.3e' % steady)
        print(indent + 'Max cross-channel spread / local c:         %.3e' % spread)
        print(indent + 'Outflow concentration: numerical %.4e, reference %.4e, c_eq %.4e'
              % (c_out, c_out_ref, c_eq))

        # The tracer advection is first-order upwind, so the profile carries
        # an O(dx / L_s) error; measured 2e-4 at dx = 2.5 m and L_s = 59 m.
        # The cross-channel spread reaches 5e-3 only in the first row of
        # cells at the inflow, where c is 2% of c_eq and the entry error of
        # the upwind scheme dominates; from the second row on it is 1e-3.
        assert edepth < 1.0e-3, 'the depth is not the normal depth: %.2e' % edepth
        assert eq < 1.0e-3, 'the discharge is not the normal discharge: %.2e' % eq
        assert ebed == 0.0, 'the fixed bed moved: %.2e' % ebed
        assert ec < 0.01, 'steady concentration off by %.2e' % ec
        assert steady < 1.0e-4, 'not steady: %.2e' % steady
        assert spread < 1.0e-2, 'concentration varies across the channel: %.2e' % spread
        assert eeq < 1.0e-3, 'outflow concentration off by %.2e' % eeq


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
