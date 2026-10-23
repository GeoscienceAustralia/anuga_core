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

    def test_bed_hump(self):
        if verbose:
            print()
            print(indent + 'Running simulation script')

        s = 'numerical_bed_hump.py'
        res = anuga.run_anuga_script(s, args=args)

        # Test that script runs ok
        assert res == 0

        if verbose:
            print(indent + 'Testing accuracy')

        import analytical_bed_hump as analytic

        with open('bed_hump_parameters.json') as f:
            prm = json.load(f)

        t, z, w, uh, cc = analytic.read_centroid_series('bed_hump.sww')
        x = cc[:, 0]
        order = numpy.argsort(x)
        w0, q, A_g, m, lam = prm['w0'], prm['q'], prm['A_g'], prm['m'], prm['porosity']
        hump = prm['hump']

        # The hump at the end against the characteristic solution, as an L1
        # error relative to the hump volume
        z_ref = analytic.bed(x, t[-1], w0, q, A_g, m, lam, **hump)
        hump_vol = numpy.sum(z[0] - hump['z_flat'])
        ez = numpy.sum(numpy.abs(z[-1] - z_ref)) / hump_vol
        ezmax = numpy.abs(z[-1] - z_ref).max() / hump['height']

        # It migrated: the crest moved by c(crest) t, within a cell
        x_crest_num = x[numpy.argmax(z[-1])]
        x_crest_ref = x[numpy.argmax(z_ref)]
        crest_shift = analytic.wave_speed(hump['z_flat'] + hump['height'],
                                          w0, q, A_g, m, lam) * t[-1]
        ecrest = abs(x_crest_num - x_crest_ref)

        # Bed volume: what enters at the inflow leaves at the outflow
        evol = abs((z[-1] - z[0]).sum()) / hump_vol

        # The flat reaches upstream and downstream of the hump stay put
        far = (x < 150.0) | (x > 900.0)
        eflat = numpy.abs(z[-1][far] - hump['z_flat']).max()

        # The flow stays close to the uniform state: the stage dips by
        # about u^2/g times the relative depth change over the hump
        estage = numpy.abs(w[-1] - w0).max()
        eq = numpy.abs(uh[-1] - q).max() / q

        print()
        print(indent + 'Shock at %.0f s, run to %.0f s; crest moves %.1f m' % (
            prm['t_shock'], t[-1], crest_shift))
        print(indent + 'Relative L1 error in the bed:               %.3e' % ez)
        print(indent + 'Max bed error / hump height:                %.3e' % ezmax)
        print(indent + 'Crest position error (m), cell size %.1f:   %.1f' % (prm['dx'], ecrest))
        print(indent + 'Bed volume change / hump volume:            %.3e' % evol)
        print(indent + 'Flat reaches moved by (m):                  %.3e' % eflat)
        print(indent + 'Max stage deviation (m) / discharge (rel):  %.3e / %.3e' % (estage, eq))

        # The bedload flux is first order (centred with Rusanov dissipation),
        # so the steepening front is smoothed: the L1 error is 9.9% of the
        # hump volume at dx = 10 m and 5.5% at dx = 5 m. The bed volume
        # drifts by 1e-4 of the hump volume because the inflow and outflow
        # cells sit at very slightly different depths under the Dirichlet
        # boundaries, so their bedload fluxes differ by that much.
        assert ez < 0.10, 'bed off by %.2e' % ez
        assert ezmax < 0.15, 'bed off by %.2e of the hump height' % ezmax
        assert ecrest <= 2.0 * prm['dx'], 'crest misplaced by %.1f m' % ecrest
        assert evol < 1.0e-3, 'bed volume not conserved: %.2e' % evol
        assert eflat < 1.0e-6, 'the flat reaches moved: %.2e' % eflat
        assert estage < 0.05, 'the stage moved: %.2e' % estage
        assert eq < 0.05, 'the discharge moved: %.2e' % eq


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
