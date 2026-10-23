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

    def test_sediment_settling(self):
        if verbose:
            print()
            print(indent + 'Running simulation script')

        s = 'numerical_sediment_settling.py'
        res = anuga.run_anuga_script(s, args=args)

        # Test that script runs ok
        assert res == 0

        if verbose:
            print(indent + 'Testing accuracy')

        import analytical_sediment_settling as analytic

        with open('settling_parameters.json') as f:
            prm = json.load(f)

        t, c_num, z_num, w_num = analytic.read_centroid_series('sediment_settling.sww')

        c_ref, z_ref, h_ref = analytic.settling(t, prm['h0'], prm['c0'], prm['v_s'],
                                                prm['d_star'], prm['porosity'])

        # The tank is uniform, so every cell should follow the same curve. Use
        # the cell mean for the comparison and the spread as a uniformity check.
        c_mean = c_num.mean(axis=1)
        dz_mean = (z_num - z_num[0]).mean(axis=1)
        dz_ref = z_ref - z_ref[0]

        # Relative L1 errors over the time series (skip t = 0: dz_ref = 0)
        ec = numpy.sum(numpy.abs(c_mean - c_ref)) / numpy.sum(numpy.abs(c_ref))
        ez = numpy.sum(numpy.abs(dz_mean[1:] - dz_ref[1:])) / numpy.sum(numpy.abs(dz_ref[1:]))
        # Uniformity: no cell may drift from the tank mean
        spread_c = numpy.abs(c_num - c_mean[:, None]).max() / prm['c0']
        # Still water: the free surface must not move
        dw = numpy.abs(w_num - w_num[0]).max()
        # Sediment mass: water column + bed = initial, at every time
        mass = (w_num - z_num) * c_num + (1.0 - prm['porosity']) * (z_num - z_num[0])
        emass = numpy.abs(mass.mean(axis=1) - prm['h0'] * prm['c0']).max() / (prm['h0'] * prm['c0'])

        print()
        print(indent + 'Relative L1 error in concentration decay: %.3e' % ec)
        print(indent + 'Relative L1 error in bed rise:            %.3e' % ez)
        print(indent + 'Max cell-to-cell spread / c0:             %.3e' % spread_c)
        print(indent + 'Max free-surface movement (m):            %.3e' % dw)
        print(indent + 'Max sediment mass imbalance (relative):   %.3e' % emass)
        print(indent + 'Concentration at the end: numerical %.4e, reference %.4e'
              % (c_mean[-1], c_ref[-1]))
        print(indent + 'Bed rise at the end:      numerical %.4e, reference %.4e'
              % (dz_mean[-1], dz_ref[-1]))

        # The deposition source is a first-order fractional step, so the
        # error scales with dt / (h0 / v_s): a fraction of a percent here.
        assert ec < 0.01, 'concentration decay off by %.2e' % ec
        assert ez < 0.01, 'bed rise off by %.2e' % ez
        assert spread_c < 1.0e-6, 'tank not uniform: %.2e' % spread_c
        assert dw < 1.0e-6, 'free surface moved by %.2e m' % dw
        assert emass < 1.0e-3, 'sediment mass not conserved: %.2e' % emass


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
