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

    def test_sediment_erosion(self):
        if verbose:
            print()
            print(indent + 'Running simulation script')

        s = 'numerical_sediment_erosion.py'
        res = anuga.run_anuga_script(s, args=args)

        # Test that script runs ok
        assert res == 0

        if verbose:
            print(indent + 'Testing accuracy')

        import analytical_sediment_erosion as analytic

        with open('erosion_parameters.json') as f:
            prm = json.load(f)

        t, c_num, z_num, w_num, uh, vh = analytic.read_centroid_series('sediment_erosion.sww')
        _, c_ctl, _, _, _, _ = analytic.read_centroid_series('sediment_erosion.sww',
                                                             name='control')

        h_cell = numpy.array(prm['h_cell'])
        law = (prm['v_s'], prm['R'], prm['diameter'], prm['tau_c_star'],
               prm['gamma0'], prm['d_star'])

        # Every cell against its own curve, at the one slope the bed was
        # built with (frozen for the run): the kernel's slope estimate must
        # reproduce it in every cell, the ones along the walls and in the
        # corners included, and the bed must lower by the entrained volume.
        c_ref, z_ref, h_ref = analytic.erosion_with_feedback(
            t, h_cell, prm['S_bed'], *law, porosity=prm['porosity'], z0=z_num[0])
        c_fixed = analytic.erosion(t, h_cell, prm['S_bed'], *law)

        # Relative L1 errors over the time series, all cells
        ec = numpy.sum(numpy.abs(c_num - c_ref)) / numpy.sum(numpy.abs(c_ref))
        # ... and the worst single cell, so a wall cell cannot hide in the mean
        ec_cell = (numpy.abs(c_num - c_ref).sum(axis=0)
                   / numpy.abs(c_ref).sum(axis=0)).max()
        dz_num = z_num - z_num[0]
        dz_ref = z_ref - z_ref[0]
        ez = numpy.sum(numpy.abs(dz_num[1:] - dz_ref[1:])) / numpy.sum(numpy.abs(dz_ref[1:]))
        # Below threshold: nothing may be entrained into the control fraction
        c_control = numpy.abs(c_ctl).max()
        # The bed feedback is real: the fixed-bed curve must be distinguishable
        feedback = numpy.abs(c_ref[-1] - c_fixed[-1]).max() / c_fixed[-1].max()
        # Sediment mass: water column + bed = 0 (nothing to start with)
        mass = (w_num - z_num) * c_num + (1.0 - prm['porosity']) * dz_num
        emass = numpy.abs(mass).max() / ((w_num - z_num) * c_num).max()
        # Still water: nothing else may move
        dw = numpy.abs(w_num - w_num[0]).max()
        mom = max(numpy.abs(uh).max(), numpy.abs(vh).max())

        print()
        print(indent + 'Relative L1 error in concentration, all cells:  %.3e' % ec)
        print(indent + 'Relative L1 error in concentration, worst cell: %.3e' % ec_cell)
        print(indent + 'Relative L1 error in bed lowering:              %.3e' % ez)
        print(indent + 'Max concentration of the control fraction:      %.3e' % c_control)
        print(indent + 'Bed feedback on c at the end (vs fixed bed):    %.3e' % feedback)
        print(indent + 'Max sediment mass imbalance (relative):         %.3e' % emass)
        print(indent + 'Max free-surface movement (m):                  %.3e' % dw)
        print(indent + 'Max |momentum|:                                 %.3e' % mom)
        print(indent + 'Concentration at the end: numerical %.4e, reference %.4e'
              % (c_num[-1].mean(), c_ref[-1].mean()))
        print(indent + 'Bed lowering at the end:  numerical %.4e, reference %.4e'
              % (dz_num[-1].mean(), dz_ref[-1].mean()))

        # The source is a first-order fractional step, so the error scales
        # with dt / (h / v_s): a few 1e-4 here.
        assert ec < 0.01, 'concentration relaxation off by %.2e' % ec
        assert ec_cell < 0.01, 'a cell is off by %.2e' % ec_cell
        assert ez < 0.01, 'bed lowering off by %.2e' % ez
        assert c_control < 1.0e-12, 'entrainment below threshold: %.2e' % c_control
        assert feedback > 1.0e-3, 'the bed did not feed back on the concentration'
        assert emass < 1.0e-3, 'sediment mass not conserved: %.2e' % emass
        assert dw < 1.0e-6, 'free surface moved by %.2e m' % dw
        assert mom < 1.0e-8, 'still water acquired momentum %.2e' % mom


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
