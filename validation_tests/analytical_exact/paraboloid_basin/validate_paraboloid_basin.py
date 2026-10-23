"""Validate paraboloid basin oscillation (Thacker) against analytical solution."""

import sys
import unittest
import os
import numpy
import anuga

args = anuga.get_args()
indent = anuga.indent
verbose = args.verbose


class Test_results(unittest.TestCase):
    def setUp(self):
        for f in os.listdir('.'):
            if f.endswith(('.sww', '.msh', '.stdout', '.png')):
                os.remove(f)

    def test_paraboloid_basin(self):
        if verbose:
            print()
            print(indent + 'Running simulation script')

        res = anuga.run_anuga_script('numerical_paraboloid_basin.py', args=args)
        assert res == 0, 'numerical_paraboloid_basin.py failed with return code %d' % res

        if verbose:
            print(indent + 'Testing accuracy against analytical solution')

        import anuga.utilities.plot_utils as util
        from analytical_paraboloid_basin import analytic_sol
        from numpy import arange

        p_st = util.get_output('paraboloid.sww')
        p2_st = util.get_centroids(p_st)

        # Use all centroids with both x and y coordinates
        v2 = arange(len(p2_st.y))
        time_level = -1
        w, u, h = analytic_sol(p2_st.x[v2], p2_st.y[v2], p2_st.time[time_level])

        stage_num = p2_st.stage[time_level, v2]
        stage_ana = w

        # Compute error only over wet cells
        wet = h > 1e-6
        if numpy.any(wet):
            err = numpy.sum(numpy.abs(stage_num[wet] - stage_ana[wet])) / numpy.sum(numpy.abs(stage_ana[wet]))
        else:
            err = 0.0

        print()
        print(indent + 'Stage L1 relative error: %.2e' % err)

        assert err < 0.05, 'Stage L1 relative error %.2e exceeds 0.05' % err


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
