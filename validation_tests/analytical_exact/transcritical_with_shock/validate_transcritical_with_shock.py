"""Validate transcritical flow with shock against analytical solution."""

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

    def test_transcritical_with_shock(self):
        if verbose:
            print()
            print(indent + 'Running simulation script')

        res = anuga.run_anuga_script('numerical_transcritical.py', args=args)
        assert res == 0, 'numerical_transcritical.py failed with return code %d' % res

        if verbose:
            print(indent + 'Testing accuracy against analytical solution')

        import anuga.utilities.plot_utils as util
        import analytical_with_shock as analytic
        from numpy import arange

        p_st = util.get_output('transcritical.sww')
        p2_st = util.get_centroids(p_st)

        v2 = arange(len(p2_st.y))
        h, z = analytic.analytic_sol(p2_st.x[v2])

        stage_num = p2_st.stage[-1, v2]
        stage_ana = h + z
        denom = numpy.sum(numpy.abs(stage_ana))
        assert denom > 0.0, 'Analytical stage sum is zero — check domain/solution mismatch'
        err = numpy.sum(numpy.abs(stage_num - stage_ana)) / denom

        print()
        print(indent + 'Stage L1 relative error: %.2e' % err)

        assert err < 0.05, 'Stage L1 relative error %.2e exceeds 0.05' % err


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
