"""Validate runup on sinusoidal beach simulation."""

import sys
import unittest
import os
import anuga

args = anuga.get_args()
indent = anuga.indent
verbose = args.verbose


class Test_results(unittest.TestCase):
    def setUp(self):
        for f in os.listdir('.'):
            if f.endswith(('.sww', '.msh', '.stdout', '.png')):
                os.remove(f)

    def test_runup_on_sinusoid_beach(self):
        if verbose:
            print()
            print(indent + 'Running simulation script')

        res = anuga.run_anuga_script('numerical_runup_sinusoid.py', args=args)
        assert res == 0, 'numerical_runup_sinusoid.py failed with return code %d' % res


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
