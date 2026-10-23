"""Validate Merewether flood scenario simulation runs without error."""

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

    def test_merewether(self):
        if not os.path.exists('topography1.zip'):
            self.skipTest('Merewether topography data not available')

        if verbose:
            print()
            print(indent + 'Running Merewether simulation')

        res = anuga.run_anuga_script('runMerewether.py', args=args)
        assert res == 0, 'runMerewether.py failed with return code %d' % res


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
