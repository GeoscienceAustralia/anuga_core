"""Validate Towradgi catchment flood simulation runs without error.

Requires data downloaded via data_download.py.  The test is skipped if the
required DEM data directory is not present.
"""

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

    def test_towradgi(self):
        if not os.path.isdir('DEM_bridges'):
            self.skipTest('Towradgi data not available — run data_download.py first')

        if verbose:
            print()
            print(indent + 'Running Towradgi simulation')

        res = anuga.run_anuga_script('run_towradgi.py', args=args)
        assert res == 0, 'run_towradgi.py failed with return code %d' % res


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
