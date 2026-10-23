"""Validate Patong beach tsunami simulation.

This test requires large data files (~300 MB) downloaded via data_download.py.
The test is skipped if the required data is not present.

Note: the full simulation (finaltime=10000 s) takes many hours.  Add
'patong' to dirs_to_skip in run_auto_validation_tests.py if you want to
exclude it from automated runs.
"""

import sys
import unittest
import os
import anuga

args = anuga.get_args()
indent = anuga.indent
verbose = args.verbose


def _data_available():
    """Return True if the required topography data is present."""
    try:
        import project
        return os.path.isdir(project.topographies_folder)
    except Exception:
        return False


class Test_results(unittest.TestCase):
    def setUp(self):
        for f in os.listdir('.'):
            if f.endswith(('.sww', '.msh', '.stdout', '.png')):
                os.remove(f)

    def test_patong(self):
        if not _data_available():
            self.skipTest('Patong data not available — run data_download.py first')

        if verbose:
            print()
            print(indent + 'Running patong simulation (this may take many hours)')

        res = anuga.run_anuga_script('run_model.py', args=args)
        assert res == 0, 'run_model.py failed with return code %d' % res


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
