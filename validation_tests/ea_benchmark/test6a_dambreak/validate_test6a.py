"""Validate EA benchmark Test 6A (dam break past a building, laboratory scale).

Runs run_test6a.py and asserts the checks in checks.py, which the report
(plot_results.py) shows as a table. Set EA_REGENERATE_BASELINE=1 to write
or rewrite the stored baseline for the algorithm being run.
"""

import os
import sys
import unittest

import anuga

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import checks  # noqa: E402

args = anuga.get_args()
indent = anuga.indent


class Test_results(unittest.TestCase):
    def setUp(self):
        # Only the default run's outputs: variants (EA_6A_MAX_AREA, ...) are kept.
        for f in ('test6a.sww', 'test6a_gauges.csv', 'test6a_volume.txt',
                  'test6a_geometry.json'):
            if os.path.exists(f):
                os.remove(f)

    def test_test6a_dambreak(self):
        res = anuga.run_anuga_script('run_test6a.py', args=args, allow_parallel=False)
        assert res == 0, 'run_test6a.py failed with return code %d' % res

        c = checks.evaluate(args.alg)
        print()
        c.report(indent)
        assert not c.failed(), 'Failed: %s' % [r[0] for r in c.failed()]


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_results)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
