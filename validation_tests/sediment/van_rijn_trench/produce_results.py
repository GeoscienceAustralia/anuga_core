import os
import sys

import anuga
from anuga import run_anuga_script
from anuga.validation_utilities import produce_report

args = anuga.get_args()

# The report is written around test 3 (side slopes 1:10); tests 1 and 2,
# the steeper trenches, join it under -l/--long, as in validate_trench.py.
# Test 3 runs last so that its parameters.tex is the one left for the report.
if getattr(args, 'long', False) and not getattr(args, 'no_run', False):
    for test in (1, 2):
        os.environ['VAN_RIJN_TRENCH_TEST'] = str(test)
        res = run_anuga_script('numerical_trench.py', args=args)
        if res != 0:
            sys.exit('numerical_trench.py (test %d) failed with return code %s'
                     % (test, res))

os.environ['VAN_RIJN_TRENCH_TEST'] = '3'
produce_report('numerical_trench.py', args=args)
