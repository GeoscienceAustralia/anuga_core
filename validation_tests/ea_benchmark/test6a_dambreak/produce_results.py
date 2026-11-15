import anuga
from anuga import run_anuga_script
from anuga.validation_utilities import produce_report

args = anuga.get_args()

# The limiter study first: each run writes parameters.tex, and the main
# run's must be the one left for the report.
if not getattr(args, 'no_run', False):
    res = run_anuga_script('run_limiter_study.py', args=args, allow_parallel=False)
    if res != 0:
        raise SystemExit('run_limiter_study.py failed with return code %s' % res)

produce_report('run_test6a.py', args=args)
