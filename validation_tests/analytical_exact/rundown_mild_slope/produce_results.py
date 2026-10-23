"""
Simple water flow example using ANUGA: Water flowing down a channel.
It was called "steep_slope" in an old validation test.
"""

import copy

import anuga
from anuga.validation_utilities import produce_report

import compare_solvers

args = anuga.get_args()

# The solver comparison runs the roll-wave case with every flow algorithm
# (or, under -nr, reuses those runs), so it has to come before the report is
# typeset.
if not getattr(args, 'no_run', False):
    compare_solvers.run_solvers(copy.copy(args))
compare_solvers.plot_solvers()

produce_report('numerical_rundown_channel.py', args=args)
