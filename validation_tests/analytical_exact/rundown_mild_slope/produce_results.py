"""
Simple water flow example using ANUGA: Water flowing down a channel.
It was called "steep_slope" in an old validation test.
"""

import anuga
from anuga.validation_utilities import produce_report

args = anuga.get_args()

produce_report('numerical_rundown_channel.py', args=args)



