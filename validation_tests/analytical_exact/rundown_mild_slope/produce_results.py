"""
Simple water flow example using ANUGA: Water flowing down a channel.
It was called "steep_slope" in an old validation test.
"""
#--------------------------------
# import modules
#--------------------------------

from anuga_validation_tests.utilities.fabricate import *
from anuga_validation_tests.utilities import run_validation_script
from anuga_validation_tests.utilities import typeset_report

# Setup the python scripts which produce the output for this
# validation test
def build():
    run_validation_script('numerical_rundown_channel.py')
    run_validation_script('plot_results.py')
    typeset_report()

def clean():
    autoclean()

main()
