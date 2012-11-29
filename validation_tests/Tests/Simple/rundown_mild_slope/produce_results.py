"""
Simple water flow example using ANUGA: Water flowing down a channel.
It was called "steep_slope" in an old validation test.
"""
#--------------------------------
# import modules
#--------------------------------

from fabricate import *
from validation_tests.utilities import run_validation_script

# Setup the python scripts which produce the output for this
# validation test
def build():
    run_validation_script('run_channel.py')
    run_validation_script('plot_results.py')

def clean():
    autoclean()

main()
