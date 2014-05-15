"""
Simple water flow example using ANUGA: Water flowing down a channel.
It was called "steep_slope" in an old validation test.
"""
#--------------------------------
# import modules
#--------------------------------

from anuga.validation_utilities.fabricate import *
from anuga.validation_utilities import run_validation_script
from anuga.validation_utilities import typeset_report

# Setup the python scripts which produce the output for this
# validation test
def build():
    run_validation_script('numerical_rundown_channel_coarse.py')
    run_validation_script('plot_results.py')
    typeset_report()

def clean():
    autoclean()

main()
