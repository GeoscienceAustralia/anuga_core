#--------------------------------
# import modules
#--------------------------------
from fabricate import *
from validation_tests.utilities import run_validation_script

# Setup the python scripts which produce the output for this
# validation test

def build():
    run_validation_script('create_okushiri.py')
    run_validation_script('run_okushiri.py')
    run_validation_script('compare_timeseries_with_measures.py')

def clean():
    autoclean()

main()
