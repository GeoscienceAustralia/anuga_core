#--------------------------------
# import modules
#--------------------------------
from fabricate import *
from validation_tests.utilities import run_validation_script


# Setup the python scripts which produce the output for this
# validation test
def build():
    run_validation_script('numerical_parabolic_basin.py')
    run_validation_script('plot_results_cross_section.py')
    run_validation_script('plot_results_origin_wrt_time.py')    

def clean():
    autoclean()

main()
