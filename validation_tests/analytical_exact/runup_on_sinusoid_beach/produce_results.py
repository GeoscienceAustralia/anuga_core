#--------------------------------
# import modules
#--------------------------------
from anuga_validation_tests.utilities.fabricate import *
from anuga_validation_tests.utilities import run_validation_script
from anuga_validation_tests.utilities import typeset_report

# Setup the python scripts which produce the output for this
# validation test
def build():
    run_validation_script('numerical_runup_sinusoid.py')
    run_validation_script('plot_results.py')
    typeset_report()

def clean():
    autoclean()

main()

