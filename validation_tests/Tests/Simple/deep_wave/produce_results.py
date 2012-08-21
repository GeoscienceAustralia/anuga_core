from fabricate import *
#--------------------------------
# import modules
#--------------------------------
from fabricate import *
from validation_tests.utilities import run_validation_script

# Setup the python scripts which produce the output for this
# validation test
def build():
    run_validation_script('run_wave.py')
    run_validation_script('plotme.py')

def clean():
    autoclean()

main()
