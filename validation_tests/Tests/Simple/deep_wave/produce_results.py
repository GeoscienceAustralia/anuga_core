

from fabricate import *

# Get the values of alg and cfl and any other
# paarameters set for these validation tests
from validation_tests.parameters import *


# Setup the python scripts which produce the output for this
# validation test
def build():
    run('python', 'run_wave.py',  '-alg', alg, '-cfl', cfl)
    run('python', 'plotme.py', '-alg', alg, '-cfl', cfl )

def clean():
    autoclean()

main()
