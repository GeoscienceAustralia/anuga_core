#--------------------------------
# import modules
#--------------------------------
from fabricate import *
from validation_tests.utilities import run_validation_script


# Setup the python scripts which produce the output for this
# validation test
def build():
    run_validation_script('numerical_dam_break_dry.py')
    run_validation_script('plot_results.py')
    run('pdflatex', 'report.tex')

def clean():
    autoclean()

main()

