"""Make directory available as a Python package
"""

# Add path of package to PYTHONPATH to allow C-extensions to be loaded
import sys
sys.path += __path__

# Make selected classes available directly
from swb_domain import Domain,\
     Transmissive_boundary, Reflective_boundary,\
     Dirichlet_boundary, Time_boundary, File_boundary,\
     Transmissive_momentum_set_stage_boundary,\
     Dirichlet_discharge_boundary,\
     Field_boundary,\
     Transmissive_stage_zero_momentum_boundary,\
     Transmissive_n_momentum_zero_t_momentum_set_stage_boundary,\
     AWI_boundary



 
