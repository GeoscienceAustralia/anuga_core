"""Make directory available as a Python package
"""

# Add path of package to PYTHONPATH to allow C-extensions to be loaded
import sys
sys.path += __path__

# Make selected classes available directly
from boundaries import Reflective_boundary,\
     Transmissive_momentum_set_stage_boundary,\
     Dirichlet_discharge_boundary,\
     Field_boundary,\
     Transmissive_stage_zero_momentum_boundary,\
     Transmissive_n_momentum_zero_t_momentum_set_stage_boundary

from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary, Dirichlet_boundary, \
            Time_boundary, File_boundary, AWI_boundary

from shallow_water_domain import Domain


#from shallow_water_balanced_domain import Swb_domain


 
