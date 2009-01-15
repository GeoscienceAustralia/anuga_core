"""Make directory available as a Python package
"""

# Add path of package to PYTHONPATH to allow C-extensions to be loaded
import sys
sys.path += __path__

# Make selected classes available directly
from shallow_water_domain import Domain,\
     create_domain_from_regions,\
     Transmissive_boundary, Reflective_boundary,\
     Dirichlet_boundary, Time_boundary, File_boundary,\
     Transmissive_momentum_set_stage_boundary,\
     Dirichlet_discharge_boundary,\
     Field_boundary


# FIXME (Ole): Deprecate
from shallow_water_domain import Transmissive_Momentum_Set_Stage_boundary
from shallow_water_domain import Dirichlet_Discharge_boundary
 
