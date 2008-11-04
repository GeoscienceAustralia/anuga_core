"""Make directory available as a Python package
"""

# Add path of package to PYTHONPATH to allow C-extensions to be loaded
import sys
sys.path += __path__

# Make selected classes available directly
from shallow_water_domain import Domain,\
    Transmissive_boundary, Reflective_boundary,\
    Dirichlet_boundary, Time_boundary, File_boundary,\
    Transmissive_Momentum_Set_Stage_boundary,\
    Dirichlet_Discharge_boundary,\
    Field_boundary


