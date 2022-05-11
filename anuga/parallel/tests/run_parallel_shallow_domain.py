import os
        
#------------------------------------------
# anuga imports
#------------------------------------------
import anuga 

from anuga.utilities.system_tools import get_pathname_from_package

from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary
from anuga import Transmissive_boundary

from anuga import rectangular_cross
from anuga import create_domain_from_file

from anuga import distribute, myid, numprocs, finalize

#----------------------------------
# set up MPI to abort on error
#----------------------------------
from anuga.utilities.parallel_abstraction import global_except_hook
import sys
sys.excepthook = global_except_hook

# ----------------
# Setup parameters
# ----------------

mod_path = get_pathname_from_package('anuga.parallel')
mesh_filename = os.path.join(mod_path, 'data', 'merimbula_10785_1.tsh')
verbose = False

# ----------------
# Setup procedures
# ----------------
class Set_Stage(object):
    """Set an initial condition with constant water height, for x<x0
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h = h

    def __call__(self, x, y):
        return self.h * ((x > self.x0) & (x < self.x1))


domain = create_domain_from_file(mesh_filename)
domain.set_quantity('stage', Set_Stage(756000.0, 756500.0, 2.0))

if numprocs > 1:
    if myid == 0 and verbose: print('PARALLEL EVOLVE')
    domain.set_name('shallow_water_parallel')        
else:
    if verbose: print('SEQUENTIAL EVOLVE')
    domain.set_name('shallow_water_sequential')        


# -----------------------------------
# Create parallel domain if requested
# -----------------------------------
if numprocs > 1:
    if myid == 0 and verbose: print('DISTRIBUTING PARALLEL DOMAIN')
    domain = distribute(domain)

# --------------------------------------------------------------
# Setup boundary conditions
# This must currently happen *after* domain has been distributed
# --------------------------------------------------------------
Br = Reflective_boundary(domain)      # Solid reflective wall - no movement

domain.set_boundary({'exterior': Br, 'open': Br})

# ---------
# Evolution
# ---------

for t in domain.evolve(yieldstep=1.0, finaltime=2.0):
    if myid == 0 and verbose: domain.write_time()

# ------------------------------------
# Wrap up parallel matters if required
# ------------------------------------
domain.sww_merge(delete_old=True)
finalize()



