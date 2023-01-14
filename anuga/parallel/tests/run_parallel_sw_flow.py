"""Run sw_flow simulation (sequentially or in parallel) to support 
   test_parallel_sw_flow.py
"""

# ------------------------
# Import necessary modules
# ------------------------
from math import exp
import anuga
from anuga import rectangular_cross_domain
from anuga import Reflective_boundary, Dirichlet_boundary
from anuga import myid, distribute, barrier, numprocs, finalize
import numpy as num


#-----------------
# Setup parameters
#-----------------
verbose = True

#-------------------------------------
# Setup function for initial condition
#-------------------------------------
def topography(x, y): 
    return -x / 2    

#------------------------------------------
# Setup computational domain and quantities
#------------------------------------------
domain = rectangular_cross_domain(29, 29)
domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.0)         # Constant friction 
domain.set_quantity('stage', expression='elevation') # Dry initial stage

#------------------------
# Setup domain parameters
#------------------------
domain.set_datadir('.')          # Set output dir

# ----------------------------------------------
# Decide if this is a sequential or parallel run
# ----------------------------------------------
if numprocs == 1:
    # This is a sequential run
    domain.set_name('sw_flow_sequential')
else:
    # This is a parallel run
    domain = distribute(domain, verbose=verbose)
    domain.set_name('sw_flow_parallel')

#---------------------------------------------------------------
# Setup boundary conditions
# This must currently happen *AFTER* domain has been distributed
#---------------------------------------------------------------
Br = Reflective_boundary(domain)         # Solid reflective wall
Bd = Dirichlet_boundary([-0.2, 0., 0.])  # Constant boundary values

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})
        
#---------------------------
# Evolve system through time
#---------------------------
for t in domain.evolve(yieldstep=0.25, finaltime=1.0):
    if myid == 0 and verbose: domain.print_timestepping_statistics()

#-------------------------------------
# Wrap up parallel matters if required
#-------------------------------------
domain.sww_merge(delete_old=True)
finalize()

