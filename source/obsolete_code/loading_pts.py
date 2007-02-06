"""Validation of the AnuGA implementation of the shallow water wave equation.

This script sets up Okushiri Island benchmark as published at the

THE THIRD INTERNATIONAL WORKSHOP ON LONG-WAVE RUNUP MODELS
June 17-18 2004
Wrigley Marine Science Center
Catalina Island, California
http://www.cee.cornell.edu/longwave/


The validation data was downloaded and made available in this directory
for convenience but the original data is available at
http://www.cee.cornell.edu/longwave/index.cfm?page=benchmark&problem=2
where a detailed description of the problem is also available.


Run create_okushiri.py to process the boundary condition and build a the
mesh before running this script.

"""

# Module imports
from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.shallow_water import Transmissive_Momentum_Set_Stage_boundary
from anuga.abstract_2d_finite_volumes.util import file_function

from anuga_parallel.parallel_api import myid, numprocs, distribute

import project


#-------------------------
# Create Domain from mesh
#-------------------------
domain = Domain(project.mesh_filename, use_cache=True, verbose=True)
print domain.statistics()


#-------------------------
# Initial Conditions
#-------------------------
domain.set_quantity('friction', 0.0)
domain.set_quantity('stage', 0.0)

import time
t0 = time.time()
bathymetry_filename=project.bathymetry_filename
print 'Starting domain.set_quantity.  Loading ', bathymetry_filename
domain.set_quantity('elevation',
                    filename=bathymetry_filename,
                    alpha=0.02,                    
                    verbose=True,
                    use_cache=False)

print 'Set_quantity elevation took %.2f seconds' %(time.time()-t0)
import sys; sys.exit() 

#-------------------------
# Distribute domain if run in parallel
#-------------------------
if numprocs > 1:
    domain = distribute(domain)


#-------------------------
# Set simulation parameters
#-------------------------
domain.set_name(project.output_filename)  # Name of output sww file
domain.set_default_order(2)               # Apply second order scheme 
domain.set_all_limiters(0.9)              # Max second order scheme (old lim)
domain.set_minimum_storable_height(0.001) # Don't store w < 0.001m
domain.set_maximum_allowed_speed(0.1)     # Allow a little runoff (0.1 is OK)


#-------------------------
# Boundary Conditions
#-------------------------

# Create boundary function from timeseries provided in file
function = file_function(project.boundary_filename,
                         domain, verbose=True)

# Create and assign boundary objects
Bts = Transmissive_Momentum_Set_Stage_boundary(domain, function)
Br = Reflective_boundary(domain)
domain.set_boundary({'wave': Bts, 'wall': Br})


#-------------------------
# Evolve through time
#-------------------------
import time
t0 = time.time()

for t in domain.evolve(yieldstep = 0.05, finaltime = 22.5):
    domain.write_time()

print 'That took %.2f seconds' %(time.time()-t0)
