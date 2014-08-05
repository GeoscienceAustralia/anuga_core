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
bathymetry_filename='Benchmark_2_Bathymetry_very_thin.pts'
print 'Starting domain.set_quantity.  Loading ', bathymetry_filename
s = "domain.set_quantity('elevation',filename=bathymetry_filename,alpha=0.02,verbose=True,use_cache=True)"


import profile, pstats
FN = 'profile.dat'

profile.run(s, FN)

print 'Set_quantity elevation took %.2f seconds' %(time.time()-t0)

S = pstats.Stats(FN)
#S.sort_stats('time').print_stats(20)
s = S.sort_stats('cumulative').print_stats(60)

#print s
