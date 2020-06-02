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


"""

# Module imports
import anuga

# Parallel routines
from anuga import distribute, myid, numprocs, finalize, barrier


from project import *
import create_okushiri

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

if hasattr(args, 'finaltime'):
    finalTime = args.finaltime

if hasattr(args, 'yieldstep'):
    yieldStep = args.yieldstep

if myid == 0 and verbose:
    print(finalTime)
    print(yieldStep)

if verbose: print('create mesh')
elevation_in_mesh = False
if myid == 0:
    create_okushiri.create_mesh(elevation_in_mesh=elevation_in_mesh, verbose=verbose)

barrier()



# configure my logging
import anuga.utilities.log as log
#log.console_logging_level = log.INFO
#log.log_logging_level = log.DEBUG
log.log_filename = './run_okushiri.log'
   

#-------------------------
# Create Domain from mesh
#-------------------------
if myid == 0:
    try:
        domain = anuga.Domain(mesh_filename, use_cache=False, verbose=verbose)
    except:
        msg = 'ERROR reading in mesh file. Have you run create_okushiri.py?'
        raise Exception(msg)
     
    if verbose: print(domain.statistics())


    #-------------------------
    # Initial Conditions
    #-------------------------
    domain.set_quantity('friction', 0.0025)
    domain.set_quantity('stage', 0.0)
    if verbose: print('set stage')
    if elevation_in_mesh is False:
#         domain.set_quantity('elevation',
#                         filename=bathymetry_filename_stem+'.pts', 
#                         alpha=0.02,                    
#                         verbose=verbose,
#                         use_cache=False)
        domain.set_quantity('elevation',
                        filename=bathymetry_filename_stem+'.asc',                   
                        verbose=verbose)

    #-------------------------
    # Set simulation parameters
    #-------------------------
    domain.set_name(output_filename)  # Name of output sww file 
    domain.set_minimum_storable_height(0.001) # Don't store w < 0.01m
    domain.set_store_vertices_smoothly(True)
    domain.set_flow_algorithm(alg)
else:
    
    domain = None

domain = distribute(domain)

#-------------------------
# Boundary Conditions
#-------------------------

# Create boundary function from timeseries provided in file
wave_function = anuga.file_function(boundary_filename,
                         domain, verbose=verbose)

# Create and assign boundary objects
Bts = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain, wave_function)
Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'wave': Bts, 'wall': Br})

#-------------------------------------------------------------------------
# Produce a documentation of parameters
#-------------------------------------------------------------------------
from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)


#-------------------------
# Evolve through time
#-------------------------
import time
t0 = time.time()

for t in domain.evolve(yieldstep = yieldStep, finaltime = finalTime):
    if myid == 0 and verbose: domain.write_time()

domain.sww_merge(delete_old=True)

if myid == 0 and verbose: print('That took %.2f seconds' %(time.time()-t0))

finalize()


