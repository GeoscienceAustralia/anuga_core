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
import anuga

import project

# configure my logging
import anuga.utilities.log as log
#log.console_logging_level = log.INFO
#log.log_logging_level = log.DEBUG
log.log_filename = './run_okushiri.log'



def main(elevation_in_mesh=False):
    #-------------------------
    # Create Domain from mesh
    #-------------------------
    domain = anuga.Domain(project.mesh_filename, use_cache=False, verbose=True)
    print domain.statistics()


    #-------------------------
    # Initial Conditions
    #-------------------------
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', 0.0)
    if elevation_in_mesh is False:
        domain.set_quantity('elevation',
                            filename=project.bathymetry_filename,
                            alpha=0.02,                    
                            verbose=True,
                            use_cache=False)


    #-------------------------
    # Set simulation parameters
    #-------------------------
    domain.set_name(project.output_filename)  # Name of output sww file
    domain.set_default_order(2)               # Apply second order scheme 
    domain.set_minimum_storable_height(0.001) # Don't store w < 0.001m
    domain.set_quantities_to_be_monitored('stage')


    #------------------------------------------------------------------------------
    # Setup Algorithm, either using command line arguments
    # or override manually yourself
    #------------------------------------------------------------------------------
    from anuga.utilities.argparsing import parse_standard_args
    alg, cfl = parse_standard_args()
    domain.set_flow_algorithm(alg)
    domain.set_CFL(cfl)

    #-------------------------
    # Boundary Conditions
    #-------------------------

    # Create boundary function from timeseries provided in file
    function = anuga.file_function(project.boundary_filename,
                             domain, verbose=True)

    # Create and assign boundary objects
    Bts = anuga.Transmissive_momentum_set_stage_boundary(domain, function)
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({'wave': Bts, 'wall': Br})


    #-------------------------
    # Evolve through time
    #-------------------------
    import time
    t0 = time.time()

    for t in domain.evolve(yieldstep = 0.05, finaltime = 22.5):
        domain.write_time()
        print domain.quantity_statistics(precision='%.12f')

    print 'That took %.2f seconds' %(time.time()-t0)

#-------------------------------------------------------------
if __name__ == "__main__":
    main()
