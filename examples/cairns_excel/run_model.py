#!/usr/bin/python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# Import necessary modules
# ------------------------------------------------------------------------------
# Standard modules

"""

Script for running a semi-generic anuga model

It's supposed to be organised in a fairly conceptual manner -- a collection of
clearly named single-line function calls which do the main part of the work. To
customize it, most edits will be required in the associated setup_* scripts

If you want to add/adjust things , consider putting detailed code as functions
in user_functions.py , and calling them from here -- to keep it clean

Gareth Davies
Geoscience Australia, 2014-present

"""


import os
import time
import sys

# ANUGA modules

import anuga
from anuga.parallel import myid, numprocs, finalize, barrier
from anuga.operators.collect_max_quantities_operator import \
    collect_max_quantities_operator

# Routines to setup the domain

from setup import setup_boundary_conditions
from setup import setup_rainfall
from setup import setup_inlets
from setup import setup_bridges
from setup import setup_pumping_stations
from setup import setup_mesh
from setup import setup_initial_conditions
from setup import setup_riverwalls
from setup import raster_outputs
from setup.prepare_data import PrepareData

# Routines defined by the user
import user_functions

# Record the time so we can report how long the simulation takes
t0 = time.time()

# Get key data for the simulation
if len(sys.argv) > 1:
    # Config file
    input_file = sys.argv[1]
else:
    input_file = 'ANUGA_setup.xls'

project = PrepareData(input_file, output_log='Simulation_logfile.log')

###########################################################################
#
# SETUP DOMAIN AND INITIAL CONDITIONS 
#
###########################################################################

set_initial_conditions_in_parallel = False

if not set_initial_conditions_in_parallel:
    print 'Making domain and initial conditions in serial'
    domain = setup_mesh.setup_mesh(project, 
        setup_initial_conditions=setup_initial_conditions)
else:
    print 'Making domain in serial'
    domain = setup_mesh.setup_mesh(project)
    print 'Making initial conditions in parallel'
    setup_initial_conditions.setup_initial_conditions(domain, project)
    

# Riverwalls must be added AFTER any distribute step
print 'Adding riverwalls'
setup_riverwalls.setup_riverwalls(domain, project)

##########################################################################
#
# SETUP FORCING TERMS
#
##########################################################################

print 'Making rainfall '
setup_rainfall.setup_rainfall(domain, project)

print 'Making inlets '
setup_inlets.setup_inlets(domain, project)

print 'Making bridges '
setup_bridges.setup_bridges(domain, project)

print 'Making pumping stations '
setup_pumping_stations.setup_pumping_stations(domain, project)

##########################################################################
#
# SETUP BOUNDARY CONDITIONS
#
##########################################################################

print 'Making boundary conditions '
setup_boundary_conditions.setup_boundary_conditions(domain, project)

##########################################################################
#
# STORE MAX QUANTITIES
#
##########################################################################

max_quantities = collect_max_quantities_operator(
    domain,
    update_frequency=project.max_quantity_update_frequency,
    collection_start_time=project.max_quantity_collection_start_time,
    # Set stored velocities to zero if height<velocity_zero_height
    velocity_zero_height=1.0e-03)

##########################################################################
#
# EVOLVE IN TIME
#
##########################################################################

print 'Evolving'

barrier()
for t in domain.evolve(yieldstep=project.yieldstep,
                       finaltime=project.finaltime):
    if myid == 0:
        domain.print_timestepping_statistics()

    if project.report_mass_conservation_statistics:
        domain.report_water_volume_statistics()

    if project.report_peak_velocity_statistics:
        user_functions.print_velocity_statistics(domain, max_quantities)

    if project.report_smallest_edge_timestep_statistics:
        domain.report_cells_with_small_local_timestep()

    # Print instantaneous operator info
    user_functions.print_operator_inputs(domain)

barrier()

###############################################################################
#
# POST-SIMULATION PROCESSING
#
###############################################################################

# Write max quantities to a file

max_quantity_file_start = domain.get_datadir() + '/Max_quantities_'
max_quantities.export_max_quantities_to_csv(max_quantity_file_start)

# Merge the parallel SWW files together

os.chdir(project.output_dir)
if myid == 0 and numprocs > 1:
    print 'Number of processors %g ' % numprocs
    print 'That took %.2f seconds' % (time.time() - t0)
    print 'Communication time %.2f seconds' % domain.communication_time
    print 'Reduction Communication time %.2f seconds' \
        % domain.communication_reduce_time
    print 'Broadcast time %.2f seconds' \
        % domain.communication_broadcast_time

    anuga.utilities.sww_merge.sww_merge_parallel(
        project.scenario,
        np=numprocs, verbose=True, delete_old=True)


# Make Geotif raster files
if myid == 0:
    try:
        raster_outputs.make_me_some_tifs(
            sww_file='./' + project.scenario + '.sww',
            bounding_polygon=project.bounding_polygon,
            proj4string=project.proj4string,
            cell_size=project.output_tif_cellsize)
    except:
        print 'GeoTif creation failed, you can try manually using' + \
              ' raster_outputs.py or anuga.utilities.plot_utils.Make_Geotif'

barrier()
finalize()
