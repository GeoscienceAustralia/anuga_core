#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
anuga_run_toml — run an ANUGA scenario from a TOML configuration file.

Usage (serial):
    anuga_run_toml  path/to/scenario.toml

Usage (parallel, N processes):
    mpirun -np N anuga_run_toml  path/to/scenario.toml

The working directory is changed to the directory that contains the TOML
file before the simulation starts, so all relative paths inside the TOML
are resolved relative to that directory.

A 'user_functions.py' module may be placed alongside the TOML file to
provide custom callbacks (print_velocity_statistics, print_operator_inputs).
If the file is absent those hooks are silently skipped.
"""

import argparse
import os
import sys
import time

parser = argparse.ArgumentParser(
    description='Run an ANUGA scenario defined by a TOML configuration file.')
parser.add_argument(
    'config',
    metavar='CONFIG.toml',
    help='Path to the TOML scenario configuration file.')
args = parser.parse_args()

config_path = os.path.abspath(args.config)
if not os.path.exists(config_path):
    sys.exit(f'ERROR: config file not found: {config_path}')

# Change to the scenario directory so all relative paths in the TOML resolve
# correctly.  Do this before any anuga imports that may write files.
scenario_dir = os.path.dirname(config_path)
os.chdir(scenario_dir)
config_basename = os.path.basename(config_path)

# ---------------------------------------------------------------------------
# ANUGA imports (after chdir so parallel init finds the right cwd)
# ---------------------------------------------------------------------------

import anuga
from anuga.parallel import myid, numprocs, finalize, barrier
from anuga.operators.collect_max_quantities_operator import \
    Collect_max_quantities_operator

from anuga.scenario import (
    setup_boundary_conditions,
    setup_rainfall,
    setup_inlets,
    setup_bridges,
    setup_pumping_stations,
    setup_mesh,
    setup_initial_conditions,
    setup_riverwalls,
    raster_outputs,
)
from anuga.scenario.prepare_data import PrepareData

# ---------------------------------------------------------------------------
# Optional user_functions module (lives alongside the TOML file)
# ---------------------------------------------------------------------------

sys.path.insert(0, scenario_dir)
try:
    import user_functions
    _have_user_functions = True
except ImportError:
    _have_user_functions = False

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

t0 = time.time()


def progress(msg):
    """Print a setup milestone to the terminal regardless of log redirection."""
    sys.__stdout__.write(msg + '\n')
    sys.__stdout__.flush()


# ---------------------------------------------------------------------------
# Load configuration
# ---------------------------------------------------------------------------

project = PrepareData(config_basename, output_log='Simulation_logfile.log')

# ---------------------------------------------------------------------------
# Build mesh and set initial conditions
# ---------------------------------------------------------------------------

progress('Building mesh')
domain = setup_mesh.setup_mesh(project)

progress('Setting initial conditions')
setup_initial_conditions.setup_initial_conditions(domain, project)

# Riverwalls must be added AFTER any distribute step
progress('Adding riverwalls')
setup_riverwalls.setup_riverwalls(domain, project)

# ---------------------------------------------------------------------------
# Forcing terms
# ---------------------------------------------------------------------------

progress('Making rainfall')
setup_rainfall.setup_rainfall(domain, project)

progress('Making inlets')
setup_inlets.setup_inlets(domain, project)

progress('Making bridges')
setup_bridges.setup_bridges(domain, project)

progress('Making pumping stations')
setup_pumping_stations.setup_pumping_stations(domain, project)

# ---------------------------------------------------------------------------
# Boundary conditions
# ---------------------------------------------------------------------------

progress('Making boundary conditions')
setup_boundary_conditions.setup_boundary_conditions(domain, project)

# ---------------------------------------------------------------------------
# Track maximum quantities
# ---------------------------------------------------------------------------

max_quantities = Collect_max_quantities_operator(
    domain,
    update_frequency=project.max_quantity_update_frequency,
    collection_start_time=project.max_quantity_collection_start_time,
    velocity_zero_height=1.0e-03)

# ---------------------------------------------------------------------------
# Evolve
# ---------------------------------------------------------------------------

domain.set_multiprocessor_mode(project.multiprocessor_mode)
domain.set_omp_num_threads(project.omp_num_threads)

progress('Evolving')

barrier()
for t in domain.evolve(yieldstep=project.yieldstep,
                       finaltime=project.finaltime,
                       outputstep=project.outputstep):
    if myid == 0:
        domain.print_timestepping_statistics()

    if project.report_mass_conservation_statistics:
        domain.report_water_volume_statistics()

    if project.report_peak_velocity_statistics and _have_user_functions:
        user_functions.print_velocity_statistics(domain, max_quantities)

    if project.report_smallest_edge_timestep_statistics:
        domain.report_cells_with_small_local_timestep()

    if project.report_operator_statistics and _have_user_functions:
        user_functions.print_operator_inputs(domain)

barrier()

# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

max_quantity_file_start = domain.get_datadir() + '/Max_quantities_'
max_quantities.export_max_quantities_to_csv(max_quantity_file_start)

os.chdir(project.output_dir)
if myid == 0 and numprocs > 1:
    print('Number of processors %g ' % numprocs)
    print('That took %.2f seconds' % (time.time() - t0))
    print('Communication time %.2f seconds' % domain.communication_time)
    print('Reduction Communication time %.2f seconds'
          % domain.communication_reduce_time)
    print('Broadcast time %.2f seconds'
          % domain.communication_broadcast_time)

    anuga.utilities.sww_merge.sww_merge_parallel(
        project.scenario,
        np=numprocs, verbose=True, delete_old=True)

if myid == 0:
    try:
        raster_outputs.make_me_some_tifs(
            sww_file='./' + project.scenario + '.sww',
            bounding_polygon=project.bounding_polygon,
            proj4string=project.proj4string,
            cell_size=project.output_tif_cellsize)
    except Exception as e:
        print('GeoTif creation failed: ' + str(e))
        print('You can try manually using raster_outputs.py or '
              'anuga.utilities.plot_utils.Make_Geotif')

barrier()
finalize()
