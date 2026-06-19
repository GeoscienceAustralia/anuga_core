#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Make the mesh (using Steve's pickling approach to distribute, which reduces memory demands)

Gareth Davies, Geoscience Australia 2014 +

"""


import glob
import os
from os.path import join
import gc

import anuga
from anuga.parallel import myid, numprocs, barrier
from anuga.utilities import spatialInputUtil as su
import anuga.utilities.log as log
from anuga.parallel.parallel_api import pypar_available
if pypar_available:
    from anuga import sequential_distribute_load
    from anuga import sequential_distribute_dump
    from anuga.parallel.sequential_distribute import \
        sequential_distribute_load_pickle_file


def build_mesh(project):
    """

    This is executed by processor 0 to build the mesh.

    """

    # Ensure mesh_breaklines include riverwalls and breaklines

    mesh_breaklines = \
        su.combine_breakLines_and_riverWalls_for_mesh(project.breaklines,
                                                      project.riverwalls)

    # Make the mesh — verbose output captured to file only
    with log.file_only():
        anuga.create_pmesh_from_regions(
            project.bounding_polygon,
            boundary_tags=project.boundary_tags,
            maximum_triangle_area=project.default_res,
            filename=project.meshname,
            interior_regions=project.interior_regions,
            use_cache=False,
            verbose=True,
            breaklines=mesh_breaklines,
            regionPtArea=project.region_point_areas,
        )

    # Make the domain using the mesh

    domain = anuga.create_domain_from_file(project.meshname)

    # Key mesh stats go to screen (info level)
    log.info('Number of triangles = %d' % len(domain))
    log.info('The extent is %s' % str(domain.get_extent()))

    # Detailed stats go to file only (verbose level)
    log.verbose(domain.statistics())

    small_areas = domain.areas.argsort()
    log.verbose('')
    log.verbose('LOCATIONS OF TRIANGLES WITH SMALLEST AREAS')
    for i in range(10):
        j = small_areas[i]
        x = domain.centroid_coordinates[j, 0] \
            + domain.geo_reference.xllcorner
        y = domain.centroid_coordinates[j, 1] \
            + domain.geo_reference.yllcorner
        log.verbose('  Area %s location: %s,%s'
                    % (domain.areas[j], round(x, 1), round(y, 1)))
    log.verbose('')

    return domain


##########################################################################

def setup_mesh(project, setup_initial_conditions=None):
    """
    Code to make the mesh (initial domain)

    The geometry is made on processor 0, then dumped and reloaded
    This reduces the memory demands

    INPUT: project == the project module

    OUTPUT: domain
    """

    # ------------------------------------------------------------------
    # Serial shortcut: build the mesh directly, no pickle/partition cycle
    # ------------------------------------------------------------------
    if numprocs == 1:
        domain = build_mesh(project)

        if setup_initial_conditions is not None:
            setup_initial_conditions.setup_initial_conditions(domain, project)

    # ------------------------------------------------------------------
    # Parallel path: partition via pickle so each rank loads only its
    # portion (reduces peak memory)
    # ------------------------------------------------------------------
    else:
        if myid == 0:

            log.verbose('Hello from processor %d' % myid)

            pickle_name = 'domain' + '_P%g_%g.pickle' % (1, 0)
            pickle_name = join(project.partition_dir, pickle_name)

            if os.path.exists(pickle_name):
                log.verbose('Saved domain seems to already exist')
            else:
                log.info('Creating partitioned domain')
                domain = build_mesh(project)

                if setup_initial_conditions is not None:
                    setup_initial_conditions.setup_initial_conditions(
                        domain, project)

                log.verbose('Saving domain')
                with log.file_only():
                    sequential_distribute_dump(domain, 1,
                                               partition_dir=project.partition_dir,
                                               verbose=True)

            par_pickle_name = 'domain' + '_P%g_%g.pickle' % (numprocs, 0)
            par_pickle_name = join(project.partition_dir, par_pickle_name)

            if os.path.exists(par_pickle_name):
                log.verbose('Saved partitioned domain seems to already exist')
            else:
                log.verbose('Load in saved sequential pickled domain')
                with log.file_only():
                    domain = sequential_distribute_load_pickle_file(
                        pickle_name, np=1, verbose=True)

                log.verbose('Dump partitioned domains')
                with log.file_only():
                    sequential_distribute_dump(
                        domain, numprocs,
                        partition_dir=project.partition_dir, verbose=True)

            domain = None
            gc.collect()

        else:
            domain = None
            log.verbose('Hello from processor %d' % myid)

        barrier()

        log.info('Loading partitioned domain')

        with log.file_only():
            domain = sequential_distribute_load(
                filename=join(project.partition_dir, 'domain'),
                verbose=True)

    # #########################################################################
    # Set output directories
    # #########################################################################

    domain.set_name(project.scenario)  # Name of sww file
    domain.set_datadir(project.output_dir)  # Store sww output here

    # Needs more changes for this to work
    # domain.set_checkpointing(checkpoint_time=project.checkpoint_time)

    # #########################################################################
    # Miscellanious numerics
    # #########################################################################

    domain.set_flow_algorithm(project.flow_algorithm)

    # Force zero beta values [hopefully integrated into source]
    # print 'Warning: Forcing everything to first order'
    # domain.beta_w=0.
    # domain.beta_uh=0.
    # domain.beta_vh=0.
    # domain.beta_w_dry=0.
    # domain.beta_uh_dry=0.
    # domain.beta_vh_dry=0.

    # Adjust velocity computation for max quantities
    # domain.velocity_protection=1.0e-05

    # Adjust CFL
    # domain.set_CFL(0.9)

    # Optionally store vertex values uniquely (large file sizes!)

    domain.set_store_vertices_uniquely(project.store_vertices_uniquely)

    if project.use_local_extrapolation_and_flux_updating:
        domain.set_local_extrapolation_and_flux_updating()

    if project.store_elevation_every_timestep:
        domain.quantities_to_be_stored['elevation'] = 2
    else:
        domain.quantities_to_be_stored['elevation'] = 1

    return domain
