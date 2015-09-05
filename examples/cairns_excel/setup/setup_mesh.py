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
from anuga.parallel.parallel_api import pypar_available
if pypar_available:
    from anuga import sequential_distribute_load
    from anuga import sequential_distribute_dump
    from anuga.parallel.sequential_distribute import \
        sequential_distribute_load_pickle_file

verbose = True

def build_mesh(project):
    """

    This is executed by processor 0 to build the mesh.

    """

    # Ensure mesh_breaklines include riverwalls and breaklines

    mesh_breaklines = \
        su.combine_breakLines_and_riverWalls_for_mesh(project.breaklines,
                                                      project.riverwalls)

    # Make the mesh

    anuga.create_mesh_from_regions(
        project.bounding_polygon,
        boundary_tags=project.boundary_tags,
        maximum_triangle_area=project.default_res,
        filename=project.meshname,
        interior_regions=project.interior_regions,
        use_cache=False,
        verbose=verbose,
        breaklines=mesh_breaklines,
        regionPtArea=project.region_point_areas,
    )

    # Make the domain using the mesh

    domain = anuga.create_domain_from_file(project.meshname)

    # Print some stats about mesh and domain

    print 'Number of triangles = ', len(domain)
    print 'The extent is ', domain.get_extent()
    print domain.statistics()

    # Print info on the smallest triangles

    small_areas = domain.areas.argsort()
    print ''
    print 'LOCATIONS OF TRIANGLES WITH SMALLEST AREAS'
    for i in range(10):
        j = small_areas[i]
        x = domain.centroid_coordinates[j, 0] \
            + domain.geo_reference.xllcorner
        y = domain.centroid_coordinates[j, 1] \
            + domain.geo_reference.yllcorner
        print '  Area ' + str(domain.areas[j]) + ' location: ' \
            + str(round(x, 1)) + ',' + str(round(y, 1))
    print ''

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

    if myid == 0:

        if verbose:
            print 'Hello from processor ', myid

        #
        # HERE, WE MAKE/PARTITION/READ THE MESH
        # This can lead to memory savings in parallel runs
        # (possibly because of limitations of python memory management)
        #

        # Let's see if we have already pickled this domain

        pickle_name = 'domain' + '_P%g_%g.pickle' % (1, 0)
        pickle_name = join(project.partition_dir, pickle_name)

        if os.path.exists(pickle_name):
            if verbose:
                print 'Saved domain seems to already exist'
        else:
            if verbose:
                print 'CREATING PARTITIONED DOMAIN'
            domain = build_mesh(project)

            if setup_initial_conditions is not None:
                # Set the initial conditions in serial
                setup_initial_conditions.setup_initial_conditions(domain,
                    project)

            # If pypar is not available don't use sequential_distribute stuff
            # (it will fail)

            if pypar_available:
                if verbose:
                    print 'Saving Domain'
                sequential_distribute_dump(domain, 1,
                                           partition_dir=project.partition_dir,
                                           verbose=verbose)

        # If pypar is not available don't use sequential_distribute stuff (it
        # will fail)

        if pypar_available:

            # Now partition the domain

            par_pickle_name = 'domain' + '_P%g_%g.pickle' % (numprocs,
                                                             0)
            par_pickle_name = join(project.partition_dir,
                                   par_pickle_name)
            if os.path.exists(par_pickle_name):
                if verbose:
                    print 'Saved partitioned domain seems to already exist'
            else:
                if verbose:
                    print 'Load in saved sequential pickled domain'
                domain = \
                    sequential_distribute_load_pickle_file(
                        pickle_name, np=1, verbose=verbose)

                if verbose:
                    print 'Dump partitioned domains'
                sequential_distribute_dump(
                    domain, numprocs,
                    partition_dir=project.partition_dir, verbose=verbose)

            # This can reduce the memory demands if the domain was made above
            # Even better to have an existing partition (i.e. stop here and
            # rerun)

            domain = None
            gc.collect()
    else:

        domain = None
        if verbose:
            print 'Hello from processor ', myid

    barrier()

    # print 'Distributing domain'
    # domain=distribute(domain)
    # barrier()

    # If pypar is not available don't use sequential_distribute stuff (it will
    # fail)

    if pypar_available:
        if myid == 0:
            print 'LOADING PARTITIONED DOMAIN'

        domain = \
            sequential_distribute_load(
                filename=join(project.partition_dir, 'domain'),
                verbose=verbose)

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
