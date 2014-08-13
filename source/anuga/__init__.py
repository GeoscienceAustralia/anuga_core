""" ANUGA models the effect of tsunamis and flooding upon a terrain mesh.
    In typical usage, a Domain class is created for a particular piece of
    terrain. Boundary conditions are specified for the domain, such as inflow
    and outflow, and then the simulation is run.

    This is the public API to ANUGA. It provides a toolkit of often-used
    modules, which can be used directly by including the following line in
    the user's code:

    import anuga
        
    This usage pattern abstracts away the internal heirarchy of the ANUGA
    system, allowing the user to concentrate on writing simulations without
    searching through the ANUGA source tree for the functions that they need.
    
    Also, it isolates the user from "under-the-hood" refactorings.
"""

pass

#Add path of package to PYTHONPATH to allow C-extensions to be loaded
import sys
sys.path += __path__


#-----------------------------------------------------
# Make selected classes available directly
#-----------------------------------------------------


from anuga.__metadata__ import __version__, __date__, __author__

#--------------------------------
# Important basic classes
#--------------------------------
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.abstract_2d_finite_volumes.region import Region
from anuga.operators.base_operator import Operator
from anuga.structures.structure_operator import Structure_operator

from anuga.abstract_2d_finite_volumes.generic_domain import Generic_Domain
from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
#------------------------------------------------------------------------------ 
# Miscellaneous
#------------------------------------------------------------------------------ 
from anuga.abstract_2d_finite_volumes.util import file_function, \
                                        sww2timeseries, sww2csv_gauges, \
                                        csv2timeseries_graphs

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross, \
                                                    rectangular

from anuga.file.csv_file import load_csv_as_building_polygons,  \
                                load_csv_as_polygons

from anuga.file.sts import create_sts_boundary

from anuga.file.ungenerate import load_ungenerate

from anuga.geometry.polygon import read_polygon
from anuga.geometry.polygon import plot_polygons
from anuga.geometry.polygon import inside_polygon
from anuga.geometry.polygon import polygon_area
from anuga.geometry.polygon_function import Polygon_function

from anuga.abstract_2d_finite_volumes.pmesh2domain import \
                                            pmesh_to_domain_instance

from anuga.utilities.system_tools import file_length
from anuga.utilities.sww_merge import sww_merge_parallel as sww_merge
from anuga.utilities.file_utils import copy_code_files
from anuga.utilities.numerical_tools import safe_acos as acos


from anuga.caching import cache
from os.path import join
from anuga.config import indent



#----------------------------
# Parallel api 
#----------------------------
from anuga_parallel.parallel_api import distribute
from anuga_parallel.parallel_api import myid, numprocs, get_processor_name
from anuga_parallel.parallel_api import send, receive
from anuga_parallel.parallel_api import pypar_available, barrier, finalize

if pypar_available:
    from anuga_parallel.parallel_api import sequential_distribute_dump
    from anuga_parallel.parallel_api import sequential_distribute_load


#-----------------------------
# Checkpointing
#-----------------------------
from anuga.shallow_water.checkpoint import load_checkpoint_file


#-----------------------------
# SwW Standard Boundaries
#-----------------------------
from anuga.shallow_water.boundaries import File_boundary
from anuga.shallow_water.boundaries import Reflective_boundary
from anuga.shallow_water.boundaries import Field_boundary
from anuga.shallow_water.boundaries import \
                    Time_stage_zero_momentum_boundary
from anuga.shallow_water.boundaries import \
                    Transmissive_stage_zero_momentum_boundary
from anuga.shallow_water.boundaries import \
                    Transmissive_momentum_set_stage_boundary
from anuga.shallow_water.boundaries import \
                    Transmissive_n_momentum_zero_t_momentum_set_stage_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import \
                    Compute_fluxes_boundary
        

#-----------------------------
# General Boundaries
#-----------------------------
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
                            import Dirichlet_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
                            import Time_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
                            import Time_space_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
                            import Transmissive_boundary



#-----------------------------
# Shallow Water Tsunamis
#-----------------------------

from anuga.shallow_water.smf import slide_tsunami, slump_tsunami


#-----------------------------
# Forcing
# These are old, should use operators
#-----------------------------
from anuga.shallow_water.forcing import Inflow, Rainfall, Wind_stress



#-----------------------------
# File conversion utilities
#-----------------------------
from anuga.file_conversion.file_conversion import sww2obj, \
                    timefile2netcdf, tsh2sww
from anuga.file_conversion.urs2nc import urs2nc
from anuga.file_conversion.urs2sww import urs2sww  
from anuga.file_conversion.urs2sts import urs2sts
from anuga.file_conversion.dem2pts import dem2pts                    
from anuga.file_conversion.esri2sww import esri2sww   
from anuga.file_conversion.sww2dem import sww2dem, sww2dem_batch 
from anuga.file_conversion.asc2dem import asc2dem     
from anuga.file_conversion.ferret2sww import ferret2sww     
from anuga.file_conversion.dem2dem import dem2dem
from anuga.file_conversion.sww2array import sww2array

#-----------------------------
# Parsing arguments
#-----------------------------
from anuga.utilities.argparsing import create_standard_parser
from anuga.utilities.argparsing import parse_standard_args


def get_args():
    """ Explicitly parse the argument list using standard anuga arguments
    
    Don't use this if you want to setup your own parser
    """
    parser = create_standard_parser()
    return parser.parse_args()


#-----------------------------
# Running Script
#-----------------------------
from anuga.utilities.run_anuga_script import run_script as run_anuga_script


#-----------------------------
# Mesh API
#-----------------------------
from anuga.pmesh.mesh_interface import create_mesh_from_regions

#-----------------------------
# SWW file access
#-----------------------------
from anuga.shallow_water.sww_interrogate import get_flow_through_cross_section
    
#---------------------------
# Operators
#---------------------------
from anuga.operators.kinematic_viscosity_operator import Kinematic_viscosity_operator

from anuga.operators.rate_operators import Rate_operator
from anuga.operators.set_friction_operators import Depth_friction_operator 

from anuga.operators.set_elevation_operator import Set_elevation_operator
from anuga.operators.set_quantity_operator import Set_quantity_operator
from anuga.operators.set_stage_operator import Set_stage_operator

from anuga.operators.set_elevation import Set_elevation
from anuga.operators.set_quantity import Set_quantity

from anuga.operators.erosion_operators import Bed_shear_erosion_operator
from anuga.operators.erosion_operators import Flat_slice_erosion_operator
from anuga.operators.erosion_operators import Flat_fill_slice_erosion_operator

#---------------------------
# Structure Operators
#---------------------------


if pypar_available:
    from anuga_parallel.parallel_operator_factory import Inlet_operator
    from anuga_parallel.parallel_operator_factory import Boyd_box_operator
    from anuga_parallel.parallel_operator_factory import Boyd_pipe_operator
    from anuga_parallel.parallel_operator_factory import Weir_orifice_trapezoid_operator
else:
    from anuga.structures.inlet_operator import Inlet_operator
    from anuga.structures.boyd_box_operator import Boyd_box_operator
    from anuga.structures.boyd_pipe_operator import Boyd_pipe_operator
    from anuga.structures.weir_orifice_trapezoid_operator import Weir_orifice_trapezoid_operator


#----------------------------
# Parallel distribute
#----------------------------


#----------------------------
# 
#Added by Petar Milevski 10/09/2013
#import time, os

from anuga.utilities.model_tools import get_polygon_from_single_file
from anuga.utilities.model_tools import get_polygons_from_Mid_Mif
from anuga.utilities.model_tools import get_polygon_list_from_files
from anuga.utilities.model_tools import get_polygon_dictionary
from anuga.utilities.model_tools import get_polygon_value_list
from anuga.utilities.model_tools import read_polygon_dir
from anuga.utilities.model_tools import read_hole_dir_multi_files_with_single_poly
from anuga.utilities.model_tools import read_multi_poly_file
from anuga.utilities.model_tools import read_hole_dir_single_file_with_multi_poly
from anuga.utilities.model_tools import read_multi_poly_file_value
from anuga.utilities.model_tools import Create_culvert_bridge_Operator


#---------------------------
# User Access Functions
#---------------------------

from anuga.utilities.system_tools import get_user_name, get_host_name, \
    get_revision_number
from anuga.utilities.mem_time_equation import estimate_time_mem




#-----------------------------
# rectangular domains
#-----------------------------
def rectangular_cross_domain(*args, **kwargs):
    """
    Create a rectangular domain with triangulation made
    up of m+1 by n+1 uniform rectangular cells divided
    into 4 triangles in a cross pattern

    Arguments
    m:      number of cells in x direction
    n:      number of cells in y direction
    len1:   length of domain in x direction (left to right)
            (default 1.0)
    len2:   length of domain in y direction (bottom to top)
            (default 1.0)
    origin: tuple (x,y) specifying location of lower left corner
            of domain (default (0,0))
    """

    try:
        verbose = kwargs.pop('verbose')
    except:
        verbose = False


    points, vertices, boundary = rectangular_cross(*args, **kwargs)
    return Domain(points, vertices, boundary, verbose= verbose)

#----------------------------
# Create domain from file
#----------------------------
def create_domain_from_file(file, DomainClass=Domain):
    """
    Create a domain from a file
    """
    return pmesh_to_domain_instance(file,DomainClass=DomainClass)

#---------------------------
# Create domain from regions
#---------------------------

def create_domain_from_regions(bounding_polygon,
                               boundary_tags,
                               maximum_triangle_area=None,
                               mesh_filename=None,
                               interior_regions=None,
                               interior_holes=None,
                               hole_tags=None,
                               poly_geo_reference=None,
                               mesh_geo_reference=None,
                               minimum_triangle_angle=28.0,
                               fail_if_polygons_outside=True,
                               use_cache=False,
                               verbose=True):
    """Create domain from bounding polygons and resolutions.

    bounding_polygon is a list of points in Eastings and Northings,
    relative to the zone stated in poly_geo_reference if specified.
    Otherwise points are just x, y coordinates with no particular 
    association to any location.

    boundary_tags is a dictionary of symbolic tags. For every tag there
    is a list of indices referring to segments associated with that tag.
    If a segment is omitted it will be assigned the default tag ''.

    maximum_triangle_area is the maximal area per triangle
    for the bounding polygon, excluding the  interior regions.

    Interior_regions is a list of tuples consisting of (polygon,
    resolution) for each region to be separately refined. Do not have
    polygon lines cross or be on-top of each other.  Also do not have
    polygon close to each other.
    
    NOTE: If a interior_region is outside the bounding_polygon it should 
    throw an error
    
    interior_holes is a list of ploygons for each hole. These polygons do not
    need to be closed, but their points must be specified in a counter-clockwise
    order.

    hole_tags  is a list of tag segment dictionaries.

    This function does not allow segments to share points - use underlying
    pmesh functionality for that

    poly_geo_reference is the geo_reference of the bounding polygon and
    the interior polygons.
    If none, assume absolute.  Please pass one though, since absolute
    references have a zone.
    
    mesh_geo_reference is the geo_reference of the mesh to be created.
    If none is given one will be automatically generated.  It was use
    the lower left hand corner of  bounding_polygon (absolute)
    as the x and y values for the geo_ref.
    
    Returns the shallow water domain instance

    Note, interior regions should be fully nested, as overlaps may cause
    unintended resolutions. 

    fail_if_polygons_outside: If True (the default) Exception in thrown
    where interior polygons fall outside bounding polygon. If False, these
    will be ignored and execution continued.
        
    
    """


    # Build arguments and keyword arguments for use with caching or apply.
    args = (bounding_polygon,
            boundary_tags)
    
    kwargs = {'maximum_triangle_area': maximum_triangle_area,
              'mesh_filename': mesh_filename,
              'interior_regions': interior_regions,
              'interior_holes': interior_holes,
              'hole_tags': hole_tags,
              'poly_geo_reference': poly_geo_reference,
              'mesh_geo_reference': mesh_geo_reference,
              'minimum_triangle_angle': minimum_triangle_angle,
              'fail_if_polygons_outside': fail_if_polygons_outside,
              'verbose': verbose} #FIXME (Ole): See ticket:14

    # Call underlying engine with or without caching
    if use_cache is True:
        try:
            from anuga.caching import cache
        except:
            msg = 'Caching was requested, but caching module'+\
                  'could not be imported'
            raise (msg)


        domain = cache(_create_domain_from_regions,
                       args, kwargs,
                       verbose=verbose,
                       compression=False)
    else:
        domain = apply(_create_domain_from_regions,
                       args, kwargs)

    return domain

        
def _create_domain_from_regions(bounding_polygon,
                                boundary_tags,
                                maximum_triangle_area=None,
                                mesh_filename=None,                           
                                interior_regions=None,
                                interior_holes=None,
                                hole_tags=None,
                                poly_geo_reference=None,
                                mesh_geo_reference=None,
                                minimum_triangle_angle=28.0,
                                fail_if_polygons_outside=True,
                                verbose=True):
    """_create_domain_from_regions - internal function.

    See create_domain_from_regions for documentation.
    """

    #from anuga.shallow_water.shallow_water_domain import Domain
    from anuga.pmesh.mesh_interface import create_mesh_from_regions
    
    create_mesh_from_regions(bounding_polygon,
                             boundary_tags,
                             maximum_triangle_area=maximum_triangle_area,
                             interior_regions=interior_regions,
                             filename=mesh_filename,
                             interior_holes=interior_holes,
                             hole_tags=hole_tags,
                             poly_geo_reference=poly_geo_reference,
                             mesh_geo_reference=mesh_geo_reference,
                             minimum_triangle_angle=minimum_triangle_angle,
                             fail_if_polygons_outside=fail_if_polygons_outside,
                             use_cache=False,
                             verbose=verbose)

    domain = Domain(mesh_filename, use_cache=False, verbose=verbose)


    return domain
    
import logging as log

from anuga.config import use_psyco, g, velocity_protection
if use_psyco:
    # try using psyco if available
    try:
        import psyco
    except:
        import os
        import sys
        if os.name == 'posix' and os.uname()[4] in ['x86_64', 'ia64']:
            pass
            # Psyco isn't supported on 64 bit systems, but it doesn't matter
        elif sys.version[:3] == '2.7' :
            pass
            # Psyco isn't available for python 2.7 (16/05/2011)
        else:
            log.critical('WARNING: psyco (speedup) could not be imported, '
                         'you may want to consider installing it')
    else:
        psyco.full() # aggressively compile everything
        #psyco.background() # attempt to profile code - only compile most used
        




