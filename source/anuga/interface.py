""" This is the public API to ANUGA. It provides a toolkit of often-used
    modules, which can be used directly by importing this module.
    
    It abstracts away the internal heirarchy of the ANUGA system, allowing the
    user to concentrate on writing simulations without searching through the
    ANUGA source tree for the functions that they need.
    
    Also, it isolates the user from "under-the-hood" refactorings.
"""

# Make selected classes available directly
from anuga.shallow_water import Domain

from anuga.abstract_2d_finite_volumes.util import file_function

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

from anuga.shallow_water.data_manager import export_grid
from anuga.shallow_water.data_manager import csv2building_polygons

from anuga.file.sts import create_sts_boundary

from anuga.geometry.polygon import read_polygon, plot_polygons, polygon_area
from anuga.geometry.polygon import Polygon_function

from anuga.abstract_2d_finite_volumes.pmesh2domain import pmesh_to_domain_instance


#-----------------------------
# Standard Boundaries
#-----------------------------
from anuga.shallow_water.boundaries import File_boundary
from anuga.shallow_water.boundaries import Reflective_boundary
from anuga.shallow_water.boundaries import Field_boundary
from anuga.shallow_water.boundaries import \
                    Transmissive_stage_zero_momentum_boundary
from anuga.shallow_water.boundaries import \
                    Transmissive_momentum_set_stage_boundary
from anuga.shallow_water.boundaries import \
                    Transmissive_n_momentum_zero_t_momentum_set_stage_boundary


#-----------------------------
# SWW-specific Boundaries
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
# Forcing
#-----------------------------
from anuga.shallow_water.forcing import Inflow

#-----------------------------
# File conversion utilities
#-----------------------------
from anuga.file_conversion.file_conversion import sww2obj, dat2obj, \
                    timefile2netcdf, tsh2sww, urs2sww
from anuga.file_conversion.urs2nc import urs2nc 
from anuga.file_conversion.dem2pts import dem2pts                    
from anuga.file_conversion.esri2sww import esri2sww   
from anuga.file_conversion.sww2dem import sww2dem     
from anuga.file_conversion.asc2dem import asc2dem     
from anuga.file_conversion.ferret2sww import ferret2sww     

#-----------------------------
# SWW file access
#-----------------------------
from anuga.shallow_water.sww_interrogate import get_flow_through_cross_section

#-----------------------------
# rectangular domains
#-----------------------------
def rectangular_cross_domain(*args, **kwargs):
    points, vertices, boundary = rectangular_cross(*args, **kwargs)
    return Domain(points, vertices, boundary)

#----------------------------
# Create domain from file
#----------------------------
def create_domain_from_file(file):
    return pmesh_to_domain_instance(file,Domain)

#---------------------------
# Create domain from regions
#---------------------------

def create_domain_from_regions(bounding_polygon,
                               boundary_tags,
                               maximum_triangle_area=None,
                               mesh_filename=None,
                               interior_regions=None,
                               interior_holes=None,
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
    
    Interior_holes is a list of ploygons for each hole.

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
            raise msg


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
                                poly_geo_reference=None,
                                mesh_geo_reference=None,
                                minimum_triangle_angle=28.0,
                                fail_if_polygons_outside=True,
                                verbose=True):
    """_create_domain_from_regions - internal function.

    See create_domain_from_regions for documentation.
    """

    from anuga.shallow_water import Domain
    from anuga.pmesh.mesh_interface import create_mesh_from_regions
    
    create_mesh_from_regions(bounding_polygon,
                             boundary_tags,
                             maximum_triangle_area=maximum_triangle_area,
                             interior_regions=interior_regions,
                             filename=mesh_filename,
                             interior_holes=interior_holes,
                             poly_geo_reference=poly_geo_reference,
                             mesh_geo_reference=mesh_geo_reference,
                             minimum_triangle_angle=minimum_triangle_angle,
                             fail_if_polygons_outside=fail_if_polygons_outside,
                             use_cache=False,
                             verbose=verbose)
    domain = Domain(mesh_filename, use_cache=False, verbose=verbose)


    return domain
    




