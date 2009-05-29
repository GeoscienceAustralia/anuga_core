"""This is the public API to ANUGA.

Ideally, all tools needed to run simulations should be 
imported from this module
"""

# FIXME(Ole): This is one step towards the API envisioned in ticket:308


from anuga.shallow_water import Domain
from anuga.shallow_water import Dirichlet_boundary
from anuga.shallow_water import File_boundary
from anuga.shallow_water import Reflective_boundary
from anuga.shallow_water import Field_boundary
from anuga.shallow_water import Transmissive_stage_zero_momentum_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import Time_boundary
from anuga.abstract_2d_finite_volumes.util import file_function

from anuga.shallow_water.data_manager import export_grid, create_sts_boundary
from anuga.shallow_water.data_manager import csv2building_polygons

from anuga.utilities.polygon import read_polygon, plot_polygons, polygon_area
from anuga.utilities.polygon import Polygon_function


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
    




