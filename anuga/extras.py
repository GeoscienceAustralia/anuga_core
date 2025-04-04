
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross, \
                                                        rectangular

from anuga.abstract_2d_finite_volumes.pmesh2domain import pmesh_to_domain_instance



#-----------------------------
# rectangular domains
#-----------------------------
def rectangular_cross_domain(*args, **kwargs):
    """
    Create a rectangular domain.
    
    The triangular mesh is made
    up of m by n uniform rectangular cells divided
    into 4 triangles in a cross pattern


    :param m:      number of cells in x direction
    :param n:      number of cells in y direction
    :param len1:   length of domain in x direction (left to right)
            (default 1.0)
    :param len2:   length of domain in y direction (bottom to top)
            (default 1.0)
    :param origin: tuple (x,y) specifying location of lower left corner
            of domain (default (0,0))
    :param verbose: boolean flag to output information (default False)

    
    :return: shallow water domain instance

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
def create_domain_from_file(filename, DomainClass=Domain):
    """
    Create a domain from a file
    """
    return pmesh_to_domain_instance(filename,DomainClass=DomainClass)

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
                               breaklines=None,
                               regionPtArea=None,
                               minimum_triangle_angle=28.0,
                               fail_if_polygons_outside=True,
                               use_cache=False,
                               verbose=False):
    

    """Create domain from bounding polygons and resolutions.

    :param bounding_polygon: list of points in Eastings and Northings,
        relative to the zone stated in poly_geo_reference if specified.
        Otherwise points are just x, y coordinates with no particular 
        association to any location.

    :param boundary_tags: dictionary of symbolic tags. For every tag there
        is a list of indices referring to segments associated with that tag.
        If a segment is omitted it will be assigned the default tag ''.

    :param maximum_triangle_area: maximal area per triangle
        for the bounding polygon, excluding the  interior regions.

    :param Interior_regions: list of tuples consisting of (polygon,
        resolution) for each region to be separately refined. Do not have
        polygon lines cross or be on-top of each other.  Also do not have
        polygon close to each other.
        NOTE: If a interior_region is outside the bounding_polygon it should 
        throw an error
    
    :param interior_holes: list of polygons for each hole. These polygons do not
        need to be closed, but their points must be specified in a counter-clockwise
        order.

    :param hole_tags: list of tag segment dictionaries. 
        This function does not allow segments to share points - use underlying
        pmesh functionality for that

    :param poly_geo_reference: geo_reference of the bounding polygon and
        the interior polygons. If none, assume absolute.  
        Please pass one though, since absolute references have a zone.
    
    :param mesh_geo_reference: geo_reference of the mesh to be created.
        If none is given one will be automatically generated.  It was use
        the lower left hand corner of  bounding_polygon (absolute)
        as the x and y values for the geo_ref.
    
    :param breaklines: list of polygons. These lines will be preserved by the
        triangulation algorithm - useful for coastlines, walls, etc.
        The polygons are not closed.    
               
    :param regionPtArea: list of 3-tuples specifing a point with max area for region
        containing point 

    :param fail_if_polygons_outside: If True (the default) Exception in thrown
        where interior polygons fall outside bounding polygon. If False, these
        will be ignored and execution continued.
    
    
    :return: shallow water domain instance

    .. note:: 
    
        Interior regions should be fully nested, as overlaps may cause
        unintended resolutions. 
  
    
    """


    # Build arguments and keyword arguments for use with caching or apply.
    args = (bounding_polygon,
            boundary_tags)
    
    if mesh_filename is None:
        import tempfile
        import time
        mesh_filename = 'mesh_%d.msh'%int(time.time())
    
    kwargs = {'maximum_triangle_area': maximum_triangle_area,
              'mesh_filename': mesh_filename,
              'interior_regions': interior_regions,
              'interior_holes': interior_holes,
              'hole_tags': hole_tags,
              'poly_geo_reference': poly_geo_reference,
              'mesh_geo_reference': mesh_geo_reference,
              'breaklines' : breaklines,
              'regionPtArea' : regionPtArea,
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
        domain = _create_domain_from_regions(*args, **kwargs)

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
                                breaklines=None,
                                regionPtArea=None,
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
                             breaklines=breaklines,
                             regionPtArea=regionPtArea,
                             minimum_triangle_angle=minimum_triangle_angle,
                             fail_if_polygons_outside=fail_if_polygons_outside,
                             use_cache=False,
                             verbose=verbose)

    domain = Domain(mesh_filename, use_cache=False, verbose=verbose)


    return domain

