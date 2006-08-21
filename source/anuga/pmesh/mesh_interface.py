
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Assume that the root AnuGA dir (inundation) is included in your pythonpath
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from coordinate_transforms.geo_reference import Geo_reference,DEFAULT_ZONE
from utilities.polygon import  point_in_polygon ,populate_polygon
from utilities.numerical_tools import ensure_numeric
from Numeric import Float
from utilities.polygon import inside_polygon

# This is due to pmesh being a package and a module and
# the current dir being unknown 
try:
    from pmesh.mesh import Mesh
except ImportError:  
    from mesh import Mesh

import exceptions
class PolygonError(exceptions.Exception): pass

def create_mesh_from_regions(bounding_polygon,
                             boundary_tags,
                             maximum_triangle_area=None,
                             filename=None,
                             interior_regions=None,
                             poly_geo_reference=None,
                             mesh_geo_reference=None,
                             minimum_triangle_angle=28.0,
                             verbose=True):
    """Create mesh from bounding polygons, and resolutions.

    bounding_polygon is a list of points in Eastings and Northings,
    relative to the poly_geo_reference.

    Boundary tags is a dictionary of symbolic tags. For every tag there
    is a list of indices referring to segments associated with that tag

    maximum_triangle_area is the maximal area per triangle
    for the bounding polygon, excluding the  interior regions.

    Interior_regions is a list of tuples consisting of (polygon, resolution)
    for each region to be separately refined.

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
    
    Returns the mesh instance if no finename is given

    Note, interior regions should be fully nested, as overlaps may cause
    unintended resolutions.
    
    """
    #FIXME (OLE-DSG)
    # check the segment indexes - throw an error if they are out of bounds
    #(DSG) Yes!
    
    #In addition I reckon the polygons could be of class Geospatial_data 
    #(DSG) Yes!

    # First check that interior polygons are fully contained in bounding polygon
    #Note, Both poly's have the same geo_ref, therefore don't take into account
    # geo_ref
    if interior_regions is not None:        
        # Test that all the interior polygons are inside the bounding_poly
        for interior_polygon, res in interior_regions:
            indices = inside_polygon(interior_polygon, bounding_polygon,
                                     closed = True, verbose = False)
    
            if len(indices) <> len(interior_polygon): 
                msg = 'Interior polygon %s is outside bounding polygon: %s'\
                      %(str(interior_polygon), str(bounding_polygon))
                raise PolygonError, msg



    # Resolve geo referencing        
    if mesh_geo_reference is None:
        bounding_polygon = ensure_numeric(bounding_polygon, Float)
        xllcorner = min(bounding_polygon[:,0])
        yllcorner = min(bounding_polygon[:,1])    
        #
        if poly_geo_reference is None:
            zone = DEFAULT_ZONE
        else:
            zone = poly_geo_reference.get_zone()
            [(xllcorner,yllcorner)] = poly_geo_reference.get_absolute( \
            [(xllcorner,yllcorner)])
        # create a geo_ref, based on the llc of the bounding_polygon
        mesh_geo_reference = Geo_reference(xllcorner = xllcorner,
                                           yllcorner = yllcorner,
                                           zone = zone)

    m = Mesh(geo_reference=mesh_geo_reference)

    # Do bounding polygon
    m.add_region_from_polygon(bounding_polygon,
                              tags=boundary_tags,
                              geo_reference=poly_geo_reference)

    # Find one point inside region automatically
    if interior_regions is not None:
        excluded_polygons = []        
        for polygon, res in interior_regions:
            excluded_polygons.append( polygon )
    else:
        excluded_polygons = None

    # Convert bounding poly to absolute values
    # this sort of thing can be fixed with the geo_points class
    if poly_geo_reference is not None:
        bounding_polygon_absolute = \
            poly_geo_reference.get_absolute(bounding_polygon)
    else:
        bounding_polygon_absolute = bounding_polygon
    
    inner_point = point_in_polygon(bounding_polygon_absolute)
    inner = m.add_region(inner_point[0], inner_point[1])
    inner.setMaxArea(maximum_triangle_area)

    # Do interior regions
    if interior_regions is not None:    
        for polygon, res in interior_regions:
            m.add_region_from_polygon(polygon,
                                      geo_reference=poly_geo_reference)
            # convert bounding poly to absolute values
            if poly_geo_reference is not None:
                polygon_absolute = \
                    poly_geo_reference.get_absolute(polygon)
            else:
                polygon_absolute = polygon
            inner_point = point_in_polygon(polygon_absolute)
            region = m.add_region(inner_point[0], inner_point[1])
            region.setMaxArea(res)
            

    if filename is None:
        return m
    else:
        m.generate_mesh(minimum_triangle_angle=minimum_triangle_angle,
                             verbose=verbose)
        m.export_mesh_file(filename)


