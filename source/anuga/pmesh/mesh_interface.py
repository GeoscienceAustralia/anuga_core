
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Assume that the directory anuga_core/source is included in your pythonpath
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from anuga.coordinate_transforms.geo_reference import Geo_reference,DEFAULT_ZONE
from anuga.utilities.polygon import  point_in_polygon ,populate_polygon
from anuga.utilities.numerical_tools import ensure_numeric
from Numeric import Float
from anuga.utilities.polygon import inside_polygon

# This is due to pmesh being a package and a module and
# the current dir being unknown 
try:
    from anuga.pmesh.mesh import Mesh
except ImportError:  
    from mesh import Mesh

import exceptions
class PolygonError(exceptions.Exception): pass

def create_mesh_from_regions(bounding_polygon,
                             boundary_tags,
                             maximum_triangle_area=None,
                             filename=None,
                             interior_regions=None,
                             interior_holes=None,
                             poly_geo_reference=None,
                             mesh_geo_reference=None,
                             minimum_triangle_angle=28.0,
                             fail_if_polygons_outside=True,
                             use_cache=False,
                             verbose=True):
    """Create mesh from bounding polygons, and resolutions.

    bounding_polygon is a list of points in Eastings and Northings,
    relative to the poly_geo_reference.

    Boundary tags is a dictionary of symbolic tags. For every tag there
    is a list of indices referring to segments associated with that tag.
    If a segment is omitted it will be assigned the default tag ''.

    maximum_triangle_area is the maximal area per triangle
    for the bounding polygon, excluding the  interior regions.

    Interior_regions is a list of tuples consisting of (polygon, resolution)
    for each region to be separately refined. 
    
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
    
    Returns the mesh instance if no filename is given

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
              'filename': filename,
              'interior_regions': interior_regions,
              'interior_holes': interior_holes,
              'poly_geo_reference': poly_geo_reference,
              'mesh_geo_reference': mesh_geo_reference,
              'minimum_triangle_angle': minimum_triangle_angle,
              'fail_if_polygons_outside': fail_if_polygons_outside,
              'verbose': verbose}   # FIXME (Ole): Should be bypassed one day
                                    # What should be bypassed? Verbose?
    
    #print 'kwargs', kwargs

    # Call underlying engine with or without caching
    if use_cache is True:
        try:
            from anuga.caching import cache
        except:
            msg = 'Caching was requested, but caching module'+\
                  'could not be imported'
            raise msg


        m = cache(_create_mesh_from_regions,
                  args, kwargs,
                  verbose=verbose,
                  compression=False)
    else:
        m = apply(_create_mesh_from_regions,
                  args, kwargs)

    return m    
        


def _create_mesh_from_regions(bounding_polygon,
                              boundary_tags,
                              maximum_triangle_area=None,
                              filename=None,                              
                              interior_regions=None,
                              interior_holes=None,
                              poly_geo_reference=None,
                              mesh_geo_reference=None,
                              minimum_triangle_angle=28.0,
                              fail_if_polygons_outside=True,
                              verbose=True):
    """_create_mesh_from_regions - internal function.

    See create_mesh_from_regions for documentation.
    """
    #FIXME (OLE-DSG)
    # check the segment indexes - throw an error if they are out of bounds
    #(DSG) Yes!

    
    #In addition I reckon the polygons could be of class Geospatial_data 
    #(DSG) If polygons were classes caching would break in places.

    # First check that interior polygons are fully contained in bounding
    # polygon
    #Note, Both poly's have the same geo_ref, therefore don't take into account
    # geo_ref

    # Simple check
    bounding_polygon = ensure_numeric(bounding_polygon, Float)
    msg = 'Bounding polygon must be a list of points or an Nx2 array'
    assert len(bounding_polygon.shape) == 2, msg
    assert bounding_polygon.shape[1] == 2, msg    

    # 
    if interior_regions is not None:        
        # Test that all the interior polygons are inside the bounding_poly
        # and throw out those that aren't fully included.

        polygons_inside_boundary = []
        for interior_polygon, res in interior_regions:
            indices = inside_polygon(interior_polygon, bounding_polygon,
                                     closed = True, verbose = False)
    
            if len(indices) <> len(interior_polygon): 
                msg = 'Interior polygon %s is not fully inside'\
                      %(str(interior_polygon))
                msg += ' bounding polygon: %s.' %(str(bounding_polygon))

                if fail_if_polygons_outside is True:
                    raise PolygonError, msg                    
                else:
                    msg += ' I will ignore it.'
                    print msg

            else:
                polygons_inside_boundary.append([interior_polygon, res])
                
        # Record only those that were fully contained        
        interior_regions = polygons_inside_boundary

            
    
# the following segment of code could be used to Test that all the 
# interior polygons are inside the bounding_poly... however it might need 
# to be change a bit   
#
#count = 0
#for i in range(len(interior_regions)):
#    region = interior_regions[i]
#    interior_polygon = region[0]
#    if len(inside_polygon(interior_polygon, bounding_polygon,
#                   closed = True, verbose = False)) <> len(interior_polygon):
#        print 'WARNING: interior polygon %d is outside bounding polygon' %(i)
#        count += 1

#if count == 0:
#    print 'interior regions OK'
#else:
#    print 'check out your interior polygons'
#    print 'check %s in production directory' %figname
#    import sys; sys.exit()
    

    if interior_holes is not None:        
        # Test that all the interior polygons are inside the bounding_poly
        for interior_polygon in interior_holes:
            indices = inside_polygon(interior_polygon, bounding_polygon,
                                     closed = True, verbose = False)
    
            if len(indices) <> len(interior_polygon): 
                msg = 'Interior polygon %s is outside bounding polygon: %s'\
                      %(str(interior_polygon), str(bounding_polygon))
                raise PolygonError, msg

    # Resolve geo referencing        
    if mesh_geo_reference is None:
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
#    if interior_regions is not None:    
#        for polygon, res in interior_regions:
#            m.add_region_from_polygon(polygon,
#                                      geo_reference=poly_geo_reference)
#            # convert bounding poly to absolute values
#            if poly_geo_reference is not None:
#                polygon_absolute = \
#                    poly_geo_reference.get_absolute(polygon)
#            else:
#                polygon_absolute = polygon
#            inner_point = point_in_polygon(polygon_absolute)
#            region = m.add_region(inner_point[0], inner_point[1])
#            region.setMaxArea(res)
            
            
    if interior_regions is not None:    
        for polygon, res in interior_regions:
            m.add_region_from_polygon(polygon,
                                      max_triangle_area=res,
                                      geo_reference=poly_geo_reference)
    
            
    # Do interior holes
    if interior_holes is not None:    
        for polygon in interior_holes:
            m.add_hole_from_polygon(polygon,
                                    geo_reference=poly_geo_reference)
       



    # NOTE (Ole): This was moved here as it is annoying if mesh is always
    # stored irrespective of whether the computation
    # was cached or not. This caused Domain to
    # recompute as it has meshfile as a dependency

    # Decide whether to store this mesh or return it    
    if filename is None:
        return m
    else:
        if verbose: print 'Generating mesh to file "%s"' %filename
        m.generate_mesh(minimum_triangle_angle=minimum_triangle_angle,
                        verbose=verbose)
        m.export_mesh_file(filename)


