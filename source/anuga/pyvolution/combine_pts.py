"""Functionality for aggregation of points files

   Ole Nielsen, Stephen Roberts, Duncan Gray
   Geoscience Australia, 2005.   
"""

from anuga.utilities.polygon import outside_polygon, inside_polygon
from Numeric import take, concatenate
import time

import exceptions
class AttributeError(exceptions.Exception): pass

def combine_rectangular_points_files(fine_points_file,
                                     coarse_points_file,
                                     output_points_file,
                                     verbose = False):
    """Add two point files by the union of the fine and coarse points.
    
       Course points that are in the extent of the fine points are removed.

       The extent of the fine points file is assumed to be a rectangle, parallel
       to the x and y axis.
    """
    
    from load_mesh.loadASCII import import_points_file, extent, point_atts2array, export_points_file, add_point_dictionaries, take_points

    
    # load fine points file
    if verbose: print "loading Point files"
    fine = import_points_file(fine_points_file)
    
    # load course points file
    coarse = import_points_file(coarse_points_file)
        
    if verbose: print "doing other stuff"
    # convert to Numeric    
    fine = point_atts2array(fine)
    coarse = point_atts2array(coarse)
    
    #check that the attribute info is the same
    #FIXME is not sorting this a problem?
    if not fine['attributelist'].keys() == coarse['attributelist'].keys():
        raise AttributeError
    
    # set points to the same geo ref -  say the fine points geo ref
    if fine.has_key('geo_reference')and not fine['geo_reference'] is None:
        if coarse.has_key('geo_reference'):
            coarse['pointlist'] = \
                fine['geo_reference'].change_points_geo_ref(coarse['pointlist'],
                                                            points_geo_ref=coarse['geo_reference'])  
    # find extent of course points
    extent = extent(fine['pointlist'])    

    # find out_points - coarse points outside of fine points extent

    if verbose:
        print "starting outside_polygon"
        t0 = time.time()
    outside_coarse_indices = outside_polygon(coarse['pointlist'],
                                             extent, closed=True)
    if verbose:
        print "Points outside determined"
        print 'That took %.2f seconds' %(time.time()-t0)

    # Remove points from the coarse data
    coarse = take_points(coarse, outside_coarse_indices)
        
    # add fine points and out_points
    if verbose: print "Adding points"
    combined = add_point_dictionaries(fine, coarse)

    # save
    if verbose: print "writing points"
    export_points_file(output_points_file, combined)

def add_points_files(fine_points_file,
                     coarse_points_file,
                     output_points_file,
                     verbose = False):
    """
    """
    
    from load_mesh.loadASCII import import_points_file, extent, point_atts2array, export_points_file, add_point_dictionaries,take_points

    
    # load fine points file
    if verbose:print "loading Point files"
    try:
        fine = import_points_file(fine_points_file,delimiter = ',')
    except SyntaxError,e:
        fine = import_points_file(fine_points_file,delimiter = ' ')
    
    
    # load course points file
    try:
        coarse = import_points_file(coarse_points_file,
                                              delimiter = ',')
    except SyntaxError,e:
        coarse = import_points_file(coarse_points_file,
                                              delimiter = ' ')
        
    if verbose:print "doing other stuff"
    # convert to Numeric    
    fine = point_atts2array(fine)
    coarse = point_atts2array(coarse)
    
    #check that the attribute info is the same
    #FIXME is not sorting this a problem?
    if not fine['attributelist'].keys() == coarse['attributelist'].keys():
        raise AttributeError
    
    # set points to the same geo ref -  say the fine points geo ref
    if fine.has_key('geo_reference')and not fine['geo_reference'] is None:
        if coarse.has_key('geo_reference'):
            coarse['pointlist'] = \
             fine['geo_reference'].change_points_geo_ref(coarse['pointlist'],
                                                         points_geo_ref=coarse['geo_reference'])  
        
    # add fine points and out_points
    if verbose: print "Adding points"
    combined = add_point_dictionaries(fine, coarse)

    # save
    if verbose: print "writing points"
    export_points_file(output_points_file, combined)

def reduce_points_to_mesh_extent(points_file,
                                 mesh_file,
                                 output_points_file,
                                 verbose = False):
    """
    """    
    
    from load_mesh.loadASCII import import_points_file, extent, point_atts2array, export_points_file, add_point_dictionaries,take_points, import_mesh_file

    
    # load  points file
    try:
        points = import_points_file(points_file,delimiter = ',')
    except SyntaxError,e:
        points = import_points_file(points_file,delimiter = ' ')
    
    # load  mesh file
    mesh = import_mesh_file(mesh_file)

    # get extent of mesh
    extent = extent(mesh['vertices'])
    #print "extent",extent 
    # set points to the same geo ref - the points geo ref
    if points.has_key('geo_reference')and not points['geo_reference'] is None:
        if mesh.has_key('geo_reference'):
            extent = \
             points['geo_reference'].change_points_geo_ref(extent,
                                                           points_geo_ref=mesh['geo_reference'])  
    #print "extent",extent 

    inside_indices = inside_polygon(points['pointlist'],
				    extent, closed=True)
    #create a list of in points
    points = take_points(points, inside_indices)

    # save the list, along with geo-ref
    export_points_file(output_points_file,points)
    
    #-------------------------------------------------------------
if __name__ == "__main__":
    """
    Load in two data points files.
    Save a new points file.
    """
    import os, sys
    usage = "usage: %s fine_input.pts coarse_input.pts output.pts " \
            %os.path.basename(sys.argv[0])

    if len(sys.argv) < 4:
        print usage
    else:
        fine_file  = sys.argv[1]
        coarse_file = sys.argv[2]
        output_file = sys.argv[3]
        
        combine_rectangular_points_files(fine_file,
                                         coarse_file,
                                         output_file)
