#!/usr/bin/env python
"""Polygon manipulations

"""


#try:
#    from scipy import Float, Int, zeros, ones, array, concatenate, reshape, dot
#except:
#    #print 'Could not find scipy - using Numeric'

from Numeric import Float, Int, zeros, ones, array, concatenate, reshape, dot, allclose


from math import sqrt
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.geospatial_data.geospatial_data import ensure_absolute


def point_on_line(point, line, rtol=0.0, atol=0.0):
    """Determine whether a point is on a line segment

    Input: 
        point is given by [x, y]
	line is given by [x0, y0], [x1, y1]] or
        the equivalent 2x2 Numeric array with each row corresponding to a point.

    Output:

    Note: Line can be degenerate and function still works to discern coinciding points from non-coinciding.
    """

    # FIXME(Ole): Perhaps make defaults as in allclose: rtol=1.0e-5, atol=1.0e-8

    point = ensure_numeric(point)
    line = ensure_numeric(line)

    res = _point_on_line(point[0], point[1],
                         line[0,0], line[0,1],
                         line[1,0], line[1,1],
                         rtol, atol)
    
    return bool(res)





def intersection(line0, line1):
    """Returns intersecting point between two line segments or None
    (if parallel or no intersection is found).

    However, if parallel lines coincide partly (i.e. shara a common segment,
    the line segment where lines coincide is returned
    

    Inputs:
        line0, line1: Each defined by two end points as in: [[x0, y0], [x1, y1]]
                      A line can also be a 2x2 numeric array with each row
                      corresponding to a point.


    Output:
        status, value

        where status is interpreted as follows
        
        status == 0: no intersection with value set to None
        status == 1: One intersection point found and returned in value as [x,y]
        status == 2: Coinciding line segment found. Value taks the form [[x0,y0], [x1,y1]]
        status == 3: Lines would coincide but only if extended. Value set to None
        status == 4: Lines are parallel with a fixed distance apart. Value set to None.
    
    """

    # FIXME (Ole): Write this in C

    line0 = ensure_numeric(line0, Float)
    line1 = ensure_numeric(line1, Float)    

    x0 = line0[0,0]; y0 = line0[0,1]
    x1 = line0[1,0]; y1 = line0[1,1]

    x2 = line1[0,0]; y2 = line1[0,1]
    x3 = line1[1,0]; y3 = line1[1,1]

    denom = (y3-y2)*(x1-x0) - (x3-x2)*(y1-y0)
    u0 = (x3-x2)*(y0-y2) - (y3-y2)*(x0-x2)
    u1 = (x2-x0)*(y1-y0) - (y2-y0)*(x1-x0)
        
    if allclose(denom, 0.0):
        # Lines are parallel - check if they coincide on a shared a segment

        if allclose( [u0, u1], 0.0 ):
            # We now know that the lines if continued coincide
            # The remaining check will establish if the finite lines share a segment

            line0_starts_on_line1 = line0_ends_on_line1 =\
            line1_starts_on_line0 = line1_ends_on_line0 = False
                
            if point_on_line([x0, y0], line1):
                line0_starts_on_line1 = True

            if point_on_line([x1, y1], line1):
                line0_ends_on_line1 = True
 
            if point_on_line([x2, y2], line0):
                line1_starts_on_line0 = True

            if point_on_line([x3, y3], line0):
                line1_ends_on_line0 = True                               

            if not(line0_starts_on_line1 or line0_ends_on_line1\
               or line1_starts_on_line0 or line1_ends_on_line0):
                # Lines are parallel and would coincide if extended, but not as they are.
                return 3, None


            # One line fully included in the other. Use direction of included line
            if line0_starts_on_line1 and line0_ends_on_line1:
                # Shared segment is line0 fully included in line1
                segment = array([[x0, y0], [x1, y1]])                

            if line1_starts_on_line0 and line1_ends_on_line0:
                # Shared segment is line1 fully included in line0
                segment = array([[x2, y2], [x3, y3]])
            

            # Overlap with lines are oriented the same way
            if line0_starts_on_line1 and line1_ends_on_line0:
                # Shared segment from line0 start to line 1 end
                segment = array([[x0, y0], [x3, y3]])

            if line1_starts_on_line0 and line0_ends_on_line1:
                # Shared segment from line1 start to line 0 end
                segment = array([[x2, y2], [x1, y1]])                                


            # Overlap in opposite directions - use direction of line0
            if line0_starts_on_line1 and line1_starts_on_line0:
                # Shared segment from line0 start to line 1 end
                segment = array([[x0, y0], [x2, y2]])

            if line0_ends_on_line1 and line1_ends_on_line0:
                # Shared segment from line0 start to line 1 end
                segment = array([[x3, y3], [x1, y1]])                

                
            return 2, segment
        else:
            # Lines are parallel but they don't coincide
            return 4, None #FIXME (Ole): Add distance here instead of None 
            
    else:
        # Lines are not parallel or coinciding
        u0 = u0/denom
        u1 = u1/denom        

        x = x0 + u0*(x1-x0)
        y = y0 + u0*(y1-y0)

        # Sanity check - can be removed to speed up if needed
        assert allclose(x, x2 + u1*(x3-x2))
        assert allclose(y, y2 + u1*(y3-y2))        

        # Check if point found lies within given line segments
        if 0.0 <= u0 <= 1.0 and 0.0 <= u1 <= 1.0: 
            # We have intersection

            return 1, array([x, y])
        else:
            # No intersection
            return 0, None


def NEW_C_intersection(line0, line1):
    #FIXME(Ole): To write in C
    """Returns intersecting point between two line segments or None
    (if parallel or no intersection is found).

    However, if parallel lines coincide partly (i.e. shara a common segment,
    the line segment where lines coincide is returned
    

    Inputs:
        line0, line1: Each defined by two end points as in: [[x0, y0], [x1, y1]]
                      A line can also be a 2x2 numeric array with each row
                      corresponding to a point.


    Output:
        status, value

        where status is interpreted as follows
        
        status == 0: no intersection with value set to None
        status == 1: One intersection point found and returned in value as [x,y]
        status == 2: Coinciding line segment found. Value taks the form [[x0,y0], [x1,y1]]
        status == 3: Lines would coincide but only if extended. Value set to None
        status == 4: Lines are parallel with a fixed distance apart. Value set to None.
    
    """


    line0 = ensure_numeric(line0, Float)
    line1 = ensure_numeric(line1, Float)    

    status, value = _intersection(line0[0,0], line0[0,1],
                                  line0[1,0], line0[1,1],
                                  line1[0,0], line1[0,1],
                                  line1[1,0], line1[1,1])

    return status, value




def is_inside_polygon(point, polygon, closed=True, verbose=False):
    """Determine if one point is inside a polygon

    See inside_polygon for more details
    """

    indices = inside_polygon(point, polygon, closed, verbose)

    if indices.shape[0] == 1:
        return True
    elif indices.shape[0] == 0:
        return False
    else:
        msg = 'is_inside_polygon must be invoked with one point only'
        raise msg
    

def inside_polygon(points, polygon, closed=True, verbose=False):
    """Determine points inside a polygon

       Functions inside_polygon and outside_polygon have been defined in
       terms af separate_by_polygon which will put all inside indices in
       the first part of the indices array and outside indices in the last

       See separate_points_by_polygon for documentation

       points and polygon can be a geospatial instance,
       a list or a numeric array
    """

    #if verbose: print 'Checking input to inside_polygon'

    try:
        points = ensure_absolute(points)
    except NameError, e:
        raise NameError, e
    except:
        # If this fails it is going to be because the points can't be
        # converted to a numeric array.
        msg = 'Points could not be converted to Numeric array' 
	raise msg

    try:
        polygon = ensure_absolute(polygon)
    except NameError, e:
        raise NameError, e
    except:
        # If this fails it is going to be because the points can't be
        # converted to a numeric array.
        msg = 'Polygon %s could not be converted to Numeric array' %(str(polygon))
	raise msg

    if len(points.shape) == 1:
        # Only one point was passed in. Convert to array of points
    	points = reshape(points, (1,2))

    indices, count = separate_points_by_polygon(points, polygon,
                                                closed=closed,
                                                verbose=verbose)

    # Return indices of points inside polygon
    return indices[:count]



def is_outside_polygon(point, polygon, closed=True, verbose=False,
                       points_geo_ref=None, polygon_geo_ref=None):
    """Determine if one point is outside a polygon

    See outside_polygon for more details
    """

    indices = outside_polygon(point, polygon, closed, verbose)
                              #points_geo_ref, polygon_geo_ref)

    if indices.shape[0] == 1:
        return True
    elif indices.shape[0] == 0:
        return False
    else:
        msg = 'is_outside_polygon must be invoked with one point only'
        raise msg
    

def outside_polygon(points, polygon, closed = True, verbose = False):
    """Determine points outside a polygon

       Functions inside_polygon and outside_polygon have been defined in
       terms af separate_by_polygon which will put all inside indices in
       the first part of the indices array and outside indices in the last

       See separate_points_by_polygon for documentation
    """

    #if verbose: print 'Checking input to outside_polygon'
    try:
        points = ensure_numeric(points, Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Points could not be converted to Numeric array'
	raise msg

    try:
        polygon = ensure_numeric(polygon, Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Polygon could not be converted to Numeric array'
	raise msg


    if len(points.shape) == 1:
        # Only one point was passed in. Convert to array of points
    	points = reshape(points, (1,2))

    indices, count = separate_points_by_polygon(points, polygon,
                                                closed=closed,
                                                verbose=verbose)

    # Return indices of points outside polygon
    if count == len(indices):
        # No points are outside
        return array([])
    else:
        return indices[count:][::-1]  #return reversed
       

def in_and_outside_polygon(points, polygon, closed = True, verbose = False):
    """Determine points inside and outside a polygon

       See separate_points_by_polygon for documentation

       Returns an array of points inside and an array of points outside the polygon
    """

    #if verbose: print 'Checking input to outside_polygon'
    try:
        points = ensure_numeric(points, Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Points could not be converted to Numeric array'
	raise msg

    try:
        polygon = ensure_numeric(polygon, Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Polygon could not be converted to Numeric array'
	raise msg

    if len(points.shape) == 1:
        # Only one point was passed in. Convert to array of points
    	points = reshape(points, (1,2))


    indices, count = separate_points_by_polygon(points, polygon,
                                                closed=closed,
                                                verbose=verbose)
    
    # Returns indices of points inside and indices of points outside
    # the polygon

    if count == len(indices):
        # No points are outside
        return indices[:count],[]
    else:
        return  indices[:count], indices[count:][::-1]  #return reversed


def separate_points_by_polygon(points, polygon,
                               closed = True, verbose = False):
    """Determine whether points are inside or outside a polygon

    Input:
       points - Tuple of (x, y) coordinates, or list of tuples
       polygon - list of vertices of polygon
       closed - (optional) determine whether points on boundary should be
       regarded as belonging to the polygon (closed = True)
       or not (closed = False)

    Outputs:
       indices: array of same length as points with indices of points falling
       inside the polygon listed from the beginning and indices of points
       falling outside listed from the end.

       count: count of points falling inside the polygon

       The indices of points inside are obtained as indices[:count]
       The indices of points outside are obtained as indices[count:]


    Examples:
       U = [[0,0], [1,0], [1,1], [0,1]] #Unit square

       separate_points_by_polygon( [[0.5, 0.5], [1, -0.5], [0.3, 0.2]], U)
       will return the indices [0, 2, 1] and count == 2 as only the first
       and the last point are inside the unit square

    Remarks:
       The vertices may be listed clockwise or counterclockwise and
       the first point may optionally be repeated.
       Polygons do not need to be convex.
       Polygons can have holes in them and points inside a hole is
       regarded as being outside the polygon.

    Algorithm is based on work by Darel Finley,
    http://www.alienryderflex.com/polygon/

    Uses underlying C-implementation in polygon_ext.c
    """


    #if verbose: print 'Checking input to separate_points_by_polygon'


    #Input checks

    assert isinstance(closed, bool), 'Keyword argument "closed" must be boolean'
    assert isinstance(verbose, bool), 'Keyword argument "verbose" must be boolean'


    try:
        points = ensure_numeric(points, Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Points could not be converted to Numeric array'
	raise msg

    #if verbose: print 'Checking input to separate_points_by_polygon 2'
    try:
        polygon = ensure_numeric(polygon, Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Polygon could not be converted to Numeric array'
	raise msg

    #if verbose: print 'check'

    assert len(polygon.shape) == 2,\
       'Polygon array must be a 2d array of vertices'

    assert polygon.shape[1] == 2,\
       'Polygon array must have two columns'

    assert len(points.shape) == 2,\
       'Points array must be a 2d array'

    assert points.shape[1] == 2,\
       'Points array must have two columns'

    N = polygon.shape[0] #Number of vertices in polygon
    M = points.shape[0]  #Number of points


    indices = zeros( M, Int )

    count = _separate_points_by_polygon(points, polygon, indices,
                                        int(closed), int(verbose))

    if verbose: print 'Found %d points (out of %d) inside polygon'\
       %(count, M)
    return indices, count


def polygon_area(polygon):
    """ Determin area of arbitrary polygon
    Reference
    http://mathworld.wolfram.com/PolygonArea.html
    """
    
    n = len(polygon)
    poly_area = 0.0

    for i in range(n):
        pti = polygon[i]
        if i == n-1:
            pt1 = polygon[0]
        else:
            pt1 = polygon[i+1]
        xi = pti[0]
        yi1 = pt1[1]
        xi1 = pt1[0]
        yi = pti[1]
        poly_area += xi*yi1 - xi1*yi
        
    return abs(poly_area/2)

def plot_polygons(polygons_points, style=None, 
                  figname=None, label=None, verbose=False):
    
    """ Take list of polygons and plot.

    Inputs:

    polygons         - list of polygons

    style            - style list corresponding to each polygon
                     - for a polygon, use 'line'
                     - for points falling outside a polygon, use 'outside'
                        
    figname          - name to save figure to

    label            - title for plot

    Outputs:

    - list of min and max of x and y coordinates
    - plot of polygons
    """ 

    from pylab import ion, hold, plot, axis, figure, legend, savefig, xlabel, ylabel, title, close, title

    assert type(polygons_points) == list,\
               'input must be a list of polygons and/or points'
               
    ion()
    hold(True)

    minx = 1e10
    maxx = 0.0
    miny = 1e10
    maxy = 0.0

    if label is None: label = ''

    n = len(polygons_points)
    colour = []
    if style is None:
        style_type = 'line' 
        style = []
        for i in range(n):
            style.append(style_type)
            colour.append('b-')
    else:
        for s in style:
            if s == 'line': colour.append('b-')            
            if s == 'outside': colour.append('r.')
            if s <> 'line':
                if s <> 'outside':
                    colour.append('g.')
            
    for i, item in enumerate(polygons_points):
        x, y = poly_xy(item)
        if min(x) < minx: minx = min(x)
        if max(x) > maxx: maxx = max(x)
        if min(y) < miny: miny = min(y)
        if max(y) > maxy: maxy = max(y)
        plot(x,y,colour[i])
        xlabel('x')
        ylabel('y')
        title(label)

    #raw_input('wait 1')
    #FIXME(Ole): This makes for some strange scalings sometimes.
    #if minx <> 0:
    #    axis([minx*0.9,maxx*1.1,miny*0.9,maxy*1.1])
    #else:
    #    if miny == 0:
    #        axis([-maxx*.01,maxx*1.1,-maxy*0.01,maxy*1.1])
    #    else:
    #        axis([-maxx*.01,maxx*1.1,miny*0.9,maxy*1.1])

    if figname is not None:
        savefig(figname)
    else:
        savefig('test_image')

    close('all')

    vec = [minx,maxx,miny,maxy]

    return vec

def poly_xy(polygon, verbose=False):
    """ this is used within plot_polygons so need to duplicate
        the first point so can have closed polygon in plot
    """

    #if verbose: print 'Checking input to poly_xy'

    try:
        polygon = ensure_numeric(polygon, Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Polygon %s could not be converted to Numeric array' %(str(polygon))
        raise msg

    x = polygon[:,0]
    y = polygon[:,1]
    x = concatenate((x, [polygon[0,0]]), axis = 0)
    y = concatenate((y, [polygon[0,1]]), axis = 0)
    
    return x, y
    
#    x = []
#    y = []
#    n = len(poly)
#    firstpt = poly[0]
#    for i in range(n):
#        thispt = poly[i]
#        x.append(thispt[0])
#        y.append(thispt[1])

#    x.append(firstpt[0])
#    y.append(firstpt[1])
    
#    return x, y

class Polygon_function:
    """Create callable object f: x,y -> z, where a,y,z are vectors and
    where f will return different values depending on whether x,y belongs
    to specified polygons.

    To instantiate:

       Polygon_function(polygons)

    where polygons is a list of tuples of the form

      [ (P0, v0), (P1, v1), ...]

      with Pi being lists of vertices defining polygons and vi either
      constants or functions of x,y to be applied to points with the polygon.

    The function takes an optional argument, default which is the value
    (or function) to used for points not belonging to any polygon.
    For example:

       Polygon_function(polygons, default = 0.03)

    If omitted the default value will be 0.0

    Note: If two polygons overlap, the one last in the list takes precedence

    Coordinates specified in the call are assumed to be relative to the
    origin (georeference) e.g. used by domain.
    By specifying the optional argument georeference,
    all points are made relative.

    FIXME: This should really work with geo_spatial point sets.
    """

    def __init__(self, regions, default=0.0, geo_reference=None):

	try:
	    len(regions)
	except:
            msg = 'Polygon_function takes a list of pairs (polygon, value).'
            msg += 'Got %s' %polygons
            raise msg


        T = regions[0]
	try:
            a = len(T)
	except:
            msg = 'Polygon_function takes a list of pairs (polygon, value).'
            msg += 'Got %s' %polygons
            raise msg

	assert a == 2, 'Must have two component each: %s' %T


        if geo_reference is None:
            from anuga.coordinate_transforms.geo_reference import Geo_reference
            geo_reference = Geo_reference()


        self.default = default

        # Make points in polygons relative to geo_reference
        self.regions = []
        for polygon, value in regions:
            P = geo_reference.change_points_geo_ref(polygon)
            self.regions.append( (P, value) )




    def __call__(self, x, y):
	x = array(x).astype(Float)
	y = array(y).astype(Float)

	N = len(x)
	assert len(y) == N

	points = concatenate( (reshape(x, (N, 1)),
	                       reshape(y, (N, 1))), axis=1 )

	if callable(self.default):
	    z = self.default(x,y)
	else:
	    z = ones(N, Float) * self.default

	for polygon, value in self.regions:
	    indices = inside_polygon(points, polygon)

	    # FIXME: This needs to be vectorised
	    if callable(value):
	        for i in indices:
		    xx = array([x[i]])
		    yy = array([y[i]])
                    z[i] = value(xx, yy)[0]
	    else:
	        for i in indices:
            	    z[i] = value

        return z


def read_polygon(filename, split=','):
    """Read points assumed to form a polygon.
       There must be exactly two numbers in each line separated by a comma.
       No header.
    """

    #Get polygon
    fid = open(filename)
    lines = fid.readlines()
    fid.close()
    polygon = []
    for line in lines:
        fields = line.split(split)
        polygon.append( [float(fields[0]), float(fields[1])] )

    return polygon


def write_polygon(polygon, filename=None):
    """Write polygon to csv file.
       There will be exactly two numbers, easting and northing,
       in each line separated by a comma.
       
       No header.    
    """

    fid = open(filename, 'w')
    for point in polygon:
        fid.write('%f, %f\n' %point)
    fid.close()
    

def populate_polygon(polygon, number_of_points, seed=None, exclude=None):
    """Populate given polygon with uniformly distributed points.

    Input:
       polygon - list of vertices of polygon
       number_of_points - (optional) number of points
       seed - seed for random number generator (default=None)
       exclude - list of polygons (inside main polygon) from where points should be excluded

    Output:
       points - list of points inside polygon

    Examples:
       populate_polygon( [[0,0], [1,0], [1,1], [0,1]], 5 )
       will return five randomly selected points inside the unit square
    """

    from random import uniform, seed as seed_function

    seed_function(seed)

    points = []

    # Find outer extent of polygon
    max_x = min_x = polygon[0][0]
    max_y = min_y = polygon[0][1]
    for point in polygon[1:]:
        x = point[0]
        if x > max_x: max_x = x
        if x < min_x: min_x = x
        y = point[1]
        if y > max_y: max_y = y
        if y < min_y: min_y = y


    while len(points) < number_of_points:
        x = uniform(min_x, max_x)
        y = uniform(min_y, max_y)

        append = False
        if is_inside_polygon([x,y], polygon):

            append = True

            #Check exclusions
            if exclude is not None:
                for ex_poly in exclude:
                    if is_inside_polygon([x,y], ex_poly):
                        append = False


        if append is True:
            points.append([x,y])

    return points


def point_in_polygon(polygon, delta=1e-8):
    """Return a point inside a given polygon which will be close to the
    polygon edge.

    Input:
       polygon - list of vertices of polygon
       delta - the square root of 2 * delta is the maximum distance from the
       polygon points and the returned point.
    Output:
       points - a point inside polygon

       searches in all diagonals and up and down (not left and right)
    """
    import exceptions
    class Found(exceptions.Exception): pass

    point_in = False
    while not point_in:
        try:
            for poly_point in polygon: #[1:]:
                for x_mult in range (-1,2):
                    for y_mult in range (-1,2):
                        x = poly_point[0]
                        y = poly_point[1]
                        if x == 0:
                            x_delta = x_mult*delta
                        else:
                            x_delta = x+x_mult*x*delta

                        if y == 0:
                            y_delta = y_mult*delta
                        else:
                            y_delta = y+y_mult*y*delta

                        point = [x_delta, y_delta]
                        #print "point",point
                        if is_inside_polygon(point, polygon, closed=False):
                            raise Found
        except Found:
            point_in = True
        else:
            delta = delta*0.1
    return point


def number_mesh_triangles(interior_regions, bounding_poly, remainder_res):
    """Calculate the approximate number of triangles inside the
    bounding polygon and the other interior regions

    Polygon areas are converted to square Kms 

    FIXME: Add tests for this function
    """
    
    from anuga.utilities.polygon import polygon_area


    # TO DO check if any of the regions fall inside one another

    print '----------------------------------------------------------------------------'
    print 'Polygon   Max triangle area (m^2)   Total area (km^2)   Estimated #triangles'
    print '----------------------------------------------------------------------------'    
        
    no_triangles = 0.0
    area = polygon_area(bounding_poly)
    
    for poly, resolution in interior_regions:
        this_area = polygon_area(poly)
        this_triangles = this_area/resolution
        no_triangles += this_triangles
        area -= this_area
        
        print 'Interior ',
        print ('%.0f' %resolution).ljust(25),
        print ('%.2f' %(this_area/1000000)).ljust(19),
        print '%d' %(this_triangles)
        
    bound_triangles = area/remainder_res
    no_triangles += bound_triangles

    print 'Bounding ',
    print ('%.0f' %remainder_res).ljust(25),
    print ('%.2f' %(area/1000000)).ljust(19),
    print '%d' %(bound_triangles)    

    total_number_of_triangles = no_triangles/0.7

    print 'Estimated total number of triangles: %d' %total_number_of_triangles
    print 'Note: This is generally about 20% less than the final amount'    

    return int(total_number_of_triangles)


def decimate_polygon(polygon, factor=10):
    """Reduce number of points in polygon by the specified
    factor (default=10, hence the name of the function) such that
    the extrema in both axes are preserved.

    Return reduced polygon
    """

    # FIXME(Ole): This doesn't work at present,
    # but it isn't critical either

    # Find outer extent of polygon
    num_polygon = ensure_numeric(polygon)
    max_x = max(num_polygon[:,0])
    max_y = max(num_polygon[:,1])
    min_x = min(num_polygon[:,0])
    min_y = min(num_polygon[:,1])        

    # Keep only some points making sure extrema are kept
    reduced_polygon = []    
    for i, point in enumerate(polygon):
        x = point[0]
        y = point[1]        
        if x in [min_x, max_x] and y in [min_y, max_y]:
            # Keep
            reduced_polygon.append(point)
        else:
            if len(reduced_polygon)*factor < i:
                reduced_polygon.append(point)                

    return reduced_polygon

##############################################
#Initialise module

from anuga.utilities.compile import can_use_C_extension
if can_use_C_extension('polygon_ext.c'):
    # Underlying C implementations can be accessed
    from polygon_ext import _point_on_line
    from polygon_ext import _separate_points_by_polygon
    #from polygon_ext import _intersection

else:
    msg = 'C implementations could not be accessed by %s.\n ' %__file__
    msg += 'Make sure compile_all.py has been run as described in '
    msg += 'the ANUGA installation guide.'
    raise Exception, msg


if __name__ == "__main__":
    pass
