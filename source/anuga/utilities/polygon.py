#!/usr/bin/env python
"""Polygon manipulations

"""


try:
    from scipy import Float, Int, zeros, ones, array, concatenate, reshape, dot
except:
    #print 'Could not find scipy - using Numeric'
    from Numeric import Float, Int, zeros, ones, array, concatenate, reshape, dot


from math import sqrt
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.geospatial_data.geospatial_data import ensure_absolute

def point_on_line(x, y, x0, y0, x1, y1):
    """Determine whether a point is on a line segment

    Input: x, y, x0, x0, x1, y1: where
        point is given by x, y
	line is given by (x0, y0) and (x1, y1)

    """

    a = array([x - x0, y - y0])
    a_normal = array([a[1], -a[0]])

    b = array([x1 - x0, y1 - y0])

    if dot(a_normal, b) == 0:
        #Point is somewhere on the infinite extension of the line

        len_a = sqrt(sum(a**2))
        len_b = sqrt(sum(b**2))
        if dot(a, b) >= 0 and len_a <= len_b:
           return True
	else:
           return False
    else:
      return False


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

    if verbose: print 'Checking input to inside_polygon'

    try:
        points = ensure_absolute(points)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Points could not be converted to Numeric array'
	raise msg

    try:
        polygon = ensure_absolute(polygon)
    except NameError, e:
        raise NameError, e
    except:
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

    if verbose: print 'Checking input to outside_polygon'
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

    if verbose: print 'Checking input to outside_polygon'
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

    """

    #Input checks


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

    px = polygon[:,0]
    py = polygon[:,1]

    #Used for an optimisation when points are far away from polygon
    minpx = min(px); maxpx = max(px)
    minpy = min(py); maxpy = max(py)


    #Begin main loop
    indices = zeros(M, Int)

    inside_index = 0    #Keep track of points inside
    outside_index = M-1 #Keep track of points outside (starting from end)

    if verbose: print 'Separating %d points' %M        
    for k in range(M):

        if verbose:
            if k %((M+10)/10)==0: print 'Doing %d of %d' %(k, M)

        #for each point
	x = points[k, 0]
	y = points[k, 1]

        inside = False

        if not x > maxpx or x < minpx or y > maxpy or y < minpy:
            #Check polygon
            for i in range(N):
	        j = (i+1)%N

	        #Check for case where point is contained in line segment
	        if point_on_line(x, y, px[i], py[i], px[j], py[j]):
	            if closed:
	    	        inside = True
	    	    else:
	                inside = False
	    	    break
	        else:
 	            #Check if truly inside polygon
                    if py[i] < y and py[j] >= y or\
                       py[j] < y and py[i] >= y:
	    	        if px[i] + (y-py[i])/(py[j]-py[i])*(px[j]-px[i]) < x:
	        	    inside = not inside

        if inside:
	    indices[inside_index] = k
	    inside_index += 1
	else:
	    indices[outside_index] = k
	    outside_index -= 1

    return indices, inside_index


def separate_points_by_polygon_c(points, polygon,
                                 closed = True, verbose = False):
    """Determine whether points are inside or outside a polygon

    C-wrapper
    """


    if verbose: print 'Checking input to separate_points_by_polygon'


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

    from polygon_ext import separate_points_by_polygon

    if verbose: print 'Allocating array for indices'

    indices = zeros( M, Int )

    #if verbose: print 'Calling C-version of inside poly'

    count = separate_points_by_polygon(points, polygon, indices,
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

def plot_polygons(polygons, figname, verbose=False):
    
    """ Take list of polygons and plot.

    Inputs:

    polygons         - list of polygons
                        
    figname          - name to save figure to

    Outputs:

    - list of min and max of x and y coordinates
    - plot of polygons
    """ 

    from pylab import ion, hold, plot, axis, figure, legend, savefig, xlabel, ylabel, title, close

    assert type(polygons) == list,\
               'input must be a list of polygons'
               
    ion()
    hold(True)

    minx = 1e10
    maxx = 0.0
    miny = 1e10
    maxy = 0.0
    
    for polygon in polygons:
        x, y = poly_xy(polygon)  
        if min(x) < minx: minx = min(x)
        if max(x) > maxx: maxx = max(x)
        if min(y) < miny: miny = min(y)
        if max(y) > maxy: maxy = max(y)
        plot(x,y,'r-')
        xlabel('x')
        ylabel('y')

    savefig(figname)

    close('all')

    vec = [minx,maxx,miny,maxy]

    return vec

def poly_xy(polygon, verbose=False):
    """ this is used within plot_polygons so need to duplicate
        the first point so can have closed polygon in plot
    """

    if verbose: print 'Checking input to poly_xy'

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

    Coordinates specified in the call are assumed to be relative to the origin (georeference)
    e.g. used by domain. By specifying the optional argument georeference, all points are made relative.

    FIXME: This should really work with geo_spatial point sets.
    """

    def __init__(self, regions, default = 0.0, geo_reference = None):

	try:
	    len(regions)
	except:
	    msg = 'Polygon_function takes a list of pairs (polygon, value). Got %s' %polygons
	    raise msg


        T = regions[0]
	try:
            a = len(T)
	except:
	    msg = 'Polygon_function takes a list of pairs (polygon, value). Got %s' %polygons
	    raise msg

	assert a == 2, 'Must have two component each: %s' %T


        if geo_reference is None:
            from anuga.coordinate_transforms.geo_reference import Geo_reference
            geo_reference = Geo_reference()


        self.default = default

        #Make points in polygons relative to geo_reference
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

	    #FIXME: This needs to be vectorised
	    if callable(value):
	        for i in indices:
		    xx = array([x[i]])
		    yy = array([y[i]])
                    z[i] = value(xx, yy)[0]
	    else:
	        for i in indices:
            	    z[i] = value

        return z

def read_polygon(filename,split=','):
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

def populate_polygon(polygon, number_of_points, seed = None, exclude = None):
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

    #Find outer extent of polygon
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



##############################################
#Initialise module

import compile
if compile.can_use_C_extension('polygon_ext.c'):
    from polygon_ext import point_on_line
    separate_points_by_polygon = separate_points_by_polygon_c


if __name__ == "__main__":
    pass
