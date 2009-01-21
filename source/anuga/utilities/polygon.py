#!/usr/bin/env python
"""Polygon manipulations
"""

import Numeric as num

from math import sqrt
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.geospatial_data.geospatial_data import ensure_absolute, Geospatial_data


def point_on_line(point, line, rtol=1.0e-5, atol=1.0e-8):
    """Determine whether a point is on a line segment

    Input: 
        point is given by [x, y]
	line is given by [x0, y0], [x1, y1]] or
        the equivalent 2x2 Numeric array with each row corresponding to a point.

    Output:

    Note: Line can be degenerate and function still works to discern coinciding points from non-coinciding.
    """

    point = ensure_numeric(point)
    line = ensure_numeric(line)

    res = _point_on_line(point[0], point[1],
                         line[0,0], line[0,1],
                         line[1,0], line[1,1],
                         rtol, atol)
    
    return bool(res)


######
# Result functions used in intersection() below for collinear lines.
# (p0,p1) defines line 0, (p2,p3) defines line 1.
######

# result functions for possible states
def lines_dont_coincide(p0,p1,p2,p3):               return (3, None)
def lines_0_fully_included_in_1(p0,p1,p2,p3):       return (2, num.array([p0,p1]))
def lines_1_fully_included_in_0(p0,p1,p2,p3):       return (2, num.array([p2,p3]))
def lines_overlap_same_direction(p0,p1,p2,p3):      return (2, num.array([p0,p3]))
def lines_overlap_same_direction2(p0,p1,p2,p3):     return (2, num.array([p2,p1]))
def lines_overlap_opposite_direction(p0,p1,p2,p3):  return (2, num.array([p0,p2]))
def lines_overlap_opposite_direction2(p0,p1,p2,p3): return (2, num.array([p3,p1]))

# this function called when an impossible state is found
def lines_error(p1, p2, p3, p4): raise RuntimeError, ("INTERNAL ERROR: p1=%s, p2=%s, p3=%s, p4=%s" %
                                                      (str(p1), str(p2), str(p3), str(p4)))

#                     0s1    0e1    1s0    1e0   # line 0 starts on 1, 0 ends 1, 1 starts 0, 1 ends 0
collinear_result = { (False, False, False, False): lines_dont_coincide,
                     (False, False, False, True ): lines_error,
                     (False, False, True,  False): lines_error,
                     (False, False, True,  True ): lines_1_fully_included_in_0,
                     (False, True,  False, False): lines_error,
                     (False, True,  False, True ): lines_overlap_opposite_direction2,
                     (False, True,  True,  False): lines_overlap_same_direction2,
                     (False, True,  True,  True ): lines_1_fully_included_in_0,
                     (True,  False, False, False): lines_error,
                     (True,  False, False, True ): lines_overlap_same_direction,
                     (True,  False, True,  False): lines_overlap_opposite_direction,
                     (True,  False, True,  True ): lines_1_fully_included_in_0,
                     (True,  True,  False, False): lines_0_fully_included_in_1,
                     (True,  True,  False, True ): lines_0_fully_included_in_1,
                     (True,  True,  True,  False): lines_0_fully_included_in_1,
                     (True,  True,  True,  True ): lines_0_fully_included_in_1
                   }

def intersection(line0, line1, rtol=1.0e-5, atol=1.0e-8):
    """Returns intersecting point between two line segments or None
    (if parallel or no intersection is found).

    However, if parallel lines coincide partly (i.e. shara a common segment,
    the line segment where lines coincide is returned
    

    Inputs:
        line0, line1: Each defined by two end points as in: [[x0, y0], [x1, y1]]
                      A line can also be a 2x2 numpy array with each row
                      corresponding to a point.


    Output:
        status, value

        where status and value is interpreted as follows
        
        status == 0: no intersection, value set to None.
        status == 1: intersection point found and returned in value as [x,y].
        status == 2: Collinear overlapping lines found. Value takes the form [[x0,y0], [x1,y1]].
        status == 3: Collinear non-overlapping lines. Value set to None.
        status == 4: Lines are parallel with a fixed distance apart. Value set to None.
    
    """

    # FIXME (Ole): Write this in C

    line0 = ensure_numeric(line0, num.Float)
    line1 = ensure_numeric(line1, num.Float)    

    x0 = line0[0,0]; y0 = line0[0,1]
    x1 = line0[1,0]; y1 = line0[1,1]

    x2 = line1[0,0]; y2 = line1[0,1]
    x3 = line1[1,0]; y3 = line1[1,1]

    denom = (y3-y2)*(x1-x0) - (x3-x2)*(y1-y0)
    u0 = (x3-x2)*(y0-y2) - (y3-y2)*(x0-x2)
    u1 = (x2-x0)*(y1-y0) - (y2-y0)*(x1-x0)
        
    if num.allclose(denom, 0.0, rtol=rtol, atol=atol):
        # Lines are parallel - check if they are collinear
        if num.allclose([u0, u1], 0.0, rtol=rtol, atol=atol):
            # We now know that the lines are collinear
            state_tuple = (point_on_line([x0, y0], line1, rtol=rtol, atol=atol),
                           point_on_line([x1, y1], line1, rtol=rtol, atol=atol),
                           point_on_line([x2, y2], line0, rtol=rtol, atol=atol),
                           point_on_line([x3, y3], line0, rtol=rtol, atol=atol))

            return collinear_result[state_tuple]([x0,y0],[x1,y1],[x2,y2],[x3,y3])
        else:
            # Lines are parallel but aren't collinear
            return 4, None #FIXME (Ole): Add distance here instead of None 
    else:
        # Lines are not parallel, check if they intersect
        u0 = u0/denom
        u1 = u1/denom        

        x = x0 + u0*(x1-x0)
        y = y0 + u0*(y1-y0)

        # Sanity check - can be removed to speed up if needed
        assert num.allclose(x, x2 + u1*(x3-x2), rtol=rtol, atol=atol)
        assert num.allclose(y, y2 + u1*(y3-y2), rtol=rtol, atol=atol)        

        # Check if point found lies within given line segments
        if 0.0 <= u0 <= 1.0 and 0.0 <= u1 <= 1.0: 
            # We have intersection
            return 1, num.array([x, y])
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


    line0 = ensure_numeric(line0, num.Float)
    line1 = ensure_numeric(line1, num.Float)    

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
    	points = num.reshape(points, (1,2))

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
        points = ensure_numeric(points, num.Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Points could not be converted to Numeric array'
	raise msg

    try:
        polygon = ensure_numeric(polygon, num.Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Polygon could not be converted to Numeric array'
	raise msg


    if len(points.shape) == 1:
        # Only one point was passed in. Convert to array of points
    	points = num.reshape(points, (1,2))

    indices, count = separate_points_by_polygon(points, polygon,
                                                closed=closed,
                                                verbose=verbose)

    # Return indices of points outside polygon
    if count == len(indices):
        # No points are outside
        return num.array([])
    else:
        return indices[count:][::-1]  #return reversed
       

def in_and_outside_polygon(points, polygon, closed = True, verbose = False):
    """Determine points inside and outside a polygon

       See separate_points_by_polygon for documentation

       Returns an array of points inside and an array of points outside the polygon
    """

    #if verbose: print 'Checking input to outside_polygon'
    try:
        points = ensure_numeric(points, num.Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Points could not be converted to Numeric array'
	raise msg

    try:
        polygon = ensure_numeric(polygon, num.Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Polygon could not be converted to Numeric array'
	raise msg

    if len(points.shape) == 1:
        # Only one point was passed in. Convert to array of points
    	points = num.reshape(points, (1,2))


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
        points = ensure_numeric(points, num.Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Points could not be converted to Numeric array'
	raise msg

    #if verbose: print 'Checking input to separate_points_by_polygon 2'
    try:
        polygon = ensure_numeric(polygon, num.Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Polygon could not be converted to Numeric array'
	raise msg

    msg = 'Polygon array must be a 2d array of vertices'
    assert len(polygon.shape) == 2, msg

    msg = 'Polygon array must have two columns' 
    assert polygon.shape[1] == 2, msg


    msg = 'Points array must be 1 or 2 dimensional.'
    msg += ' I got %d dimensions' %len(points.shape)
    assert 0 < len(points.shape) < 3, msg


    if len(points.shape) == 1:
        # Only one point was passed in.
        # Convert to array of points
        points = num.reshape(points, (1,2))

    
    msg = 'Point array must have two columns (x,y), '
    msg += 'I got points.shape[1] == %d' %points.shape[0]
    assert points.shape[1] == 2, msg

       
    msg = 'Points array must be a 2d array. I got %s' %str(points[:30])
    assert len(points.shape) == 2, msg

    msg = 'Points array must have two columns'
    assert points.shape[1] == 2, msg


    N = polygon.shape[0] #Number of vertices in polygon
    M = points.shape[0]  #Number of points


    indices = num.zeros( M, num.Int )

    count = _separate_points_by_polygon(points, polygon, indices,
                                        int(closed), int(verbose))

    if verbose: print 'Found %d points (out of %d) inside polygon'\
       %(count, M)
    return indices, count


def polygon_area(input_polygon):
    """ Determin area of arbitrary polygon
    Reference
    http://mathworld.wolfram.com/PolygonArea.html
    """
    
    # Move polygon to origin (0,0) to avoid rounding errors
    # This makes a copy of the polygon to avoid destroying it
    input_polygon = ensure_numeric(input_polygon)
    min_x = min(input_polygon[:,0])
    min_y = min(input_polygon[:,1])    
    polygon = input_polygon - [min_x, min_y]

    # Compute area    
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
        polygon = ensure_numeric(polygon, num.Float)
    except NameError, e:
        raise NameError, e
    except:
        msg = 'Polygon %s could not be converted to Numeric array' %(str(polygon))
        raise msg

    x = polygon[:,0]
    y = polygon[:,1]
    x = num.concatenate((x, [polygon[0,0]]), axis = 0)
    y = num.concatenate((y, [polygon[0,1]]), axis = 0)
    
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

    def __init__(self,
                 regions,
                 default=0.0,
                 geo_reference=None):

	try:
	    len(regions)
	except:
            msg = 'Polygon_function takes a list of pairs (polygon, value).'
            msg += 'Got %s' %regions
            raise msg


        T = regions[0]

        if isinstance(T, basestring):
            msg = 'You passed in a list of text values into polygon_function'
            msg += ' instead of a list of pairs (polygon, value): "%s"' %T
            
            raise Exception, msg
        
	try:
            a = len(T)
	except:
            msg = 'Polygon_function takes a list of pairs (polygon, value).'
            msg += 'Got %s' %str(T)
            raise msg

        msg = 'Each entry in regions have two components: (polygon, value).'
        msg +='I got %s' %str(T)
	assert a == 2, msg


        if geo_reference is None:
            from anuga.coordinate_transforms.geo_reference import Geo_reference
            geo_reference = Geo_reference()


        self.default = default

        # Make points in polygons relative to geo_reference
        self.regions = []
        for polygon, value in regions:
            P = geo_reference.change_points_geo_ref(polygon)
            self.regions.append((P, value))




    def __call__(self, x, y):
	x = num.array(x, num.Float)
	y = num.array(y, num.Float)

	N = len(x)
	assert len(y) == N

        points = num.concatenate((num.reshape(x, (N, 1)),
                                  num.reshape(y, (N, 1))), axis=1)

	if callable(self.default):
	    z = self.default(x,y)
	else:
	    z = num.ones(N, num.Float) * self.default

	for polygon, value in self.regions:
	    indices = inside_polygon(points, polygon)

	    # FIXME: This needs to be vectorised
	    if callable(value):
	        for i in indices:
		    xx = num.array([x[i]])
		    yy = num.array([y[i]])
                    z[i] = value(xx, yy)[0]
	    else:
	        for i in indices:
            	    z[i] = value

        if len(z) == 0:
            msg = 'Warning: points provided to Polygon function did not fall within'
            msg += 'its regions'
            msg += 'x in [%.2f, %.2f], y in [%.2f, %.2f]' % (min(x), max(x),
                                                             min(y), max(y))
            print msg

            
        return z


        
# Functions to read and write polygon information        
def read_polygon(filename, split=','):
    """Read points assumed to form a polygon.
       There must be exactly two numbers in each line separated by a comma.
       No header.
    """

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
    

def read_tagged_polygons(filename):
    """
    """
    pass
    
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
    
    
        

##
# @brief Interpolate linearly from polyline nodes to midpoints of triangles.
# @param data The data on the polyline nodes.
# @param polyline_nodes ??
# @param gauge_neighbour_id ??  FIXME(Ole): I want to get rid of this
# @param point_coordinates ??
# @param verbose True if this function is to be verbose.
def interpolate_polyline(data,
                         polyline_nodes,
                         gauge_neighbour_id,
                         interpolation_points=None,
                         rtol=1.0e-6,
                         atol=1.0e-8,
                         verbose=False):
    """Interpolate linearly between values data on polyline nodes
    of a polyline to list of interpolation points. 

    data is the data on the polyline nodes.


    Inputs:
      data: Vector or array of data at the polyline nodes.
      polyline_nodes: Location of nodes where data is available.      
      gauge_neighbour_id: ?
      interpolation_points: Interpolate polyline data to these positions.
          List of coordinate pairs [x, y] of
          data points or an nx2 Numeric array or a Geospatial_data object
      rtol, atol: Used to determine whether a point is on the polyline or not. See point_on_line.

    Output:
      Interpolated values at interpolation points
    """
    
    if isinstance(interpolation_points, Geospatial_data):
        interpolation_points = interpolation_points.get_data_points(absolute=True)

    interpolated_values = num.zeros(len(interpolation_points), num.Float)

    data = ensure_numeric(data, num.Float)
    polyline_nodes = ensure_numeric(polyline_nodes, num.Float)
    interpolation_points = ensure_numeric(interpolation_points, num.Float)
    gauge_neighbour_id = ensure_numeric(gauge_neighbour_id, num.Int)

    n = polyline_nodes.shape[0] # Number of nodes in polyline        
    # Input sanity check
    msg = 'interpolation_points are not given (interpolate.py)'
    assert interpolation_points is not None, msg
    msg = 'function value must be specified at every interpolation node'
    assert data.shape[0]==polyline_nodes.shape[0], msg
    msg = 'Must define function value at one or more nodes'
    assert data.shape[0]>0, msg

    if n == 1:
        msg = 'Polyline contained only one point. I need more. ' + str(data)
        raise Exception, msg
    elif n > 1:
        _interpolate_polyline(data,
                              polyline_nodes,
                              gauge_neighbour_id,
                              interpolation_points,                               
                              interpolated_values,
                              rtol,
                              atol)
        
    return interpolated_values

        
def _interpolate_polyline(data,
                          polyline_nodes, 
                          gauge_neighbour_id, 
                          interpolation_points, 
                          interpolated_values,
                          rtol=1.0e-6,
                          atol=1.0e-8):
    """Auxiliary function used by interpolate_polyline
    
    NOTE: OBSOLETED BY C-EXTENSION
    """
    
    number_of_nodes = len(polyline_nodes)                
    number_of_points = len(interpolation_points)
    
    for j in range(number_of_nodes):                
        neighbour_id = gauge_neighbour_id[j]
        
        # FIXME(Ole): I am convinced that gauge_neighbour_id can be discarded, but need to check with John J.
        # Keep it for now (17 Jan 2009)
        # When gone, we can simply interpolate between neighbouring nodes, i.e. neighbour_id = j+1.
        # and the test below becomes something like: if j < number_of_nodes...  
        
        if neighbour_id >= 0:
            x0, y0 = polyline_nodes[j,:]
            x1, y1 = polyline_nodes[neighbour_id,:]
            
            segment_len = sqrt((x1-x0)**2 + (y1-y0)**2)
            segment_delta = data[neighbour_id] - data[j]            
            slope = segment_delta/segment_len
            
                
            for i in range(number_of_points):                
                
                x, y = interpolation_points[i,:]
                if point_on_line([x, y], 
                                 [[x0, y0], [x1, y1]], 
                                 rtol=rtol,
                                 atol=atol):
                                 

                    dist = sqrt((x-x0)**2 + (y-y0)**2)
                    interpolated_values[i] = slope*dist + data[j]
      

    

##############################################
#Initialise module

from anuga.utilities import compile
if compile.can_use_C_extension('polygon_ext.c'):
    # Underlying C implementations can be accessed
    from polygon_ext import _point_on_line
    from polygon_ext import _separate_points_by_polygon
    from polygon_ext import _interpolate_polyline    
    #from polygon_ext import _intersection

else:
    msg = 'C implementations could not be accessed by %s.\n ' %__file__
    msg += 'Make sure compile_all.py has been run as described in '
    msg += 'the ANUGA installation guide.'
    raise Exception, msg


if __name__ == "__main__":
    pass
