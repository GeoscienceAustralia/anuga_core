#!/usr/bin/env python

"""Polygon manipulations"""


from .polygon_ext import _is_inside_triangle
from .polygon_ext import _interpolate_polyline
from .polygon_ext import _line_intersect
from .polygon_ext import _polygon_overlap
from .polygon_ext import _separate_points_by_polygon
from .polygon_ext import _point_on_line

import numpy as num
import math

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.geospatial_data.geospatial_data import ensure_absolute, \
    Geospatial_data
import anuga.utilities.log as log

from .aabb import AABB


def point_on_line(point, line, rtol=1.0e-5, atol=1.0e-8):
    """Determine whether a point is on a line segment

    Input:
        point is given by [x, y]
        line is given by [x0, y0], [x1, y1]] or
        the equivalent 2x2 numeric array with each row corresponding to a point.

    Output:

    Note: Line can be degenerate and function still works to discern coinciding
          points from non-coinciding.
    """

    point = ensure_numeric(point)
    line = ensure_numeric(line)

    res = _point_on_line(point[0], point[1],
                         line[0, 0], line[0, 1],
                         line[1, 0], line[1, 1],
                         rtol, atol)

    return bool(res)


######
# Result functions used in intersection() below for collinear lines.
# (p0,p1) defines line 0, (p2,p3) defines line 1.
######

# result functions for possible states
def lines_dont_coincide(p0, p1, p2, p3):
    return (3, None)


def lines_0_fully_included_in_1(p0, p1, p2, p3):
    return (2, num.array([p0, p1]))


def lines_1_fully_included_in_0(p0, p1, p2, p3):
    return (2, num.array([p2, p3]))


def lines_overlap_same_direction(p0, p1, p2, p3):
    return (2, num.array([p0, p3]))


def lines_overlap_same_direction2(p0, p1, p2, p3):
    return (2, num.array([p2, p1]))


def lines_overlap_opposite_direction(p0, p1, p2, p3):
    return (2, num.array([p0, p2]))


def lines_overlap_opposite_direction2(p0, p1, p2, p3):
    return (2, num.array([p3, p1]))

# this function called when an impossible state is found


def lines_error(p1, p2, p3, p4):
    raise RuntimeError('INTERNAL ERROR: p1=%s, p2=%s, p3=%s, p4=%s'
                       % (str(p1), str(p2), str(p3), str(p4)))


collinear_result = {
    # line 0 starts on 1, 0 ends 1, 1 starts 0, 1 ends 0
    #       0s1    0e1    1s0    1e0
    (False, False, False, False): lines_dont_coincide,
    (False, False, False, True): lines_error,
    (False, False, True,  False): lines_error,
    (False, False, True,  True): lines_1_fully_included_in_0,
    (False, True,  False, False): lines_error,
    (False, True,  False, True): lines_overlap_opposite_direction2,
    (False, True,  True,  False): lines_overlap_same_direction2,
    (False, True,  True,  True): lines_1_fully_included_in_0,
    (True,  False, False, False): lines_error,
    (True,  False, False, True): lines_overlap_same_direction,
    (True,  False, True,  False): lines_overlap_opposite_direction,
    (True,  False, True,  True): lines_1_fully_included_in_0,
    (True,  True,  False, False): lines_0_fully_included_in_1,
    (True,  True,  False, True): lines_0_fully_included_in_1,
    (True,  True,  True,  False): lines_0_fully_included_in_1,
    (True,  True,  True,  True): lines_0_fully_included_in_1
}


def intersection(line0, line1, rtol=1.0e-5, atol=1.0e-8):
    """Returns intersecting point between two line segments.

    However, if parallel lines coincide partly (i.e. share a common segment),
    the line segment where lines coincide is returned

    Inputs:
        line0, line1: Each defined by two end points as in: [[x0, y0], [x1, y1]]
                      A line can also be a 2x2 numpy array with each row
                      corresponding to a point.

    Output:
        status, value - where status and value is interpreted as follows:
        status == 0: no intersection, value set to None.
        status == 1: intersection point found and returned in value as [x,y].
        status == 2: Collinear overlapping lines found.
                     Value takes the form [[x0,y0], [x1,y1]].
        status == 3: Collinear non-overlapping lines. Value set to None.
        status == 4: Lines are parallel. Value set to None.
    """

    # FIXME (Ole): Write this in C

    line0 = ensure_numeric(line0, float)
    line1 = ensure_numeric(line1, float)

    x0 = line0[0, 0]
    y0 = line0[0, 1]
    x1 = line0[1, 0]
    y1 = line0[1, 1]

    x2 = line1[0, 0]
    y2 = line1[0, 1]
    x3 = line1[1, 0]
    y3 = line1[1, 1]

    denom = (y3-y2)*(x1-x0) - (x3-x2)*(y1-y0)
    u0 = (x3-x2)*(y0-y2) - (y3-y2)*(x0-x2)
    u1 = (x2-x0)*(y1-y0) - (y2-y0)*(x1-x0)

    if num.allclose(denom, 0.0, rtol=rtol, atol=atol):
        # Lines are parallel - check if they are collinear
        if num.allclose([u0, u1], 0.0, rtol=rtol, atol=atol):
            # We now know that the lines are collinear
            state_tuple = (point_on_line([x0, y0], line1, rtol=rtol, atol=atol),
                           point_on_line([x1, y1], line1,
                                         rtol=rtol, atol=atol),
                           point_on_line([x2, y2], line0,
                                         rtol=rtol, atol=atol),
                           point_on_line([x3, y3], line0, rtol=rtol, atol=atol))
            # print state_tuple
            return collinear_result[state_tuple]([x0, y0], [x1, y1],
                                                 [x2, y2], [x3, y3])
        else:
            # Lines are parallel but aren't collinear
            return 4, None  # FIXME (Ole): Add distance here instead of None
    else:
        # Lines are not parallel, check if they intersect
        u0 = u0/ denom
        u1 = u1/ denom

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
    """Returns intersecting point between two line segments.

    However, if parallel lines coincide partly (i.e. share a common segment),
    the line segment where lines coincide is returned

    Inputs:
        line0, line1: Each defined by two end points as in: [[x0, y0], [x1, y1]]
                      A line can also be a 2x2 numpy array with each row
                      corresponding to a point.

    Output:
        status, value - where status and value is interpreted as follows:
        status == 0: no intersection, value set to None.
        status == 1: intersection point found and returned in value as [x,y].
        status == 2: Collinear overlapping lines found.
                     Value takes the form [[x0,y0], [x1,y1]].
        status == 3: Collinear non-overlapping lines. Value set to None.
        status == 4: Lines are parallel. Value set to None.
    """

    line0 = ensure_numeric(line0, float)
    line1 = ensure_numeric(line1, float)

    status, value = _intersection(line0[0, 0], line0[0, 1],
                                  line0[1, 0], line0[1, 1],
                                  line1[0, 0], line1[0, 1],
                                  line1[1, 0], line1[1, 1])

    return status, value


def polygon_overlap(triangles, polygon, verbose=False):
    """Determine if a polygon and triangle overlap

    """
    polygon = ensure_numeric(polygon)
    triangles = ensure_numeric(triangles)

    M = triangles.shape[0]//3 # Number of triangles

    indices = num.zeros(M, int)

    count = _polygon_overlap(polygon, triangles, indices)

    if verbose:
        log.critical('Found %d triangles (out of %d) that polygon' %
                     (count, M))

    return indices[:count]


def not_polygon_overlap(triangles, polygon, verbose=False):
    """Determine if a polygon and triangle overlap

    """
    polygon = ensure_numeric(polygon)
    triangles = ensure_numeric(triangles)

    M = triangles.shape[0]// 3  # Number of triangles

    indices = num.zeros(M, int)

    count = _polygon_overlap(polygon, triangles, indices)

    if verbose:
        log.critical('Found %d triangles (out of %d) that polygon' %
                     (count, M))

    return indices[count:]


def line_intersect(triangles, line, verbose=False):
    """Determine which of a list of triangles intersect a line

    """
    line = ensure_numeric(line)
    triangles = ensure_numeric(triangles)

    M = triangles.shape[0]// 3  # Number of triangles

    indices = num.zeros(M, int)

    count = _line_intersect(line, triangles, indices)

    if verbose:
        log.critical(
            'Found %d triangles (out of %d) that intersect line' % (count, M))

    return indices[:count]


def line_length(line):
    """Determine the length of the line
    """

    l12 = line[1]-line[0]

    return math.sqrt(num.dot(l12, l12))


def not_line_intersect(triangles, line, verbose=False):
    """Determine if a polyline and triangle overlap

    """
    line = ensure_numeric(line)
    triangles = ensure_numeric(triangles)

    M = triangles.shape[0]// 3  # Number of triangles

    indices = num.zeros(M, int)

    count = _line_intersect(line, triangles, indices)

    if verbose:
        log.critical(
            'Found %d triangles (out of %d) that intersect the line' % (count, M))

    return indices[count:]


def is_inside_triangle(point, triangle,
                       closed=True,
                       rtol=1.0e-12,
                       atol=1.0e-12,
                       check_inputs=True):
    """Determine if one point is inside a triangle

    This uses the barycentric method:

    Triangle is A, B, C
    Point P can then be written as

    P = A + alpha * (C-A) + beta * (B-A)
    or if we let 
    v=P-A, v0=C-A, v1=B-A    

    v = alpha*v0 + beta*v1 

    Dot this equation by v0 and v1 to get two:

    dot(v0, v) = alpha*dot(v0, v0) + beta*dot(v0, v1)
    dot(v1, v) = alpha*dot(v1, v0) + beta*dot(v1, v1)    

    or if a_ij = dot(v_i, v_j) and b_i = dot(v_i, v)
    the matrix equation:

    a_00 a_01   alpha     b_0
                       = 
    a_10 a_11   beta      b_1

    Solving for alpha and beta yields:

    alpha = (b_0*a_11 - b_1*a_01)/denom
    beta =  (b_1*a_00 - b_0*a_10)/denom

    with denom = a_11*a_00 - a_10*a_01

    The point is in the triangle whenever
    alpha and beta and their sums are in the unit interval.

    rtol and atol will determine how close the point has to be to the edge
    before it is deemed to be on the edge.

    """

    triangle = ensure_numeric(triangle)
    point = ensure_numeric(point, float)

    if check_inputs is True:
        msg = 'is_inside_triangle must be invoked with one point only'
        assert num.allclose(point.shape, [2]), msg

    # Use C-implementation
    return bool(_is_inside_triangle(point, triangle, int(closed), rtol, atol))


def is_complex(polygon, closed=True, verbose=False):
    """Check if a polygon is complex (self-intersecting).
       Uses a sweep algorithm that is O(n^2) in the worst case, but
       for most normal looking polygons it'll be O(n log n). 

       polygon is a list of points that define a closed polygon.
       verbose will print a list of the intersection points if true

       Return True if polygon is complex.
    """

    def key_xpos(item):
        """ Return the x coord out of the passed point for sorting key. """
        return (item[0][0])

    def segments_joined(seg0, seg1):
        """ See if there are identical segments in the 2 lists. """
        for i in seg0:
            for j in seg1:
                if i == j:
                    return True
        return False

    polygon = ensure_numeric(polygon, float)

    # build a list of discrete segments from the polygon
    unsorted_segs = []
    for i in range(0, len(polygon)-1):
        unsorted_segs.append([list(polygon[i]), list(polygon[i+1])])

    if closed:
        unsorted_segs.append([list(polygon[0]), list(polygon[-1])])

    # all segments must point in same direction
    for val in unsorted_segs:
        if val[0][0] > val[1][0]:
            val[0], val[1] = val[1], val[0]

    l_x = sorted(unsorted_segs, key=key_xpos)

    comparisons = 0

    # loop through, only comparing lines that partially overlap in x
    for index, leftmost in enumerate(l_x):
        cmp = index+1
        while cmp < len(l_x) and leftmost[1][0] > l_x[cmp][0][0]:
            if not segments_joined(leftmost, l_x[cmp]):
                (type, point) = intersection(leftmost, l_x[cmp])
                comparisons += 1
                if type != 0 and type != 4 and type != 3 or (type == 2 and list(point[0]) !=
                                                             list(point[1])):
                    if verbose:
                        print('Self-intersecting polygon found, type ', type)
                        print('point', point, end=' ')
                        print('vertices: ', leftmost, ' - ', l_x[cmp])
                    return True
            cmp += 1

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
        raise Exception(msg)


def inside_polygon(points, polygon, closed=True, verbose=False):
    """Determine points inside a polygon

       Functions inside_polygon and outside_polygon have been defined in
       terms of separate_by_polygon which will put all inside indices in
       the first part of the indices array and outside indices in the last

       See separate_points_by_polygon for documentation

       points and polygon can be a geospatial instance,
       a list or a numeric array
    """

    try:
        points = ensure_absolute(points)
    except NameError as err:
        raise NameError(err)
    except:
        # If this fails it is going to be because the points can't be
        # converted to a numeric array.
        msg = 'Points could not be converted to numeric array'
        raise Exception(msg)

    try:
        polygon = ensure_absolute(polygon)
    except NameError as e:
        raise NameError(e)
    except:
        # If this fails it is going to be because the points can't be
        # converted to a numeric array.
        msg = ('Polygon %s could not be converted to numeric array'
               % (str(polygon)))
        raise Exception(msg)

    if len(points.shape) == 1:
        # Only one point was passed in. Convert to array of points
        points = num.reshape(points, (1, 2))

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

    if indices.shape[0] == 1:
        return True
    elif indices.shape[0] == 0:
        return False
    else:
        msg = 'is_outside_polygon must be invoked with one point only'
        raise Exception(msg)


def outside_polygon(points, polygon, closed=True, verbose=False):
    """Determine points outside a polygon

       Functions inside_polygon and outside_polygon have been defined in
       terms of separate_by_polygon which will put all inside indices in
       the first part of the indices array and outside indices in the last

       See separate_points_by_polygon for documentation
    """

    try:
        points = ensure_numeric(points, float)
    except NameError as e:
        raise NameError(e)
    except:
        msg = 'Points could not be converted to numeric array'
        raise Exception(msg)

    try:
        polygon = ensure_numeric(polygon, float)
    except NameError as e:
        raise NameError(e)
    except:
        msg = 'Polygon could not be converted to numeric array'
        raise Exception(msg)

    if len(points.shape) == 1:
        # Only one point was passed in. Convert to array of points
        points = num.reshape(points, (1, 2))

    indices, count = separate_points_by_polygon(points, polygon,
                                                closed=closed,
                                                verbose=verbose)

    # Return indices of points outside polygon
    if count == len(indices):
        # No points are outside
        return num.array([])
    else:
        return indices[count:][::-1]  # return reversed


def in_and_outside_polygon(points, polygon, closed=True, verbose=False):
    """Determine points inside and outside a polygon

       See separate_points_by_polygon for documentation

       Returns an array of points inside and array of points outside the polygon
    """

    try:
        points = ensure_numeric(points, float)
    except NameError as e:
        raise NameError(e)
    except:
        msg = 'Points could not be converted to numeric array'
        raise Exception(msg)

    try:
        polygon = ensure_numeric(polygon, float)
    except NameError as e:
        raise NameError(e)
    except:
        msg = 'Polygon could not be converted to numeric array'
        raise Exception(msg)

    if len(points.shape) == 1:
        # Only one point was passed in. Convert to array of points
        points = num.reshape(points, (1, 2))

    indices, count = separate_points_by_polygon(points, polygon,
                                                closed=closed,
                                                verbose=verbose)

    # Returns indices of points inside and indices of points outside
    # the polygon
    if count == len(indices):
        # No points are outside
        return indices[:count], []
    else:
        return indices[:count], indices[count:][::-1]  # return reversed


def separate_points_by_polygon(points, polygon,
                               closed=True,
                               check_input=True,
                               verbose=False):
    """Determine whether points are inside or outside a polygon

    Input:
       points - Tuple of (x, y) coordinates, or list of tuples
       polygon - list of vertices of polygon
       closed - (optional) determine whether points on boundary should be
       regarded as belonging to the polygon (closed = True)
       or not (closed = False)
       check_input: Allows faster execution if set to False

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

    if check_input:
        # Input checks
        assert isinstance(closed, bool), \
            'Keyword argument "closed" must be boolean'
        assert isinstance(verbose, bool), \
            'Keyword argument "verbose" must be boolean'

        try:
            points = ensure_numeric(points, float)
        except NameError as e:
            raise NameError(e)
        except:
            msg = 'Points could not be converted to numeric array'
            raise Exception(msg)

        try:
            polygon = ensure_numeric(polygon, float)
        except NameError as e:
            raise NameError(e)
        except:
            msg = 'Polygon could not be converted to numeric array'
            raise Exception(msg)

        msg = 'Polygon array must be a 2d array of vertices'
        assert len(polygon.shape) == 2, msg

        msg = 'Polygon array must have two columns'
        assert polygon.shape[1] == 2, msg

        msg = ('Points array must be 1 or 2 dimensional. '
               'I got %d dimensions' % len(points.shape))
        assert 0 < len(points.shape) < 3, msg

        if len(points.shape) == 1:
            # Only one point was passed in.  Convert to array of points.
            points = num.reshape(points, (1, 2))

            msg = ('Point array must have two columns (x,y), '
                   'I got points.shape[1]=%d' % points.shape[0])
            assert points.shape[1] == 2, msg

            msg = ('Points array must be a 2d array. I got %s.'
                   % str(points[:30]))
            assert len(points.shape) == 2, msg

            msg = 'Points array must have two columns'
            assert points.shape[1] == 2, msg

    N = polygon.shape[0]  # Number of vertices in polygon
    M = points.shape[0]  # Number of points

    indices = num.zeros(M, int)

    count = _separate_points_by_polygon(points, polygon, indices,
                                        int(closed), int(verbose))

    if verbose:
        log.critical('Found %d points (out of %d) inside polygon' % (count, M))

    return indices, count


def polygon_area(input_polygon):
    """ Determine area of arbitrary polygon.

        input_polygon The polygon to get area of.

        return A scalar value for the polygon area.

        Reference:     http://mathworld.wolfram.com/PolygonArea.html
    """
    # Move polygon to origin (0,0) to avoid rounding errors
    # This makes a copy of the polygon to avoid destroying it
    input_polygon = ensure_numeric(input_polygon)
    min_x = min(input_polygon[:, 0])
    min_y = min(input_polygon[:, 1])
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

    return abs(poly_area/ 2)


def plot_polygons(polygons_points,
                  style=None,
                  figname=None,
                  label=None,
                  alpha=None):
    """ Take list of polygons and plot.

    Inputs:

    polygons         - list of polygons

    style            - style list corresponding to each polygon
                     - for a polygon, use 'line'
                     - for points falling outside a polygon, use 'outside'
                     - style can also be user defined as in normal pylab plot.

    figname          - name to save figure to

    label            - title for plotA

    alpha            - transparency of polygon fill, 0.0=none, 1.0=solid
                       if not supplied, no fill.

    Outputs:

    - plot of polygons
    """

    try:
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib.pyplot import plot, savefig, xlabel, \
            ylabel, title, close, title, fill
    except:
        return

    assert type(polygons_points) == list, \
        'input must be a list of polygons and/or points'

    if label is None:
        label = ''

    # clamp alpha to sensible range
    if alpha:
        try:
            alpha = float(alpha)
        except ValueError:
            alpha = None
        else:
            alpha = max(0.0, min(1.0, alpha))

    num_points = len(polygons_points)
    colour = []
    if style is None:
        style_type = 'line'
        style = []
        for i in range(num_points):
            style.append(style_type)
            colour.append('b-')
    else:
        for style_name in style:
            if style_name == 'line':
                colour.append('b-')
            if style_name == 'outside':
                colour.append('r.')
            if style_name == 'point':
                colour.append('g.')
            if style_name not in ['line', 'outside', 'point']:
                colour.append(style_name)

    for i, item in enumerate(polygons_points):
        pt_x, pt_y = _poly_xy(item)
        plot(pt_x, pt_y, colour[i])
        if alpha:
            fill(pt_x, pt_y, colour[i], alpha=alpha)
        xlabel('x')
        ylabel('y')
        title(label)

    if figname is not None:
        savefig(figname)
    else:
        savefig('test_image')

    close('all')


def _poly_xy(polygon):
    """ this is used within plot_polygons so need to duplicate
        the first point so can have closed polygon in plot

        polygon A set of points defining a polygon.
        verbose True if this function is to be verbose.

        Returns a tuple (x, y) of X and Y coordinates of the polygon.
        We duplicate the first point so can have closed polygon in plot.
    """

    try:
        polygon = ensure_numeric(polygon, float)
    except NameError as err:
        raise NameError(err)
    except:
        msg = ('Polygon %s could not be converted to numeric array'
               % (str(polygon)))
        raise Exception(msg)

    pts_x = num.concatenate((polygon[:, 0], [polygon[0, 0]]), axis=0)
    pts_y = num.concatenate((polygon[:, 1], [polygon[0, 1]]), axis=0)

    return pts_x, pts_y


################################################################################
# Functions to read and write polygon information
################################################################################

def read_polygon(filename, delimiter=',', closed=True, verbose=False):
    """ Read points assumed to form a (closed) polygon.
        Can also be used to read  in a polyline (closed=False)

        Also checks to make sure polygon (polyline)
        is not complex (self-intersecting).

        filename Path to file containing polygon data.
        delimiter Delimiter to split polygon data with.
        A list of point data from the polygon file.

        There must be exactly two numbers in each line
        separated by the delimiter.
        No header.
    """

    fid = open(filename)
    lines = fid.readlines()
    fid.close()
    polygon = []
    for line in lines:
        fields = line.split(delimiter)
        polygon.append([float(fields[0]), float(fields[1])])

    # check this is a valid polygon (polyline).
    if is_complex(polygon, closed, verbose=verbose):
        msg = 'ERROR: Self-intersecting polygon detected in file '
        msg += filename + '. A complex polygon will not '
        msg += 'necessarily break the algorithms within ANUGA, but it'
        msg += 'usually signifies pathological data. Please fix this file.'
        raise Exception(msg)

    return polygon


def write_polygon(polygon, filename=None):
    """Write polygon to csv file.

    There will be exactly two numbers, easting and northing, in each line
    separated by a comma.

    No header.
    """

    fid = open(filename, 'w')
    for point in polygon:
        fid.write('%f, %f\n' % point)
    fid.close()


def populate_polygon(polygon, number_of_points, seed=None, exclude=None):
    """Populate given polygon with uniformly distributed points.

    Input:
       polygon - list of vertices of polygon
       number_of_points - (optional) number of points
       seed - seed for random number generator (default=None)
       exclude - list of polygons (inside main polygon) from where points
                 should be excluded

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
    extents = AABB(polygon)

    while len(points) < number_of_points:
        rand_x = uniform(extents.xmin, extents.xmax)
        rand_y = uniform(extents.ymin, extents.ymax)

        append = False
        if is_inside_polygon([rand_x, rand_y], polygon):
            append = True

            # Check exclusions
            if exclude is not None:
                for ex_poly in exclude:
                    if is_inside_polygon([rand_x, rand_y], ex_poly):
                        append = False

        if append is True:
            points.append([rand_x, rand_y])

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

       searches in all diagonals and up and down (not left and right).
    """

    polygon = ensure_numeric(polygon)

    while True:
        for poly_point in polygon:
            for x_mult in range(-1, 2):
                for y_mult in range(-1, 2):
                    pt_x, pt_y = poly_point

                    if pt_x == 0:
                        x_delta = x_mult * delta
                    else:
                        x_delta = pt_x + x_mult*pt_x*delta

                    if pt_y == 0:
                        y_delta = y_mult * delta
                    else:
                        y_delta = pt_y + y_mult*pt_y*delta

                    point = [x_delta, y_delta]

                    if is_inside_polygon(point, polygon, closed=False):
                        return point
        delta = delta * 0.1


def number_mesh_triangles(interior_regions, bounding_poly, remainder_res):
    """Calculate the approximate number of triangles inside the
    bounding polygon and the other interior regions

    Polygon areas are converted to square Kms

    FIXME: Add tests for this function
    """

    # TO DO check if any of the regions fall inside one another

    log.info('-' * 80)
    log.info('Polygon  Max triangle area (m^2)  Total area (km^2)  '
             'Estimated #triangles')
    log.info('-' * 80)

    no_triangles = 0.0
    area = polygon_area(bounding_poly)

    for poly, resolution in interior_regions:
        this_area = polygon_area(poly)
        this_triangles = this_area/ resolution
        no_triangles += this_triangles
        area -= this_area

        log.info('Interior %s%s%d'
                 % (('%.0f' % resolution).ljust(25),
                    ('%.2f' % (this_area/ 1_000_000)).ljust(19),
                    this_triangles))

    bound_triangles = area/remainder_res
    no_triangles += bound_triangles

    log.info('Bounding %s%s%d'
             % (('%.0f' % remainder_res).ljust(25),
                ('%.2f' % (area/ 1_000_000)).ljust(19),
                bound_triangles))

    total_number_of_triangles = no_triangles/0.7

    log.info('Estimated total number of triangles: %d'
             % total_number_of_triangles)
    log.info('Note: This is generally about 20% less than the final amount')

    return int(total_number_of_triangles)


def decimate_polygon(polygon, factor=10):
    """Reduce number of points in polygon by the specified
    factor (default=10, hence the name of the function) such that
    the extrema in both axes are preserved.

    Reduce number of points in polygon by the specified factor.
    polygon The polygon to reduce.
    factor The factor to reduce polygon points by (default 10).

    The extrema of both axes are preserved.

    Return reduced polygon
    """

    # FIXME(Ole): This doesn't work at present,
    # but it isn't critical either

    # Find outer extent of polygon
    num_polygon = ensure_numeric(polygon)
    max_x = max(num_polygon[:, 0])
    max_y = max(num_polygon[:, 1])
    min_x = min(num_polygon[:, 0])
    min_y = min(num_polygon[:, 1])

    # Keep only some points making sure extrema are kept
    reduced_polygon = []
    for i, point in enumerate(polygon):
        if point[0] in [min_x, max_x] and point[1] in [min_y, max_y]:
            # Keep
            reduced_polygon.append(point)
        else:
            if len(reduced_polygon)*factor < i:
                reduced_polygon.append(point)

    return reduced_polygon


def interpolate_polyline(data,
                         polyline_nodes,
                         gauge_neighbour_id,
                         interpolation_points=None,
                         rtol=1.0e-6,
                         atol=1.0e-8):
    """Interpolate linearly between values data on polyline nodes
    of a polyline to list of interpolation points.

    data is the data on the polyline nodes.

    Inputs:
      data: Vector or array of data at the polyline nodes.
      polyline_nodes: Location of nodes where data is available.
      gauge_neighbour_id: ?
      interpolation_points: Interpolate polyline data to these positions.
          List of coordinate pairs [x, y] of
          data points or an nx2 numeric array or a Geospatial_data object
      rtol, atol: Used to determine whether a point is on the polyline or not.
                  See point_on_line.

    Output:
      Interpolated values at interpolation points
    """

    if isinstance(interpolation_points, Geospatial_data):
        interpolation_points = interpolation_points.\
            get_data_points(absolute=True)

    interpolated_values = num.zeros(len(interpolation_points), float)

    data = ensure_numeric(data, float)
    polyline_nodes = ensure_numeric(polyline_nodes, float)
    interpolation_points = ensure_numeric(interpolation_points, float)
    gauge_neighbour_id = ensure_numeric(gauge_neighbour_id, int)

    num_nodes = polyline_nodes.shape[0]    # Number of nodes in polyline

    # Input sanity check
    assert_msg = 'interpolation_points are not given (interpolate.py)'
    assert interpolation_points is not None, assert_msg

    assert_msg = 'function value must be specified at every interpolation node'
    assert data.shape[0] == polyline_nodes.shape[0], assert_msg

    assert_msg = 'Must define function value at one or more nodes'
    assert data.shape[0] > 0, assert_msg

    if num_nodes == 1:
        assert_msg = 'Polyline contained only one point. I need more. '
        assert_msg += str(data)
        raise Exception(assert_msg)
    elif num_nodes > 1:
        _interpolate_polyline(data,
                              polyline_nodes,
                              gauge_neighbour_id,
                              interpolation_points,
                              interpolated_values,
                              rtol,
                              atol)

    return interpolated_values


def polylist2points_verts(polylist):
    """ Convert a list of polygons to discrete points and vertices.
    """

    offset = 0
    points = []
    vertices = []
    for poly in polylist:
        points.extend(poly)
        vertices.extend([[i, i+1] for i in range(offset, offset+len(poly)-1)])
        offset += len(poly)

    return points, vertices


################################################################################
# Initialise module
################################################################################

#from polygon_ext import _intersection


if __name__ == "__main__":
    pass
