"""
    Callable function to determine if points lie inside or outside a polygon.
    
    As of June 2010 this module has a pylint quality rating of 8.85/10.
"""

import anuga.utilities.log as log
import numpy as num
from .polygon import inside_polygon

class Polygon_function(object):
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
        """Create instance of a polygon function.

        regions A list of (x,y) tuples defining a polygon.
        default Value or function returning value for points outside poly.
        geo_reference ??
        """

        try:
            len(regions)
        except:
            msg = ('Polygon_function takes a list of pairs (polygon, value).'
                   'Got %s' % str(regions))
            raise Exception(msg)

        first_region = regions[0]

        if isinstance(first_region, str):
            msg = ('You passed in a list of text values into polygon_function '
                   'instead of a list of pairs (polygon, value): "%s"'
                   % str(first_region))
            raise Exception(msg)

        try:
            num_region_components = len(first_region)
        except:
            msg = ('Polygon_function takes a list of pairs (polygon, value). '
                   'Got %s' % str(num_region_components))
            raise Exception(msg)

        msg = ('Each entry in regions have two components: (polygon, value). '
               'I got %s' % str(num_region_components))
        assert num_region_components == 2, msg

        if geo_reference is None:
            from anuga.coordinate_transforms.geo_reference import Geo_reference
            geo_reference = Geo_reference()

        self.default = default

        # Make points in polygons relative to geo_reference
        self.regions = []
        for polygon, value in regions:
            georeffed_poly = geo_reference.change_points_geo_ref(polygon)
            self.regions.append((georeffed_poly, value))

    def __call__(self, pts_x, pts_y):
        """Implement the 'callable' property of Polygon_function.

        x List of x coordinates of points ot interest.
        y List of y coordinates of points ot interest.
        """
        pts_x = num.array(pts_x, float)
        pts_y = num.array(pts_y, float)

        # x and y must be one-dimensional and same length
        assert len(pts_x.shape) == 1 and len(pts_y.shape) == 1
        pts_len = pts_x.shape[0]
        assert pts_y.shape[0] == pts_len, 'x and y must be same length'

        points = num.ascontiguousarray(num.concatenate((pts_x[:, num.newaxis],
                                                        pts_y[:, num.newaxis]),
                                                       axis = 1 ))

        if callable(self.default):
            result = self.default(pts_x, pts_y)
        else:
            result = num.ones(pts_len, float) * self.default

        for polygon, value in self.regions:
            indices = inside_polygon(points, polygon)

            # FIXME: This needs to be vectorised
            if callable(value):
                for i in indices:
                    xx = num.array([pts_x[i]])
                    yy = num.array([pts_y[i]])
                    result[i] = value(xx, yy)[0]
            else:
                for i in indices:
                    result[i] = value

        if len(result) == 0:
            msg = ('Warning: points provided to Polygon function did not fall '
                   'within its regions in [%.2f, %.2f], y in [%.2f, %.2f]'
                   % (min(pts_x), max(pts_x), min(pts_y), max(pts_y)))
            log.critical(msg)

        return result
