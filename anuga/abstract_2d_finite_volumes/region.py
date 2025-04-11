"""
Define region
"""

__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"

from anuga.utilities.numerical_tools import ensure_numeric

from anuga import Domain
from anuga import Quantity
import numpy as num
from pprint import pprint


from anuga.geometry.polygon import inside_polygon, line_intersect

from anuga.utilities.function_utils import determine_function_type

#from anuga import indent


class Region(object):
    """ Object which defines a region within the domain

    """

    def __init__(self,
                 domain,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 line=None,
                 poly=None,
                 expand_polygon=False,
                 verbose = False):

        """Create a Region object

        :param domain: Region must be defined wrt a domain
        :param indices: Define the region by triangle IDs
        :param polygon: List of [x,y] points to define region
        :param center: point [x,y] which defines the centre of a circle
        :param radius: radius of a circle which defines a region
        :param line: List of [x,y] points defining a polyline
        :param poly: An old argument which was used to define a polyline or polygon
        :param expand_polygon: If set true, then calculation of intersection of polygon with triangles based on vertices, otherwise based just on centroids
        :param verbose: Set to True for more verbose output 

        Setup region (defined by indices, polygon or center/radius).
        Useful in defining where to apply certain operations
        
        """


        #------------------------------------------
        # Local variables
        #------------------------------------------
        self.indices = indices
        self.center = center
        self.radius = radius
        self.polygon = polygon
        self.line = line
        self.poly = poly
        self.type = 'full_region'
        self.expand_polygon = expand_polygon
        self.verbose =  verbose

        #------------------------------------------
        #  Useful aliases
        #------------------------------------------
        self.domain = domain
        self.coord_c = self.domain.centroid_coordinates
        self.areas = self.domain.areas


        #-------------------------------------------
        # Work out indices associated with region:
        #-------------------------------------------
        if self.indices is not None:
            # This overrides polygon, center and radius, line

            assert self.radius is None
            assert self.center is None
            assert self.polygon is None
            assert self.line is None

            self.indices = num.asarray(self.indices)

            if self.indices.size == 0:
                self.indices = []
                self.type = 'empty'
            else:
                self.type = 'indices_specified'

        elif (self.center is not None) and (self.radius is not None):

            assert self.indices is None
            assert self.polygon is None
            assert self.line is None

            self._setup_indices_circle()
            self.type = 'circle'

        elif (self.polygon is not None):

            assert self.indices is None
            assert self.radius is None
            assert self.center is None
            assert self.line is None

            self._setup_indices_polygon()
            self.type = 'polygon'

        elif (self.line is not None):

            assert self.indices is None
            assert self.radius is None
            assert self.center is None
            assert self.polygon is None

            self._setup_indices_line()
            self.type = 'line'

        elif (self.poly is not None):
            # could be either a line or a polygon
            # This is essentially for backwards compatibility

            assert self.indices is None
            assert self.radius is None
            assert self.center is None
            assert self.polygon is None
            assert self.line is None

            self.poly = num.asarray(self.poly)
            if len(self.poly) > 2:
                self.polygon = self.poly
                self._setup_indices_polygon()
                self.type = 'polygon'
            else:
                self.line = self.poly
                self._setup_indices_line()
                self.type = 'line'
        else:
            assert self.indices is None or self.indices is []


        if self.indices is None:
            self.full_indices = num.where(self.domain.tri_full_flag ==1)[0]
        elif len(self.indices) == 0:
            self.full_indices = []
        else:
            self.full_indices = num.array(self.indices)[self.domain.tri_full_flag[self.indices]==1]


    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)

    def get_type(self):
        return self.type


    def plot_region(self, filename=None):

        try:
            import matplotlib
            #matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import matplotlib.tri as tri
        except:
            msg ="Couldn't import module from matplotlib, probably you need to update matplotlib"
            raise msg

        vertices = self.domain.get_vertex_coordinates()
        full_mask = num.repeat(self.domain.tri_full_flag == 1, 3)

        region_mask = num.zeros(len(self.domain),int).astype(bool)
        region_mask[self.indices] = True
        region_mask = num.repeat(region_mask,3)


        #num.repeat(self.indices, 3)

        #pprint(region_mask)

        # Gather full and region nodes
        fx = vertices[full_mask,0]
        fy = vertices[full_mask,1]
        gx = vertices[region_mask,0]
        gy = vertices[region_mask,1]


        # Plot mesh
        n = len(fx) // 3
        triang = num.array(list(range(0,3*n)))
        triang.shape = (n, 3)
        plt.triplot(fx, fy, triang, 'b-')

        # Plot region
        n = len(gx) // 3
        if n > 0:
            triang = num.array(list(range(0,3*n)))
            triang.shape = (n, 3)
            plt.triplot(gx, gy, triang, 'r-')

        # Save triangulation to location pointed by filename
        if filename is not None: plt.savefig(filename)

        plt.show()



    def _setup_indices_circle(self):

        # Determine indices in circular region
        N = self.domain.get_number_of_triangles()
        points = self.domain.get_centroid_coordinates(absolute=True)

        indices = []

        c = self.center
        r = self.radius


        intersect = False
        for k in range(N):
            x, y = points[k,:]    # Centroid

            if ((x-c[0])**2+(y-c[1])**2) < r**2:
                intersect = True
                indices.append(k)

        if len(indices) == 0:
            self.indices = indices
        else:
            self.indices = num.asarray(indices)

        if not self.domain.parallel:
            msg = 'No centroids intersect circle center'+str(c)+' radius '+str(r)
            if not intersect: raise Exception(msg)


    def _setup_indices_polygon(self):

        # Determine indices for polygonal region
        points = self.domain.get_centroid_coordinates(absolute=True)
        vertex_coordinates = self.domain.get_vertex_coordinates(absolute=True)

        indices = inside_polygon(points, self.polygon)

        if self.expand_polygon :
            n = len(self.polygon)
            for j in range(n):
                tris_0 = line_intersect(vertex_coordinates,
                                        [self.polygon[j],self.polygon[(j+1)%n]])
                indices = num.union1d(tris_0, indices)

        if len(indices) == 0:
            self.indices = indices
        else:
            self.indices = num.asarray(indices)


        if not self.domain.parallel:
            # only warn if not parallel as we should get lots of subdomains without indices
            if len(indices) == 0:
                msg = 'No centroids found for polygon %s '% str(self.polygon)
                import warnings
                warnings.warn(msg)



    def _setup_indices_line(self):

        # Determine indices for triangles intersecting a line  region

        vertex_coordinates = self.domain.get_vertex_coordinates(absolute=True)

        indices = line_intersect(vertex_coordinates, self.line)

        if len(indices) == 0:
            self.indices = indices
        else:
            self.indices = num.asarray(indices)

        if not self.domain.parallel:
            msg = 'No centroids intersecting line %s '% str(self.line)
            if len(indices) == 0: raise Exception(msg)


    def get_indices(self, full_only=True):

        if full_only:
            return self.full_indices
        else:
            return self.indices

    def set_verbose(self, verbose=True):

        self.verbose = verbose

class Centroid_field(object):

    def __init__(self, region, value, verbose=None):

        self.region = region
        self.value = value
        self.domain = self.region.domain


    def set_value(self,value):
        """Set value
        Can change value while running
        Can be a scalar, or a function of t or x,y or x,y,t or a quantity
        """

        # Test if rate is a quantity
        if isinstance(value, Quantity):
            self.value_type = 'quantity'
        else:
            # Possible types are 'scalar', 't', 'x,y' and 'x,y,t'
            from anuga.utilities.function_utils import determine_function_type
            self.value_type = determine_function_type(value)

        self.value = value

        if self.value_type == 'scalar':
            self.value_callable = False
            self.value_spatial = False
        elif self.value_type == 'quantity':
            self.value_callable = False
            self.valuespatial = False
        elif self.value_type == 't':
            self.value_callable = True
            self.value_spatial = False
        else:
            self.value_callable = True
            self.value_spatial = True
