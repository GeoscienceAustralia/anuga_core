"""
Define region
"""

__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"


from anuga import Domain
from anuga import Quantity
import numpy as num


from anuga.geometry.polygon import inside_polygon

from anuga.utilities.function_utils import determine_function_type

#from anuga import indent


class Region(object):
    """
    Helper class to setup region (defined by indices, polygon or center/radius)
    """

    def __init__(self,
                 domain,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 verbose = False):


        #------------------------------------------
        # Local variables
        #------------------------------------------
        self.indices = indices
        self.center = center
        self.radius = radius
        self.polygon = polygon

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
            # This overrides polygon, center and radius
            self.indices = num.asarray(self.indices)

            if self.indices.size == 0:
                self.indices = []

        elif (self.center is not None) and (self.radius is not None):

            assert self.indices is None
            assert self.polygon is None

            self.setup_indices_circle()

        elif (self.polygon is not None):

            assert self.indices is None

            self.setup_indices_polygon()
        else:
            assert self.indices is None or self.indices is []



        if self.indices == []:
            self.full_indices = []
        elif self.indices is None:
            self.full_indices = num.where(self.domain.tri_full_flag ==1)[0]
        else:
            self.full_indices = self.indices[self.domain.tri_full_flag[self.indices]==1]


    def setup_indices_circle(self):

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

        if indices is []:
            self.indices = indices
        else:
            self.indices = num.asarray(indices)

        msg = 'No centroids intersect circle center'+str(c)+' radius '+str(r)
        assert intersect, msg


    def setup_indices_polygon(self):

        # Determine indices for polygonal region
        points = self.domain.get_centroid_coordinates(absolute=True)

        indices = num.asarray(inside_polygon(points, self.polygon))

        if indices is []:
            self.indices = indices
        else:
            self.indices = num.asarray(indices)


    def get_indices(self, full_only=True):

        if full_only:
            return self.full_indices
        else:
            return self.indices

        


