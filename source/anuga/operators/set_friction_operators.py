"""
Set friction operators

Constraints: See GPL license in the user guide
Version: 1.0 ($Revision: 7731 $)
"""

__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"


from anuga import Domain
from anuga import Quantity
import numpy as num
import anuga.utilities.log as log

from anuga.geometry.polygon import inside_polygon

from anuga.operators.base_operator import Operator
from anuga.fit_interpolate.interpolate import Modeltime_too_early, \
                                              Modeltime_too_late
from anuga import indent


default_friction_min = 0.01
default_friction_max = 0.035

class Set_depth_friction_operator(Operator):
    """
    Set the friction in a region

    indices: None == all triangles, Empty list [] no triangles

    rate can be a function of time.

    """

    def __init__(self,
                 domain,
                 friction_min=default_friction_min,
                 friction_max=default_friction_max,
                 indices = None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        Operator.__init__(self, domain, description, label, logging, verbose)

        #------------------------------------------
        # Local variables
        #------------------------------------------

        assert friction_min >= 0.0
        assert friction_min <= friction_max

        self.friction_min = friction_min
        self.friction_max = friction_max
        self.friction_c = self.domain.get_quantity('friction').centroid_values
        self.indices = indices


    def __call__(self):
        """
        Change friction based on depth
        """

        if self.indices is []:
            return


        #-----------------------------------------
        # Here is where the important formula is applied
        #----------------------------------------
        if self.indices is None:
            height = self.stage_c - self.elev_c
            self.friction_c[:] = (self.friction_max - self.friction_min)/(1.0 - height) + self.friction_min
        else:
            ind = self.indices
            height = self.stage_c[ind] - self.elev_c[ind]
            self.friction_c[ind] = (self.friction_max - self.friction_min)/(1.0 - height) + self.friction_min



    def parallel_safe(self):
        """Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True

    def statistics(self):

        message = self.label + ': Set_depth_friction_operator'
        return message


    def timestepping_statistics(self):

        message  = indent + self.label + ': Set_depth_friction, '
        message  += ' min '+str(self.friction_min)
        message  += ' max '+str(self.friction_max)
        return message



#===============================================================================
# Specific Bed Operators for circular region.
#===============================================================================
class Circular_set_depth_friction_operator(Set_depth_friction_operator):
    """
    Set elevation over a circular region

    """

    def __init__(self, domain,
                 friction_min=default_friction_min,
                 friction_max=default_friction_max,
                 center=None,
                 radius=None,
                 verbose=False):

        assert center is not None
        assert radius is not None


        # Determine indices in update region
        N = domain.get_number_of_triangles()
        points = domain.get_centroid_coordinates(absolute=True)


        indices = []

        c = center
        r = radius

        self.center = center
        self.radius = radius

        intersect = False
        for k in range(N):
            x, y = points[k,:]    # Centroid

            if ((x-c[0])**2+(y-c[1])**2) < r**2:
                intersect = True
                indices.append(k)


        msg = 'No centroids intersect circle center'+str(center)+' radius '+str(radius)
        assert intersect, msg




        # It is possible that circle doesn't intersect with mesh (as can happen
        # for parallel runs


        Set_depth_friction_operator.__init__(self,
                                    domain,
                                    friction_min=friction_min,
                                    friction_max=friction_max,
                                    indices=indices,
                                    verbose=verbose)





#===============================================================================
# Specific Bed Operators for polygonal region.
#===============================================================================
class Polygonal_depth_friction_operator(Set_depth_friction_operator):
    """
    Add water at certain rate (ms^{-1} = vol/Area/sec) over a
    polygonal region

    rate can be a function of time.

    """

    def __init__(self, domain,
                 friction_min=default_friction_min,
                 friction_max=default_friction_max,
                 polygon=None,
                 verbose=False):


        # Determine indices in update region
        N = domain.get_number_of_triangles()
        points = domain.get_centroid_coordinates(absolute=True)


        indices = inside_polygon(points, polygon)
        self.polygon = polygon

        # It is possible that circle doesn't intersect with mesh (as can happen
        # for parallel runs


        Set_depth_friction_operator.__init__(self,
                               domain,
                               friction_min=friction_min,
                               friction_max=friction_max,
                               indices=indices,
                               verbose=verbose)



