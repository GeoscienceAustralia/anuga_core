"""
Set value operators

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
from anuga.config import indent



class Set_w_uh_vh_operator(Operator):
    """
    Set the w, uh and vh in a region

    indices: None == all triangles, Empty list [] no triangles

    rate can be a function of time.

    """

    def __init__(self,
                 domain,
                 w_uh_vh=None,
                 indices=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        Operator.__init__(self, domain, description, label, logging, verbose)

        #------------------------------------------
        # Local variables
        #------------------------------------------
        self.w_uh_vh = w_uh_vh
        self.indices = indices


    def __call__(self):
        """
        Apply w_uh_vh to those triangles defined in indices

        indices == [], then don't apply anywhere
        indices == None, then apply everywhere
        otherwise apply for the specific indices
        """


        if self.indices is []:
            return

        w_uh_vh = self.get_w_uh_vh()

        if w_uh_vh is None:
            return

        if self.verbose is True:
            log.critical('w_uh_vh of %s at time = %.2f = %f'
                         % (self.quantity_name, domain.get_time(), stage))

        if self.indices is None:
            self.stage_c[:] = w_uh_vh[0]
            self.xmom_c[:]  = w_uh_vh[1]
            self.ymom_c[:]  = w_uh_vh[2]
        else:
            self.stage_c[self.indices] = w_uh_vh[0]
            self.xmom_c[self.indices]  = w_uh_vh[1]
            self.ymom_c[self.indices]  = w_uh_vh[2]


    def get_w_uh_vh(self, t=None):
        """Get value of w_uh_vh at time t.
        If t not specified, return stage at current domain time
        """

        if t is None:
            t = self.domain.get_time()

        if callable(self.w_uh_vh):
            try:
                w_uh_vh = self.w_uh_vh(t)
            except Modeltime_too_early, e:
                raise Modeltime_too_early(e)
            except Modeltime_too_late, e:
                msg = '%s: ANUGA is trying to run longer than specified data.\n' %str(e)
                msg += 'You can specify keyword argument default_rate in the '
                msg += 'function to tell it what to do in the absence of time data.'
                raise Modeltime_too_late(msg)
        else:
            w_uh_vh = self.w_uh_vh


        if w_uh_vh  is None:
            msg = ('Attribute w_uh_vh must be specified in '+self.__name__+
                   ' before attempting to call it')
            raise Exception(msg)

        return w_uh_vh



    def parallel_safe(self):
        """Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True

    def statistics(self):

        message = self.label + ': Set_w_uh_vh_operator'
        message = message + ' on triangles '+ str(self.indices)
        return message


    def timestepping_statistics(self):

        message  = indent + self.label + ': Set_w_uh_vh = ' + str(self.get_w_uh_vh())
        return message




#===============================================================================
# Specific Stage Operators for circular region.
#===============================================================================
class Circular_set_w_uh_vh_operator(Set_w_uh_vh_operator):
    """
    Set w_uh_vh over a circular region

    """

    def __init__(self, domain,
                 w_uh_vh=[0.0, 0.0, 0.0],
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


        Set_w_uh_vh_operator.__init__(self,
                                    domain,
                                    w_uh_vh=w_uh_vh,
                                    indices=indices,
                                    verbose=verbose)





#===============================================================================
# Specific Stage Operators for polygonal region.
#===============================================================================
class Polygonal_set_w_uh_vh_operator(Set_w_uh_vh_operator):
    """
    Add water at certain rate (ms^{-1} = vol/Area/sec) over a
    polygonal region

    rate can be a function of time.

    """

    def __init__(self, domain,
                 w_uh_vh=[0.0, 0,0, 0.0],
                 polygon=None,
                 verbose=False):


        # Determine indices in update region
        N = domain.get_number_of_triangles()
        points = domain.get_centroid_coordinates(absolute=True)


        indices = inside_polygon(points, polygon)
        self.polygon = polygon

        # It is possible that circle doesn't intersect with mesh (as can happen
        # for parallel runs


        Set_w_uh_vh_operator.__init__(self,
                               domain,
                               w_uh_vh=w_uh_vh,
                               indices=indices,
                               verbose=verbose)



        




