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
from anuga import indent



class Set_stage_operator(Operator):
    """
    Set the stage in a region

    indices: None == all triangles, Empty list [] no triangles

    rate can be a function of time.

    """

    def __init__(self,
                 domain,
                 stage=None,
                 indices=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        Operator.__init__(self, domain, description, label, logging, verbose)

        #------------------------------------------
        # Local variables
        #------------------------------------------
        self.stage = stage
        self.indices = indices


    def __call__(self):
        """
        Apply rate to those triangles defined in indices

        indices == [], then don't apply anywhere
        indices == None, then apply everywhere
        otherwise apply for the specific indices
        """


        if self.indices is []:
            return

        stage = self.get_stage()

        if stage is None:
            return

        if self.verbose is True:
            log.critical('Stage of %s at time = %.2f = %f'
                         % (self.quantity_name, domain.get_time(), stage))

        if self.indices is None:
            self.stage_c[:] = stage
        else:
            self.stage_c[self.indices] = stage


    def get_stage(self, t=None):
        """Get value of stage at time t.
        If t not specified, return stage at current domain time
        """

        if t is None:
            t = self.domain.get_time()

        if callable(self.stage):
            try:
                stage = self.stage(t)
            except Modeltime_too_early, e:
                raise Modeltime_too_early(e)
            except Modeltime_too_late, e:
                msg = '%s: ANUGA is trying to run longer than specified data.\n' %str(e)
                msg += 'You can specify keyword argument default_rate in the '
                msg += 'stage function to tell it what to do in the absence of time data.'
                raise Modeltime_too_late(msg)
        else:
            stage = self.stage


        if stage is None:
            msg = ('Attribute stage must be specified in '+self.__name__+
                   ' before attempting to call it')
            raise Exception(msg)

        return stage



    def parallel_safe(self):
        """Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True

    def statistics(self):

        message = self.label + ': Set_stage_operator'
        message = message + ' on triangles '+ str(self.indices)
        return message


    def timestepping_statistics(self):

        message  = indent + self.label + ': Set_stage = ' + str(self.get_stage())
        message  += ' at center '+str(self.center)
        return message




#===============================================================================
# Specific Stage Operators for circular region.
#===============================================================================
class Circular_set_stage_operator(Set_stage_operator):
    """
    Set stage over a circular region

    """

    def __init__(self, domain,
                 stage=0.0,
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


        Set_stage_operator.__init__(self,
                                    domain,
                                    stage=stage,
                                    indices=indices,
                                    verbose=verbose)





#===============================================================================
# Specific Stage Operators for polygonal region.
#===============================================================================
class Polygonal_set_stage_operator(Set_stage_operator):
    """
    Add water at certain rate (ms^{-1} = vol/Area/sec) over a
    polygonal region

    rate can be a function of time.

    """

    def __init__(self, domain,
                 stage=0.0,
                 polygon=None,
                 verbose=False):


        # Determine indices in update region
        N = domain.get_number_of_triangles()
        points = domain.get_centroid_coordinates(absolute=True)


        indices = inside_polygon(points, polygon)
        self.polygon = polygon

        # It is possible that circle doesn't intersect with mesh (as can happen
        # for parallel runs


        Set_stage_operator.__init__(self,
                               domain,
                               stage=stage,
                               indices=indices,
                               verbose=verbose)



        




