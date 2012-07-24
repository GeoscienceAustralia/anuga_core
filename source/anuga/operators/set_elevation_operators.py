"""
Set elevation operators


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



class Set_elevation_operator(Operator):
    """
    Set the elevation in a region (careful to maintain continuitiy of elevation)

    indices: None == all triangles, Empty list [] no triangles

    rate can be a function of time.

    """

    def __init__(self,
                 domain,
                 elevation=None,
                 indices=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        Operator.__init__(self, domain, description, label, logging, verbose)

        #------------------------------------------
        # Local variables
        #------------------------------------------
        self.elevation = elevation
        self.indices = indices


    def __call__(self):
        """
        Apply rate to those triangles defined in indices

        indices == [], then don't apply anywhere
        indices == None, then apply everywhere
        otherwise apply for the specific indices
        """

        #if self.indices is []:
        #    return

        #elevation = self.get_elevation()

#        if self.verbose is True:
#            log.critical('Bed of %s at time = %.2f = %f'
#                         % (self.quantity_name, domain.get_time(), elevation))

        #if self.indices is None:
        #    self.elev_c[:] = elevation
        #else:
        #    self.elev_c[self.indices] = elevation

        t = self.get_time()
        dt = self.get_timestep()

        v_coors = self.domain.vertex_coordinates
        self.elev_v = self.domain.quantities['elevation'].vertex_values


        if self.indices is None:
            self.elev_v[:] = self.elev_v + 0.0
        else:
            self.elev_v[self.indices] += self.elevation(t)*dt

        ### make sure centroid is correct as well
        
        #self.domain.add_quantity('elevation', lambda x,y: dt*self.elevation(x,y,t))



        # clean up discontinuities for now
        #self.domain.quantities['elevation'].smooth_vertex_values()




    def get_elevation(self, t=None):
        """Get value of elevation at time t.
        If t not specified, return elevation at current domain time
        """

        if t is None:
            t = self.domain.get_time()

        if callable(self.elevation):
            try:
                elevation = self.elevation(t)
            except Modeltime_too_early, e:
                raise Modeltime_too_early(e)
            except Modeltime_too_late, e:
                msg = '%s: ANUGA is trying to run longer than specified data.\n' %str(e)
                msg += 'You can specify keyword argument default_rate in the '
                msg += 'rate operator to tell it what to do in the absence of time data.'
                raise Modeltime_too_late(msg)
        else:
            elevation = self.elevation


        if elevation is None:
            msg = ('Attribute elevation must be specified in '+self.__name__+
                   ' before attempting to call it')
            raise Exception(msg)

        return elevation



    def parallel_safe(self):
        """Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True

    def statistics(self):

        message = self.label + ': Set_elevation_operator'
        message = message + ' on triangles '+ str(self.indices)
        return message


    def timestepping_statistics(self):

        #message  = indent + self.label + ': Set_elevation = ' + str('')
        #message  += ' at center '+str(self.center)
        return 'test'




#===============================================================================
# Specific Bed Operators for circular region.
#===============================================================================
class Circular_set_elevation_operator(Set_elevation_operator):
    """
    Set elevation over a circular region

    """

    def __init__(self, domain,
                 elevation=0.0,
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


        Set_elevation_operator.__init__(self,
                                    domain,
                                    elevation=elevation,
                                    indices=indices,
                                    verbose=verbose)





#===============================================================================
# Specific Bed Operators for polygonal region.
#===============================================================================
class Polygonal_set_elevation_operator(Set_elevation_operator):
    """
    Add water at certain rate (ms^{-1} = vol/Area/sec) over a
    polygonal region

    rate can be a function of time.

    """

    def __init__(self, domain,
                 elevation=0.0,
                 polygon=None,
                 verbose=False):


        # Determine indices in update region
        N = domain.get_number_of_triangles()
        points = domain.get_centroid_coordinates(absolute=True)


        indices = inside_polygon(points, polygon)
        self.polygon = polygon

        # It is possible that circle doesn't intersect with mesh (as can happen
        # for parallel runs


        Set_elevation_operator.__init__(self,
                               domain,
                               elevation=elevation,
                               indices=indices,
                               verbose=verbose)



        




