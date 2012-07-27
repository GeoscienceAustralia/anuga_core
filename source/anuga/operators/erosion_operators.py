"""
Erosion operators


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

def lineno():
    """Returns the current line number in our program."""
    import inspect
    return inspect.currentframe().f_back.f_back.f_lineno



def stage_elev_info(self):
    print 80*"="

    print 'In Evolve: line number ', lineno()
    import inspect
    print inspect.getfile(lineno)

    print 80*"="
    ind = num.array([ 976,  977,  978,  979,  980,  981,  982,  983, 1016, 1017, 1018,
             1019, 1020, 1021, 1022, 1023])
    elev_v = self.get_quantity('elevation').vertex_values
    stage_v = self.get_quantity('stage').vertex_values
    elev_c = self.get_quantity('elevation').centroid_values
    stage_c = self.get_quantity('stage').centroid_values

    from pprint import pprint
    print 'elev_v, elev_c, elev_avg \n'
    pprint( num.concatenate( (elev_v[ind], (elev_c[ind]).reshape(16,1),
                               num.mean(elev_v[ind],axis=1).reshape(16,1)), axis = 1))
    print 'stage_v, stage_c, stage_avg \n'
    pprint( num.concatenate( (stage_v[ind], (stage_c[ind]).reshape(16,1),
                               num.mean(stage_v[ind],axis=1).reshape(16,1)), axis = 1))

    

    print 80*"="


class Erosion_operator(Operator):
    """
    Simple erosion operator in a region (careful to maintain continuitiy of elevation)

    indices: None == all triangles, Empty list [] no triangles

    rate can be a function of time.

    """

    def __init__(self,
                 domain,
                 threshold= 0.0,
                 indices=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        Operator.__init__(self, domain, description, label, logging, verbose)

        #------------------------------------------
        # Local variables
        #------------------------------------------
        self.threshold = threshold
        self.indices = indices


        #------------------------------------------
        # Find vertex nodes
        #------------------------------------------
        node_id = set()

        for ind in self.indices:
            for k in [0,1,2]:
                node_id.add(self.domain.triangles[ind,k])

        self.node_id = [ id for id in node_id ]


        node_index = num.zeros((self.domain.number_of_nodes)+1, dtype = num.int)

        k = 0
        node_index[0] = self.domain.number_of_triangles_per_node[0]
        for i in range(self.domain.number_of_nodes):
            if i == 0:
                continue

            node_index[i] = node_index[i-1] + self.domain.number_of_triangles_per_node[i]

        n0 = self.node_id[0]
        i0 = node_index[n0]                                      # index of first vertex associated to node
        i1 = i0 + self.domain.number_of_triangles_per_node[i]    # index of last vertex associated to node

        print node_index
        print n0,i0,i1
        print len(node_index),self.domain.number_of_nodes


        print self.node_id


    def __call__(self):
        """
        Apply rate to those triangles defined in indices

        indices == [], then don't apply anywhere
        indices == None, then apply everywhere
        otherwise apply for the specific indices
        """



        if self.indices is []:
            return



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


        #print self.indices
        

        self.elev_v  = self.domain.quantities['elevation'].vertex_values
        self.stage_v = self.domain.quantities['stage'].vertex_values


        # Need to store water heights before change to ensure
        # no water lost or produced
        height_c = self.stage_c - self.elev_c


        #stage_elev_info(self.domain)
        #--------------------------------------------
        # Here we do the actual erosion
        #--------------------------------------------

        if t < 15.0 and t > 10.0:
            if self.indices is None:
                self.elev_v[:] = self.elev_v + 0.0
            else:
                self.elev_v[self.indices] -= 0.1*dt


        for nid in self.node_id:
            print 'nid ', nid
            print 'vvi ', self.domain.vertex_value_indices[nid]


        #stage_elev_info(self.domain)
        

        # FIXME SR: At present need to ensure the elevation is continuous
        # In future with discontinuous bed we will not need to do this.
        self.domain.quantities['elevation'].smooth_vertex_values()

        # Need to do this faster.

        #stage_elev_info(self.domain)

        self.domain.quantities['elevation'].interpolate()

        #stage_elev_info(self.domain)

        #self.elev_c = self.domain.quantities['elevation'].centroid_values

#        # Fix up water conservation
        self.stage_c[:] = self.elev_c +  height_c

        #stage_elev_info(self.domain)


        old_flag = self.domain.optimise_dry_cells
        self.domain.optimise_dry_cells = 0
        self.domain.distribute_to_vertices_and_edges()
        self.domain.optimise_dry_cells = old_flag

        #stage_elev_info(self.domain)
        #print 'time in erosion ',self.get_time(), dt



    def parallel_safe(self):
        """Operator is applied independently on each cell and
        so is parallel safe.
        """
        return False

    def statistics(self):

        message = self.label + ': Erosion_operator'
        message = message + ' on triangles '+ str(self.indices)
        return message


    def timestepping_statistics(self):

        message  = indent + self.label + ': Erosion_operator'
        return message




#===============================================================================
# Specific Erosion Operator for circular region.
#===============================================================================
class Circular_erosion_operator(Erosion_operator):
    """
    Erosion over a circular region

    """

    def __init__(self, domain,
                 threshold=0.0,
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


        Erosion_operator.__init__(self,
                                  domain,
                                  threshold,
                                  indices=indices,
                                  verbose=verbose)





#===============================================================================
# Specific Bed Operators for polygonal region.
#===============================================================================
class Polygonal_erosion_operator(Erosion_operator):
    """
    Erosion over a ploygon

    """

    def __init__(self, domain,
                 threshold=0.0,
                 polygon=None,
                 verbose=False):


        # Determine indices in update region
        N = domain.get_number_of_triangles()
        points = domain.get_centroid_coordinates(absolute=True)


        indices = inside_polygon(points, polygon)
        self.polygon = polygon

        # It is possible that circle doesn't intersect with mesh (as can happen
        # for parallel runs


        Erosion_operator.__init__(self,
                                  domain,
                                  threshold=threshold,
                                  indices=indices,
                                  verbose=verbose)



        




