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
        
        #------------------------------------------
        # Extra aliases for changing elevation at 
        # vertices and edges
        #------------------------------------------
        self.elev_v  = self.domain.quantities['elevation'].vertex_values
        self.elev_e = self.domain.quantities['elevation'].edge_values

        #------------------------------------------
        # Need to turn off this optimization as it
        # doesn't fixup the relationship between
        # bed and stage vertex values in dry region
        #------------------------------------------
        self.domain.optimise_dry_cells = 0

        #-----------------------------------------
        # Extra structures to support maintaining
        # continuity of elevation
        #-----------------------------------------
        self.setup_node_structures()


    def __call__(self):
        """
        Apply rate to those triangles defined in indices

        indices == [], then don't apply anywhere
        indices == None, then apply everywhere
        otherwise apply for the specific indices
        """

        if self.indices is []:
            return

        #------------------------------------------
        # Apply changes to elevation vertex values
        # via the update_quantites routine
        #------------------------------------------
        if not self.update_quantities():
            return


        #------------------------------------------
        # Cleanup elevation and stage quantity values
        #-----------------------------------------
        if self.indices is None:

            #--------------------------------------
            # Make elevation continuous and clean up
            # stage values to ensure conservation
            #--------------------------------------
            height_c = self.stage_c - self.elev_c
            self.domain.quantities['elevation'].smooth_vertex_values()
            self.domain.quantities['elevation'].interpolate()
            self.stage_c[:] = self.elev_c +  height_c


        else:

            #--------------------------------------
            # Make elevation continuous and clean up
            # stage values to ensure conservation
            #--------------------------------------
            height_c = self.stage_c[self.vols] - self.elev_c[self.vols]
            for nid in self.node_ids:
                non = self.domain.number_of_triangles_per_node[nid]

                vid = num.arange(self.node_index[nid], self.node_index[nid+1],dtype=num.int)
                vidd = self.domain.vertex_value_indices[vid]

                self.elev_v[vidd/3,vidd%3] = num.sum(self.elev_v[vidd/3,vidd%3])/non


            #--------------------------------------
            # clean up the centroid values and edge values
            #--------------------------------------
            self.elev_c[self.vols] = num.mean(self.elev_v[self.vols],axis=1)

            self.elev_e[self.vols,0] = 0.5*(self.elev_v[self.vols,1]+ self.elev_v[self.vols,2])
            self.elev_e[self.vols,1] = 0.5*(self.elev_v[self.vols,2]+ self.elev_v[self.vols,0])
            self.elev_e[self.vols,2] = 0.5*(self.elev_v[self.vols,0]+ self.elev_v[self.vols,1])

            self.stage_c[self.vols] = self.elev_c[self.vols] +  height_c



    def update_quantities(self):
        """Update the vertex values of the quantities to model erosion
        """


        elevation = self.get_elevation()

        updated = True

        if self.indices is None:

            #--------------------------------------
            # Update all three vertices for each cell
            #--------------------------------------
            self.elev_v[:] = elevation

        else:

            #--------------------------------------
            # Update all three vertices for each cell
            #--------------------------------------
            self.elev_v[self.indices] = elevation


        return updated

    def get_elevation(self, t=None):
        """Get value of elevation at time t.
        If t not specified, return elevation at current domain time
        """

        from anuga.fit_interpolate.interpolate import Modeltime_too_early, \
                                                      Modeltime_too_late
                                                      
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

    def setup_node_structures(self):
        """ For setting elevation we need to
        ensure that the elevation quantity remains
        continuous (at least for version 1.3 of anuga)

        So we need to find all the vertices that need to
        update within each timestep.
        """

        node_ids = set()

        for ind in self.indices:
            for k in [0,1,2]:
                node_ids.add(self.domain.triangles[ind,k])

        self.node_ids = [ id for id in node_ids ]


        node_index = num.zeros((self.domain.number_of_nodes)+1, dtype = num.int)

        k = 0
        node_index[0] = 0
        for i in range(self.domain.number_of_nodes):
            node_index[i+1] = node_index[i] + self.domain.number_of_triangles_per_node[i]

        self.node_index = node_index

        vertex_ids =[]
        for nid in self.node_ids:
            #print nid,self.domain.number_of_triangles_per_node[nid]
            for vid in range(node_index[nid], node_index[nid+1]):
                vidd = self.domain.vertex_value_indices[vid]
                vertex_ids.append(vidd)
                #print '   ',nid, vid, vidd, vidd/3, vidd%3

        self.vol_ids  = num.array(vertex_ids,dtype=num.int)/3
        self.vols = num.array(list(set(self.vol_ids)), dtype=num.int)
        self.vert_ids = num.array(vertex_ids,dtype=num.int)%3



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



        




