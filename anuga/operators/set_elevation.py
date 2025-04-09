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

from anuga.operators.set_quantity import Set_quantity
from anuga.config import indent

class Set_elevation(Set_quantity):
    """
    Helper class to setup calculation of elevation
    associated with a region (defined by indices, polygon or center/radius
    """
    
    set_elevation = Set_quantity.set_value


    def __init__(self,
                 domain,
                 elevation=None,
                 region=None,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 line=None,
                 verbose = False):

        Set_quantity.__init__(self, domain, 'elevation',
                              value=elevation,
                              region=region,
                              indices=indices,
                              polygon=polygon,
                              center=center,
                              radius=radius,
                              line=line,
                              verbose=verbose,
                              test_elevation=False)



        #-----------------------------------------
        # Extra structures to support maintaining
        # continuity (smoothness) of elevation
        # FIXME SR: Not needed if using discontinuous elevation
        #-----------------------------------------
        self.setup_node_structures()


        #------------------------------------------
        # Extra aliases for changing elevation
        # to ensure conservation of mass and
        # continuity at vertices and edges
        #------------------------------------------
        self.elev_v  = self.domain.quantities['elevation'].vertex_values
        self.elev_e = self.domain.quantities['elevation'].edge_values
        self.stage_c = self.domain.quantities['stage'].centroid_values
        self.elev_c = self.domain.quantities['elevation'].centroid_values

        #------------------------------------------
        # x,y coordinates of vertices of cells that are
        # updated
        #------------------------------------------
        coord_v = self.domain.vertex_coordinates
        ids = self.indices
        if ids is None:
            self.v_x = coord_v[:,0].reshape((-1,3))
            self.v_y = coord_v[:,1].reshape((-1,3))
        else:
            self.v_x = coord_v[:,0].reshape((-1,3))[ids]
            self.v_y = coord_v[:,1].reshape((-1,3))[ids]

        #print self.v_x.shape
        #print self.v_y.shape

        #------------------------------------------
        # Need to turn off this optimization as it
        # doesn't fixup the relationship between
        # bed and stage vertex values in dry region
        #------------------------------------------
        self.domain.optimise_dry_cells = 0




    def __call__(self):
        """
        Apply rate to those triangles defined in indices

        indices == [], then don't apply anywhere
        indices is None, then apply everywhere
        otherwise apply for the specific indices
        """

        if self.value is None:
            return

        if self.indices is []:
            return

        #------------------------------------------
        # Apply changes to elevation vertex values
        # via the update_quantites routine
        # Assume vertex values updated and need to 
        # fix up centroid values unless 
        # domain.get_discontinuous_elevation is true
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
            if self.domain.get_using_discontinuous_elevation():
                pass
            else:
                height_c = self.stage_c - self.elev_c
                self.domain.quantities['elevation'].smooth_vertex_values()
                self.domain.quantities['elevation'].interpolate()
                self.stage_c[:] = self.elev_c +  height_c


        else:

            #--------------------------------------
            # Make elevation continuous and clean up
            # stage values to ensure conservation
            #--------------------------------------
            
            if self.domain.get_using_discontinuous_elevation():
                pass
            else:
                height_c = self.stage_c[self.vols] - self.elev_c[self.vols]
                for nid in self.node_ids:
                    non = self.domain.number_of_triangles_per_node[nid]
    
                    vid = num.arange(self.node_index[nid], self.node_index[nid+1], dtype=int)
                    vidd = self.domain.vertex_value_indices[vid]

                    # Replaced this (Ole)
                    #self.elev_v[vidd//3,vidd%3] = num.sum(self.elev_v[vidd//3,vidd%3])/non

                    # with this to get it working in both Python2 and Python3
                    res = num.sum(self.elev_v[vidd // 3, vidd % 3]) / non
                    self.elev_v[vidd // 3, vidd % 3] = res

                    
                    
                #--------------------------------------
                # clean up the centroid values and edge values
                #--------------------------------------
                self.elev_c[self.vols] = num.mean(self.elev_v[self.vols],axis=1)
    
                self.elev_e[self.vols,0] = 0.5*(self.elev_v[self.vols,1]+ self.elev_v[self.vols,2])
                self.elev_e[self.vols,1] = 0.5*(self.elev_v[self.vols,2]+ self.elev_v[self.vols,0])
                self.elev_e[self.vols,2] = 0.5*(self.elev_v[self.vols,0]+ self.elev_v[self.vols,1])
    
                self.stage_c[self.vols] = self.elev_c[self.vols] +  height_c



    def update_quantities(self):
        """Update the vertex and centroid values of the quantities to model erosion
        """

        if self.value is None:
            return False


        updated = True

        if self.indices is None:

            if self.domain.get_using_discontinuous_elevation():
                try:
                    height_c = self.stage_c - self.elev_c
                    value = self.get_value(x=self.coord_c[:,0], y=self.coord_c[:,1])
                    self.elev_c[:] = value
                    self.stage_c[:] =  self.elev_c + height_c
                except ValueError:
                    updated = False
                    pass
            else:
            #--------------------------------------
            # Update all three vertices for each cell
            #--------------------------------------            
                try:
                    value = self.get_value(self.v_x, self.v_y)
                    self.elev_v[:] = value
                except ValueError:
                    updated = False
                    pass
     
        #----------------------------------
        # Apply just to indices
        #----------------------------------
        else: 

            if self.domain.get_using_discontinuous_elevation():
                ids = self.indices
                x = self.coord_c[ids,0]
                y = self.coord_c[ids,1]
                try:
                    height_c = self.stage_c[ids] - self.elev_c[ids]
                    value = self.get_value(x=x,y=y)
                    self.elev_c[ids] = value
                    self.stage_c[ids] = self.elev_c[ids] + height_c
                except ValueError:
                    updated = False
                    pass
            else:
                ids = self.indices
                try:
                    value = self.get_value(self.v_x, self.v_y)
                    self.elev_v[ids] = value
                except ValueError:
                    updated = False
                    pass


        return updated




    def setup_node_structures(self):
        """ For setting elevation we need to
        ensure that the elevation quantity remains
        continuous (at least for version 1.3 of anuga)

        So we need to find all the vertices that need to
        update within each timestep.
        """



        if self.indices is None or self.indices is []:
            self.vol_ids  = None
            self.vols = None
            self.vert_ids = None
            return


        
        node_ids = set()

        for ind in self.indices:
            for k in [0,1,2]:
                node_ids.add(self.domain.triangles[ind,k])

        self.node_ids = [ id for id in node_ids ]



        node_index = num.zeros((self.domain.number_of_nodes)+1, dtype = int)

        # FIXME: SR Don't we calculate this for the domain already!
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

        self.vol_ids = num.array(vertex_ids, dtype=int) // 3  # FIXME(Ole): Tests past both with / and //
        self.vols = num.array(list(set(self.vol_ids)), dtype=int)
        self.vert_ids = num.array(vertex_ids,dtype=int)%3






