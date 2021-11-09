"""
Set stage

Ensures water height is non-negative
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

class Set_stage(Set_quantity):
    """
    Helper class to setup calculation of elevation
    associated with a region (defined by indices, polygon or center/radius
    """
    
    set_elevation = Set_quantity.set_value


    def __init__(self,
                 domain,
                 stage=None,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 line=None,
                 verbose = False):

        Set_quantity.__init__(self, domain, 'stage',
                              value=stage,
                              indices=indices,
                              polygon=polygon,
                              center=center,
                              radius=radius,
                              line=line,
                              verbose=verbose,
                              test_stage=False)

        #------------------------------------------
        # Extra aliases for changing stage
        #------------------------------------------
        self.stage_c = self.domain.quantities['stage'].centroid_values
        self.elev_c = self.domain.quantities['elevation'].centroid_values
        self.height_c = self.domain.quantities['height'].centroid_values

        #------------------------------------------
        # x,y coordinates of vertices of cells that are
        # updated
        #------------------------------------------
        coord_v = self.domain.vertex_coordinates

        if isinstance(self.indices, list):
            if not self.indices:
                return

        if isinstance(self.indices, num.ndarray):
            if self.indices.size == 0:
                return

        #print('init ', type(self.indices), self.indices)

        ids = self.indices
        
        if ids is None:
            self.v_x = coord_v[:,0].reshape((-1,3))
            self.v_y = coord_v[:,1].reshape((-1,3))
        else:
            self.v_x = coord_v[:,0].reshape((-1,3))[ids]
            self.v_y = coord_v[:,1].reshape((-1,3))[ids]

        #------------------------------------------
        # Need to turn off this optimization as it
        # doesn't fixup the relationship between
        # bed and stage vertex values in dry region
        #------------------------------------------
        self.domain.optimise_dry_cells = 0



    def __call__(self):
        """
        Apply value to those triangles defined by indices

        indices == [], don't apply anywhere
        indices is None, apply everywhere
        otherwise apply for the specific indices
        """

        if isinstance(self.indices, list):
            if not self.indices:
                return

        if isinstance(self.indices, num.ndarray):
            if self.indices.size == 0:
                return

        #print('call ', type(self.indices), self.indices)

        #value = self.get_value()
        
        from pprint import pprint
        #print 'value'
        #pprint(value)



        if self.indices is None:

            #--------------------------------------
            # Update centroid values
            #--------------------------------------
            try:
                value = self.get_value(x=self.coord_c[:,0], y=self.coord_c[:,1])
                #print value
                self.quantity_c[:] = num.maximum(self.elev_c, value)
            except ValueError:
                pass

        else:

            #--------------------------------------
            # Update centroid values
            #--------------------------------------
            ids = self.indices
            x = self.coord_c[ids,0]
            y = self.coord_c[ids,1]
            try:
                value = self.get_value(x=x,y=y)
                self.quantity_c[ids] = num.maximum(self.elev_c[ids], value)
            except ValueError:
                pass






