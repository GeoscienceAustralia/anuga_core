"""
Set elevation operators


"""

__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"

import numpy as num

from anuga import Domain
from anuga import Quantity

import anuga.utilities.log as log
from anuga.geometry.polygon import inside_polygon
from anuga.utilities.function_utils import determine_function_type
from anuga.operators.region import Region
from anuga import indent

class Set_quantity(Region):
    """
    Helper class to setup calculation of quantity
    associated with a region (defined by indices, polygon or center/radius
    """

    def __init__(self,
                 domain,
                 quantity,
                 value=None,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 verbose = False,
                 test_elevation=True):


        Region.__init__(self, domain,
                        indices=indices,
                        polygon=polygon,
                        center=center,
                        radius=radius,
                        verbose=verbose)
        

        self.set_value(value)

        #-------------------------------------------
        # Test quantity
        #-------------------------------------------
        self.quantity = quantity
        msg = 'quantity not found in domain'
        assert quantity in domain.quantities, msg

        if test_elevation:
            msg ='Use Set_elevation to maintain continuity'
            assert quantity is not 'elevation'
        
        #-------------------------------------------
        # Useful quantity alias
        #------------------------------------------
        self.quantity_c = self.domain.quantities[quantity].centroid_values

        self.coord_c = self.domain.centroid_coordinates
        self.areas = self.domain.areas


    def __call__(self):
        """
        Apply value to those triangles defined by indices

        indices == [], don't apply anywhere
        indices == None, apply everywhere
        otherwise apply for the specific indices
        """

        if self.indices is []:
            return

        value = self.get_value()

        if value is None:
            return


        if self.indices is None:

            #--------------------------------------
            # Update all three vertices for each cell
            #--------------------------------------
            self.quantity_c[:] = self.get_value(self.coord_c[:,0], self.coord_c[:,1])

        else:

            #--------------------------------------
            # Update all three vertices for each cell
            #--------------------------------------
            ids = self.indices
            x = self.coord_c[ids,0]
            y = self.coord_c[ids,1]
            self.quantity_c[ids] = self.get_value(x,y)



    def set_value(self, value = None):

        self.value = value
        self.value_type = determine_function_type(value)


        
    def get_value(self, x = None, y = None, t = None):
        """Get value of quantity at time t.
        If t not specified, return value at current domain time
        """

        #from anuga.fit_interpolate.interpolate import Modeltime_too_early, \
        #                                              Modeltime_too_late


        if t is None:
            t = self.domain.get_time()

        #try:
        if self.value_type == 't':
            value = self.value(t)
        elif self.value_type == 'x,y':
            value = self.value(x,y)
        elif self.value_type == 'x,y,t':
            value = self.value(x,y,t)
        else:
            value = float(self.value)

        #except Modeltime_too_early, e:
        #    raise Modeltime_too_early(e)
        #except Modeltime_too_late, e:
        #    raise Modeltime_too_late(e)


        return value







