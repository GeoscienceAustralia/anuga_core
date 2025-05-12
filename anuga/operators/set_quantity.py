"""
Set elevation operators


"""

__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"

import numpy as num

from anuga import Domain
from anuga import Quantity
from anuga import Region

import anuga.utilities.log as log
from anuga.geometry.polygon import inside_polygon
from anuga.utilities.function_utils import determine_function_type
from anuga import Region
from anuga.config import indent

class Set_quantity(object):
    """
    Helper class to setup calculation of quantity
    associated with a region (defined by indices, polygon or center/radius
    """

    def __init__(self,
                 domain,
                 quantity,
                 value=None,
                 region=None,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 line=None,
                 verbose = False,
                 test_elevation=True,
                 test_stage=True):

        #-----------------------------------------------------
        # Make sure region is actually an instance of a region
        # Otherwise create a new region based on the other 
        # input arguments
        #-----------------------------------------------------
        if isinstance(region,Region):
            region.set_verbose(verbose)
            self.region = region

        else:
            self.region = Region(domain,
                        indices=indices,
                        polygon=polygon,
                        center=center,
                        radius=radius,
                        line=line,
                        verbose=verbose)

        self.set_value(value)
        self.domain = domain
        self.indices = self.region.indices

        #-------------------------------------------
        # Test quantity
        #-------------------------------------------
        self.quantity = quantity
        msg = 'quantity not found in domain'
        assert quantity in domain.quantities, msg

        # FIXME SR: These should be dealt with in this class
        if test_elevation:
            msg ='Use Set_elevation to maintain mass continuity'
            assert quantity != 'elevation', msg
            
        if test_stage:
            msg ='Use Set_stage to maintain non-negative water depth'
            assert quantity != 'stage', msg
        
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
        indices is None, apply everywhere
        otherwise apply for the specific indices
        """

        if self.region.indices is []:
            return



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
                self.quantity_c[:] = value
            except ValueError:
                pass

        else:

            #--------------------------------------
            # Update centroid values
            #--------------------------------------
            rids = self.indices
            x = self.coord_c[rids,0]
            y = self.coord_c[rids,1]
            try:
                value = self.get_value(x=x,y=y)
                self.quantity_c[rids] = value
            except ValueError:
                pass



    def set_value(self, value = None):

        self.value = value
        self.value_type = determine_function_type(value)


        
    def get_value(self, x = None, y = None, t = None):
        """Get value of quantity at time t.
        If t not specified, return value at current domain time
        """

        #from anuga.fit_interpolate.interpolate import Modeltime_too_early, \
        #                                              Modeltime_too_late


        #print 'x,y,t'
        #print x,y,t

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
            #print self.value
            value = float(self.value)

        #except Modeltime_too_early, e:
        #    raise Modeltime_too_early(e)
        #except Modeltime_too_late, e:
        #    raise Modeltime_too_late(e)


        return value







