"""
Set w_uh_vh operator

Constraints: See GPL license in the user guide
Version: 1.0 ($Revision: 7731 $)
"""

from builtins import str
__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"


from anuga import Domain
from anuga import Quantity
import numpy as num
import anuga.utilities.log as log

from anuga.geometry.polygon import inside_polygon
from anuga.utilities.function_utils import determine_function_type
from anuga.utilities.numerical_tools import ensure_numeric

from anuga.operators.base_operator import Operator
from anuga import Region
from anuga.fit_interpolate.interpolate import Modeltime_too_early, \
                                              Modeltime_too_late
from anuga.config import indent



class Set_w_uh_vh_operator(Operator, Region):
    """
    Set the w, uh and vh in a region

    indices: None == all triangles, Empty list [] no triangles

    rate can be a function of time.

    """

    def __init__(self,
                 domain,
                 w_uh_vh=None,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        Operator.__init__(self, domain, description, label, logging, verbose)

        Region.__init__(self, domain,
                        indices=indices,
                        polygon=polygon,
                        center=center,
                        radius=radius,
                        verbose=verbose)

        #print self.indices
        self.set_w_uh_vh(w_uh_vh)
        
        #------------------------------------------
        # Local variables
        #------------------------------------------
        #self.w_uh_vh = w_uh_vh


    def __call__(self):
        """
        Apply w_uh_vh to those triangles defined in indices

        indices == [], then don't apply anywhere
        indices is None, then apply everywhere
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


    def set_w_uh_vh(self, w_uh_vh = None):

        self.w_uh_vh = w_uh_vh

        self.w_uh_vh_type = determine_function_type(w_uh_vh)
        if self.w_uh_vh_type == 'array':
            self.w_uh_vh = ensure_numeric(self.w_uh_vh)
        elif self.w_uh_vh_type == 'scalar':
            self.w_uh_vh = float(self.w_uh_vh)

        
    def get_w_uh_vh(self, x = None, y = None, t = None):
        """Get w_uh_vh value at time t.
        If t not specified, return value at current domain time
        """

        #from anuga.fit_interpolate.interpolate import Modeltime_too_early, \
        #                                              Modeltime_too_late


        if t is None:
            t = self.domain.get_time()

        #print  self.w_uh_vh_type
        
        #try:
        if self.w_uh_vh_type == 't':
            w_uh_vh = self.w_uh_vh(t)
        elif self.w_uh_vh_type == 'x,y':
            w_uh_vh = self.w_uh_vh(x,y)
        elif self.w_uh_vh_type == 'x,y,t':
            w_uh_vh = self.w_uh_vh(x,y,t)
        else:
            w_uh_vh = self.w_uh_vh

        #except Modeltime_too_early, e:
        #    raise Modeltime_too_early(e)
        #except Modeltime_too_late, e:
        #    raise Modeltime_too_late(e)


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







