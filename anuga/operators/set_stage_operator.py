"""
Set stage operator

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

from anuga.operators.set_quantity_operator import Set_quantity_operator
from anuga.config import indent



class Set_stage_operator(Set_quantity_operator):
    """
    Set the stage over a region
    """

    get_stage = Set_quantity_operator.get_value
    set_stage = Set_quantity_operator.set_value

    def __init__(self,
                 domain,
                 stage=None,
                 region=None,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 line=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        Set_quantity_operator.__init__(self, domain,
                                       quantity = 'stage',
                                       value = stage,
                                       region = region,
                                       indices = indices,
                                       polygon = polygon,
                                       center = center,
                                       radius = radius,
                                       line = line,
                                       description = description,
                                       label = label,
                                       logging = logging,
                                       verbose = verbose,
                                       test_stage=False)



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

        Set_stage_operator.__init__(self,
                                    domain,
                                    stage=stage,
                                    center=center,
                                    radius=radius,
                                    verbose=verbose)

#===============================================================================
# Specific Stage Operators for polygonal region.
#===============================================================================
class Polygonal_set_stage_operator(Set_stage_operator):
    """
    Set stage over a polygonal region
    """

    def __init__(self, domain,
                 stage=0.0,
                 polygon=None,
                 verbose=False):



        Set_stage_operator.__init__(self,
                               domain,
                               stage=stage,
                               polygon=polygon,
                               verbose=verbose)



        




