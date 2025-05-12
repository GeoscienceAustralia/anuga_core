"""
Set quantity operators


"""

from builtins import str
__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"




from anuga.operators.base_operator import Operator
from anuga.operators.set_quantity import Set_quantity
from anuga.config import indent


#===============================================================================
# General Set Quantity  Operator
#===============================================================================

class Set_quantity_operator(Operator, Set_quantity):
    """
    Set the elevation in a region (careful to maintain continuitiy of elevation)

    indices: None == all triangles, Empty list [] no triangles

    elevation can be a function of time.

    """


    __call__ = Set_quantity.__call__

    
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
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False,
                 test_stage=True,
                 test_elevation=True):


 

        Set_quantity.__init__(self, 
                              domain, 
                              quantity, 
                              value,
                              region=region,   
                              indices=indices, 
                              polygon=polygon, 
                              center=center, 
                              radius=radius, 
                              line=line,
                              test_stage=test_stage,
                              test_elevation=test_elevation)

        Operator.__init__(self, domain, description, label, logging, verbose)




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

        return 'Set_quantity_operator active at time '+str(self.domain.get_time())



#===============================================================================
# Specific Bed Operators for circular region.
#===============================================================================
class Circular_set_quantity_operator(Set_quantity_operator):
    """
    Set elevation over a circular region

    """

    def __init__(self, domain, quantity,
                 value=0.0,
                 center=None,
                 radius=None,
                 verbose=False):


        assert center is not None
        assert radius is not None


        Set_quantity_operator.__init__(self,
                                       domain,
                                       quantity,
                                       value=value,
                                       center=center,
                                       radius=radius,
                                       verbose=verbose)





#===============================================================================
# Specific Bed Operators for polygonal region.
#===============================================================================
class Polygonal_set_quantity_operator(Set_quantity_operator):
    """
    Add water at certain rate (ms^{-1} = vol/Area/sec) over a
    polygonal region

    rate can be a function of time.

    """

    def __init__(self, domain, quantity,
                 value=0.0,
                 polygon=None,
                 verbose=False):



        Set_elevation_operator.__init__(self,
                               domain,
                               quantity,
                               value=value,
                               polygon=polygon,
                               verbose=verbose)



        




