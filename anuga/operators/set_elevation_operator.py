"""
Set elevation operators


"""

from builtins import str
__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"



from anuga.operators.base_operator import Operator
from anuga.operators.set_elevation import Set_elevation

from anuga.config import indent


#===============================================================================
# General Set Elevation Operator
#===============================================================================

class Set_elevation_operator(Operator, Set_elevation):

    
    def __init__(self,
                 domain,
                 elevation=None,
                 region=None,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
        """
Set the elevation in a region (careful to maintain continuitiy of elevation)

indices: None == all triangles, Empty list [] no triangles

elevation can be a function of time.

        """


 

        Set_elevation.__init__(self, domain, elevation, region, indices, polygon, center, radius)

        Operator.__init__(self, domain, description, label, logging, verbose)


    # Use the __call__ method from Set_elevation
    # to set the elevation in this operator
    __call__ = Set_elevation.__call__
    



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

        return 'Set_elevation_operator active at time '+str(self.domain.get_time())



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




        #print "Depreciating this operator, just use Set_elevation_operator with center and radius set"


        assert center is not None
        assert radius is not None



        Set_elevation_operator.__init__(self,
                                    domain,
                                    elevation=elevation,
                                    center=center,
                                    radius=radius,
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



        Set_elevation_operator.__init__(self,
                               domain,
                               elevation=elevation,
                               polygon=polygon,
                               verbose=verbose)



        




