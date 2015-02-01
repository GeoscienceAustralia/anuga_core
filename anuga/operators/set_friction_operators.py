"""
Set friction operators

Constraints: See GPL license in the user guide
Version: 1.0 ($Revision: 7731 $)
"""

__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"



import numpy

from anuga.operators.base_operator import Operator
from anuga import Region
from anuga.config import indent


default_friction_min = 0.01
default_friction_max = 0.035

class Depth_friction_operator(Operator, Region):
    """
    Set the friction in a region
    """

    def __init__(self,
                 domain,
                 friction=lambda h: 0.03,
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

        #------------------------------------------
        # Local variables
        #------------------------------------------

        self.friction = friction
        self.friction_c = self.domain.get_quantity('friction').centroid_values



    def __call__(self):
        """
        Change friction based on depth
        """

        if self.indices is []:
            return


        #-----------------------------------------
        # Here is where the important formula is applied
        #----------------------------------------
        if self.indices is None:
            height = self.stage_c - self.elev_c
            self.friction_c[:] = self.friction(height)
        else:
            ind = self.indices
            height = self.stage_c[ind] - self.elev_c[ind]
            self.friction_c[ind] = self.friction(height)



    def parallel_safe(self):
        """Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True

    def statistics(self):

        message = self.label + ': Set_depth_friction_operator'
        return message


    def timestepping_statistics(self):

        message  = indent + self.label + ': Set_depth_friction, '

        if self.indices is None:
            message  += str(numpy.min(self.friction_c) ) + ' '
            message  += str(numpy.max(self.friction_c) )
        else:
            ind = self.indices
            message  += str(numpy.min(self.friction_c[ind]) ) + ' '
            message  += str(numpy.max(self.friction_c[ind]) )


        return message



#===============================================================================
# Specific Operator for circular region.
#===============================================================================
class Circular_depth_friction_operator(Depth_friction_operator):
    """
    Set friction over a circular region

    """

    def __init__(self, domain,
                 friction,
                 center=None,
                 radius=None,
                 verbose=False):




        Depth_friction_operator.__init__(self,
                                    domain,
                                    friction=friction,
                                    center=center,
                                    radius=radius,
                                    verbose=verbose)


#===============================================================================
# Specific Operator for polygonal region.
#===============================================================================
class Polygonal_depth_friction_operator(Depth_friction_operator):
    """
    Set Friction over a polygon

    """

    def __init__(self, domain,
                 friction,
                 polygon=None,
                 verbose=False):



        Depth_friction_operator.__init__(self,
                               domain,
                               friction=friction,
                               polygon=polygon,
                               verbose=verbose)



