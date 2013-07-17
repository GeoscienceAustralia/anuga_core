"""
Set friction operators

Constraints: See GPL license in the user guide
Version: 1.0 ($Revision: 7731 $)
"""

__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"




from anuga.operators.base_operator import Operator
from anuga.operators.region import Region


from anuga import indent


default_friction_min = 0.01
default_friction_max = 0.035

class Depth_friction_operator(Operator, Region):
    """
    Set the friction in a region
    """

    def __init__(self,
                 domain,
                 friction_min=default_friction_min,
                 friction_max=default_friction_max,
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

        assert friction_min >= 0.0
        assert friction_min <= friction_max

        self.friction_min = friction_min
        self.friction_max = friction_max
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
            self.friction_c[:] = (self.friction_max - self.friction_min)/(1.0 - height) + self.friction_min
        else:
            ind = self.indices
            height = self.stage_c[ind] - self.elev_c[ind]
            self.friction_c[ind] = (self.friction_max - self.friction_min)/(1.0 - height) + self.friction_min



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
        message  += ' min '+str(self.friction_min)
        message  += ' max '+str(self.friction_max)
        return message



#===============================================================================
# Specific Operator for circular region.
#===============================================================================
class Circular_depth_friction_operator(Depth_friction_operator):
    """
    Set friction over a circular region

    """

    def __init__(self, domain,
                 friction_min=default_friction_min,
                 friction_max=default_friction_max,
                 center=None,
                 radius=None,
                 verbose=False):




        Dpth_friction_operator.__init__(self,
                                    domain,
                                    friction_min=friction_min,
                                    friction_max=friction_max,
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
                 friction_min=default_friction_min,
                 friction_max=default_friction_max,
                 polygon=None,
                 verbose=False):



        Depth_friction_operator.__init__(self,
                               domain,
                               friction_min=friction_min,
                               friction_max=friction_max,
                               polygon=polygon,
                               verbose=verbose)



