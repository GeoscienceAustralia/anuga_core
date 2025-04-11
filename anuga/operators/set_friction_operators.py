"""
Set friction operators


"""

from builtins import str
__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"



import numpy as np

from anuga.operators.base_operator import Operator
from anuga import Region
from anuga.config import indent


default_friction_min = 0.01
default_friction_max = 0.035

class Set_depth_friction_operator(Operator):
    """
    Set the friction in a region as a function of water depth
    """

    def __init__(self,
                 domain,
                 friction=lambda h: 0.03,
                 region=None,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        Operator.__init__(self, domain, description, label, logging, verbose)

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
                        verbose=verbose)

        # Region.__init__(self, domain,
        #         indices=indices,
        #         polygon=polygon,
        #         center=center,
        #         radius=radius,
        #         verbose=verbose)

        #------------------------------------------
        # Local variables
        #------------------------------------------

        if friction is None:
            msg = 'Friction function not specified'
            raise ValueError(msg)

        if not callable(friction):
            msg = 'Friction function not callable'
            if type(friction) is float or type(friction) is int:
                friction_fun = lambda h: friction
            else:
                msg += ' (maybe you forgot to use lambda?)'
                raise ValueError(msg)
        else:
            friction_fun = friction

        # vectorize friction_fun
        friction_fun = np.vectorize(friction_fun)

        self.friction = friction_fun
        self.friction_c = self.domain.get_quantity('friction').centroid_values



    def __call__(self):
        """
        Change friction based on depth
        """

        if self.region.indices is []:
            return


        #-----------------------------------------
        # Here is where the important formula is applied
        #----------------------------------------
        if self.region.indices is None:
            height = self.stage_c - self.elev_c
            self.friction_c[:] = self.friction(height)
        else:
            ind = self.region.indices
            height = self.stage_c[ind] - self.elev_c[ind]
            self.friction_c[ind] = self.friction(height)



    def parallel_safe(self):
        """Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True

    def statistics(self):

        message = self.label + ': Set_friction_operator'
        return message


    def timestepping_statistics(self):

        message  = indent + self.label + ': Set_friction_operator, '

        if self.indices is None:
            message  += str(np.min(self.friction_c) ) + ' '
            message  += str(np.max(self.friction_c) )
        else:
            ind = self.indices
            message  += str(np.min(self.friction_c[ind]) ) + ' '
            message  += str(np.max(self.friction_c[ind]) )


        return message


