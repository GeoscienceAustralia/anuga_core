
__author__="steve"
__date__ ="$11/11/2011 1:52:17 PM$"

from anuga.operators.base_operator import Operator
import numpy as num


class Mannings_operator(Operator):
    """
    Class for setting up a mannings opeator to apply Mannings fricition

    Applying

    d/dt uh =  -g nu^2 uh  sqrt( uh^2 + vh^2 ) / h^{7/3}

    d/dt vh =  -g nu^2 vh  sqrt( uh^2 + vh^2 ) / h^{7/3}


    """

    def __init__(self, domain, verbose=False):
        if verbose: log.critical('Mannings Operator: Beginning Initialisation')


        Operator.__init__(self,domain)

        self.gamma_c  = num.zeros_like(self.stage_c)
        self.height_c = num.zeros_like(self.stage_c)
        self.friction_c = self.domain.quantities['friction'].centroid_values
        self.g = self.domain.g


        self.exp_gamma_max = 0.0
        self.exp_gamma_min = 1.0
        
        if verbose: log.critical('Mannings Operator: Initialisation Done')


    def __call__(self):

        timestep = self.domain.get_timestep()

        self.height_c[:] = self.stage_c - self.elev_c

        self.gamma_c[:] = -self.g * self.friction_c**2 * num.sqrt( self.xmom_c**2 + self.ymom_c**2 )
        self.gamma_c[:] = num.where(self.height_c > 0.0, self.gamma_c/num.power(self.height_c,7.0/3.0),-100.0)

        exp_gamma = num.exp(self.gamma_c*timestep)

        self.exp_gamma_max = max(self.exp_gamma_max,num.max(exp_gamma))
        self.exp_gamma_min = min(self.exp_gamma_min,num.min(exp_gamma))

        #print "Mannings: ",exp_gamma_max,exp_gamma_min


        self.xmom_c[:] = exp_gamma*self.xmom_c
        self.ymom_c[:] = exp_gamma*self.ymom_c



    def parallel_safe(self):
        """
        This operator only works with centroid values of quantities
        and so is parallel safe.
        """
        return True


    def statistics(self):

        message = 'Manning Operator'
        return message


    def timestepping_statistics(self):

        message = '   Manning Operator: Max and Min factors %f %f ' % (self.exp_gamma_max ,self.exp_gamma_min)
        self.exp_gamma_max = 0.0
        self.exp_gamma_min = 1.0

        return message


