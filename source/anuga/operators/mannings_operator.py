# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$11/11/2011 1:52:17 PM$"

from anuga.operators.base_operator import Operator



class Mannings_operator(Operator):
    """
    Class for setting up a mannings opeator to apply Mannings fricition

    Applying

    d/dt uh =  uh  sqrt( uh^2 + vh^2 ) / h^{7/3}

    d/dt vh =  vh  sqrt( uh^2 + vh^2 ) / h^{7/3}


    """

    def __init__(self, domain, verbose=False):
        if verbose: log.critical('Mannings Operator: Beginning Initialisation')


        Operator.__init__(self,domain)



        if verbose: log.critical('Mannings Operator: Initialisation Done')


    def __call__(self):

        #timestep = self.domain.get_timestep()
        raise Exception('Need to implement __call__ for your operator')


    def __parallel_safe(self):
        """
        This operator only works with centroid values of quantities
        and so is parallel safe.
        """
        return True


    def statistics(self):

        message = 'Manning Operator '
        return message


    def timestepping_statistics(self):

        message = '    Manning Operator: '
        return message


