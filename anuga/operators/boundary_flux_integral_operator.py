"""

Integrate the fluxes through the boundaries

Relies on boundary_flux_sum being computed in
compute_fluxes.



"""

__author__="gareth"
__date__ ="$10/07/2014 $"



import numpy as num
from anuga.operators.base_operator import Operator

class boundary_flux_integral_operator(Operator):
    """
    Simple operator to collect the integral of the boundary fluxes during a run
    """

    def __init__(self,
                 domain,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        Operator.__init__(self, domain, description, label, logging, verbose)

        #------------------------------------------
        # Setup a quantity to store the boundary flux integral 
        #------------------------------------------
        self.boundary_flux_integral=num.array([0.])

        # Alias for domain
        self.domain=domain
        

    def __call__(self):
        """
        Accumulate boundary flux for each timestep
        """

        dt=self.domain.timestep
        ts_method=self.domain.timestepping_method        
        
        if(ts_method=='euler'): 
            self.boundary_flux_integral = self.boundary_flux_integral + dt*self.domain.boundary_flux_sum[0]
        elif(ts_method=='rk2'):
            self.boundary_flux_integral = self.boundary_flux_integral + 0.5*dt*self.domain.boundary_flux_sum[0:2].sum()
        elif(ts_method=='rk3'):
            self.boundary_flux_integral = self.boundary_flux_integral + 1.0/6.0*dt*(self.domain.boundary_flux_sum[0] + self.domain.boundary_flux_sum[1] + 4.0*self.domain.boundary_flux_sum[2])
        else:
            raise Exception('Cannot compute boundary flux integral with this timestepping method')
     
        # Zero the boundary_flux_sum 
        self.domain.boundary_flux_sum[:]=0.

    def parallel_safe(self):
        """Operator is applied independently on each parallel domain

        """
        return True

    def statistics(self):

        message = self.label + ': Boundary_flux_integral operator'
        return message


    def timestepping_statistics(self):
        from anuga import indent

        message  = indent + self.label + ': Integrating the boundary flux'
        return message

