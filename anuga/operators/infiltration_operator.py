"""
Infiltration operator

"""

import numpy as num

from anuga import Domain
from anuga import Quantity
from anuga.operators.base_operator import Operator


class Infiltration_operator(Operator, object):
    """
    Infiltration operator implements the Green-Ampt equation to calculate
    flow loss to infiltration
    """

    def __init__(self, domain,
                 indices=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
                
        Operator.__init__(self, domain, description, label, logging, verbose)
            
            
        try:
            wet_depth = self.domain.get_quantity('wetting_front_depth')
            
        except:
            Quantity(domain, name='wetting_front_depth', register=True)
            
        # initialize as 0.5 cm to avoid too-fast rates
        self.domain.set_quantity('wetting_front_depth', 0.0005)
        self.wet_depth = self.domain.quantities['wetting_front_depth'].centroid_values 
            
        self.depth = self.domain.quantities['height'].centroid_values   
        
        # values for sand from Rawls and Brakensiek (1993) Hydrology Handbook
        
        self.Ks = 0.00006544444 # Saturated hydraulic conductivity (meters/second)
        self.psi = 0.0495 # Wetting front soil suction head (meters)
        theta_s = 0.5
        theta_i = 0.25
        self.dtheta = theta_s - theta_i


    def __call__(self):
        """
        Apply rate to those triangles defined in indices

        indices == [], then don't apply anywhere
        indices == None, then apply everywhere
        otherwise apply for the specific indices
        """
        
        self.dt = self.get_timestep()
        self.wet_cells = self.depth > 0.001 # mask for flow depth > 1 cm
        
        if len(self.depth[self.wet_cells]) > 1:
                
            self.calculate_infiltration_rate()
            self.calculate_wetting_front_speed()
                

    def calculate_infiltration_rate(self):
        '''
        Calculates the infiltration rate using the Green-Ampt equation
        for all cells with a flow depth > 1 cm
        
        Updates the flow depth by subtracting the loss of flow to infiltration
        '''
        
        self.f = - self.Ks * (self.psi - self.depth[self.wet_cells]
                 - self.wet_depth[self.wet_cells]) / self.wet_depth[self.wet_cells]
        
        
        
        d_depth = self.f * self.dt
    
        self.f[d_depth < 0] = self.depth[self.wet_cells][d_depth < 0] / self.dt
        d_depth[d_depth < 0] = 0

        
        self.depth[self.wet_cells] = self.depth[self.wet_cells] - d_depth
        
        self.domain.quantities['height'].\
            set_values(self.depth, location = 'centroids')
            
      
            
            
    def calculate_wetting_front_speed(self):
        '''
        Calculates the speed of the wetting front due to infiltration
        
        Updates the quantity "wetting_front_depth"
        '''
        
        self.wet_speed = self.f / self.dtheta
        
        d_wet_depth = self.wet_speed * self.dt
        assert num.min(d_wet_depth) >= 0, "Wetting front moving the wrong way!"
        
        
        self.wet_depth[self.wet_cells] = self.wet_depth[self.wet_cells] + d_wet_depth
        
        self.domain.quantities['wetting_front_depth'].\
            set_values(self.wet_depth, location = 'centroids')


    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        return False
        

    def statistics(self):

        message = self.label + ': Infiltration_operator'
        message = message + ' on triangles '+ str(self.indices)
        return message


    def timestepping_statistics(self):
        from anuga import indent

        message  = indent + self.label + ': Infiltration_operator, time '
        message += str(self.get_time())
        return message