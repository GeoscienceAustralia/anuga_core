"""
Vegetation operators

"""

import numpy as num

from anuga import Domain
from anuga import Quantity
from anuga.operators.base_operator import Operator

from math import sqrt, log
from anuga.config import epsilon, g

#===============================================================================
# Specific Erosion operator trying to implement bed shear
#===============================================================================
class Vegetation_operator(Operator, object):
    """
    Vegetation operator that applies a drag on the flow due to the presence of veg
    """

    def __init__(self, domain,
                 indices=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
                
        Operator.__init__(self, domain, description, label, logging, verbose)
        
        self.Cd = None
        self.veg = None
        
        try:
            # the value in quantity 'veg_diameter' should be the stem diameter in meters
            self.veg_diameter = self.domain.quantities['veg_diameter'].centroid_values
        except:
            self.veg_diameter = None
            
        try:
            # the value in quantity 'veg_spacing' should be the stem spacing in meters
            self.veg_spacing = self.domain.quantities['veg_spacing'].centroid_values
        except:
            self.veg_spacing = None
            
            
        try:
            diff = self.domain.get_quantity('diffusivity')
        except:
            Quantity(domain, name='diffusivity', register=True)
            
        self.domain.set_use_kinematic_viscosity(True)
            
            
        self.xmom = self.domain.quantities['xmomentum'].centroid_values
        self.ymom = self.domain.quantities['ymomentum'].centroid_values
        self.depth = self.domain.quantities['height'].centroid_values   


    def __call__(self):
        """
        Apply rate to those triangles defined in indices

        indices == [], then don't apply anywhere
        indices == None, then apply everywhere
        otherwise apply for the specific indices
        """
        
        self.dt = self.get_timestep()
        
        if self.veg_diameter is None:
            self.veg_diameter = self.domain.quantities['veg_diameter'].centroid_values
            
        if self.veg_spacing is None:
            self.veg_spacing = self.domain.quantities['veg_spacing'].centroid_values
            
        if self.veg is None:
        
            self.veg = self.veg_diameter / self.veg_spacing**2
            self.ad = self.veg * self.veg_diameter
            
            
        if self.Cd is None:
            self.calculate_drag_coefficient()
            
        self.calculate_diffusivity()
            
    
        xvel = self.xmom / (self.depth + epsilon)
        yvel = self.ymom / (self.depth + epsilon)  
        
        Fd_x = 0.5 * self.Cd * self.veg * xvel**2
        Fd_y = 0.5 * self.Cd * self.veg * yvel**2 
        
        xvel_v = xvel - Fd_x * self.dt
        yvel_v = yvel - Fd_y * self.dt 
        
        self.domain.quantities['xmomentum'].\
            set_values(xvel_v * self.depth, location = 'centroids')
            
        self.domain.quantities['ymomentum'].\
            set_values(yvel_v * self.depth, location = 'centroids')
                


    
    def calculate_drag_coefficient(self):
        '''
        Calculate the drag coefficient Cd as a function of ad using
        the curve fitted to Figure 6 in Nepf (1999)
        '''
        
        self.Cd = (56.11 * self.ad**2
                   - 15.28 * self.ad
                   + 1.3
                   - 0.0005465 / self.ad)
        self.Cd[self.ad < 0.006] = 1.2
        

        
    def calculate_diffusivity(self):
        
        self.momentum = num.sqrt(self.xmom**2 + self.ymom**2)
        self.velocity = self.momentum / (self.depth + epsilon)
    
        # mixing length
        mix_length = num.zeros_like(self.veg_diameter)
        
        # for all cells, the transition from mix_length = depth to a linear
        # form happens at veg_spacing = depth, which is ad = diameter**2 / depth**2
        
        ad_deltaS = self.veg_diameter**2 / (self.depth + epsilon)**2
        
        mix_length_slope = (self.veg_diameter - self.depth) / (0.01 - ad_deltaS)
        
        mix_length = (mix_length_slope * self.ad +
                    (self.depth - mix_length_slope * ad_deltaS))
        
        mix_length[self.veg_spacing > self.depth] = \
                   self.depth[self.veg_spacing > self.depth]
                   
        mix_length[self.ad > 0.01] = self.veg_diameter[self.ad > 0.01]
                        
        # turbulent kinetic energy
        Cb = 0.001
        k = ((1 - self.ad) * Cb + (self.Cd * self.ad)**0.66) * self.velocity**2
        
        # total diffusivity
        diffusivity = (num.sqrt(k)**0.5 * mix_length +
                    self.ad * self.velocity * self.veg_diameter)
        
        self.domain.quantities['diffusivity'].\
                set_values(diffusivity, location = 'centroids')
        



    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        return False
        

    def statistics(self):

        message = self.label + ': Veg_operator'
        message = message + ' on triangles '+ str(self.indices)
        return message


    def timestepping_statistics(self):
        from anuga import indent

        message  = indent + self.label + ': Veg_operator, time '
        message += str(self.get_time())
        return message