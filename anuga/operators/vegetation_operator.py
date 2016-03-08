"""
Erosion operators


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
            self.Cd = self.domain.quantities['drag_coefficient'].centroid_values
        except:
            self.Cd = 1.2 # drag coefficient of a cylinder
            print 'Drag coefficient set to default value Cd = 1.2'
            
            
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
            
        self.veg = self.veg_diameter / self.veg_spacing**2
        self.veg[self.depth == 0] = 0
        
        xvel = self.xmom / (self.depth + epsilon)
        
        Fd_x = 0.5 * self.Cd * self.veg * xvel**2
        
        xvel_v = xvel - Fd_x * self.dt
        
        self.domain.quantities['xmomentum'].\
                set_values(xvel_v * self.depth, location = 'centroids')
           
           
                
        yvel = self.ymom / (self.depth + epsilon)        
        
        Fd_y = 0.5 * self.Cd * self.veg * yvel**2        
        
        yvel_v = yvel - Fd_y * self.dt        
        
        self.domain.quantities['ymomentum'].\
                set_values(yvel_v * self.depth, location = 'centroids')
                
            
        self.calculate_diffusivity()
        
        
        
    def calculate_diffusivity(self):
    
        Cb = 0.001
        Cd = 1.2
    
        self.momentum = num.sqrt(self.xmom**2 + self.ymom**2)
        self.velocity = self.momentum / (self.depth + epsilon)
        
        ad = self.veg * self.veg_diameter
    
        # mixing length
        mix_length = num.minimum(((self.veg_diameter - self.depth) / 0.01) * ad +
                        self.depth, self.depth)
        mix_length[ad >= 0.01] = self.veg_diameter[ad >= 0.01]
                        
        # turbulent kinetic energy
        k = ((1 - ad) * Cb + (Cd * ad)**0.66) * self.velocity**2
        
        # total diffusivity
        diffusivity = num.sqrt(k) * mix_length + ad * self.velocity * self.veg_diameter
        
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

        message  = indent + self.label + ': Sed_operator, time '
        message += str(self.get_time())
        return message