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
            # the values in quantity 'vegetation' should be alpha, the ratio
            # of stem diameters and the square of their spacing
            self.veg = self.domain.quantities['vegetation'].centroid_values
        except:
            self.veg = None
            
        try:
            self.Cd = self.domain.quantities['drag_coefficient'].centroid_values
        except:
            self.Cd = 1.2 # drag coefficient of a cylinder
            print 'Drag coefficient set to default value Cd = 1.2'
            
            
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
        
        if self.veg is None:
            self.veg = self.domain.quantities['vegetation'].centroid_values
            
            
        
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