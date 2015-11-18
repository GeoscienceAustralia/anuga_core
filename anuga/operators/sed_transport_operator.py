"""
Erosion operators


"""

__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"



import numpy as num


from anuga import Domain
from anuga import Quantity
from anuga.operators.base_operator import Operator
from anuga import Region

from erosion_operators import Erosion_operator
from math import sqrt

import sed_transport_operator_config as st

from anuga.config import epsilon, g



#===============================================================================
# Specific Erosion operator trying to implement bed shear
#===============================================================================
class Sed_transport_operator(Erosion_operator):
    """
    Sed transport operator based on Bed shear erosion operator


    """

    def __init__(self, domain,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
                 

        Erosion_operator.__init__(self,
                 domain,
                 indices=indices,
                 polygon=polygon,
                 center=center,
                 radius=radius,
                 description=description,
                 label=label,
                 logging=logging,
                 verbose=verbose)
                 
                 
        self.D50 = st.D50
        self.porosity = st.porosity

        self.Ke_star = st.Ke_star
        self.R = st.R

        self.criticalshear_star = st.criticalshear_star

        self.c1 = st.c1
        self.c2 = st.c2
        self.mu = st.mu
        self.rousecoeff = st.rousecoeff
        
        self.kappa = st.kappa
        
        
        self.Ke = self.Ke_star * self.D50 * sqrt(self.R * g * self.D50)
        
        self.settlingvelocity = ((self.R * g * self.D50**2) /
                            ((self.c1 * self.mu) +
                            (0.75 * self.c2 * (self.R * g * self.D50**3)**0.5)))
        


    @property
    def grain_size(self):
        """Grain size in m"""
        
        return self.D50


    @grain_size.setter
    def grain_size(self, new_D50):
        """Set the grain size, update the quantities that use it.

        Parameters
        ----------
        new_D50 : float
            New grain size in m.
        """
        self.D50 = new_D50
        
        self.Ke = self.Ke_star * self.D50 * sqrt(self.R * g * self.D50)
        
        self.settlingvelocity = ((self.R * g * self.D50**2) /
                            ((self.c1 * self.mu) +
                            (0.75 * self.c2 * (self.R * g * self.D50**3)**0.5)))
    
        
        



    def update_quantities(self):
        """Update the vertex values of the quantities to model erosion
        """

        
        t = self.get_time()
        self.dt = self.get_timestep()
        
        self.depth = self.stage_c - self.elev_c
        
        self.ind = self.depth > 0.1 # 10 cm
        
        if len(self.ind)>0:
        
            self.conc = self.domain.quantities['concentration'].centroid_values
            
            
            self.momentum = num.sqrt(self.xmom_c[self.ind]**2 + self.ymom_c[self.ind]**2)
        
            self.velocity = self.momentum / (self.depth[self.ind] + epsilon)


            edot = self.erosion()
            ddot = self.deposition()

            dzdt = (ddot - edot) / (1 - self.porosity)
            self.update_bed(dzdt)
        
            dChdt = (edot - ddot)
            self.update_concentration(dChdt)

        updated = True

        return updated
        
        
    def update_concentration(self, dChdt):
    
        # sediment vol already in the water column
        Ch_o = self.conc[self.ind] * self.depth[self.ind]
        
        # sediment vol added or removed from the water column
        Ch_i = dChdt * self.dt
        
        # change in sed vol
        Ch_new = Ch_o + Ch_i
        
        # protect against negative concentrations
        Ch_new[Ch_new < 0] = 0.
        
        new_conc = num.zeros_like(self.conc)
        new_conc[self.ind] = Ch_new / (self.depth[self.ind] + epsilon)
        
        self.domain.quantities['concentration'].\
                set_values(new_conc, location = 'centroids')
        
        
    
        
    def update_bed(self, dzdt):
    
        dz = dzdt * self.dt
        water_loss = self.porosity * dz # water lost to pores in bed
        
        if self.domain.get_using_discontinuous_elevation():
        
        
            self.elev_c[self.ind] = self.elev_c[self.ind] + dz
            self.stage_c[self.ind] = (self.elev_c[self.ind] +
                                        self.depth[self.ind] - water_loss)

        else:
            
            dz = num.vstack((dz,dz,dz)).T
            self.elev_v[self.ind] = self.elev_v[self.ind] + dz
        

        
    def erosion(self):
        
        u_star = (self.velocity * self.kappa /
                    num.log(self.depth[self.ind] / self.D50))

        shear_stress_star = u_star**2 / (g * self.R * self.D50)

        edot = self.Ke * (shear_stress_star - self.criticalshear_star)
        edot[edot<0.0] = 0.0
        
        return edot        
        

    def deposition(self):
        
        ddot = self.rousecoeff * self.conc[self.ind] * self.settlingvelocity
        ddot[ddot<0.0] = 0.0
    
        return ddot

