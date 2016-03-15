"""
Erosion operators


"""


import numpy as num


from anuga import Domain
from anuga import Quantity
from anuga.operators.base_operator import Operator

from math import sqrt, log

import sed_transport_operator_config as st

from anuga.config import epsilon, g

from anuga import Dirichlet_boundary


#===============================================================================
# Specific Erosion operator trying to implement bed shear
#===============================================================================
class Sed_transport_operator(Operator, object):
    """
    Sed transport operator based on Bed shear erosion operator


    """

    def __init__(self, domain,
                 indices=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
                
        Operator.__init__(self, domain, description, label, logging, verbose)
        
                
        self.normals = self.domain.normals  
        self.neighbours = self.domain.neighbours
        self.edgelengths = self.domain.edgelengths
        
        try:
            self.conc = self.domain.quantities['concentration'].centroid_values
        except:
            self.conc = None
            
        self.depth = self.domain.quantities['height'].centroid_values
        
        self.depth_e = self.domain.quantities['height'].edge_values  
        self.xmom_e = self.domain.quantities['xmomentum'].edge_values
        self.ymom_e = self.domain.quantities['ymomentum'].edge_values
        
        quant = self.domain.quantities['elevation']
        
        self.porosity = st.porosity
        self.Ke_star = st.Ke_star
        self.R = st.R
        self.criticalshear_star = st.criticalshear_star
        self.c1 = st.c1
        self.c2 = st.c2
        self.mu = st.mu
        self.kappa = st.kappa
        self.grain_size = st.D50
        
        self.bdry_indices = None
        self.inflow_concentration = None
        if self.domain.boundary_map:
            self.initialize_inflow_boundary()





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
                            
                            
        self.z_o = self.D50 / 30
        
        self.prepare_d_star()
        
    
    
        
            
    def set_inflow_concentration(self, bdry_conc):
        """ Set the concentration of sediment in the flow
        along Dirichlet boundaries.
        
        Default value for inflow_concentration is the maximum of the
        quantity "concentration" when this operator is instantiated
        """

        self.inflow_concentration = bdry_conc
        
        
        
    def initialize_inflow_boundary(self):
        """Identify the Dirichlet boundaries indices
        """
        
        self.bdry_indices = []
    
        for tag in self.domain.tag_boundary_cells:
            
            B = self.domain.boundary_map[tag]
            
            if isinstance(B, Dirichlet_boundary):
            
                bdry_indices = self.domain.tag_boundary_cells[tag]
                
                for i in bdry_indices:
            
                    self.bdry_indices.append(self.domain.boundary_cells[i])
                
        

    def __call__(self):
        """
        Apply rate to those triangles defined in indices

        indices == [], then don't apply anywhere
        indices == None, then apply everywhere
        otherwise apply for the specific indices
        """
        
        if self.conc is None:
             self.conc = self.domain.quantities['concentration'].centroid_values
        
        if not self.inflow_concentration:
            self.inflow_concentration = self.conc.max()
        
        if not self.bdry_indices:
            self.initialize_inflow_boundary()

        self.ind = (self.depth > 0.05) * (self.xmom_c > 0) # 5 cm (and moving)
        self.update_quantities()    
        

        


    def update_quantities(self):
        """Update the vertex values of the quantities to model erosion
        """

        t = self.get_time()
        self.dt = self.get_timestep()
        
        
        if sum(self.ind) > 0:
            
#             self.velocity = (num.sqrt(self.xmom_c**2 + self.ymom_c**2)
#                              / (self.depth + epsilon))
                             
            
            quant = self.domain.quantities['elevation']
            quant.compute_gradients()

            S = num.maximum(num.abs(quant.x_gradient), num.abs(quant.y_gradient))
            
            self.u_star = num.zeros_like(self.depth)
            
            self.u_star[self.ind] = num.sqrt(g * S[self.ind] * self.depth[self.ind])
            
            
#             self.u_star = (self.velocity * self.kappa /
#                     (num.log(self.depth / self.z_o) - 1))


            edot = self.erosion()
            ddot = self.deposition()

            dzdt = (ddot - edot) / (1 - self.porosity)
        
            dChdt = (edot - ddot)
            
            self.update_concentration(dChdt)
            self.sediment_flux()
            
            self.update_bed(dzdt)



    def update_concentration(self, dChdt):
    
        # sediment vol already in the water column
        sed_vol_in_cell = (self.conc[self.ind] *
                           self.depth[self.ind] * self.areas[self.ind])
        
        # sediment vol added or removed from the water column
        change_sed_vol = dChdt * self.areas[self.ind] * self.dt
        
        # change in sed vol
        new_sed_vol = num.maximum(sed_vol_in_cell + change_sed_vol, 0)
        
        new_conc = num.zeros_like(self.conc)
        new_conc[self.ind] = new_sed_vol / (self.depth[self.ind] * self.areas[self.ind])
        
        self.domain.quantities['concentration'].\
                set_values(new_conc, location = 'centroids') 
                

        

    def sediment_flux(self):
    
        normal_vels = num.zeros_like(self.depth_e)
        
        normal_vels[self.ind,:] = ((self.normals[self.ind,0::2] *
                                    self.xmom_e[self.ind,:] +
                                    self.normals[self.ind, 1::2] *
                                    self.ymom_e[self.ind,:]) /
                                    self.depth_e[self.ind,:])

        
        edge_flux = (self.depth_e * self.edgelengths * normal_vels * self.dt)
        
        # negative fluxes are inwards and must use the concentration of the neighbour
        # positive fluxes are outwards and use the concentration of this cell
        
        neighbour_conc = self.conc[self.neighbours]
        
        sed_flux = edge_flux * self.conc[:,num.newaxis]    
        sed_flux[edge_flux < 0] = (
            edge_flux[edge_flux < 0] * neighbour_conc[edge_flux < 0])
            
        for k in self.bdry_indices:
        
            for i in range(3):
                n = self.neighbours[k,i]
                
                if n < 0:
                    sed_flux[k,i] = edge_flux[k,i] * self.inflow_concentration
        
        sed_vol_change = num.sum(-sed_flux, axis=1)
        
        sed_vol_in_cell = self.conc * self.depth * self.areas
        new_sed_vol_in_cell = num.maximum(sed_vol_in_cell + sed_vol_change, 0)
        
        new_conc = num.zeros_like(self.conc)
        new_conc[self.ind] = (new_sed_vol_in_cell[self.ind] /
                             (self.depth[self.ind] * self.areas[self.ind]))
        
        
        self.domain.quantities['concentration'].\
                set_values(new_conc, location = 'centroids') 

               
  
  
        
    def update_bed(self, dzdt):
    
        dz = num.zeros_like(self.elev_c)
        dz[self.ind] = dzdt * self.dt
        
        new_bed = self.elev_c + dz
        
        self.domain.set_quantity('elevation', new_bed, location='centroids')
        
        

        
    def erosion(self):
    

        shear_stress_star = self.u_star[self.ind]**2 / (g * self.R * self.D50)

        edot = self.Ke * (shear_stress_star - self.criticalshear_star)
        edot[edot<0.0] = 0.0
        
        return edot        
        

    def deposition(self):
    
        self.calculate_d_star()
        
        ddot = (self.d_star[self.ind] * self.conc[self.ind] * self.settlingvelocity)
        ddot[ddot<0.0] = 0.0
    
        return ddot        
        
        

    def prepare_d_star(self):
        """
        Calculates part of the values needed to obtain d* (for deposition)
        following the formulation of Davy and Lague (2009). The equations are
        not very sensitive to flow depth so we save time by calculating things once.
        """
        
        self.d_star_counter = 0
        
        H = 1. # generic depth - not sensitive
    
        a = 0.05 * H
        self.z = num.arange(a, H+a, a)
        
        self.diff_z = num.diff(self.z)

        self.d_dz_u = map(log, self.z / self.z_o)
        self.integral_u = (num.sum((self.d_dz_u[1:] + num.diff(self.d_dz_u)) *
                            self.diff_z)) # numerator of eq 13 in equations.pdf

        self.d_dz_rouse_partial = ((self.z - a) / (H - a)) * (a / self.z)
        
        self.d_star = num.zeros_like(self.depth)
        
        

    def calculate_d_star(self):
        
        if self.d_star_counter % 10 == 0:
        
            wet = num.compress(self.ind, num.arange(len(self.depth)))
            
            for i in wet:
            
                rouse_number = self.settlingvelocity / (self.kappa * self.u_star[i])
        
                d_dz_rouse = (self.d_dz_rouse_partial ** rouse_number) * self.d_dz_u
                # interior of integral in denominator
            
                integral_rouse = (num.sum((d_dz_rouse[1:] + num.diff(d_dz_rouse)) *
                                self.diff_z))
        
                self.d_star[i] = self.integral_u / integral_rouse

            

                
        self.d_star_counter += 1
        
    
        
        


    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        return False
        

    def statistics(self):

        message = self.label + ': Sed_operator'
        message = message + ' on triangles '+ str(self.indices)
        return message


    def timestepping_statistics(self):
        from anuga import indent

        message  = indent + self.label + ': Sed_operator, time '
        message += str(self.get_time())
        return message