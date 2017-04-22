"""
Erosion operators
"""

import numpy as num
from scipy.integrate import quad
from scipy import interpolate

from anuga import Domain
from anuga import Quantity
from anuga.operators.base_operator import Operator

from math import sqrt, log
from anuga.config import epsilon, g

from anuga import Dirichlet_boundary, Reflective_boundary
# 
# import warnings
# 
# warnings.filterwarnings('error')


#===============================================================================
# Erosion and sed transport operator
#===============================================================================
class Sed_transport_operator(Operator, object):
    """
    Sed transport operator
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
        self.x = self.domain.quantities['x'].centroid_values
        self.y = self.domain.quantities['y'].centroid_values
        
        self.depth_e = self.domain.quantities['height'].edge_values  
        self.xmom_e = self.domain.quantities['xmomentum'].edge_values
        self.ymom_e = self.domain.quantities['ymomentum'].edge_values

        self.porosity = 0.3
        self.c1 = 18.
        self.c2 = 0.4
        self.nu = 1.0e-6
        self.kappa = 0.408
        self.rho_s = 2650.
        self.rho_w = 1000.
        self.R = (self.rho_s - self.rho_w) / self.rho_w
        
        self.grain_size = 0.00013 # from Griffin et al 2010
        
        self.num_cells = len(self.depth)
        
        self.bdry_indices = None
        self.inflow_concentration = None
        
        if self.domain.boundary_map:
            self.initialize_inflow_boundary()    
            
        self.dx = self.get_dx()
        
        self.initialize_arrays()
        
        
        
    def get_dx(self):
    
        edges = self.domain.get_edge_midpoint_coordinates()
        centroids = self.domain.centroid_coordinates
        
        dx = num.zeros((self.num_cells,))
        dy = num.zeros((self.num_cells,))
        
        for j in range(len(centroids)):
        
            edges_j = edges[3*j:3*j+3,:]
            dx[j] = num.max((num.abs(edges_j[:,0] - centroids[j,0]))*2.)
            dy[j] = num.max((num.abs(edges_j[:,1] - centroids[j,1]))*2.)
        
        return (dx + dy).mean()
        
        
        
    def initialize_arrays(self):
        """
        Initialize sed transport arrays and energy quantity
        """
    
        Quantity(self.domain, name='energy', register=True)
        
        self.edot = num.zeros((self.num_cells,))
        self.ddot = num.zeros((self.num_cells,))
        self.u_star = num.zeros((self.num_cells,))
        self.dzdt = num.zeros((self.num_cells,))
        self.dChdt = num.zeros((self.num_cells,))
        self.S = num.zeros((self.num_cells,))


    @property
    def grain_size(self):
        """Grain size in m"""
        return self.D50


    @grain_size.setter
    def grain_size(self, new_D50):
        """
        Set the grain size, update the quantities that use it.

        Parameters
        ----------
        new_D50 : float
            New grain size in m.
        """
        
        self.D50 = new_D50
        
        self.settlingvelocity = ((self.R * g * self.D50**2) /
                            ((self.c1 * self.nu) +
                            (0.75 * self.c2 * (self.R * g * self.D50**3)**0.5)))
        
        self.A = 5.7e-7
        self.Re = (sqrt(self.R * g * self.grain_size) * self.grain_size) / self.nu
        
            
            
    def set_inflow_concentration(self, bdry_conc):
        """
        Set the concentration of sediment in the flow
        along Dirichlet boundaries.
        
        Default value for inflow_concentration is the maximum of the
        quantity "concentration" when this operator is instantiated
        """

        self.inflow_concentration = bdry_conc
        
        
        
    def initialize_inflow_boundary(self):
        """
        Identify the Dirichlet boundaries indices
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
            
        self.edot[:] = 0
        self.ddot[:] = 0
        self.u_star[:] = 0
        self.dzdt[:] = 0
        self.dChdt[:] = 0
        self.conc[self.depth <= 0.01] = 0.

        self.ind = (self.depth > 0.1)
        
        self.update_quantities()    
        
        

    def calculate_energy_slope(self):
        """
        Calculate the energy slope of the flow and store
        as quantity "energy"
        """
    
        e = self.elev_c + self.depth + self.U**2 / (2. * g)
    
        self.domain.set_quantity('energy', e, location='centroids')
    
        quant = self.domain.quantities['energy']
        quant.compute_gradients()

        self.S = num.sqrt(quant.x_gradient**2 + quant.y_gradient**2) / 10
                 
        self.S[self.S <= 0.000001] = 0.000001
        self.S[self.depth < 0.01] = 0.000001
        
        return self.S
    
    
    
    def calculate_velocity(self):
        """
        Calcultes the magnitude of the flow velocity using flux limiting, as described
        in section 16.5 of the ANUGA docs
        """
        
        ho = 1e-6
        
        xvel = (self.xmom_c * self.depth) / (self.depth**2 + ho)
        yvel = (self.ymom_c * self.depth) / (self.depth**2 + ho)
        
        U = num.sqrt(xvel**2. + yvel**2)
        
        return U



    def update_quantities(self):
        """
        Calculates erosion and deposition, moves sediment, and updates the centroid values
        of elevation and concentration to model topographic change
        """

        t = self.get_time()
        self.dt = self.get_timestep()
        
        if sum(self.ind) > 0: 
        
            # self.S = self.calculate_slope()
            self.U = self.calculate_velocity()
            self.S = self.calculate_energy_slope()
            
            self.u_star[self.ind] = num.sqrt(g * self.S[self.ind] * self.depth[self.ind])
            
            self.edot[self.ind] = self.erosion()
            self.ddot[self.ind] = self.deposition()
            
            self.edot[self.edot > 0.001] = 0.001
            self.ddot[self.ddot > 0.001] = 0.001
            
        self.dChdt = (self.edot - self.ddot)
        self.dzdt = -1. * self.dChdt / (1 - self.porosity)

        self.update_concentration(self.dChdt) 

        self.sediment_flux()
        
        self.update_bed(self.dzdt)


    def update_bed(self, dzdt):
        """
        Updates the elevation of the bed at centroids based on erosion and deposition
        """
        
        self.elev_c[:] = self.elev_c + dzdt * self.dt
        self.domain.set_quantity('elevation', self.elev_c, location='centroids')
        
        self.depth[:] = self.depth - dzdt * self.dt
        self.domain.set_quantity('height', self.depth, location='centroids')



    def erosion(self):
        """
        Calculates the erosion rate as an excess shear stress
        """

        tau_crit = 0.103
        self.Ke = 0.2e-6 / tau_crit**(0.5) * 10.
        
        shear_stress = self.rho_w * self.u_star[self.ind]**2
        
        edot = self.Ke * (shear_stress - tau_crit)
        edot[edot<0.0] = 0.0
        
        return edot     


    def deposition(self):
        """
        Calculate the deposition rate following Davy and Lague (2004)
        """
        
        Z = self.settlingvelocity / (self.kappa * self.u_star[self.ind])
        
        z = num.array([ -2.83615213e+00,   5.74715420e+01,  -3.91878963e+02,
         1.14269529e+03,  -1.45387728e+03,   1.15735710e+03,
        -2.97363091e+02,   3.78936399e+01,   5.46546143e-01])
        p = num.poly1d(z)
        
        self.d_star = p(Z)
        self.d_star[Z>4] = -20000 + 10000 * Z[Z > 4]
        
        ddot = self.d_star * self.conc[self.ind] * self.settlingvelocity
        ddot[ddot<0.0] = 0.0
    
        return ddot           


    def update_concentration(self, dChdt):
        """
        Updates the centroid values of concentration based on the volume of material
        that was added or removed from the water column through erosion and deposition
        """
    
        # sediment vol already in the water column
        sed_vol_in_cell = (self.conc *
                           self.depth * self.areas)
        
        # sediment vol added or removed from the water column
        change_sed_vol = dChdt * self.areas * self.dt
        
        # change in sed vol
        new_sed_vol = num.maximum(sed_vol_in_cell + change_sed_vol, 0)
        
        self.conc[:] = 0.
        self.conc[self.ind] = new_sed_vol[self.ind] / (self.depth[self.ind] * self.areas[self.ind])
        

    def sediment_flux(self):
        """
        Calculates the flux of sediment between cells based on the flux of water
        to calculate the new concentration at centroids
        
        Assumes that sediment moves at the same speed as the flow
        
        Negative edge fluxes are inwards and must use the concentration of the neighbour
        Positive edge fluxes are outwards and use the concentration of the cell
        """
    
        normal_vels = num.zeros((self.num_cells,3))
        
        normal_vels[self.ind,:] = ((self.normals[self.ind,0::2] *
                                    self.xmom_e[self.ind,:] +
                                    self.normals[self.ind, 1::2] *
                                    self.ymom_e[self.ind,:]) /
                                    self.depth_e[self.ind,:])

        
        edge_flux = (self.depth_e * self.edgelengths * normal_vels * self.dt)
        
        neighbour_conc = self.conc[self.neighbours]
        
        sed_flux = edge_flux * self.conc[:,num.newaxis]
        
        sed_flux[edge_flux < 0] = (
            edge_flux[edge_flux < 0] * neighbour_conc[edge_flux < 0])
            
            
        for k in self.bdry_indices:
            for i in range(3):

                n = self.neighbours[k,i]
    
                if n < 0:
                    sed_flux[k,i] =  edge_flux[k,i] * self.inflow_concentration
        
        
        sed_vol_change = num.sum(-sed_flux, axis=1)
        
        sed_vol_in_cell = self.conc * self.depth * self.areas
        new_sed_vol_in_cell = num.maximum(sed_vol_in_cell + sed_vol_change, 0)
        
        self.conc[:] = 0.
        self.conc[self.ind] = (new_sed_vol_in_cell[self.ind] /
                             (self.depth[self.ind] * self.areas[self.ind]))

        self.conc[self.conc > 0.5] = 0.5
        
        self.domain.quantities['concentration'].\
                set_values(self.conc, location = 'centroids') 
     



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