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

from anuga import Dirichlet_boundary

import warnings

warnings.filterwarnings('error')


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
        
        self.I_bounds = {'Hmin' : 0.2, 'Hmax' : 2., 'dH' : 0.25,
                 'Umin' : 0.1, 'Umax' : 2.5, 'dU' : 0.25,
                 'Smin' : 0.00001, 'Smax' : 0.0001, 'dS' : 0.00005}
        self.I_flag = True


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
            
        self.dx, self.dy = self.get_dx()
        
        self.initialize_arrays()
        
        self.x = self.domain.quantities['x'].centroid_values
        self.y = self.domain.quantities['y'].centroid_values
        
        
        
    def get_dx(self):
    
        edges = self.domain.get_edge_midpoint_coordinates()
        centroids = self.domain.centroid_coordinates
        
        dx = num.zeros((self.num_cells,))
        dy = num.zeros((self.num_cells,))
        
        for j in range(len(centroids)):
        
            edges_j = edges[3*j:3*j+3,:]
            dx[j] = num.max((num.abs(edges_j[:,0] - centroids[j,0]))*2.)
            dy[j] = num.max((num.abs(edges_j[:,1] - centroids[j,1]))*2.)
        
        return dx, dy
        
        
    def initialize_arrays(self):
    
        self.edot = num.zeros((self.num_cells,))
        self.ddot = num.zeros((self.num_cells,))
        self.u_star = num.zeros((self.num_cells,))
        self.dzdt = num.zeros((self.num_cells,))
        self.dChdt = num.zeros((self.num_cells,))
            
                         
                         


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
        
        self.settlingvelocity = ((self.R * g * self.D50**2) /
                            ((self.c1 * self.nu) +
                            (0.75 * self.c2 * (self.R * g * self.D50**3)**0.5)))
        
        self.A = 5.7e-7
        self.Re = (sqrt(self.R * g * self.grain_size) * self.grain_size) / self.nu
        
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
            
        self.edot[:] = 0
        self.ddot[:] = 0
        self.u_star[:] = 0
        self.dzdt[:] = 0
        self.dChdt[:] = 0
        
        self.conc[self.depth <= 0.01] = 0.

        self.ind = (self.depth > 0.01) & ((self.xmom_c != 0) | (self.ymom_c != 0)) # 5 cm (and moving)
#         self.ind_s = (self.depth > 0.01) & ~self.ind
        self.ind_a = (self.depth > 0.01)
        
#         self.ind_b = (self.x > 90) & (self.x < 95) & (self.y > 555) & (self.y < 560)
        
        self.update_quantities()    
        
        

    def calculate_slope(self):
    
        quant = self.domain.quantities['elevation']
        quant.compute_gradients()
        
        S = num.maximum(num.abs(quant.x_gradient[self.ind]) / self.dx[self.ind],
                 num.abs(quant.y_gradient[self.ind]) / self.dy[self.ind])
        S[S <= 0.000001] = 0.000001
        
        return S
    
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
        
            self.S = self.calculate_slope()
            self.U = self.calculate_velocity()
            
            self.u_star[self.ind] = num.sqrt(g * self.S * self.depth[self.ind])
            
            self.edot[self.ind] = self.erosion()
            self.ddot[self.ind] = self.deposition()
            
            
#         if sum(self.ind_s) > 0:
#         
#             self.ddot[self.ind_s] = self.settlingvelocity * self.conc[self.ind_s]

        self.dzdt = (self.ddot - self.edot) / (1 - self.porosity)
#         self.dzdt[self.bdry_indices] = 0
        self.dChdt = (self.edot - self.ddot)
        
#         if sum(self.ind) > 0: 
#             print self.depth[self.bdry_indices]
#             print self.dzdt[self.bdry_indices]
#             print '-' * 10
        
        
        self.update_concentration(self.dChdt)
            
        self.sediment_flux()
        
        self.update_bed(self.dzdt)



    def update_bed(self, dzdt):
        """
        Updates the elevation of the bed at centroids based on erosion and deposition
        """
        
        self.elev_c[:] = self.elev_c + dzdt * self.dt
        self.domain.set_quantity('elevation', self.elev_c, location='centroids')


        
    def erosion(self):
        """
        Calculates an erosion rate from an excess shear stress formulation
        """
        
        self.calculate_d_star()
        
        Z5 = ((self.u_star[self.ind] / self.settlingvelocity) * self.Re**0.6 * self.S**0.07)**5.
        
        edot = self.settlingvelocity * (self.A * Z5 / (1 + (self.A / 0.3) * Z5))
        edot[edot<0.0] = 0.0
        
        return edot        
        

    def deposition(self):
        """
        Calculates a rate of deposition from the sediment concentration and the
        settling velocity of the particles
        """
       
        ddot = self.d_star * self.conc[self.ind]
        
        
        ddot[ddot<0.0] = 0.0
        
#         print d.max(), ddot.max()
    
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
        self.conc[self.ind_a] = new_sed_vol[self.ind_a] / (self.depth[self.ind_a] * self.areas[self.ind_a])
        
        self.domain.quantities['concentration'].\
                set_values(self.conc, location = 'centroids') 
                

        

    def sediment_flux(self):
        """
        Calculates the flux of sediment between cells based on the flux of water
        to calculate the new concentration at centroids
        
        Assumes that sediment moves at the same speed as the flow
        """
    
        normal_vels = num.zeros((self.num_cells,3))
        
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
                    sed_flux[k,i] =  -1. * num.abs(edge_flux[k,i] * self.inflow_concentration)
#                     print sed_flux[k,:], self.neighbours[k,:], edge_flux[k,:], self.depth_e[k,:]
                    #MP################ removed /2
        
        sed_vol_change = num.sum(-sed_flux, axis=1)
        
       
        sed_vol_in_cell = self.conc * self.depth * self.areas
        new_sed_vol_in_cell = num.maximum(sed_vol_in_cell + sed_vol_change, 0)
        
#         self.conc[:] = 0
        self.conc[self.ind] = (new_sed_vol_in_cell[self.ind] /
                             (self.depth[self.ind] * self.areas[self.ind]))
                             
                            
        assert self.conc.max() < 0.5, 'Max concentration is %d' % self.conc.max()
        
        
        self.domain.quantities['concentration'].\
                set_values(self.conc, location = 'centroids') 
     



    def integrand(self, x, H, U, S):
    
        db = 0.05

        u_star = num.sqrt(g * H * S)
        
#         Cz = U / u_star
#         if Cz > 500:
#             Cz = 500.
# 
#         Kc = 11. * H / num.exp(self.kappa * Cz)

        zo = self.nu / (9. * u_star)
            
        ln = num.log(H * x / zo)

        intv = self.settlingvelocity / (self.kappa * u_star)
        
        integ = (((1 - x) / x)/((1 - db) / db))**intv * ln

        return integ


    def integrand_top(self, x, H, U, S):

        u_star = num.sqrt(g * H * S)

        zo = self.nu / (9. * u_star)
            
        ln = num.log(H * x / zo)
        
        integ = ln

        return integ                                 
                             


    def prepare_d_star(self):
    
        if self.I_flag:
        
            self._H = num.linspace(self.I_bounds['Hmin'],
                                 self.I_bounds['Hmax'],
                                 5)
                                 
            self._U = num.linspace(self.I_bounds['Umin'],
                                 self.I_bounds['Umax'],
                                 5)
                                 
            self._S = num.linspace(self.I_bounds['Smin'],
                                 self.I_bounds['Smax'],
                                 10)
                                 
            self.I_flag = False

        Ix = num.zeros((len(self._H),len(self._U),len(self._S)))

        for nh, h in enumerate(self._H):
            for nu, u in enumerate(self._U):
                for ns, s in enumerate(self._S):

                    I_base = quad(self.integrand, 0.05 , 0.95, args=(h,u,s))
                    I_top = quad(self.integrand_top, 0.05 , 0.95, args=(h,u,s))
                    Ix[nh,nu,ns] = I_top[0] / I_base[0]

        Ix[Ix < 0] = 0
        Ix[Ix > 1] = 1
        
        Hx, Ux = num.meshgrid(self._H,self._U)
        Hs = Hx.flatten()
        Us = Ux.flatten()

        self.ff = []

        for i in range(len(self._S)):
    
            Is = Ix[:,:,i].flatten()
            f = interpolate.SmoothBivariateSpline(Hs,Us,Is)
    
            self.ff.append(f)
                


    def check_d_star_bounds(self, H, U, S):
        
        if (((H.min() < self.I_bounds['Hmin']) and
                (self.I_bounds['Hmin'] != 0.1)) or
           (H.max() >= self.I_bounds['Hmax']) or
           ((U.min() < self.I_bounds['Umin']) and
                (self.I_bounds['Umin'] != 0.05)) or
           (U.max() >= self.I_bounds['Umax']) or
           ((S.min() < self.I_bounds['Smin']) and
                (self.I_bounds['Smin'] != 0.000001)) or
           (S.max() > self.I_bounds['Smax'])):
           
            self.I_flag = True
            
        
        if self.I_flag:
        
            self.I_bounds['Hmin'] = max(0.1, H.min() - self.I_bounds['dH'])
            self.I_bounds['Hmax'] = H.max() + 2. * self.I_bounds['dH']
            self.I_bounds['Umin'] = max(0.05, U.min() - 2. * self.I_bounds['dU'])
            self.I_bounds['Umax'] = U.max() + 2. * self.I_bounds['dU']
            self.I_bounds['Smin'] = max(0.000001, S.min() - 2. * self.I_bounds['dS'])
            self.I_bounds['Smax'] = S.max() + 2. * self.I_bounds['dS']
        
            self.prepare_d_star()
            
                
                

    def calculate_d_star(self):
    
        S_pts = self.S
        U_pts = self.U[self.ind]
        H_pts = self.depth[self.ind]
        
        
        self.check_d_star_bounds(H_pts, U_pts, S_pts)
        
        slope_slice = map(lambda i: num.where(self._S <= i)[0][-1], S_pts.flatten())

        
        self.d_star = num.array(
                      map(lambda f,h,u: f(h,u)[0][0],
                            map(lambda s: self.ff[s], slope_slice), H_pts, U_pts))
                            
        self.d_star[self.d_star > 1] = 1.
                            


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