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
                            
        self.normals = self.domain.normals
        
        


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
            self.x = self.domain.quantities['x'].centroid_values
            self.y = self.domain.quantities['y'].centroid_values
            
            self.conc[1250] = 0.1
            
#             print 'conc', self.conc[1248:1253], self.conc[1328]
#             print 'depth', self.depth[5:10]
#             print 'z', self.elev_c[5:10]
            
            
            self.momentum = num.sqrt(self.xmom_c[self.ind]**2 + self.ymom_c[self.ind]**2)
        
            self.velocity = self.momentum / (self.depth[self.ind] + epsilon)


            edot = self.erosion()
            ddot = self.deposition()

            dzdt = (ddot - edot) / (1 - self.porosity)
            self.update_bed(dzdt)
        
            dChdt = (edot - ddot)
            
#             self.update_concentration(dChdt)
            
#             print 'conc', self.conc[1248:1253], self.conc[1328]
            
#             print 'Edot', edot[5:10]
#             print 'Ddot', ddot[5:10]
#             print 'dChdt', dChdt[5:10]
#             print 'dzdt', dzdt[5:10]
#             print 'z', self.elev_c[5:10]
            
            self.sediment_flux()
            
#             print 'conc', self.conc[1248:1253], self.conc[1328]
#             print 'x', self.x[1248:1253], self.x[1328]
#             print 'y', self.y[1248:1253], self.y[1328]
#             print self.conc.max()
#             print self.depth[self.conc == self.conc.max()]
#             print self.elev_c[self.conc == self.conc.max()]
#             print (self.x[self.conc == self.conc.max()], self.y[self.conc == self.conc.max()])
            print '-'*10
            

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
                
#         print 'Newconc', new_conc[5:10]




                
    def sediment_flux_bad(self):
    
        self.depth = self.stage_c - self.elev_c
    
        # Average the velocities and depths at all edges

        N = len(self.domain)

        neighbors = self.domain.neighbours
        neighbor_edges = self.domain.neighbour_edges
        
        areas = self.domain.areas
        edgelengths = self.domain.edgelengths

        depth_neighbors = self.depth[neighbors]
        depth_edges = (self.depth[:,num.newaxis] + depth_neighbors)/2
        
        xmom_edges = (self.xmom_c[:,num.newaxis] + self.xmom_c[neighbors])/2
        ymom_edges = (self.ymom_c[:,num.newaxis] + self.ymom_c[neighbors])/2
        
        
        sed_volume = self.conc * self.depth * areas
        
        
        
        
        n1 = self.domain.normals[:,[0,2,4]]
        n2 = self.domain.normals[:,[1,3,5]]
        mom_rotated = n1 * xmom_edges + n2 * ymom_edges
        
        mom_rotated[depth_neighbors < 0.1] = 0.0
        
        
        edge_areas = edgelengths * depth_edges
        
        flux = num.zeros(1, num.float)
        
        
        for k in range(N):
        
            if self.depth[k] > 0.1:
        
                q_edges = mom_rotated[k] * edgelengths[k] * self.dt
                conc_neigh = self.conc[neighbors[k]]
            
                qs_edges = num.where(q_edges>0, self.conc[k], conc_neigh) * q_edges
                
                flux = -num.sum(qs_edges)
                
                sed_volume[k] = num.max([sed_volume[k] + flux, 0.0])
                
        
        new_conc = sed_volume / (self.depth + epsilon) / areas

        self.domain.quantities['concentration'].\
                set_values(new_conc, location = 'centroids')   

        

    def sediment_flux(self):
    
        self.depth = self.stage_c - self.elev_c
        
        
        N = len(self.domain)   
        neighbours = self.domain.neighbours
        
        neighbour_edges = self.domain.neighbour_edges
        normals = self.domain.normals   
        
        areas = self.domain.areas
        edgelengths = self.domain.edgelengths

        xvel = self.xmom_c / (self.depth + epsilon)
        yvel = self.ymom_c / (self.depth + epsilon)
        
        sed_vol_in_cell = self.conc * self.depth * areas
        new_conc = num.zeros_like(self.conc)
        
        stage_bdry = self.domain.quantities['stage'].boundary_values
        elev_bdry = self.domain.quantities['elevation'].boundary_values
        xmom_bdry = self.domain.quantities['xmomentum'].boundary_values
        ymom_bdry = self.domain.quantities['ymomentum'].boundary_values
        
        depth_bdry = stage_bdry - elev_bdry
        xvel_bdry = xmom_bdry / (depth_bdry + epsilon)
        yvel_bdry = ymom_bdry / (depth_bdry + epsilon)
        
        conc_bdry = self.domain.quantities['concentration'].boundary_values
        
        sed_vol_in_cell = self.conc * self.depth * areas
        sed_vol_in_cell_diff = num.zeros_like(sed_vol_in_cell)
        print sed_vol_in_cell.max()
        print sed_vol_in_cell[1250]
        
        for k in range(N):

            depth_l = self.depth[k]
            
            flux_sed = 0.
            
            if depth_l > 0.1:
            
                xvel_l = xvel[k]
                yvel_l = yvel[k]
                conc_l = self.conc[k]

                for i in range(3):

                    n = neighbours[k,i]
                    
                    if n < 0:
                    
                        m = -n-1 #Convert neg flag to index
                        depth_r = depth_bdry[m]
                        xvel_r = xvel_bdry[m]
                        yvel_r = yvel_bdry[m]
                        conc_r = conc_bdry[m]
                        
                        
                    else:
                        depth_r = self.depth[n]
                        xvel_r = xvel[n]
                        yvel_r = yvel[n]
                        conc_r = self.conc[n]
                        


                    depth_edge = (depth_l + depth_r) / 2
                    xvel_edge = (xvel_l + xvel_r) / 2
                    yvel_edge = (yvel_l + yvel_r) / 2
                    conc_edge = (conc_l + conc_r) / 2
                
                
                    normal = self.normals[k, 2*i:2*i+2]
                    edge_vel = normal[0] * xvel_edge + normal[1] * yvel_edge
                
                    edge_area = edgelengths[k,i] * depth_edge
                    edge_discharge = edge_vel * edge_area
                    edge_volume_water = edge_discharge * self.dt
                    
                    print depth_l, depth_r
                    print conc_l, conc_r
                    print edge_volume_water, depth_l * areas[k]
                    print '-' * 4
        
         
        
        
        
                
    def sediment_flux_bad_2(self):
    
        self.depth = self.stage_c - self.elev_c
    
        # Average the velocities and depths at all edges

        N = len(self.domain)
        

        neighbours = self.domain.neighbours
        
        neighbour_edges = self.domain.neighbour_edges
        normals = self.domain.normals
        

        areas = self.domain.areas
        edgelengths = self.domain.edgelengths

        xvel = self.xmom_c / (self.depth + epsilon)
        yvel = self.ymom_c / (self.depth + epsilon)

        stage_bdry = self.domain.quantities['stage'].boundary_values
        elev_bdry = self.domain.quantities['elevation'].boundary_values
        xmom_bdry = self.domain.quantities['xmomentum'].boundary_values
        ymom_bdry = self.domain.quantities['ymomentum'].boundary_values
        
        depth_bdry = stage_bdry - elev_bdry
        xvel_bdry = xmom_bdry / (depth_bdry + epsilon)
        yvel_bdry = ymom_bdry / (depth_bdry + epsilon)
        
        conc_bdry = self.domain.quantities['concentration'].boundary_values
        
        bdry_indices = self.domain.boundary_cells
        
        
        
#         conc_bdry[bdry_indices < 100] = 0.01
# #         conc_bdry[bdry_indices >= 60] = 0.0

        if self.conc[1250]>0.0:
            print self.conc[1250], self.depth[1250]


        
        sed_vol_in_cell = self.conc * self.depth * areas
        new_conc = num.zeros(N, num.float)
        

        #Loop
        for k in range(N):

            depth_l = self.depth[k]
            xvel_l = xvel[k]
            yvel_l = yvel[k]
            conc_l = self.conc[k]
            
            flux_sed = 0.
            
            if depth_l > 0.1:

                for i in range(3):

                    n = neighbours[k,i]
                    if n < 0:
                    
                        m = -n-1 #Convert neg flag to index
                        depth_r = depth_bdry[m]
                        xvel_r = xvel_bdry[m]
                        yvel_r = yvel_bdry[m]
                        conc_r = conc_bdry[m]
                        
                        
                    else:
                        depth_r = self.depth[n]
                        xvel_r = xvel[n]
                        yvel_r = yvel[n]
                        conc_r = self.conc[n]
                        
                    if depth_r >= 0.1:
                    
                    
                        depth_edge = (depth_l + depth_r) / 2
                        xvel_edge = (xvel_l + xvel_r) / 2
                        yvel_edge = (yvel_l + yvel_r) / 2
                        conc_edge = (conc_l + conc_r) / 2
                    
                    
                        normal = self.normals[k, 2*i:2*i+2]
                        edge_vel = normal[0] * xvel_edge + normal[1] * yvel_edge
                    
                    
                    
                        edge_area = edgelengths[k,i] * depth_edge
                        edge_discharge = edge_vel * edge_area / 2
                    
                        edge_volume_water = edge_discharge * self.dt
                    
                        if edge_vel > 0:
                    
                            edge_volume_sed = edge_volume_water * conc_l
                        
                        if edge_vel <= 0:
                    
                            edge_volume_sed = edge_volume_water * conc_r
                    
                    
                        if k == 1328:
                    
                            print depth_l, conc_l
                            print n
                            print depth_r, conc_r
    #                         print depth_edge, conc_edge
    #                         print edge_vel, edge_discharge
                            print edge_volume_sed
                        
                        
                        flux_sed -= edge_volume_sed

                        
                sed_vol_in_cell[k] += flux_sed

                new_conc[k] = sed_vol_in_cell[k] / (depth_l * areas[k])
                
                if k == 1328:
                    print '*'*5
                    print flux_sed, sed_vol_in_cell[k]
                    print new_conc[k]
                
                
                
        new_conc[new_conc < 0.] = 0.
        
        self.domain.quantities['concentration'].\
                set_values(new_conc, location = 'centroids')
             
#         conc_bdry = self.domain.quantities['concentration'].boundary_values  
#         print 'Max conc_bdry', conc_bdry.max()

        print self.conc[1250], self.conc[1328]

        
        
        
        
    
        
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

