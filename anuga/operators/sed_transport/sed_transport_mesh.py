"""
Mesh handlers for sediment transport operator

"""

import numpy as np
from shallow_water_ext import rotate
from anuga.config import epsilon, g, minimum_allowed_height
import sed_transport_config as st
from anuga.operators.set_quantity_operator import Set_quantity_operator

"""
Updaters
"""

def update_quantity_nonconserved(self,expression,newvals):

    valsInst = self.domain.quantities[expression]
    vals = np.zeros((len(self.elev_v),3))
    
    vals[self.indices] = newvals

    valsInst.set_values(numeric = vals, location = 'vertices')

def update_elevation_stage(self, porosity):
    """Update the vertex values of the quantities to model erosion
    Modified from erosion_operators.Erosion_Operator
    """

    de = self.dzdt * self.domain.timestep
    pe = porosity * de # water lost to pores in bed
    
    elev = np.copy(self.elev_v)
    elev[self.indices] = elev[self.indices] + de
    
    self.q_elev.set_values(elev, location = 'vertices')
    self.q_elev.smooth_vertex_values()
    
    self.depth_v = self.stage_v - self.elev_v
    self.depth = self.depth_v[self.indices]

    self.max_change = np.max(np.absolute(de))
    
    de_c = np.zeros((len(self.stage_v),3))
    de_c[self.indices] = de + pe
    stage_de = np.mean(de_c, axis=1)
    
    stage = self.stage_c + stage_de
    
    set_quantity(self, 'stage', stage)
    
    
def update_momentum(self):
    """Update the vertex values of momentum to model erosion
    Modified from erosion_operators.Erosion_Operator
    """

    ''' the centroids have higher momenti than the vertices - need to update the
    centroids and not the vertices, but don't want to just assign the centroids 
    the average of the centroids around them (which is what happens with when it
    interpolates, reducing the centroid velocities too much)
    
    Try this:
    - Calculate the changes in momentum at the centroids
    - Extrapolate to vertices and edges (can't do this! using:
            self.q_xmom.extrapolate_second_order_and_limit_by_vertex()
            self.q_xmom.interpolate_from_vertices_to_edges()
        increases the velocity at the vertices and edges even when the centroids are decreasing!
    - Need to update vertices and edges by hand
    - The method of averaging values at centroids around each vertex reduces the values too quickly because it averages with dry cells early on!
    
    Try this:
    - calculate and set_value at the vertices
    - this messes up the centroid values, so then set centroid values from pre-calculated numbers
    
    - Check that it's not making the vertices drop too far
    '''
    
    uh = np.zeros(len(self.xmom_c))
    uh = np.sign(self.xmom_c) * \
                       np.absolute(self.duhdt) * self.domain.timestep
    uh[np.absolute(self.xmom_c) < np.absolute(uh)] = 0
    xv = self.xmom_c - uh
    
    vh = np.zeros(len(self.ymom_c))
    vh = np.sign(self.ymom_c) * \
                       np.absolute(self.dvhdt) * self.domain.timestep
    vh[np.absolute(self.ymom_c) < np.absolute(vh)] = 0
    yv = self.ymom_c - vh
    
#     print self.xmom_v[90,:]
#     print self.xmom_c[90]
#     print self.domain.quantities['xmomentum'].edge_values[90,:]
#     print '---'
    
    set_quantity(self, 'xmomentum', xv)
    set_quantity(self, 'ymomentum', yv)

#     print self.xmom_v[90,:]
#     print self.xmom_c[90]
#     print self.domain.quantities['xmomentum'].edge_values[90,:]
#     print 10*'-'
#     print 10*'-'


def set_quantity(self, quantity_name, Q_c):

    Q_v = np.zeros((len(self.domain),3))
    for i in range(self.domain.number_of_nodes):
        L = self.domain.get_triangles_and_vertices_per_node(node=i)
        val = np.mean(Q_c[L[:,0]])
        Q_v[L[:,0],L[:,1]] = val
    
    # Automatically calculates quantity at centroids and edges
    Q = self.domain.quantities[quantity_name]
    
    Q.set_values(Q_v, location = 'vertices')
    Q.centroid_values[:] = Q_c
    

def protect_concentration(self):
    ''' prevents the concentration from going negative'''
    
    newCh = self.q_conc.vertex_values[self.indices] + \
            self.dChdt * self.domain.timestep
    
    newdChdt = self.dChdt
    newdChdt[newCh<0.0] = self.dChdt[newCh<0.0] \
                        - (newCh[newCh<0.0] / self.domain.timestep)
    
    # Alert if concentrations get too high
#     newC = newCh / (self.depth_v[self.indices] + epsilon)
# 
#     newC[self.depth_v[self.indices]<st.min_depth] = 0
# 
#     msg = 'Concentrations are too high! '
#     msg += ' Concentration: '
#     msg += str(newC.max()*100)
#     msg += ' %. Might be a problem with small depths... '
#     msg += str(newC[newC==newC.max()])
#     msg += ' '
#     msg += str(self.depth[newC==newC.max()])
#     assert newC.max() < st.max_conc, msg
    
                        
    self.dChdt = newdChdt
    self.ddot[newCh<0.0] = self.edot[newCh<0.0] - self.dChdt[newCh<0.0]


# def sediment_dispersion(self):
# 
#     """
#     Calculate the volume of sediment that crosses each edge during the
#     timestep due to turbulent diffusion
#     Sediment does not diffuse across the boundaries
#     """
#     if self.verbose:
#         print 'Including sediment dispersion'
# 
#     self.sdiff = np.absolute(0.1 * self.ustar * self.depth)
# #     self.sdiff[self.depth<=st.min_depth] = 0.0
#     
#     # concentration gradients across each edge ( 1/length )
#     neigh = self.domain.neighbours
#     neigh_conc = np.zeros((len(self.q_conc.edge_values),3))
#     neigh_conc[neigh>=0] = self.q_conc.centroid_values[neigh[neigh>=0]] / \
#                            (self.depth_c[neigh[neigh>=0]] + epsilon)
#     
#     conc_c = self.q_conc.centroid_values / \
#              (self.depth_c + epsilon)
#     conc_c = conc_c[np.newaxis,:].transpose()
#     
#     conc_gradient = np.absolute(np.subtract(neigh_conc,conc_c)) / \
#                     (self.neighbour_distance + epsilon)
#     conc_gradient[neigh<0] = 0.0
#     
#     flux = np.zeros((len(self.q_conc.edge_values),3))
#     flux[self.indices] = self.sdiff * conc_gradient[self.indices]
#                          
#     self.sed_vol_e = flux * self.depth_e * \
#                      self.domain.timestep * self.domain.edgelengths


def compute_sed_flux(self):
    
    """
    Change concentrations from topo change
    """
    protect_concentration(self)
    self.conc_v[self.indices] += self.dChdt * self.domain.timestep
    
    self.conc_v[self.conc_v<0.0] = 0.0
    self.conc_v[self.depth_v<=st.min_depth] = 0.0
    
    self.q_conc.set_values(self.conc_v, location = 'vertices')
    
    """
    Calculate sed flux
    """
    
    # Velocity normal to the edges
    
    u_e = self.q_xmom.edge_values / (self.depth_e + epsilon)
    v_e = self.q_ymom.edge_values / (self.depth_e + epsilon)
    
    n1 = self.domain.normals[:,[0,2,4]]
    n2 = self.domain.normals[:,[1,3,5]]
    u_rotated = n1*u_e + n2*v_e
    
    # Volume of sediment in each cell
    sed_vol_c = self.q_conc.centroid_values * self.domain.areas
    
    # Ch at the edges, replaced by boundary values (so Dirichlet_Sed comes in)
    b_indices = - self.domain.neighbours - 1     # indices of boundary edges
    edge_conc = np.copy(self.q_conc.edge_values)
    edge_conc[b_indices>=0] = \
                self.q_conc.boundary_values[b_indices[b_indices>=0]]
                
    self.sed_vol_e = np.zeros((len(self.conc_v),3))
    
    ###
    """
    Extra: Sediment dispersion
    """
#     if self.use_sed_dispersion:
#         sediment_dispersion(self)    
#     ####        
    
    # Volume of sediment going through each edge
    self.sed_vol_e += edge_conc * u_rotated * self.domain.timestep * \
                self.domain.edgelengths

    # Net difference in sediment flux through edges
    vol_diff = np.sum(self.sed_vol_e,axis=1)
    
    # Update volume of sediment in each cell
    sed_vol_c = sed_vol_c - vol_diff
    sed_vol_c[sed_vol_c<0.0] = 0.0

    Ch_c = sed_vol_c / self.domain.areas
    C_c = Ch_c / (self.depth_c + epsilon)
    
    ''' can't use set_quantity(self,...) because need to average the
    fractional concentrations, not Ch'''
    
    # Calculate sediment concentration at vertices
    # as mean of the centroids of the triangles around it
    Ch_v = np.zeros((len(self.domain),3))
    for i in range(self.domain.number_of_nodes):
        L = self.domain.get_triangles_and_vertices_per_node(node=i)
        val = np.mean(C_c[L[:,0]])
        Ch_v[L[:,0],L[:,1]] = val
    Ch_v = Ch_v * self.depth_v
    Ch_v[Ch_v<0] = 0
    
    # Set the new values for concentration at vertices
    # Automatically calculates quantity at centroids and edges
    self.q_conc.set_values(Ch_v, location = 'vertices')
    
    # Update sediment concentration at centroids
    self.q_conc.centroid_values[:] = Ch_c
    
    """
    Recalculate variables in use
    """
    self.conc = self.conc_v[self.indices] / (self.depth_v[self.indices] + epsilon)
    
    
    
    
    








