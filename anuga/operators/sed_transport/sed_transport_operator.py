"""
Sediment transport and vegetation operators
"""

__author__="mariela"
__date__ ="$23/05/2014 01:33:00 PM$"


import sys
from anuga import Domain
from anuga import Quantity
from anuga.operators.base_operator import Operator
# from anuga import Region
import numpy as np
from anuga.config import epsilon, g
from anuga.abstract_2d_finite_volumes.general_mesh import General_mesh

from anuga.shallow_water.shallow_water_domain import manning_friction_implicit

import sed_transport_mesh as st_mesh
import sed_transport_utils as st_util
import sed_transport_config as st

class sv_common(Operator):
    """
    Common functions for sed transport and vegetation drag operators.
    Separating them allows for easy switching between processes
    """

    def __init__(self,
                 domain,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):

        Operator.__init__(self, domain, description, label, logging, verbose)

        """
        Switches for optional processes
        """     
        if self.use_turbulence:
            
            try:
                diff = self.domain.get_quantity('diffusivity')
                self.domain.set_use_kinematic_viscosity(True)
            except KeyError:
                msg = 'Quantity diffusivity not yet created. '
                msg += 'Initialize when creating domain or set '
                msg += 'turbulence as False. '
                msg += 'Continuing without turbulence.'
                print msg
                self.use_turbulence = False
            

        """
        Call quantities
        """
        self.q_elev = self.domain.quantities['elevation']
        self.q_stage = self.domain.quantities['stage']
        self.q_xmom=self.domain.quantities['xmomentum']
        self.q_ymom=self.domain.quantities['ymomentum']

        self.elev_v  = self.q_elev.vertex_values

        self.stage_v  = self.q_stage.vertex_values
        
        self.xmom_v = self.q_xmom.vertex_values
        self.ymom_v = self.q_ymom.vertex_values
        
        
        self.elev_c  = self.q_elev.centroid_values

        self.stage_c  = self.q_stage.centroid_values
        
        self.xmom_c = self.q_xmom.centroid_values
        self.ymom_c = self.q_ymom.centroid_values

        #-----------------------------------------
        # Some extras
        #-----------------------------------------
        self.max_change = 0
        self.domain.optimise_dry_cells = 0
        
        # to allow for one-equation turbulence calculations
        self.stem_diameter = 0 * self.elev_v
        self.stem_spacing = 0 * self.elev_v
        self.ad = 0 * self.elev_v
    



    def init_sed(self):

        self.sd_operator = None

        """
        Create quantities
        """
        # Create the quantity for concentration, set to 0.0 everywhere.
        # Can be overwritten by a value set in the input file
        
        
        try:
            self.domain.get_quantity('concentration')
        except KeyError:
            msg = 'Error! Quantity concentration not yet created. '
            msg += 'Initialize when creating domain '
            msg += 'as an evolved quantity. Cannot continue!'
            sys.exit(msg)

        """
        Switches for optional processes
        """    
#         if self.use_sed_dispersion:
#             self.neighbour_distance = st_util.get_distance_centroids(self)
            
           
            
#             try:
# #                 sv_common.set_use_sed_dispersion(self,True)
#                 diff = self.domain.get_quantity('sed_diffusivity')
# #                 st.operators_added += 1
#                 
#                 self.neighbour_distance = st_util.get_distance_centroids(self)
# 
#             except KeyError:
#                 msg = 'Error! Quantity sed_diffusivity not yet created. '
#                 msg += 'Initialize when creating domain or turn '
#                 msg += 'use_sed_diffusivity off in sed_transport_config. '
#                 msg += 'Continuing without sediment diffusivity.'
#                 print msg
#                 self.use_sed_dispersion = False

        """
        Call quantities
        """
        self.q_conc=self.domain.quantities['concentration']
        self.conc_v=self.q_conc.vertex_values

        """
        Initiate sources and sinks of momentum
        """
        if self.use_momentum_sinks:
        
            if self.verbose:
                print 'Using momentum sources and sinks'
                
            self.domain.forcing_terms.append(manning_friction_implicit)
            st.momentum_sink_set.add(sv_common.mom_Concentration)
            st.momentum_sink_set.add(sv_common.mom_BedExchange)


    def init_veg(self):

        """
        Initiate vegetation quantities, get values from files
        """
        try:
            self.vegtype = \
                self.domain.quantities['vegetation'].vertex_values
        except KeyError:
            msg = 'Error! Quantity vegetation not yet created. '
            msg += 'Initialize when creating domain or use '
            msg += 'Sed_transport_operator only. Cannot continue!'
            sys.exit(msg)
        
        if self.vegtype.max() > 0:
        
            from anuga.file.csv_file import load_csv_as_array as load_array
            
            try:
                vegcodes = load_array(self.vegfile)
            except:
                msg = 'Error! Text file of stem spacing '
                msg += 'and diameters missing. Cannot continue!'
                sys.exit(msg)
        
            for t in range(len(vegcodes['vegcode'])):

                code=vegcodes['vegcode'][t]
                st_d=vegcodes['stem_diameter'][t]
                st_s=vegcodes['stem_spacing'][t]
            
                self.stem_diameter[self.vegtype==code] = st_d
                self.stem_spacing[self.vegtype==code] = st_s
                

        self.xarea = 0 * self.vegtype
        self.xarea[self.stem_spacing > 0] = \
            self.stem_diameter[self.stem_spacing > 0] / \
            self.stem_spacing[self.stem_spacing > 0]**2
        self.ad = self.xarea * self.stem_diameter
        
        self.xarea_c = np.mean(self.xarea, axis=1)

        """
        Initiate sources and sinks of momentum
        """
        if self.use_momentum_sinks:

            if self.verbose:
                print 'Using momentum sources and sinks'

            self.domain.forcing_terms.append(manning_friction_implicit)
        
        # Vegetation drag is always used as a momentum sink
        st.momentum_sink_set.add(sv_common.mom_VegDrag)
            
        """
        Other variables
        """
        self.dt = self.domain.flux_timestep



    def call_init(self):
        
        self.depth_e = self.q_stage.edge_values - self.q_elev.edge_values
        
        self.depth_c = self.stage_c - self.elev_c
        self.u_c = self.xmom_c / (self.depth_c + epsilon)
        self.v_c = self.ymom_c / (self.depth_c + epsilon)
        
        self.depth = self.depth_v[self.indices]
        
        self.u_v = self.xmom_v[self.indices] / (self.depth + epsilon)
        self.v_v = self.ymom_v[self.indices] / (self.depth + epsilon)

        self.Uref = np.sqrt(self.u_v**2 + self.v_v**2)
        


    def veg_drag(self):

        # drag around a cylinder
        self.drag_x = 0.5 * st.Cd * self.xarea[self.indices] * self.u_v**2
        self.drag_y = 0.5 * st.Cd * self.xarea[self.indices] * self.v_v**2
        wetveg = self.vegtype[self.indices]

        # reduce the flow velocity
        u_veg = np.sign(self.u_v) * \
                    np.maximum(abs(self.u_v) - self.drag_x * self.dt, 0.0)
        u_veg[wetveg==0] = self.u_v[wetveg==0]

        v_veg = np.sign(self.v_v) * \
                    np.maximum(abs(self.v_v) - self.drag_y * self.dt, 0.0)
        v_veg[wetveg==0] = self.v_v[wetveg==0]

        Uref_v = np.sqrt( u_veg**2 + v_veg**2 )
        Uref_v[wetveg==0] = self.Uref[wetveg==0]

        # replace the non-veg velocities with the veg-altered velocities
        self.u_v = u_veg
        self.v_v = v_veg
        self.Uref = Uref_v


    def sed_processes(self):

        #-----------------------------------------
        # Extra structures to support maintaining
        # continuity of elevation
        #----------------------------------------- 

        self.conc = self.conc_v[self.indices] / (self.depth + epsilon)
#         self.conc = np.clip(self.conc, 0, 1)
        
        """
        Calculate erosion + deposition
        """

        self.edot = self.erosion()
        self.ddot = self.deposition()
        
        """
        Changes to concentration
        """
        self.dChdt = (self.edot - self.ddot)
        
        st_mesh.compute_sed_flux(self)
        
        """
        Changes to the bed
        """
        self.dzdt = (self.ddot - self.edot) / (1 - st.porosity)
        

        st_mesh.update_elevation_stage(self, st.porosity)
        
        


    def erosion(self):

        self.ustar = np.sqrt(st.frict/8) * self.Uref
        
        shear_stress_star = st.ss_coeff * self.ustar**2
        edot = st.ero_coeff * (shear_stress_star - st.criticalshear_star)
        edot[edot<0.0] = 0.0
        
        maxe = st.max_rate_e * self.domain.timestep
        edot[edot>maxe] = maxe
        
        if not self.e:
            edot[:]=0.0
    
        return edot

    def deposition(self):

        ddot = st.rousecoeff * self.conc * st.settlingvelocity
        ddot[ddot<0.0] = 0.0
        
        maxd = st.max_rate_d * self.domain.timestep
        ddot[ddot>maxd] = maxd
        
        if not self.d:
            ddot[:]=0.0
    
        return ddot

#     def set_use_sed_dispersion(self, flag=False):
# 
#         from sed_dispersion_operator import Sed_dispersion_operator
# 
#         if flag :
#             # Create Operator if necessary
#             if self.sd_operator is None:
#                 self.sd_operator = Sed_dispersion_operator(self.domain,\
#                 'sed_diffusivity',use_triangle_areas=True)
#         else:
#             if self.sd_operator is None:
#                 return
#             else:
#                 # Remove operator from fractional_step_operators
#                 self.domain.fractional_step_operators.remove(self.sd_operator)
#                 self.sd_operator = None     



    def change_momentum(self):

        self.duhdt = 0
        self.dvhdt = 0
        
        for sink in st.momentum_sink_set:

            duhdt, dvhdt = sink(self)
            
            self.duhdt += duhdt
            self.dvhdt += dvhdt

        st_mesh.update_momentum(self)


    def mom_Concentration(self):

        c = self.q_conc.centroid_values / (self.depth_c + epsilon)

        mom_cg = st.cg_coeff * self.depth_c**2 / \
                (st.rho_w * (1-c) + st.rho_s * c)
                
        mom_cg = mom_cg.transpose()

        self.q_conc.compute_local_gradients()
        dCdx = np.array(self.q_conc.x_gradient).transpose()
        dCdy = np.array(self.q_conc.y_gradient).transpose()
        
        ax = dCdx * mom_cg
        ay = dCdy * mom_cg

        return np.absolute(ax).transpose(), np.absolute(ay).transpose()


    def mom_BedExchange(self):
    
        c = self.q_conc.centroid_values / (self.depth_c + epsilon)
        dz = np.zeros((len(self.stage_v),3))
        dz[self.indices] = self.dzdt
        dz_c = np.mean(dz, axis=1)

        mom_be = dz_c * ((st.rhoo / (st.rho_w * (1-c) + \
                st.rho_s * c)) - 1)
#         mom_be = - np.absolute(mom_be) # for opposite sign to velocity
    
        return np.absolute(mom_be * self.u_c), np.absolute(mom_be * self.v_c)


    def mom_VegDrag(self):
        
        vd_coeff = 0.5 * st.Cd * self.xarea_c * self.depth_c

        duhdx_c =  - np.sign(self.u_c) * vd_coeff * self.u_c**2
        dvhdy_c = - np.sign(self.v_c) * vd_coeff * self.v_c**2

        return duhdx_c, dvhdy_c




    def turbulence(self):
    
        if self.verbose:
            print 'Calculating turbulence'
            
        self.mix_length = self.depth - (self.depth - \
                  self.stem_diameter[self.indices]) * \
                 (self.ad[self.indices] - 0.005) / 0.005

        self.ke = st.Cb * self.Uref**2 + (st.Cd * \
                  self.ad[self.indices])**(0.66) * self.Uref**2

        self.diffusivity = np.sqrt( self.ke ) * self.mix_length + \
                           self.ad[self.indices] * \
                           abs(self.Uref) * self.stem_diameter[self.indices]

        st_mesh.update_quantity_nonconserved(self,\
                         'diffusivity',self.diffusivity)
 
 
    def verbose_output_sed(self):
    
        c = self.conc_v[self.indices]
            
        print 'edot:', self.edot.max()
        print 'ddot:', self.ddot.max()
        print 'dzdt:', self.dzdt.max(), ', ', self.dzdt.min()
        print 'dChdt:', self.dChdt.max()
        print 'conc:', c.max(), c.min()

        
    def verbose_output_common(self):
    
        d = self.depth_v[self.indices]
            
        print 'u_c:', np.absolute(self.u_v).max()
        print 'v_c:', np.absolute(self.v_v).max()
        print 'depth:', d.max(), d.min()
        print 60*'-'

        
    def get_wet_elements(self):
        
        self.depth_v = self.stage_v - self.elev_v

        indices = np.where(np.min(self.depth_v,axis=1) > st.min_depth)
        
        return indices[0]
   
   
    def clear_operators(self, domain, name):

        for i, o in enumerate(domain.fractional_step_operators):
            if repr(o) == name:
                del domain.fractional_step_operators[i]
                
    def set_arguments(self,
                      turbulence,
                      momentum_sinks,
                      verbose):
        
        # change the flags only if the new ones are True
        if turbulence:
            st.use_turbulence = True
        if momentum_sinks:
            st.use_momentum_sinks = True
        if verbose:
            st.verbose = True




######### Publicly callable operators #########
class Sed_transport_operator(sv_common):

    def __init__(self,
                 domain,
                 erosion = True,
                 deposition = True,
                 sed_dispersion = False,
                 turbulence = False,
                 momentum_sinks = False,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
        
        st.erosion = erosion
        st.deposition = deposition
        st.sed_dispersion = False
        
        instance_check = any(repr(x)=='Veg' for x in \
                         domain.fractional_step_operators)

        if not instance_check:
            # no operator has been called
            # flag and run init
            
            st.turbulence = turbulence
            st.momentum_sinks = momentum_sinks
            st.verbose = verbose
            
            Sed_transport_only(domain,
                               description = description,
                               label = label,
                               logging = logging)
            
        else:
            # if the other operator has been turned on
            # clean up and call Sed+Veg operator
            # want to clean up so Veg+Sed can also be called
            
            self.clean_init_veg(domain)
            
            sv_common.set_arguments(self,
                                    turbulence = turbulence,
                                    momentum_sinks = momentum_sinks,
                                    verbose = verbose)
            
            Vegetation_and_Sed(domain,
                               description = description,
                               label = label,
                               logging = logging)
        
        

    def clean_init_veg(self,domain):
    
        if st.verbose:
            print 'Cleaning up Vegetation Only init'
    
        try:
            st.momentum_sink_set.remove(sv_common.mom_VegDrag)
            self.domain.forcing_terms.remove(manning_friction_implicit)
        except:
            pass
            
        domain.set_use_kinematic_viscosity(False)
        sv_common.clear_operators(self, domain, 'Veg')



class Vegetation_operator(sv_common):
    """
    Operator for vegetation drag
    For now, it assumes that you also want sediment and did NOT call the sediment component too
    """
    
    def __init__(self,
                 domain,
                 vegfile = st.vegfile,
                 turbulence = st.turbulence,
                 momentum_sinks = st.momentum_sinks,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = st.verbose):
                 
        st.vegfile = vegfile

        instance_check = any(repr(x)=='Sed' for x in \
                         domain.fractional_step_operators)
        

        if not instance_check:
            # no operator has been called
            # flag and run init
            
            st.use_turbulence = turbulence
            st.use_momentum_sinks = momentum_sinks
            st.verbose = verbose
            
            Vegetation_only(domain,
                            description = description,
                            label = label,
                            logging = logging)
            
        else:
            # if the other operator has been turned on
            # clean up and call Sed+Veg operator
            # want to clean up so Veg+Sed can also be called
            self.clean_init_sed(domain)
            
            sv_common.set_arguments(self,
                                    turbulence = turbulence,
                                    momentum_sinks = momentum_sinks,
                                    verbose = verbose)
            
            Vegetation_and_Sed(domain,
                               description = description,
                               label = label,
                               logging = logging)      
          
    def clean_init_sed(self,domain):
    
        if st.verbose:
            print 'Cleaning up Sed Transport Only init'
    
        try:
            self.domain.forcing_terms.remove(manning_friction_implicit)
            st.momentum_sink_set.remove(sv_common.mom_Concentration)
            st.momentum_sink_set.remove(sv_common.mom_BedExchange)
        except:
            pass

        domain.set_use_kinematic_viscosity(False)
        sv_common.clear_operators(self,domain,'Sed')
        
        
############ Second-level operators #############

class Sed_transport_only(sv_common):
    """
    Operator for sediment entrainment, deposition, and transport
    """

    def __init__(self,
                 domain,
                 description,
                 label,
                 logging):
                 
                 
        self.e = st.erosion
        self.d = st.deposition
        self.use_sed_dispersion = st.sed_dispersion
        self.use_turbulence = st.turbulence
        self.use_momentum_sinks = st.momentum_sinks
        self.verbose = st.verbose

        sv_common.__init__(self, domain, \
                        description, label, logging, self.verbose)

        sv_common.init_sed(self)
        
        if self.verbose:
            print 'Using Sed Transport only'
            
    def __repr__(self):
        
        name = 'Sed'
        
        return name

    def __call__(self):
    
        
        self.indices = sv_common.get_wet_elements(self)
        
        if len(self.indices)>0:
        
            """ operate only on wet triangles """

            sv_common.call_init(self)

            sv_common.sed_processes(self)

            """
            Sources and sinks of Momentum
            """
            if self.use_momentum_sinks:
                sv_common.change_momentum(self)

            """
            Turbulence
            """
            if self.use_turbulence:
                sv_common.turbulence(self)
            
            if self.verbose:
                sv_common.verbose_output_sed(self)
                sv_common.verbose_output_common(self)
    
    
class Vegetation_and_Sed(sv_common):

    def __init__(self,
                 domain,
                 description,
                 label,
                 logging):
                 
        self.e = st.erosion
        self.d = st.deposition
        self.use_sed_dispersion = st.sed_dispersion
        self.use_turbulence = st.turbulence
        self.use_momentum_sinks = st.momentum_sinks
        self.verbose = st.verbose
        self.vegfile = st.vegfile

        sv_common.__init__(self, domain, \
                        description, label, logging, self.verbose)
        sv_common.init_sed(self)
        sv_common.init_veg(self)
        
        if self.verbose:
            print 'Using Sediment and Vegetation'
            
    def __repr__(self):
        
        name = 'Sed_Veg'
        
        return name

    def __call__(self):
    
        
        self.indices = sv_common.get_wet_elements(self)
        
        if len(self.indices)>0:
        
            """ operate only on wet triangles """

            sv_common.call_init(self)
            
            sv_common.veg_drag(self)

            sv_common.sed_processes(self)

            """
            Sources and sinks of Momentum
            Always compute sink to veg drag, so must always do this
            """
            sv_common.change_momentum(self)

            """
            Turbulence
            """
            if self.use_turbulence:
                sv_common.turbulence(self)

            
            if self.verbose:
                sv_common.verbose_output_sed(self)
                sv_common.verbose_output_common(self)
                
                
class Vegetation_only(sv_common):

    def __init__(self,
                 domain,
                 description = None,
                 label = None,
                 logging = False):
                                
                 
        self.use_turbulence = st.turbulence
        self.use_momentum_sinks = st.momentum_sinks
        self.vegfile = st.vegfile
        self.verbose = st.verbose

        sv_common.__init__(self, domain, \
                        description, label, logging, self.verbose)

        sv_common.init_veg(self)
        
        if self.verbose:
            print 'Using Vegetation only'
    
    
    def __repr__(self):
        
        name = 'Veg'
        
        return name
    

    def __call__(self):
        
        self.indices = sv_common.get_wet_elements(self)
        
        if len(self.indices)>0:
            """ operate only on wet triangles """

            sv_common.call_init(self)

            sv_common.veg_drag(self)
            

            """
            Sources and sinks of Momentum
            Always compute sink to veg drag, so must always do this
            """
            sv_common.change_momentum(self)

            """
            Turbulence
            """
            if self.use_turbulence:
                sv_common.turbulence(self)
                
            if self.verbose:
                sv_common.verbose_output_common(self)
           
    
################################################################################

"""
Common methods for using Operator
"""

def parallel_safe(self):
    """If Operator is applied independently on each cell and
    so is parallel safe.
    """
    
    return False



def statistics(self):

    message = self.label + ': Sed_transport_operator'
    message = message + ' on triangles '+ str(self.indices)

    return message


def timestepping_statistics(self):

    from anuga import indent

    message  = indent + self.label + ': Sed_transport_operator, time '
    message += str(self.get_time())+ ' max(Delta Elev) '+ str(self.max_change)
    return message
    
