#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$23/08/2010 5:18:51 PM$"


import culvert_routine
from anuga.config import velocity_protection
from anuga.utilities.numerical_tools import safe_acos as acos

from math import pi, sqrt, sin, cos
from anuga.config import g


class Boyd_box_routine(culvert_routine.Culvert_routine):
    """Boyd's generalisation of the US department of transportation culvert methods

        WARNING THIS IS A SIMPLISTIC APPROACH and OUTLET VELOCITIES ARE LIMITED TO EITHER
        FULL PIPE OR CRITICAL DEPTH ONLY
        For Supercritical flow this is UNDERESTIMATING the Outlet Velocity
        The obtain the CORRECT velocity requires an iteration of Depth to Establish
        the Normal Depth of flow in the pipe.

        It is proposed to provide this in a seperate routine called
        boyd_generalised_culvert_model_complex

        The Boyd Method is based on methods described by the following:
        1.
        US Dept. Transportation Federal Highway Administration (1965)
        Hydraulic Chart for Selection of Highway Culverts.
        Hydraulic Engineering Circular No. 5 US Government Printing
        2.
        US Dept. Transportation Federal Highway Administration (1972)
        Capacity charts for the Hydraulic design of highway culverts.
        Hydraulic Engineering Circular No. 10 US Government Printing
        These documents provide around 60 charts for various configurations of culverts and inlets.

        Note these documents have been superceded by:
        2005 Hydraulic Design of Highway Culverts, Hydraulic Design Series No. 5 (HDS-5),
        Which combines culvert design information previously contained in Hydraulic Engineering Circulars
        (HEC) No. 5, No. 10, and No. 13 with hydrologic, storage routing, and special culvert design information.
        HEC-5 provides 20 Charts
        HEC-10 Provides an additional 36 Charts
        HEC-13 Discusses the Design of improved more efficient inlets
        HDS-5 Provides 60 sets of Charts

        In 1985 Professor Michael Boyd Published "Head-Discharge Relations for Culverts", and in
        1987 published "Generalised Head Discharge Equations for Culverts".
        These papers reviewed the previous work by the US DOT and provided a simplistic approach for 3 configurations.

        It may be possible to extend the same approach for additional charts in the original work, but to date this has not been done.
        The additional charts cover a range of culvert shapes and inlet configurations


        """

    def __init__(self, culvert):

        culvert_routine.Culvert_routine.__init__(self, culvert)
        
        self.manning = culvert.manning



    def __call__(self):

        self.determine_inflow()

        local_debug ='false'
        
        if self.inflow.get_enquiry_height() > 0.01: #this value was 0.01:
            if local_debug =='true':
                log.critical('Specific E & Deltat Tot E = %s, %s'
                             % (str(self.inflow.get_enquiry_specific_energy()),
                                str(self.delta_total_energy)))
                log.critical('culvert type = %s' % str(culvert_type))
            # Water has risen above inlet

            if self.log_filename is not None:
                s = 'Specific energy  = %f m' % self.inflow.get_enquiry_specific_energy()
                log_to_file(self.log_filename, s)

            msg = 'Specific energy at inlet is negative'
            assert self.inflow.get_enquiry_specific_energy() >= 0.0, msg

            if self.use_velocity_head :
                driving_energy = self.inflow.get_enquiry_specific_energy()
            else:
                driving_energy = self.get_enquiry_height

            height = self.culvert_height
            width = self.culvert_width
            flow_width = self.culvert_width

            Q_inlet_unsubmerged = 0.540*g**0.5*width*driving_energy**1.50 # Flow based on Inlet Ctrl Inlet Unsubmerged
            Q_inlet_submerged = 0.702*g**0.5*width*height**0.89*driving_energy**0.61  # Flow based on Inlet Ctrl Inlet Submerged

            # FIXME(Ole): Are these functions really for inlet control?
            if Q_inlet_unsubmerged < Q_inlet_submerged:
                Q = Q_inlet_unsubmerged
                dcrit = (Q**2/g/width**2)**0.333333
                if dcrit > height:
                    dcrit = height
                flow_area = width*dcrit
                outlet_culvert_depth = dcrit
                case = 'Inlet unsubmerged Box Acts as Weir'
            else:
                Q = Q_inlet_submerged
                flow_area = width*height
                outlet_culvert_depth = height
                case = 'Inlet submerged Box Acts as Orifice'

            dcrit = (Q**2/g/width**2)**0.333333

            outlet_culvert_depth = dcrit
            if outlet_culvert_depth > height:
                outlet_culvert_depth = height  # Once again the pipe is flowing full not partfull
                flow_area = width*height  # Cross sectional area of flow in the culvert
                perimeter = 2*(width+height)
                case = 'Inlet CTRL Outlet unsubmerged PIPE PART FULL'
            else:
                flow_area = width * outlet_culvert_depth
                perimeter = width+2*outlet_culvert_depth
                case = 'INLET CTRL Culvert is open channel flow we will for now assume critical depth'

            if self.delta_total_energy < driving_energy:
                # Calculate flows for outlet control

                # Determine the depth at the outlet relative to the depth of flow in the Culvert
                if self.outflow.get_enquiry_height() > height:        # The Outlet is Submerged
                    outlet_culvert_depth=height
                    flow_area=width*height       # Cross sectional area of flow in the culvert
                    perimeter=2.0*(width+height)
                    case = 'Outlet submerged'
                else:   # Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
                    dcrit = (Q**2/g/width**2)**0.333333
                    outlet_culvert_depth=dcrit   # For purpose of calculation assume the outlet depth = Critical Depth
                    if outlet_culvert_depth > height:
                        outlet_culvert_depth=height
                        flow_area=width*height
                        perimeter=2.0*(width+height)
                        case = 'Outlet is Flowing Full'
                    else:
                        flow_area=width*outlet_culvert_depth
                        perimeter=(width+2.0*outlet_culvert_depth)
                        case = 'Outlet is open channel flow'

                hyd_rad = flow_area/perimeter

                if self.log_filename is not None:
                    s = 'hydraulic radius at outlet = %f' % hyd_rad
                    log_to_file(self.log_filename, s)

                # Outlet control velocity using tail water
                culvert_velocity = sqrt(self.delta_total_energy/((self.sum_loss/2/g)+(self.manning**2*self.culvert_length)/hyd_rad**1.33333))
                Q_outlet_tailwater = flow_area * culvert_velocity

                if self.log_filename is not None:
                    s = 'Q_outlet_tailwater = %.6f' % Q_outlet_tailwater
                    log_to_file(self.log_filename, s)
                Q = min(Q, Q_outlet_tailwater)
            else:
                pass
                #FIXME(Ole): What about inlet control?

            culv_froude=sqrt(Q**2*flow_width/(g*flow_area**3))
            if local_debug =='true':
                log.critical('FLOW AREA = %s' % str(flow_area))
                log.critical('PERIMETER = %s' % str(perimeter))
                log.critical('Q final = %s' % str(Q))
                log.critical('FROUDE = %s' % str(culv_froude))

            # Determine momentum at the outlet
            barrel_velocity = Q/(flow_area + velocity_protection/flow_area)

        # END CODE BLOCK for DEPTH  > Required depth for CULVERT Flow

        else: # self.inflow.get_enquiry_height() < 0.01:
            Q = barrel_velocity = outlet_culvert_depth = 0.0

        # Temporary flow limit
        if barrel_velocity > self.max_velocity:
            barrel_velocity = self.max_velocity
            Q = flow_area * barrel_velocity



        

        return Q, barrel_velocity, outlet_culvert_depth



