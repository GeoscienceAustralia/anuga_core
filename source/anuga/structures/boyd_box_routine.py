#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$23/08/2010 5:18:51 PM$"



def boyd_box(height, width, flow_width, inflow_specific_energy):
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

    # Calculate flows for inflow control

    Q_inflow_unsubmerged = 0.540*g**0.5*width*inflow_specific_energy**1.50 # Flow based on inflow Ctrl inflow Unsubmerged
    Q_inflow_submerged = 0.702*g**0.5*width*height**0.89*inflow_specific_energy**0.61  # Flow based on inflow Ctrl inflow Submerged

    if log_filename is not None:
        s = 'Q_inflow_unsubmerged = %.6f, Q_inflow_submerged = %.6f' %(Q_inflow_unsubmerged, Q_inflow_submerged)
        log_to_file(log_filename, s)

    # FIXME(Ole): Are these functions really for inflow control?
    if Q_inflow_unsubmerged < Q_inflow_submerged:
        Q = Q_inflow_unsubmerged
        dcrit = (Q**2/g/width**2)**0.333333
        if dcrit > height:
            dcrit = height
        flow_area = width*dcrit
        outflow_culvert_depth = dcrit
        case = 'inflow unsubmerged Box Acts as Weir'
    else:
        Q = Q_inflow_submerged
        flow_area = width*height
        outflow_culvert_depth = height
        case = 'inflow submerged Box Acts as Orifice'

    dcrit = (Q**2/g/width**2)**0.333333

    outflow_culvert_depth = dcrit
    if outflow_culvert_depth > height:
        outflow_culvert_depth = height  # Once again the pipe is flowing full not partfull
        flow_area = width*height  # Cross sectional area of flow in the culvert
        perimeter = 2*(width+height)
        case = 'inflow CTRL outflow unsubmerged PIPE PART FULL'
    else:
        flow_area = width * outflow_culvert_depth
        perimeter = width+2*outflow_culvert_depth
        case = 'inflow CTRL Culvert is open channel flow we will for now assume critical depth'

    if delta_total_energy < inflow_specific_energy:
        # Calculate flows for outflow control

        # Determine the depth at the outflow relative to the depth of flow in the Culvert
        if outflow_depth > height:        # The outflow is Submerged
            outflow_culvert_depth=height
            flow_area=width*height       # Cross sectional area of flow in the culvert
            perimeter=2.0*(width+height)
            case = 'outflow submerged'
        else:   # Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
            dcrit = (Q**2/g/width**2)**0.333333
            outflow_culvert_depth=dcrit   # For purpose of calculation assume the outflow depth = Critical Depth
            if outflow_culvert_depth > height:
                outflow_culvert_depth=height
                flow_area=width*height
                perimeter=2.0*(width+height)
                case = 'outflow is Flowing Full'
            else:
                flow_area=width*outflow_culvert_depth
                perimeter=(width+2.0*outflow_culvert_depth)
                case = 'outflow is open channel flow'

        hyd_rad = flow_area/perimeter

        if log_filename is not None:
            s = 'hydraulic radius at outflow = %f' % hyd_rad
            log_to_file(log_filename, s)

        # outflow control velocity using tail water
        culvert_velocity = sqrt(delta_total_energy/((sum_loss/2/g)+(manning**2*culvert_length)/hyd_rad**1.33333))
        Q_outflow_tailwater = flow_area * culvert_velocity

        if log_filename is not None:
            s = 'Q_outflow_tailwater = %.6f' % Q_outflow_tailwater
            log_to_file(log_filename, s)
        Q = min(Q, Q_outflow_tailwater)

    return Q


if __name__ == "__main__":


    g=9.81
    culvert_slope=0.1  # Downward

    inlet_depth=2.0
    outlet_depth=0.0

    inlet_velocity=0.0,
    outlet_velocity=0.0,

    culvert_length=4.0
    culvert_width=1.2
    culvert_height=0.75

    culvert_type='box'
    manning=0.013
    sum_loss=0.0

    inlet_specific_energy=inlet_depth #+0.5*v**2/g
    z_in = 0.0
    z_out = -culvert_length*culvert_slope/100
    E_in = z_in+inlet_depth # +
    E_out = z_out+outlet_depth # +
    delta_total_energy = E_in-E_out

    Q = boyd_box(culvert_height, culvert_width, culvert_width, inlet_specific_energy)

    print 'Q ',Q
