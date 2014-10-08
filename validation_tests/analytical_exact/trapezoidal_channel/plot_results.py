from anuga.utilities import plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
import numpy

from project import *

filename = 'channel_floodplain.sww'

# Time-index to plot outputs from
p2 = util.get_output(filename)
p=util.get_centroids(p2, velocity_extrapolation=True)
v = (p.x>6.0)*(p.x<8.0)

print numpy.any(v)
# Numerical results along a central channel 'slice'
index= -1
V1 = p.stage[index,v] - p.elev[v]
V2 = p.yvel[index,v]
V3 = p.xvel[index,v]

##########################################################################
# Analytical solution of steady uniform 2D flow in a trapezoidal channel.
##########################################################################

Qin=4.6932   # Inflow discharge
slp=1./300.  # Floodplain slope (= water slope for steady uniform flow)
man_n=0.03   # Manning's n
Bcentral=6.0 #Flat bed width of the trapezoidal channel
alpha=0.5    # Side slope of the trapezoidal banks

k = (slp*(1./man_n)**2)**0.5 # At any point, the analytical solution says U = k*d^(2/3)

# Function to calculate the discharge, given the channel centre depth dc, assuming
# steady uniform flow
def discharge_su(dc):
    if(alpha>0.):
        out = 2*k*( 3./(8.*alpha)*(dc)**(8./3.)) +Bcentral*k*(dc)**(5./3.)
    else:
        out = Bcentral*k*(dc)**(5./3.)    
    return out 

# Function that will be minimized to find the depth associated with discharge Qin
def minme(dc):
    q1 = discharge_su(dc)
    return (q1-Qin)**2.

# Minimise the function mimne, to find the centre depth.
import scipy.optimize
dc_analytical = scipy.optimize.fmin(minme, x0=1.0)[0]

print 'dc_analytic ',dc_analytical

##################################
# Plots
##################################

# Analytical solution has U*abs(U)*n^2 / D^(4./3.) = Sf = bed slope
# Hence, the following two variables should be equal -- I have checked that x velocities are fairly small
pyplot.clf()
pyplot.figure(figsize=(12.,8.))
pyplot.plot(p.y[v], V2*0+k*dc_analytical**(2./3.), 'o', label='analytical velocity')
pyplot.plot(p.y[v], (V2**2)**0.5,'o', label='numerical velocity')
pyplot.plot(p.y[v], V1, 'o',label='numerical depth')
pyplot.plot(p.y[v], V1*0. + dc_analytical, 'o', label='analytical depth')
pyplot.title('Mid channel numerical velocities and depths, vs analytical velocities and depths')
pyplot.legend(loc='best')
pyplot.xlabel('Down-channel distance (m)')
pyplot.ylabel('Generic scale (m or m/s)')
pyplot.savefig('fig1mid_channel.png')



# Plot velocity over the cross-section
pyplot.clf()
v1 = (p.y<105.0)&(p.y>95.0)


analytical_stage = min(p.elev[v1]) + dc_analytical
analytic_vel = ( (1./300.)*numpy.maximum(analytical_stage-p.elev[v1],0.0)**(4./3.)*(1./0.03)**2.)**0.5
analytic_vel = analytic_vel*(analytical_stage>p.elev[v1])

temp0 = p.stage[index,v1]*0. + analytical_stage
temp1 = (temp0) * (temp0 > p.elev[v1])
temp2 = (p.elev[v1]) * (temp0 < p.elev[v1])
Analytic_Stage = temp1 + temp2
pyplot.figure(figsize=(12.,8.))
pyplot.plot(p.x[v1], analytic_vel,'o', label='analytical velocity')
pyplot.plot(p.x[v1], p.yvel[index,v1],'o', label='numerical velocity')
#pyplot.plot(p.x[v1],p.stage[index,v1]-p.elev[v1],'ko', label='numerical height')
pyplot.plot(p.x[v1],p.stage[index,v1],'o', label='numerical stage')
pyplot.plot(p.x[v1],Analytic_Stage,'o', label='analytical stage')
pyplot.plot(p.x[v1],p.elev[v1],'o', label='bed elevation')
pyplot.ylim([-4,2])
pyplot.legend(loc=8)
pyplot.title('Velocity (analytical and numerical) and Stage:' + '\n' +'Downstream channel regions (95 to 105m)' +'\n')
pyplot.xlabel('Cross-channel distance (m)')
pyplot.ylabel('Generic scale (m or m/s)')
pyplot.savefig('fig2upstream_channel.png')



# Plot velocity over the cross-section
pyplot.clf()
v1 = (p.y<505.0)&(p.y>495.0)



analytical_stage = min(p.elev[v1]) + dc_analytical
analytic_vel = ( (1./300.)*numpy.maximum(analytical_stage-p.elev[v1],0.0)**(4./3.)*(1./0.03)**2.)**0.5
analytic_vel = analytic_vel*(analytical_stage>p.elev[v1])

temp0 = p.stage[index,v1]*0. + analytical_stage
temp1 = (temp0) * (temp0 > p.elev[v1])
temp2 = (p.elev[v1]) * (temp0 < p.elev[v1])
Analytic_Stage = temp1 + temp2
pyplot.figure(figsize=(12.,8.))
pyplot.plot(p.x[v1], analytic_vel,'o', label='analytical velocity')
pyplot.plot(p.x[v1], p.yvel[index,v1],'o', label='numerical velocity')
#pyplot.plot(p.x[v1],p.stage[index,v1]-p.elev[v1],'ko', label='numerical height')
pyplot.plot(p.x[v1],p.stage[index,v1],'o', label='numerical stage')
pyplot.plot(p.x[v1],Analytic_Stage,'o', label='analytical stage')
pyplot.plot(p.x[v1],p.elev[v1],'o', label='bed elevation')
pyplot.ylim([-4,2])
pyplot.legend(loc=10)
pyplot.title('Velocity (analytical and numerical) and Stage:' + '\n' +'Central channel regions (495 to 505m)' +'\n')
pyplot.xlabel('Cross-channel distance (m)')
pyplot.ylabel('Generic scale (m or m/s)')
pyplot.savefig('fig3central_channel.png') 



# Plot velocity over the cross-section
pyplot.clf()
v1 = (p.y<705.0)&(p.y>695.0)




analytical_stage = min(p.elev[v1]) + dc_analytical
analytic_vel = ( (1./300.)*numpy.maximum(analytical_stage-p.elev[v1],0.0)**(4./3.)*(1./0.03)**2.)**0.5
analytic_vel = analytic_vel*(analytical_stage>p.elev[v1])

temp0 = p.stage[index,v1]*0. + analytical_stage
temp1 = (temp0) * (temp0 > p.elev[v1])
temp2 = (p.elev[v1]) * (temp0 < p.elev[v1])
Analytic_Stage = temp1 + temp2
pyplot.figure(figsize=(12.,8.))
pyplot.plot(p.x[v1], analytic_vel,'o', label='analytical velocity')
pyplot.plot(p.x[v1], p.yvel[index,v1],'o', label='numerical velocity')
#pyplot.plot(p.x[v1],p.stage[index,v1]-p.elev[v1],'ko', label='numerical height')
pyplot.plot(p.x[v1],p.stage[index,v1],'o', label='numerical stage')
pyplot.plot(p.x[v1],Analytic_Stage,'o', label='analytical stage')
pyplot.plot(p.x[v1],p.elev[v1],'o', label='bed elevation')
pyplot.ylim([-4,2])
pyplot.legend(loc=10)
pyplot.title('Velocity (analytical and numerical) and Stage:' + '\n' +'Downstream channel regions (695 to 705m)' +'\n')
pyplot.xlabel('Cross-channel distance (m)')
pyplot.ylabel('Generic scale (m or m/s)')
pyplot.savefig('fig4downstream_channel.png')





print '#======================================================================'
print '# Extract some cross section info'
print '#======================================================================'

from anuga.shallow_water.sww_interrogate import get_flow_through_multiple_cross_sections

polyline0 = [ [floodplain_width, 10.0], [0., 10.0]]
polyline1 = [[floodplain_width, floodplain_length-300.0], [0., floodplain_length-300.0]]
polyline2 = [[floodplain_width, floodplain_length-1.0], [0., floodplain_length-1.0]]
        
polylines= [polyline0, polyline1, polyline2]        

time, [Q0,Q1,Q2]  = get_flow_through_multiple_cross_sections(filename, polylines, verbose=True)

pyplot.figure(figsize=(12.,8.))
pyplot.plot(time, Q0, label='10m')
pyplot.plot(time, Q1, label='500m')
pyplot.plot(time, Q2, label='799m')
pyplot.plot([0,time[-1]], [Qin,Qin], label='Input Q')
pyplot.ylim([0,7])
pyplot.legend(loc=10)
pyplot.title(' (Approximate) Cross sectional flow across transect at 10m, 500m and 799m')
pyplot.xlabel('Time (sec)')
pyplot.ylabel('Discharge (m^3/sec)')
pyplot.savefig('cross_section_10_500_790.png')




