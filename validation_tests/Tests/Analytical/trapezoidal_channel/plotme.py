from anuga.utilities import plot_utils as util
from matplotlib import pyplot as pyplot
#import pylab

# Time-index to plot outputs from
index=600

#p2 = util.get_output('channel_floodplain1_bal_dev.sww', minimum_allowed_height=0.01)
#p=p2
#p=util.get_centroids(p2, velocity_extrapolation=True)

p2 = util.get_output('channel_floodplain1.sww')
p=util.get_centroids(p2, velocity_extrapolation=True)


#p2 = util.get_output('channel_floodplain1_test_edge_lim.sww', minimum_allowed_height=0.01)
#p2 = util.get_output('channel_floodplain1_test.sww', minimum_allowed_height=0.01)
#p=util.get_centroids(p2, velocity_extrapolation=False)

#p2 = util.get_output('channel_floodplain1_standard.sww', minimum_allowed_height=0.01)
#p2 = util.get_output('channel_floodplain1_balanced_basic.sww', minimum_allowed_height=0.01)
#p2 = util.get_output('channel_floodplain1_dev.sww')
#p=util.get_centroids(p2, velocity_extrapolation=True)
#p=p2
v = (p.x>6.0)*(p.x<8.0)
#util.animate_1D(p.time, p.stage[:, v], p.y[v])

# Numerical results along a central channel 'slice'
V1 = p.stage[index,v] - p.elev[v]
V2 = p.yvel[index,v]
V3 = p.xvel[index,v]

##########################################################################
# Analytical solution of steady uniform 2D flow in a trapezoidal channel.
##########################################################################

#Qin=6.6339  # Inflow discharge
#Qin=5.2
Qin=4.6932
#Qin=4.693286 # Inflow discharge
slp=1./300. # Floodplain slope (= water slope for steady uniform flow)
man_n=0.03  # Manning's n
Bcentral=6.0 #Flat bed width of the trapezoidal channel
alpha=0.5  # Side slope of the trapezoidal banks

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



##################################
# Plots
##################################

# Analytical solution has U*abs(U)*n^2 / D^(4./3.) = Sf = bed slope
# Hence, the following two variables should be equal -- I have checked that x velocities are fairly small
pyplot.clf()
pyplot.figure(figsize=(12.,8.))
pyplot.plot(p.y[v], (V2**2)**0.5,'o', label='computed velocity')
pyplot.plot(p.y[v], V2*0+k*dc_analytical**(2./3.), 'o', label='Analytical velocity')
pyplot.plot(p.y[v], V1, 'o',label='computed depth')
pyplot.plot(p.y[v], V1*0. + dc_analytical, 'o', label='analytical depth')
#pyplot.plot( ( (1./300.)*V1**(4./3.)*(1./0.03)**2.)**0.5,'o', label='Analytical velocity based on computed depth')
pyplot.title('Mid channel velocities and depths, vs analytical velocities and depths')
# But in my tests, they are not equal
pyplot.legend( ('computed velocity', 'Analytical velocity', 'computed depth', 'analytical depth'), loc=4)
pyplot.savefig('trapz_velocity_downstream_l0_eq_1_EL.png')

# Plot velocity over the cross-section
v1 = (p.y<500.0)&(p.y>490.0)

pyplot.clf()
analytical_stage = min(p.elev[v1]) + dc_analytical
analytic_vel = ( (1./300.)*(analytical_stage-p.elev[v1])**(4./3.)*(1./0.03)**2.)**0.5
analytic_vel = analytic_vel*(analytical_stage>p.elev[v1])
pyplot.figure(figsize=(12.,8.))
pyplot.plot(p.x[v1], p.yvel[index,v1],'o', label='computed velocity (m/s)')
pyplot.plot(p.x[v1], analytic_vel,'o', label='analytical velocity (m/s)')
pyplot.plot(p.x[v1],p.elev[v1],'o', label='bed elevation (m)')
pyplot.plot(p.x[v1],p.stage[index,v1],'o', label='computed stage (m)')
pyplot.plot(p.x[v1],p.stage[index,v1]*0. + analytical_stage,'o', label='analytical stage (m)')

pyplot.legend( ('computed velocity (m/s)', 'analytical velocity (m/s)', 'bed elevation (m)', 'computed stage (m)', 'analytical_stage (m)') ,loc=10)
pyplot.title('Velocity (analytical and numerical) and Stage:' + '\n' +'Central channel regions (470 to 500m)' +'\n')
pyplot.savefig('trapz_velocity_cross_channel_l0_eq_1_EL.png') 


# Plot velocity over the cross-section
v1 = (p.y<800.0)&(p.y>790.0)

pyplot.clf()
analytical_stage = min(p.elev[v1]) + dc_analytical
analytic_vel = ( (1./300.)*(analytical_stage-p.elev[v1])**(4./3.)*(1./0.03)**2.)**0.5
analytic_vel = analytic_vel*(analytical_stage>p.elev[v1])
pyplot.figure(figsize=(12.,8.))
pyplot.plot(p.x[v1], p.yvel[index,v1],'o', label='computed velocity (m/s)')
pyplot.plot(p.x[v1], analytic_vel,'o', label='analytical velocity (m/s)')
pyplot.plot(p.x[v1],p.elev[v1],'o', label='bed elevation (m)')
pyplot.plot(p.x[v1],p.stage[index,v1],'o', label='computed stage (m)')
pyplot.plot(p.x[v1],p.stage[index,v1]*0. + analytical_stage,'o', label='analytical stage (m)')

pyplot.legend( ('computed velocity (m/s)', 'analytical velocity (m/s)', 'bed elevation (m)', 'computed stage (m)', 'analytical_stage (m)') ,loc=10)
pyplot.title('Velocity (analytical and numerical) and Stage:' + '\n' +'Central channel regions (470 to 500m)' +'\n')
pyplot.savefig('trapz_velocity_cross_channel_l0_eq_1b_EL.png') 

# Plot velocity over the cross-section
v1 = (p.y<900.0)&(p.y>890.0)

pyplot.clf()
analytical_stage = min(p.elev[v1]) + dc_analytical
analytic_vel = ( (1./300.)*(analytical_stage-p.elev[v1])**(4./3.)*(1./0.03)**2.)**0.5
analytic_vel = analytic_vel*(analytical_stage>p.elev[v1])
pyplot.figure(figsize=(12.,8.))
pyplot.plot(p.x[v1], p.yvel[index,v1],'o', label='computed velocity (m/s)')
pyplot.plot(p.x[v1], analytic_vel,'o', label='analytical velocity (m/s)')
pyplot.plot(p.x[v1],p.elev[v1],'o', label='bed elevation (m)')
pyplot.plot(p.x[v1],p.stage[index,v1],'o', label='computed stage (m)')
pyplot.plot(p.x[v1],p.stage[index,v1]*0. + analytical_stage,'o', label='analytical stage (m)')

pyplot.legend( ('computed velocity (m/s)', 'analytical velocity (m/s)', 'bed elevation (m)', 'computed stage (m)', 'analytical_stage (m)') , loc=10)
pyplot.title('Velocity (analytical and numerical) and Stage:' + '\n' +'Central channel regions (870 to 900m)' +'\n')
pyplot.savefig('trapz_velocity_cross_channel_l0_eq_1c_EL.png') 
