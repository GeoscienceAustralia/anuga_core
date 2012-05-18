"""View results of channel_SU.py
"""
#---------------
# Import Modules
#---------------
import time
import anuga
import numpy
import scipy
from matplotlib import pyplot as pyplot
from anuga.utilities import plot_utils as util
#--------------
# Get variables
#--------------
p=util.get_output('channel_SU_2_v2.sww', 0.001)
p2=util.get_centroids(p,velocity_extrapolation=True)

#------------------
# Select line
#------------------
#v=(y==2.5)
v=(p2.y==p2.y[3])

#-------------------------------
# Define variables of case study
#-------------------------------
mann=0.03 # Manning's coef
bedslope=-0.1
fluxin=20./100. #The momentum flux at the upstream boundary ( = discharge / width)

#---------------------------------------------
# Analytical solution for steady, uniform flow
#---------------------------------------------
# This comes from solving 
# 1) Friction_slope=bedslope , and
# 2) flux = depth*velocity = flux_in
uana= ( mann**(-2.)*abs(bedslope)*fluxin**(4./3.) )**(3./10.) # Velocity
dana= fluxin/uana # Depth

#--------------------
# Make plot animation
#--------------------
pyplot.close() #If the plot is open, there will be problems
pyplot.ion()

line, = pyplot.plot( (p2.x[v].min(),p2.x[v].max()) ,( (p2.stage[:,v]-p2.elev[:,v]).max(),(p2.stage[:,v]-p2.elev[v]).min() ) )
line.set_label('Numerical')
pyplot.plot( (0,100),(dana,dana), 'r',label='Analytical' )
pyplot.legend()
#pyplot.plot(x[v],elev[v],'green')
for i in range(p2.xmom.shape[0]):
    line.set_xdata(p2.x[v])
    line.set_ydata(p2.stage[i,v]-p2.elev[v])
    pyplot.draw()
    pyplot.title('Flow depth: should be converging to steady uniform flow ' )
    pyplot.xlabel('Distance down-slope (m)')
    pyplot.ylabel('depth (m)')

pyplot.title('Flow depth at 400s -- there should be a central region with steady uniform flow ' )
pyplot.savefig('final_depth_v2.png')

#--------------------------------------------
# Compare x-velocity with analytical solution
#--------------------------------------------
#pyplot.clf()
#pyplot.plot(x[v],xvel[200,v], label='Numerical')
#pyplot.plot((0,100),(uana,uana),'r',label='Analytical')
#pyplot.legend()
#pyplot.xlabel('x')
#pyplot.ylabel('Velocity m/s')
#plottitle = "Final velocity in the x-direction -- some wiggles. \n Can this be improved with\
# a different downstream boundary condition?" #\n Setting depth only at the downstream boundary does not help "
#pyplot.title(plottitle)
#pyplot.savefig('SU_x_velocity.png')


# Find an x value close to x==50
tmp=(abs(p2.x-50.)).argmin()
v=(p2.x==p2.x[tmp])

#--------------------------------------------------
# Check if y velocity is zero everywhere (it should be)
#--------------------------------------------------
pyplot.clf()
pyplot.plot(p2.y[v],p2.yvel[100,v],'o', label='Numerical')
pyplot.plot((0,100),(0.0, 0.0),label='Analytical')
pyplot.xlabel('Distance along the line x==50 (across the slope) (m)')
pyplot.ylabel('Velocity m/s')
pyplot.title('Final y-velocity near x=50.')
pyplot.legend()
pyplot.savefig('y_velocity_v2.png')


pyplot.clf()
pyplot.plot(p2.y[v],p2.xvel[100,v], 'o', label='Numerical')
pyplot.plot((0,100),(uana,uana),label='Analytical')
pyplot.xlabel('Distance along the line x==50 (across the slope) (m)')
pyplot.ylabel('Velocity m/s')
pyplot.title('Final velocity in the x-direction (downstream) \n at x=50. (similar for other x) ')
pyplot.legend()
pyplot.savefig('x_velocity_v2.png')
pyplot.clf()

