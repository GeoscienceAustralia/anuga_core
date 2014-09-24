"""
View results of run_channel.py

Simple water flow example using ANUGA: Water flowing down a channel.
It was called "steep_slope" in an old validation test.
"""
#---------------
# Import Modules
#---------------
import time
import anuga
import numpy
import scipy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from anuga.utilities import plot_utils as util
#--------------
# Get variables
#--------------
p=util.get_output('channel.sww', 0.001)
p2=util.get_centroids(p,velocity_extrapolation=True)

#------------------
# Select line
#------------------
#v=(y==2.5)
#v=(p2.y==p2.y[3])

#-------------------------------
# Define variables of case study
#-------------------------------
mann=0.03 # Manning's coef
bedslope=-0.1
fluxin=0.1/100. #The momentum flux at the upstream boundary ( = discharge / width)

#---------------------------------------------
# Analytical solution for steady, uniform flow
#---------------------------------------------
# This comes from solving 
# 1) Friction_slope=bedslope , and
# 2) flux = depth*velocity = flux_in
uana= ( mann**(-2.)*abs(bedslope)*fluxin**(4./3.) )**(3./10.) # Velocity
dana= fluxin/uana # Depth

#--------------------
# Make plots
#--------------------
# Find an y value close to y==50
#tmp=(abs(p2.y-50.)).argmin()
#vx=(p2.y==p2.y[tmp])
vx=(abs(p2.y - 50.0)<10.0)

pyplot.clf()
pyplot.plot(p2.x[vx],p2.stage[-1,vx]-p2.elev[vx], 'o', label='numerical')
pyplot.plot((0,400),(dana,dana),label='analytical')
pyplot.ylim([0.0,0.01])
pyplot.xlabel('Xposition m')
pyplot.ylabel('Depth m')
pyplot.title('Depth down the slope (along y=50.)')
pyplot.legend(loc='best')
pyplot.savefig('depth_x.png')



pyplot.clf()
pyplot.plot(p2.x[vx],p2.xvel[-1,vx], 'o', label='numerical')
pyplot.plot((0,400),(uana,uana),label='analytical')
#pyplot.ylim([0.0,0.05])
pyplot.xlabel('Xposition m')
pyplot.ylabel('Velocity m/s')
pyplot.title('X Velocity down the slope (along y=50.)')
pyplot.legend(loc='best')
pyplot.savefig('xvelocity_x.png')




#--------------------------------------------
# Compare velocity with analytical solution
#--------------------------------------------
# Find an x value close to x==50
#tmp=(abs(p2.x-50.)).argmin()
v=(abs(p2.x - 50.0)<10.0)

pyplot.clf()
pyplot.plot(p2.y[v],p2.stage[-1,v]-p2.elev[v], 'o', label='numerical')
pyplot.plot((0,100),(dana,dana),label='analytical')
pyplot.ylim([0.0,0.005])
pyplot.xlabel('Yposition m')
pyplot.ylabel('Depth m')
pyplot.title('Depth across the slope (x=50.)')
pyplot.legend(loc='best')
pyplot.savefig('depth_y.png')


pyplot.clf()
pyplot.plot(p2.y[v],p2.xvel[-1,v], 'o', label='numerical')
pyplot.plot((0,100),(uana,uana),label='analytical')
pyplot.ylim([0.0,0.35])
pyplot.xlabel('Yposition along the line x=50')
pyplot.ylabel('Velocity m/s')
pyplot.title('X Velocity around the line x=50.')
pyplot.legend(loc='best')
pyplot.savefig('x_velocity.png')

pyplot.clf()
pyplot.plot(p2.y[v],p2.yvel[-1,v],'o', label='numerical')
pyplot.plot((0,100),(0.0, 0.0),label='analytical')
#pyplot.ylim([0.0,0.3])
pyplot.xlabel('Yposition along the line x=50')
pyplot.ylabel('Velocity m/s')
pyplot.title('Y Velocity around the line x=50.')
pyplot.legend(loc='best')
pyplot.savefig('y_velocity.png')




