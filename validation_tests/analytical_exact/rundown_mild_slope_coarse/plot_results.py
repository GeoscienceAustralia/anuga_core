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

tmp = abs(p2.y-50.)
n = tmp.argmin()
v=(p2.y==p2.y[n])

ll=p2.vel.shape[0]-1

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
# Make plot animation
#--------------------
pyplot.clf()


print (p2.x[v].min(),p2.x[v].max())
print (p2.stage[:,v]-p2.elev[v]).max()
print (p2.stage[:,v]-p2.elev[v]).min()

line, = pyplot.plot( (p2.x[v].min(),p2.x[v].max()) ,( (p2.stage[:,v]-p2.elev[v]).max(),(p2.stage[:,v]-p2.elev[v]).min() ) )
line.set_label('numerical')
pyplot.plot( (0,140),(dana,dana), 'r',label='analytical' )
pyplot.yscale('log')
pyplot.legend()
#pyplot.plot(x[v],elev[v],'green')
for i in range(p2.xmom.shape[0]):
    line.set_xdata(p2.x[v])
    line.set_ydata(p2.stage[i,v]-p2.elev[v])
    pyplot.draw()
    pyplot.title('Flow depth: should be converging to steady uniform flow ' )
    pyplot.xlabel('Xposition')
    pyplot.ylabel('Depth')
#pyplot.ylim([0.085,0.095])
pyplot.legend(loc='best')
pyplot.title('Flow depth at 800s')# -- there should be a central region with steady uniform flow ' )
pyplot.savefig('final_depth.png')

#--------------------------------------------
# Compare velocity with analytical solution
#--------------------------------------------
# Find an x value close to x==50
tmp=(abs(p2.x-50.)).argmin()
v=(p2.x==p2.x[tmp])

pyplot.clf()
pyplot.plot(p2.y[v],p2.xvel[ll,v], 'o', label='numerical')
pyplot.plot((0,100),(uana,uana),label='analytical')
#pyplot.ylim([1.0,3.0])
pyplot.xlabel('Yposition along the line x=50')
pyplot.ylabel('Xvelocity m/s')
pyplot.title('Final Xvelocity around the line x=50.')
pyplot.legend(loc='best')
pyplot.ylim((0,1.5*uana))
pyplot.savefig('x_velocity.png')

pyplot.clf()
pyplot.plot(p2.y[v],p2.yvel[ll,v],'o', label='numerical')
pyplot.plot((0,100),(0.0, 0.0),label='analytical')
pyplot.xlabel('Yposition along the line x=50')
pyplot.ylabel('Yvelocity')
pyplot.title('Final Yvelocity around the line x=50.')
pyplot.legend(loc='best')
pyplot.savefig('y_velocity.png')




