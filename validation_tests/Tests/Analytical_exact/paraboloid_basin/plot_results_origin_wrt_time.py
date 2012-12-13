"""View results of runup.py
"""
#---------------
# Import Modules
#---------------
#import anuga
#import struct
import numpy
from matplotlib import pyplot as pyplot
from anuga.utilities import plot_utils as util
from numpy import sqrt, sin, cos, pi

# Parameters
D0 = 1000.
L = 2500.
R0 = 2000.
g = 9.81

A = (L**4 - R0**4)/(L**4 + R0**4)
omega = 2./L*sqrt(2.*g*D0)
T = pi/omega
index_origin = 11248#4898


# Get the sww file
p = util.get_output('paraboloid.sww')
p2=util.get_centroids(p)


# Free surface over time at the origin of spatial domain.
w = D0*((sqrt(1.-A*A))/(1.-A*cos(omega*p2.time)) -1.)

#Plot stages
pyplot.clf()
pyplot.plot(p2.time,w, color='blue', label='analytical')
pyplot.plot(p2.time,p2.stage[:,index_origin], color='green', label='numerical')
pyplot.legend(loc='best')
pyplot.title('Stage at the centre of the paraboloid basin over time')
pyplot.xlabel('Time')
pyplot.ylabel('Stage')
pyplot.savefig('Stage_origin.png')

#Plot xmomentums
pyplot.clf()
pyplot.plot(p2.time,0.0*p2.time, color='blue', label='analytical') # At the origin of space, w=h as always wet.
pyplot.plot(p2.time,p2.xmom[:,index_origin], color='green', label='numerical')
pyplot.legend(loc='best')
pyplot.title('Xmomentum at the centre of the paraboloid basin over time')
pyplot.xlabel('Time')
pyplot.ylabel('Xmomentum')
pyplot.savefig('Xmom_origin.png')

#Plot velocities
pyplot.clf()
pyplot.plot(p2.time,0.0*p2.time, color='blue', label='analytical')
pyplot.plot(p2.time,p2.xvel[:,index_origin], color='green', label='numerical')
pyplot.legend(loc='best')
pyplot.title('Xvelocity at the centre of the paraboloid basin over time')
pyplot.xlabel('Time')
pyplot.ylabel('Xvelocity')
pyplot.savefig('Xvel_origin.png')

