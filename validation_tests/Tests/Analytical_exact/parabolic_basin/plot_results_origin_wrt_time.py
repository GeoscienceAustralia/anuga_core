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
#from analytical_parabolic_basin import analytic_cannal

p=util.get_output('parabola.sww', 0.01)
p2=util.get_centroids(p, velocity_extrapolation=True)

### Define some 'lines' along which to plot
##v=(p2.y==p2.y[0])

# Compute the analytical solution
D0=4.0
L=10.0
A=2.0
g=9.81

omega=numpy.sqrt(2*D0*g)/L
T= 2*numpy.pi/omega
index_origin = 3978

# Free surface over time at the origin of spatial domain.
w = D0 + 2*A*D0/(L**2)*numpy.cos(omega*p2.time)*( p2.x[index_origin] -A/2.*numpy.cos(omega*p2.time))
#w = w*(w>p2.elev[index_origin])+p2.elev[index_origin]*(w<=p2.elev[index_origin]) #not needed

# Velocity at the origin of the spatial domain.
u = -A*omega*numpy.sin(omega*p2.time)

#Plot stages
pyplot.clf()
pyplot.plot(p2.time,w, color='blue', label='analytical')
pyplot.plot(p2.time,p2.stage[:,index_origin], color='green', label='numerical')
pyplot.legend(loc='best')
pyplot.title('Stage at the centre of the parabolic basin over time')
pyplot.xlabel('Time')
pyplot.ylabel('Stage')
pyplot.savefig('Stage_centre.png')

#Plot xmomentums
pyplot.clf()
pyplot.plot(p2.time,u*w, color='blue', label='analytical') # At the origin of space, w=h as always wet.
pyplot.plot(p2.time,p2.xmom[:,index_origin], color='green', label='numerical')
pyplot.legend(loc='best')
pyplot.title('Xmomentum at the centre of the parabolic basin over time')
pyplot.xlabel('Time')
pyplot.ylabel('Xmomentum')
pyplot.savefig('Xmom_centre.png')

#Plot velocities
pyplot.clf()
pyplot.plot(p2.time,u, color='blue', label='analytical')
pyplot.plot(p2.time,p2.xvel[:,index_origin], color='green', label='numerical')
pyplot.legend(loc='best')
pyplot.title('Xvelocity at the centre of the parabolic basin over time')
pyplot.xlabel('Time')
pyplot.ylabel('Xvelocity')
pyplot.savefig('Xvel_centre.png')

