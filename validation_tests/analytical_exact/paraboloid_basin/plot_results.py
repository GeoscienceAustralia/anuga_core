"""
    Quick plot of the dam break outputs
    Note: v = p2.y[?], where ? is the cell position of the origin.

"""
import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
#from math import sqrt, pi,sin,cos
from numpy import sqrt, sin, cos, pi
from numpy import ones
from analytical_paraboloid_basin import analytic_sol

# Get the sww file
p = util.get_output('paraboloid.sww')
p2=util.get_centroids(p)

v = p2.y[11248]#[4898]
v2=(p2.y==v)
time_level = 50
w2,u2,h2 = analytic_sol(p2.x[v2], v2, p2.time[time_level])


#Plot stages
pyplot.clf()
pyplot.plot(p2.x[v2], p2.stage[time_level,v2],'b.', label='numerical stage')
pyplot.plot(p2.x[v2], p2.elev[v2],'k-', label='bed elevation')
pyplot.plot(p2.x[v2], w2,'r-', label='analytical stage')
pyplot.title('Stage at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('cross_section_stage.png')


#Plot xmomentum
pyplot.clf()
pyplot.plot(p2.x[v2], p2.xmom[time_level,v2], 'b.', label='numerical')
pyplot.plot(p2.x[v2], u2*h2,'r-', label='analytical')
pyplot.title('Xmomentum at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('cross_section_xmom.png')


#Plot velocities
pyplot.clf()
pyplot.plot(p2.x[v2], p2.xvel[time_level,v2], 'b.', label='numerical')
pyplot.plot(p2.x[v2], u2,'r-', label='analytical')
pyplot.title('Xvelocity at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('cross_section_xvel.png')


# Parameters
D0 = 1000.
L = 2500.
R0 = 2000.
g = 9.81

A = (L**4 - R0**4)/(L**4 + R0**4)
omega = 2./L*sqrt(2.*g*D0)
T = pi/omega
index_origin = 11248#4898





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



