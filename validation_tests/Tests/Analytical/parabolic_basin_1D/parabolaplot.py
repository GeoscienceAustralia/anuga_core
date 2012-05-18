"""View results of runup.py
"""
#---------------
# Import Modules
#---------------
#import anuga
#import struct
import numpy
#import scipy
from matplotlib import pyplot as pyplot
from anuga.utilities import plot_utils as util

p=util.get_output('parabola_v2.sww', 0.01)
p2=util.get_centroids(p, velocity_extrapolation=True)

# Define some 'lines' along which to plot
v=(p2.y==p2.y[0])

#--------------------
# Make plot animation
#--------------------
pyplot.close() #If the plot is open, there will be problems
pyplot.ion()

if False:
    line, = pyplot.plot( (p2.x[v].min(),p2.x[v].max()) ,(p2.stage[:,v].min(),p2.stage[:,v].max() ) )
    for i in range(700,p2.xmom.shape[0]):
        line.set_xdata(p2.x[v])
        line.set_ydata(p2.stage[i,v])
        pyplot.draw()
        #pyplot.plot( (0,1),(0,0), 'r' )
        pyplot.plot(p2.x[v],p2.elev[v]) 
        pyplot.title(str(i)+'/200') # : velocity does not converge to zero' )
        pyplot.xlabel('x')
        pyplot.ylabel('stage (m/s)')

if False:
    pyplot.clf()
 
    line, = pyplot.plot( (p2.x[v].min(),p2.x[v].max()) ,(p2.xvel[:,v].min(),p2.xvel[:,v].max() ), 'ro' )
    for i in range(700,p2.xmom.shape[0]):
        line.set_xdata(p2.x[v])
        line.set_ydata(p2.xvel[i,v])
        pyplot.draw()
        pyplot.plot( (p2.x[v].min(),p2.x[v].max()),(0,0), 'ro' )
        #pyplot.plot(x[v],elev[v]) 
        pyplot.title(str(i)+'/200') # : velocity does not converge to zero' )
        pyplot.xlabel('x')
        pyplot.ylabel('xvel (m/s)')


# Compute the analytical solution
D0=4.0
L=10.0
A=2.0
g=9.8
#t=time[30]

omega=numpy.sqrt(2*D0*g)/L
T= 2*numpy.pi/omega
ppp= (abs(p2.x-4.*L/2.)).argmin()
ppp2= (abs(p2.x-4.*L/4.)).argmin()

# Free surface over time in the channel centre
w = D0 + 2*A*D0/(L**2)*numpy.cos(omega*p2.time)*( (p2.x[ppp]-4.*L/2.) -A/2.*numpy.cos(omega*p2.time))
# Free surface over time at a wet/dry point on the parabola.
w2 = D0 + 2*A*D0/(L**2)*numpy.cos(omega*p2.time)*( (p2.x[ppp2]-4.*L/2.) -A/2.*numpy.cos(omega*p2.time))
w2 = w2*(w2>p2.elev[ppp2])+p2.elev[ppp2]*(w2<=p2.elev[ppp2])


pyplot.clf()
pyplot.plot(p2.time,w, color='blue', label='analytical')
pyplot.plot(p2.time,p2.stage[:,ppp], color='green', label='numerical')
pyplot.legend()
pyplot.xlabel('time (s)')
pyplot.ylabel('Stage (m)')
pyplot.savefig('Stage_centre_v2.png')
#       pyplot.savefig('runup_x_velocities.png')
pyplot.clf()
pyplot.plot(p2.time,w2, color='blue', label='analytical')
pyplot.plot(p2.time,p2.stage[:,ppp2], color='green', label='numerical')
pyplot.legend(loc=10)
pyplot.xlabel('time (s)')
pyplot.ylabel('Stage (m)')
pyplot.savefig('Stage_centre_v3.png')


pltind=numpy.argmin(abs(p2.time-T*3))

# Free surface at time p2.time[pltind]
w = D0 + 2*A*D0/(L**2)*numpy.cos(omega*p2.time[pltind])*( (p2.x[v]-4.*L/2.) -A/2.*numpy.cos(omega*p2.time[pltind]))
pyplot.clf()
pyplot.plot(p2.x[v],p2.xvel[pltind,v], label='Numerical')
pyplot.plot(p2.x[v], -omega*A*numpy.sin(omega*p2.time[pltind])*(w>p2.elev[v]), label='Analytical')
pyplot.legend()
pyplot.xlabel('x (m)')
pyplot.ylabel('x velocity (m/s)')
pyplot.title('Velocity near to time=3T')
pyplot.savefig('Vel_3T_v2.png')

pltind=numpy.argmin(abs(p2.time-T*3.25))
# Free surface at time p2.time[pltind]
w = D0 + 2*A*D0/(L**2)*numpy.cos(omega*p2.time[pltind])*( (p2.x[v]-4.*L/2.) -A/2.*numpy.cos(omega*p2.time[pltind]))

pyplot.clf()
pyplot.plot(p2.x[v],p2.xvel[pltind,v], label='Numerical')
pyplot.plot(p2.x[v], -omega*A*numpy.sin(omega*p2.time[pltind])*(w>p2.elev[v]), label='Analytical')
pyplot.legend()
pyplot.xlabel('x (m)')
pyplot.ylabel('x velocity (m/s)')
pyplot.title('Velocity near to time=3.25T')
pyplot.savefig('Vel_3_5T_v2.png')
