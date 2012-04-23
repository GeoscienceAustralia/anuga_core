"""View results of runup.py
"""
#---------------
# Import Modules
#---------------
import anuga
from anuga.utilities import plot_utils as util
import numpy
import scipy
from matplotlib import pyplot as pyplot
#import util # Routines to read in and work with ANUGA output

p2=util.get_output('runup_v2.sww', minimum_allowed_height=1.0e-03)
p=util.get_centroids(p2, velocity_extrapolation=True)

#p=util.get_output('runup_v2.sww', minimum_allowed_height=1.0e-03)

#------------------
# Select line
#------------------
py_central=p.y[scipy.argmin(abs(p.y-0.5))]
v=(p.y==p.y[py_central])

#--------------------
# Make plot animation
#--------------------
pyplot.close() #If the plot is open, there will be problems
pyplot.ion()

if False:
    line, = pyplot.plot( (p.x[v].min(),p.x[v].max()) ,(p.xvel[:,v].min(),p.xvel[:,v].max() ) )
    for i in range(p.xmom.shape[0]):
        line.set_xdata(p.x[v])
        line.set_ydata(p.xvel[i,v])
        pyplot.draw()
        pyplot.plot( (0,1),(0,0), 'r' )
        pyplot.title(str(i)+'/200') # : velocity does not converge to zero' )
        pyplot.xlabel('x')
        pyplot.ylabel('Velocity (m/s)')

    pyplot.savefig('runup_x_velocities.png')

#pyplot.clf()
#pyplot.close()

#------------------------------------------------
# Maximum y velocities -- occurs in output step 3
#------------------------------------------------
#print yvel[3,:].max(), yvel[3,:].min()
#highx=yvel[3,:].argmax()
#v=(x==x[highx])
#pyplot.plot(yvel[3,v])
#pyplot.title('y-velocity is not always zero at the boundaries, e.g. x='+str(x[highx])+' , t=0.3s')
#pyplot.xlabel('y')
#pyplot.ylabel('Velocity (m/s)')
#pyplot.savefig('runup_y_velocities.png')

pyplot.clf()
pyplot.plot(p.x[v],p.stage[5,v])
pyplot.plot(p.x[v],p.stage[5,v],'o')
pyplot.plot(p.x[v],p.elev[v])
pyplot.xlabel('x (m)')
pyplot.ylabel('z (m)')
pyplot.title('Free surface and bed at y==0.5, time = 1.0 second')
pyplot.savefig('elev_1s_v2.png')

pyplot.clf()
pyplot.plot(p.x[v],p.xvel[5,v])
pyplot.plot(p.x[v],p.xvel[5,v],'o')
pyplot.title('Velocity at y==0.500, time = 1.0 second')
pyplot.xlabel('x (m)')
pyplot.ylabel('Velocity (m/s)')
pyplot.savefig('vel1d_1s_v2.png')

#pyplot.clf()
#pyplot.plot(x_cent[y_cent==y_cent[80]],xvel_cent[15,y_cent==y_cent[80]])
#pyplot.plot(x_cent[y_cent==y_cent[80]],xvel_cent[15,y_cent==y_cent[80]],'o')
#pyplot.title('Velocity at y==0.500, time = 3.0 second')
#pyplot.xlabel('x (m)')
#pyplot.ylabel('Velocity (m/s)')
#pyplot.savefig('vel1d_3s_v2.png')
#-------------------------------------
# Final velocities plot
#-------------------------------------
#pyplot.clf()
#pyplot.quiver(x,y,xvel[200,:],yvel[200,:])
#pyplot.xlabel('x')
#pyplot.ylabel('y')
#pyplot.title('The maximum speed is '+ str(vel[200,:].max()) + ' m/s')
#pyplot.savefig('final_vel_field.png')
#print vel[200,:].max()

pyplot.clf()
pyplot.scatter(p.x,p.y,c=p.elev,edgecolors='none', s=25)
pyplot.colorbar()
#pyplot.quiver(x_cent,y_cent,xvel_cent[15,:],yvel_cent[15,:])
#pyplot.title('The maximum speed is '+ str(vel_cent[15,:].max()) + ' m/s at time 3.0s')
pyplot.quiver(p.x,p.y,p.xvel[5,:],p.yvel[5,:])
pyplot.title('The maximum speed is '+ str(p.vel[5,:].max()) + ' m/s at time 1.0s')
pyplot.savefig('vel_1s_v2.png')

pyplot.clf()
pyplot.scatter(p.x,p.y,c=p.elev,edgecolors='none', s=25)
pyplot.colorbar()
#pyplot.quiver(x_cent,y_cent,xvel_cent[150,:],yvel_cent[150,:])
#pyplot.title('The maximum speed is '+ str(vel_cent[150,:].max()) + ' m/s at time 30.0s')
pyplot.quiver(p.x,p.y,p.xvel[150,:],p.yvel[150,:])
pyplot.title('The maximum speed is '+ str(p.vel[150,:].max()) + ' m/s at time 30.0s')
pyplot.savefig('vel_30s_v2.png')

pyplot.clf()
#pyplot.plot(p.x[v],p.stage[150,v])
pyplot.plot(p.x[v],p.stage[150,v],'o-',label='numerical')
pyplot.plot(p.x[v],p.elev[v])
pyplot.plot([0.,1.], [-0.1, -0.1], label='Analytical (wet regions)')
pyplot.xlabel('x (m)')
pyplot.ylabel('z (m)')
pyplot.title('Free surface and bed at y==0.5, time = 30.0 second')
pyplot.legend()
pyplot.savefig('elev_30s_v2.png')
pyplot.clf()
#pyplot.plot(p.x[v],p.xvel[150,v])
pyplot.plot(p.x[v],p.xvel[150,v],'-o', label='numerical')
pyplot.plot([0.,1.], [-0.0, -0.0], label='Analytical')
pyplot.legend()
pyplot.title('Velocity at y==0.500, time = 30.0 second')
pyplot.xlabel('x (m)')
pyplot.ylabel('Velocity (m/s)')
pyplot.savefig('vel1d_30s_v2.png')
#pyplot.clf()
#pyplot.plot(x_cent[y_cent==y_cent[80]],xvel_cent[150,y_cent==y_cent[80]])
#pyplot.plot(x_cent[y_cent==y_cent[80]],xvel_cent[150,y_cent==y_cent[80]],'o')
#pyplot.title('Velocity at y==0.500, time = 30.0 second')
#pyplot.xlabel('x (m)')
#pyplot.ylabel('Velocity (m/s)')
#pyplot.savefig('vel1d_30s_v2.png')
