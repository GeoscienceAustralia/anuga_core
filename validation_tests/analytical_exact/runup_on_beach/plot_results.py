"""View results of runup.py
"""
#---------------
# Import Modules
#---------------
import scipy
import anuga
from anuga.utilities import plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from numpy import zeros

p2=util.get_output('runup.sww', minimum_allowed_height=1.0e-03)
p=util.get_centroids(p2, velocity_extrapolation=True)

#------------------
# Select line
#------------------
py_central=p.y[scipy.argmin(abs(p.y-0.5))]
v=(p.y==p.y[py_central])

#--------------------
# Make plot animation
#--------------------
##pyplot.close() #If the plot is open, there will be problems
##if False:
##    line, = pyplot.plot( (p.x[v].min(),p.x[v].max()) ,(p.xvel[:,v].min(),p.xvel[:,v].max() ) )
##    for i in range(p.xmom.shape[0]):
##        line.set_xdata(p.x[v])
##        line.set_ydata(p.xvel[i,v])
##        pyplot.draw()
##        pyplot.plot( (0,1),(0,0), 'r' )
##        pyplot.title(str(i)+'/200') # : velocity does not converge to zero' )
##        pyplot.xlabel('x')
##        pyplot.ylabel('Velocity (m/s)')
##
##    pyplot.savefig('runup_x_velocities.png')

def analytic_sol(x):
    y = zeros(len(x))
    for i in range(len(x)):
        if x[i] < 0.8:
            y[i] = -0.5*x[i]
        else:
            y[i] = -0.4
    return y

#------------------------------------------------
# Maximum y velocities -- occurs in output step 3
#------------------------------------------------
pyplot.clf()
pyplot.plot(p.x[v],p.stage[5,v],'o',label='numerical stage')
pyplot.plot(p.x[v],p.elev[v],'k-',label='bed elevation')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.title('Free surface and bed at y = 0.5, time = 1.0 second')
pyplot.legend(loc='best')
pyplot.savefig('stage_1s.png')

pyplot.clf()
pyplot.plot(p.x[v],p.xvel[5,v],'-o', label='numerical')
pyplot.title('Xvelocity at y = 0.5, time = 1.0 second')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.legend(loc='best')
pyplot.savefig('xvel_1s.png')



pyplot.clf()
pyplot.plot(p.x[v],p.stage[150,v],'o-',label='numerical stage')
pyplot.plot(p.x[v], analytic_sol(p.x[v]), label='analytical stage')
pyplot.plot(p.x[v],p.elev[v],'k-',label='bed elevation')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.title('Free surface and bed at y = 0.5, time = 30.0 second')
pyplot.legend(loc='best')
pyplot.savefig('stage_30s.png')

pyplot.clf()
pyplot.plot(p.x[v],p.xvel[150,v],'-o', label='numerical')
pyplot.plot(p.x[v], zeros(len(p.x[v])), label='analytical')
pyplot.legend(loc='best')
pyplot.title('Xvelocity at y = 0.5, time = 30.0 second')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_30s.png')
