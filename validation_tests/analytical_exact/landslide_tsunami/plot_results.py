#from anuga.utilities import plot_utils as util
from anuga.utilities import plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
import scipy, numpy

p= util.get_output('runup_v2.sww')
p2=util.get_centroids(p,velocity_extrapolation=True)

# Find points with x<1km
v = p2.x<1000.

# Get analytical solution snapshots
t160=numpy.genfromtxt('./DATA/t160.csv',delimiter=',',skip_header=1)
t175=numpy.genfromtxt('./DATA/t175.csv',delimiter=',',skip_header=1)
t220=numpy.genfromtxt('./DATA/t220.csv',delimiter=',',skip_header=1)
shoreline=numpy.genfromtxt('./DATA/Shoreline.csv',delimiter=',',skip_header=1)

# Offset x by 200m, to reflect how the mesh is set up
t160[:,0] = t160[:,0]+200.
t175[:,0] = t175[:,0]+200.
t220[:,0] = t220[:,0]+200.

# Make plots
t_wanted=[160., 175., 220.]

pyplot.clf()
for i,t in enumerate(t_wanted):
    tmp=abs(p2.time-t)
    p2_ind=tmp.argmin()
    pyplot.plot(p2.x[v], p2.stage[p2_ind,v], 'o',label='Numerical')
    if(i==0):
        pyplot.plot(t160[:,0], t160[:,1],'-', color='red', label='Analytical')
    if(i==1):
        pyplot.plot(t175[:,0], t175[:,1],'-', color='red', label='Analytical' )
    if(i==2):
        pyplot.plot(t220[:,0], t220[:,1],'-',color='red', label='Analytical' )

    pyplot.plot(p2.x[v], p2.elev[v], '--', color='green')
    pyplot.title('Stage, Time ' + str(round(t,2)))
    pyplot.legend(loc='best')
    pyplot.ylim((-25,20))
    pyplot.savefig('figure_stage_'+str(i)+'.png')
    pyplot.clf()

    pyplot.plot(p2.x[v], p2.xvel[p2_ind,v], 'o', label='Numerical')
    if(i==0):
        pyplot.plot(t160[:,0], t160[:,2],'-',color='red', label='Analytical')
    if(i==1):
        pyplot.plot(t175[:,0], t175[:,2],'-',color='red', label='Analytical')
    if(i==2):
        pyplot.plot(t220[:,0], t220[:,2],'-',color='red', label='Analytical')
    pyplot.ylim((-25,20))
    pyplot.title('Velocity, Time ' + str(round(t,2)))
    pyplot.legend(loc='best')
    pyplot.savefig('figure_vel_'+str(i)+'.png')
    pyplot.clf()

# Plot shoreline position
model_shore_x=p2.time*0.
model_shore_u=p2.time*0.

# Hacks to identify the location of the shoreline
shoreline_depth_thresh1=0.01
shoreline_depth_thresh2=4.0

for i in range(len(model_shore_x)):
    # Compute index of shoreline
    vtmp = p2.stage[i,:]>p2.elev+shoreline_depth_thresh1

    model_shore_x[i]=p2.x[vtmp].min()

    # Extract shoreline velocity. It is tricky, to avoid getting
    # zero velocity points -- might be a better way to do it?
    vtmp2= (p2.stage[i,:]>p2.elev+shoreline_depth_thresh1)*\
           (p2.stage[i,:]<p2.elev+shoreline_depth_thresh2)

    #print 'vtmp2 ', numpy.any(vtmp2)
    mloc=abs(p2.xvel[i,vtmp2]).argmax()
    model_shore_u[i]=p2.xvel[i,vtmp2][mloc]

pyplot.plot(p2.time, model_shore_x-200.,'o', label='Numerical')
pyplot.plot(shoreline[:,0], shoreline[:,1],'-', color='red', label='Analytical')
pyplot.legend(loc='best')
pyplot.title('Shoreline position (where depth >'+str(shoreline_depth_thresh1)+')')
pyplot.savefig('Shoreline_position.png')

# Can plot velocity as well 
pyplot.clf()
pyplot.plot(p2.time, model_shore_u,'o',label='Numerical')
pyplot.plot(shoreline[:,0], shoreline[:,2],'-',color='red',label='Analytical')
pyplot.title('Shoreline velocity: Peak speed where depth >'\
             + str(shoreline_depth_thresh1) + ' and depth < '\
             + str(shoreline_depth_thresh2))
pyplot.legend(loc='best')
pyplot.savefig('Shoreline_velocity.png')
