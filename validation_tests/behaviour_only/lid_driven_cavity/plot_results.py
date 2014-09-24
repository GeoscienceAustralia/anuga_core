"""View results of runup_sinusoid.py
"""
#---------------
# Import Modules
#---------------
from anuga.utilities import plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
import numpy

p1=util.get_output('dimensional_lid_driven.sww')#, 0.001)
#p1=util.get_output('runup_sinusoid.sww')
p2=util.get_centroids(p1)#, velocity_extrapolation=True)
#p2=util.get_centroids(p1)

#------------------
# Select line
#------------------
#v=(p.y==0.5)
v=(p2.y==p2.y[80])

#--------------------
# Make plot animation
#--------------------
pyplot.close() #If the plot is open, there will be problems

# Plot vertex values
pyplot.clf()
t1=int(len(p1.time)/2)
t2=-1
#pyplot.scatter(p1.x,p1.y,c=p1.elev,edgecolors='none', s=25)
#pyplot.colorbar()
pyplot.quiver(p1.x,p1.y,p1.xvel[t1,:],p1.yvel[t1,:])
pyplot.title('The maximum VERTEX speed is '+ str(p1.vel[t1,:].max()) + ' m/s at time '+ str(p1.time[t1])+' s')
pyplot.xlabel('Xposition')
pyplot.ylabel('Yposition')
pyplot.savefig('vel_t1_vertex.png')

# Plot vertex values
pyplot.clf()
#pyplot.scatter(p2.x,p2.y,c=p2.elev,edgecolors='none', s=25)
#pyplot.colorbar()
pyplot.quiver(p2.x,p2.y,p2.xvel[t1,:],p2.yvel[t1,:])
pyplot.title('The maximum CENTROID speed is '+ str(p2.vel[t1,:].max()) + ' m/s at time ' + str(p1.time[t1]) + ' s')
pyplot.xlabel('Xposition')
pyplot.ylabel('Yposition')
pyplot.savefig('vel_t1_centroid.png')

# Plot vertex values
pyplot.clf()
#pyplot.scatter(p1.x,p1.y,c=p1.elev,edgecolors='none', s=25)
#pyplot.colorbar()
pyplot.quiver(p1.x,p1.y,p1.xvel[t2,:],p1.yvel[t2,:])
pyplot.title('The maximum VERTEX speed is '+ str(p1.vel[t2,:].max()) + ' m/s at time ' + str(p1.time[t2]) +' s')
pyplot.xlabel('Xposition')
pyplot.ylabel('Yposition')
pyplot.savefig('vel_t2_vertex.png')

# Plot vertex values
pyplot.clf()
#pyplot.scatter(p2.x,p2.y,c=p2.elev,edgecolors='none', s=25)
#pyplot.colorbar()
pyplot.quiver(p2.x,p2.y,p2.xvel[t2,:],p2.yvel[t2,:])
pyplot.title('The maximum CENTROID speed is '+ str(p2.vel[t2,:].max()) + ' m/s at time ' + str(p1.time[t2])+ ' s')
pyplot.xlabel('Xposition')
pyplot.ylabel('Yposition')
pyplot.savefig('vel_t2_centroid.png')


v=(p2.y>0.49)*(p2.y<0.51)
pyplot.clf()
pyplot.plot(p2.x[v], p2.xvel[t2,v],'b.', label='Xvelocity')
pyplot.title('Xvelocity at x=0.5')
pyplot.legend(loc='best')
pyplot.xlabel('Yposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvelocity_at_y05.png')
#pyplot.show()
pyplot.close()
