from anuga.utilities import plot_utils as util
from matplotlib import pyplot as pyplot
import numpy


verbose= True

swwfile = 'merewether_1m.sww'
p=util.get_output(swwfile)
p2=util.get_centroids(p)
# Time index at last time
tindex = len(p2.time)-1



if verbose: print 'calculating experimental transect'

x_data =    [ 0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0]
#vel =   [ 0.0, 0.0, 1.1, 3.2,  3.4, 2.4,  3.2,  3.2,  3.7,  3.1,  0.4,  0.0]
vel_data =   [ 0.0, 0.4, 3.1, 3.7,  3.2, 3.2,  2.4,  3.4,  3.2,  1.1,  0.0,  0.0]
#depth = [ 0.0, 0.0, 0.1, 0.5,  0.45, 0.4, 0.55, 0.1, 0.1,  0.05,  0.04, 0.0]
depth_data = [ 0.0, 0.04, 0.05, 0.1,  0.1, 0.55, 0.4, 0.45, 0.5,  0.1,  0.0, 0.0]

from scipy import interpolate

fvel = interpolate.interp1d(x_data, vel_data)
fdepth = interpolate.interp1d(x_data, depth_data)





if verbose: print 'calculating model heights at observation points'
# Get nearest wet points to 'point observations'
point_observations = numpy.genfromtxt(
    'Observations/ObservationPoints.csv',
    delimiter=",",skip_header=1)

nearest_points = point_observations[:,0]*0. - 1
for i in range(len(nearest_points)):
    # Compute distance of ANUGA points to observation, and
    # if the ANUGA point is dry then add a large value
    # Then find index of minimum
    n = ( (p2.x+p2.xllcorner-point_observations[i,0])**2 + \
          (p2.y+p2.yllcorner-point_observations[i,1])**2 + \
          (p2.stage[tindex,:] <= p2.elev)*1.0e+06).argmin()
    nearest_points[i] = n

f = open('Stage_point_comparison.csv','w')
f.writelines( 'Field, ANUGA, TUFLOW, ANUGA minus Field, ANUGA minus TUFLOW \n' )
for i in range(len(nearest_points)):
    po = point_observations[i,-2]
    tu = point_observations[i,-1]
    anuga_data = p2.stage[tindex, nearest_points.tolist()[i]]
    newline = str(round(po,2)) + ', ' + str(round(anuga_data,2)) + ', ' + str(tu) + ', ' + \
          str(round(anuga_data - po,2)) + ', ' + str(round(anuga_data - tu,2)) + '\n'
    f.writelines(newline)

f.flush()
f.close()


if verbose: print 'Plot transect'
## Plot transect 1 [need to guess appropriate end points as these are not so
## clear from the report]
xx=util.near_transect(p2,[103, 100.], [130.,80.],tol=0.5)
xx2=xx[0]

pyplot.clf()
pyplot.figure(figsize=(16,10.5))
pyplot.subplot(121)
pyplot.scatter(p2.x, p2.y, c=p2.elev,edgecolors='none')
# Add nice elevation data
colVals = numpy.maximum(numpy.minimum(p2.elev, 25.), 19.)
util.plot_triangles(p, values = colVals, edgecolors='none')

pyplot.gca().set_aspect('equal')
pyplot.scatter(p2.x[xx2],p2.y[xx2],color='green')
pyplot.xlim( (40., 160.))
pyplot.ylim( (0.,140.))
pyplot.title('Transect points in green')

pyplot.subplot(222)
pyplot.scatter(xx[1],p2.vel[tindex,xx[0]],color='green',label='model')
pyplot.scatter(xx[1],fvel(xx[1]),color='blue',label='data')
pyplot.legend(loc='upper left')
#pyplot.xlim(0,25)
pyplot.title('Final flow speed along the transect')

pyplot.subplot(224)
pyplot.scatter(xx[1],p2.stage[tindex,xx[0]]-p2.elev[xx[0]],color='green',label='model')
pyplot.scatter(xx[1],fdepth(xx[1]),color='blue',label='data')
pyplot.legend(loc='upper left')
#pyplot.xlim(0,25)
pyplot.title('Final depth along the transect')
pyplot.savefig('Transect1.png', bbox_inches='tight')


if verbose: print 'Plot velocity field'
pyplot.clf()

# Velocity vector plot
pyplot.figure(figsize=(16,22))
pyplot.scatter(p2.x,p2.y,c=(p2.elev>24.),edgecolors='none', s=0.2)
pyplot.gca().set_aspect('equal')
pyplot.xlim((100,180))
pyplot.ylim((100,210))
#k=range(0,len(p2.x),2) # Thin out the vectors for easier viewing
colVals = numpy.maximum(numpy.minimum(p2.elev, 25.), 19.)
util.plot_triangles(p, values = colVals, edgecolors='white')
k = range(len(p2.x))
# Thin out the triangles
#k = (((10.*(p2.x - p2.x.round())).round()%2 == 0.0)*((10.*(p2.y - p2.y.round())).round()%2 == 0.0)).nonzero()[0]
pyplot.quiver(p2.x[k],p2.y[k],p2.xvel[tindex,k], p2.yvel[tindex,k],
              scale_units='xy',units='xy',width=0.1,
              color='black',scale=1.0)
pyplot.savefig('velocity_stationary.png',dpi=100, bbox_inches='tight')



## Froude number plot
if verbose: print 'Plot Froude number plot'
pyplot.clf()
pyplot.figure(figsize=(6,8))
froude_number = p2.vel[tindex]/(numpy.maximum(p2.height[tindex], 1.0e-03)*9.8)**0.5
froude_category = (froude_number>1.).astype(float) + (froude_number > 0.).astype(float)
pyplot.scatter(p2.x,p2.y,edgecolors='none', s=0.2)

## Fake additions to plot to hack matplotlib legend
pyplot.scatter(0.,0., color='FireBrick',label='>1', marker='s')
pyplot.scatter(0.,0., color='PaleGreen',label='0-1', marker='s')
pyplot.scatter(0.,0., color='blue',label='0',marker='s')

pyplot.gca().set_aspect('equal')
util.plot_triangles(p, values = froude_category, edgecolors='none')
pyplot.xlim((p.x.min(), p.x.max()))
pyplot.ylim((p.y.min(), p.y.max()))
pyplot.title("Froude Number zones: 0, (0,1], or >1")

import matplotlib.patches as mpatches

#red_patch = mpatches.Patch(color='red', label='>1')
#green_patch = mpatches.Patch(color='green', label='(0-1]')
#blue_patch = mpatches.Patch(color='blue', label='0.')
#pyplot.legend(handles=[red_patch, green_patch, blue_patch], labels=['>1', '(0-1]', '0.'], loc='best')
pyplot.legend(loc='upper left')
pyplot.savefig('froudeNumber.png',dpi=100,bbox_inches='tight')





