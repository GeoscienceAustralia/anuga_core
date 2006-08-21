"""Read in sww file, interpolate at specified locations and plot time series

"""

from anuga.pyvolution.util import file_function
from anuga.coordinate_transforms.redfearn import degminsec2decimal_degrees, redfearn
from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


swwfile = 'island.sww' 
#gauges = [[56, 48], [58, 48], [60, 48], [62, 48], [64, 48]]
#gauges = [[10, 10], [60, 48], [62, 48], [64, 48]]    
gauges = [[10, 10], [40, 10], [70, 10]]    


store = True

#Read model output
quantities = ['stage', 'elevation', 'xmomentum', 'ymomentum']
f = file_function(swwfile,
                  quantities = quantities,
                  interpolation_points = gauges,
                  verbose = True,
                  use_cache = True)


T = f.get_time()
               
from math import sqrt, atan, degrees
from Numeric import ones
N = len(gauges)
for k, g in enumerate(gauges):
    if k%((N+10)/10)==0: # diagnostics - print 10 lines
        print 'Doing row %d of %d' %(k, N)

    model_time = []
    stages = []
    elevations = []
    momenta = []
    velocity = []
    xmom = []
    ymom = []
    bearings = []
    depths = []

    max_depth = 0
    max_momentum = 0
    max_velocity = 0

    due_east = 90.0*ones([len(T)])
    due_west = 270.0*ones([len(T)])
    maxT = max(T)
    tstep = maxT/(len(T)-1)


    myloc = 'G%d' %k #locations[k]

    for i, t in enumerate(T): # T is a list of times
        #if tmin < t < tmax:
        w = f(t, point_id = k)[0]
        z = f(t, point_id = k)[1]
        uh = f(t, point_id = k)[2]
        vh = f(t, point_id = k)[3]
        depth = w-z

        m = sqrt(uh*uh + vh*vh)   #Absolute momentum
        vel = sqrt(uh*uh + vh*vh) / (w-z + 1.e-30) #Absolute velocity
        angle = degrees(atan(vh/(uh+1.e-15)))
        if (0 < angle < 90.0):
            if vh > 0:
                bearing = 90.0 - abs(angle)
            if vh < 0:
                bearing = 270.0 - abs(angle)
        if (-90 < angle < 0):
            if vh < 0:
                bearing = 90.0 - abs(angle)
            if vh > 0:
                bearing = 270.0 - abs(angle)
        if angle == 0:
            bearing = 0.0
                
        model_time.append(t)        
        stages.append(w)
        elevations.append(z)  #Should be constant over time
        momenta.append(m)
        velocity.append(vel)
        xmom.append(uh)
        ymom.append(vh)
        bearings.append(bearing)
        depths.append(depth)

        if w-z > max_depth:
            max_depth = w-z
        if m > max_momentum:            
            max_momentum = m
        if vel > max_velocity:
            max_velocity = vel

    #Plot only those gauges that have been inundated by more than a threshold
    #if max_depth < 0.2:
    #    print 'Skipping gauge %d' %k
    #    continue

    ion()
    hold(False)

    if elevations[0] <= 0:
        plot(model_time, stages, '-b')
    else:    
        plot(model_time, stages, '-b',
             model_time, elevations, '-k')
    #axis([0, 100, z-0.01, z+0.005])
        
    #name = 'Gauge_%d: (%.1f, %.1f)' %(k, g[0], g[1])
    name = 'Gauge_%d: (%.1f, %.1f) Location: %s' %(k, g[0], g[1], myloc)
    title(name)

    title('%s (stage)' %name)
    xlabel('time [s]')
    ylabel('elevation [m]')    
    legend(('Stage', 'Bed = %.1f' %elevations[0]),
           shadow=True,
           loc='upper right')

    if store is True: savefig('Gauge_%s_stage' %myloc)

    raw_input('Next')

    """
    #Momentum plot
    ion()
    hold(False)
    plot(model_time, momenta, '-r')
    title(name)

    title('%s (momentum)' %name)
    xlabel('time [s]')
    ylabel('sqrt( uh^2 + vh^2 ) [m^2/s]')    
    #savefig('Gauge_%d_momentum' %k)
    if store is True: savefig('Gauge_%s_momentum' %myloc)
    
    raw_input('Next')
    """

    """
    #Bearing plot
    ion()
    hold(False)
    ax = plot(model_time, bearings, '-b', model_time, due_west, '-.b',
         model_time, due_east, '-.b')
    title(name)
    ax = axis([0, maxT, 0, 360])
    text(maxT+tstep, 90, 'East')
    text(maxT+tstep, 270, 'West')
    #majorLocator = MultipleLocator(3600)
    #print 'major', majorLocator[1]
    #ax.xaxis.set_major_locator(majorLocator) #'yticklabels', range(30,390,30))
    # set(labels,color='g',rotation=45)

    title('%s (bearing)' %name)
    xlabel('time [s]')
    ylabel(' atan(vh/uh) [degrees from North]')    
    #savefig('Gauge_%d_bearing' %k)
    if store is True: savefig('Gauge_%s_bearing' %myloc)
    
    raw_input('Next')
    """

    """
    #Speed plot
    ion()
    hold(False)
    plot(model_time, velocity, '-r')
    title(name)

    title('%s (velocity)' %name)
    xlabel('time [s]')
    ylabel('sqrt( uh^2 + vh^2 ) / depth [m/s]')    
    #savefig('Gauge_%d_speed' %k)
    if store is True: savefig('Gauge_%s_speed' %myloc)
    
    raw_input('Next')
    """
    

    whichone = '_%s' %myloc
    thisfile = 'island_gauge'+whichone+'.csv'
    fid = open(thisfile, 'w')
    s = 'Time, Depth, Momentum, Velocity \n'
    fid.write(s)
    for i_t, i_d, i_m, i_vel in zip(model_time, depths, momenta, velocity):
        s = '%.2f, %.2f, %.2f, %.2f\n' %(i_t, i_d, i_m, i_vel)
        fid.write(s)
        
show()


