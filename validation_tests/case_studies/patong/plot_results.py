"""Compare selected timeseries from model outputs

Usage:
python plot_results.py 

"""

from anuga import file_function
from anuga.utilities.numerical_tools import cov, get_machine_precision

import project

from os.path import join

import pylab
import numpy as num
import sys, os

try:
    import pylab
    pylab.hold(False)  # Check if this command can be issued
except:
    print 'Could not import pylab'
    plotting = False
else:
    # Create plots as png files
    plotting = True



# Model specific data
gauges = {'g10a': [422233.4, 874380.2],
          'g10b': [422107.9, 873556.8],
          'g10c': [421966.8, 872890.3],
          'g10d': [421606.1, 872106.2],
          'g11a': [422628.2, 873933.2],
          #'g11b': [422716.2, 873420.6],
          'g11c': [422689.1, 872859.8],
          'g11d': [422408.7, 871940.3],
          'north': [422572.4, 873992.6],
          'south': [422004.4, 872014.4]}


#---------------------------------------

def get_timeseries(sww_filename,gauges):
    """Generate time series for sww file based on gauges
    """

    gauge_locations = gauges.values()
    gauge_names = gauges.keys()

    #tempfile = 'xyz1234tempfile.sww' # Has to end with sww

    #os.system('cp %s %s' % (sww_filename, tempfile))

    f = file_function(sww_filename, 
                      quantities='stage',
                      interpolation_points=gauge_locations,
                      use_cache=True,
                      verbose=True)

    timevector = f.get_time()

    timeseries = {}
    for k, name in enumerate(gauge_names):
        model = timeseries[name] = []
        
        for t in timevector:
            model.append(f(t, point_id=k)[0])


    return num.array(timevector), timeseries


    
def plot_timeseries(timevector, 
                       timeseries,
                       name=''):
    """Compare and plot two timeseries
    """

    
    timeseries = num.array(timeseries)
    N = timevector.shape[0]
    assert timeseries.shape[0] == N    

    print 'Plotting gauge "%s"' % name

    if True:
        # Generate plots
        #pylab.ion() # No plotting on screen
        pylab.hold(False)
    
        pylab.plot(timevector, timeseries, 'r-')

        pylab.title('Gauge %s' % name)
        pylab.xlabel('time(s)')
        pylab.ylabel('stage (m)')    
        pylab.savefig(name, dpi = 300)        




timevector, timeseries = get_timeseries(join('outputs','patong.sww'), gauges )

print timevector.shape
print timeseries.keys()

for gauge, ts in timeseries.iteritems():
    plot_timeseries(timevector, ts, name=join('outputs',gauge+'.png'))


