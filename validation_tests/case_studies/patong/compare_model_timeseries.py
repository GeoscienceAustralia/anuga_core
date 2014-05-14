"""Compare selected timeseries from model outputs

Usage:
python compare_model_timeseries.py swwfile1 swwfile2

Return 0 if timeseries are 'close enough', otherwise return non zero error code

Typically swwfile1 would be model output from unit test and swwfile2 would
be a reference model.
"""

from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.utilities.numerical_tools import cov, get_machine_precision


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

def usage():
    print 'Usage:'
    print 'python compare_model_timeseries.py swwfile1 swwfile2 epsilon'



def get_timeseries(sww_filename, gauges):
    """Generate time series for sww file based on gauges
    """

    gauge_locations = gauges.values()
    gauge_names = gauges.keys()

    tempfile = 'xyz1234tempfile.sww' # Has to end with sww

    os.system('cp %s %s' % (sww_filename, tempfile))
    f = file_function(tempfile, 
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


    
def compare_timeseries(timevector, 
                       timeseries1,
                       timeseries2,
                       name='',
                       legend = ('Time series 1', 'Time series 2'),
                       eps=1.0e-6):
    """Compare and plot two timeseries
    """

    
    timeseries1 = num.array(timeseries1)
    timeseries2 = num.array(timeseries2) 
    N = timevector.shape[0]
    assert timeseries1.shape[0] == N    
    assert timeseries2.shape[0] == N        

    print 'Testing gauge "%s"' % name
    print 'epsilon', eps
    # Covariance measure    
    res = cov(timeseries1-timeseries2)
    print '2 norm diff', res        
    msg = '2-norm of timeseries difference was too large: %e' % res
    assert res < eps, msg
     
    # Maximum norm
    nominator = max(abs(timeseries1-timeseries2))    
    denominator = max(abs(timeseries1))
    if denominator > 0:
        # Relative measure
        res = nominator/denominator
    else:
        # Absolute measure
        res = nominator/N

    print 'Max abs diff', res    
    msg = '%s: Difference between timeseries was too large: %e' % (name, res)
    assert res < eps, msg
    
    
    nominator = sum(abs(timeseries1-timeseries2))    
    denominator = sum(abs(timeseries1))
    if denominator > 0:
        # Relative measure
        res = nominator/denominator
    else:
        # Absolute measure
        res = nominator/N

    print 'Sum abs diff', res    
    msg = '%s: Difference between timeseries was too large: %e' % (name, res)
    assert res < eps, msg

    # Extrema
    max1 = max(timeseries1)
    max2 = max(timeseries2)
    res = abs(max1-max2)
    print 'Max diff', res
    msg = '%s: Difference between maxima was too large: %e' % (name, res)
    assert res < eps, msg    

    min1 = min(timeseries1)
    min2 = min(timeseries2)
    res = abs(min1-min2)    
    print 'Min diff', res
    msg = '%s: Difference between minima was too large: %e' % (name, res)
    assert res < eps, msg

    
    # Locations of extrema
    i1 = num.argmax(timeseries1)
    i2 = num.argmax(timeseries2)
    
    res = abs(timevector[i1]-timevector[i2])
    print 'Max loc diff', res
    msg = '%s: Difference between location of maxima was too large: %e' % (name, res)
    assert res < eps, msg

    i1 = num.argmin(timeseries1)
    i2 = num.argmin(timeseries2)
    
    res = abs(timevector[i1]-timevector[i2])    
    print 'Min loc diff', res    
    msg = '%s: Difference between location of minima was too large: %e' % (name, res)
    assert res < eps, msg
    

    if plotting:
        # Generate plots
        #pylab.ion() # No plotting on screen
        pylab.hold(False)
    
        pylab.plot(timevector, timeseries1, 'r-',
                   timevector, timeseries2, 'k-')

        pylab.title('Gauge %s' % name)
        pylab.xlabel('time(s)')
        pylab.ylabel('stage (m)')    
        pylab.legend(legend, shadow=True, loc='upper left')
        pylab.savefig(name, dpi = 300)        

        # Error vector
        #pylab.ion() # No plotting on screen
        pylab.hold(False)
        pylab.plot(timevector, timeseries1-timeseries2, 'b-')
        pylab.title('Gauge %s (difference)' % name)    
        pylab.xlabel('time(s)')
        pylab.ylabel('stage difference (m)')  
        pylab.savefig(name + '_diff', dpi = 300)                  


    print 'Gauge "%s" OK' % name
    print
    

def compare_all_timeseries(swwfile1, swwfile2, eps=1.0e-6):
    
    timevector1, timeseries1 = get_timeseries(swwfile1, gauges)
    timevector2, timeseries2 = get_timeseries(swwfile2, gauges)    

    msg = 'Time vectors were different in models %s and %s' %\
          (swwfile1, swwfile2)
    assert num.allclose(timevector1, timevector2), msg
    timevector = timevector1
    
    # Check that both timeseries exist for all gauges
    for name in timeseries2:
        assert name in timeseries1
    
    for name in timeseries1:
        assert name in timeseries2
            
    # Compare all timeseries data individually
    for name in timeseries1:
        compare_timeseries(timevector, 
                           timeseries1[name],
                           timeseries2[name],
                           name=name,
                           eps=eps)
        
        
        
            
    
    print 'All gauges OK'
    
    

if __name__ == '__main__':

    if len(sys.argv) != 4:
        usage()
        sys.exit(-1) 

    eps = float(sys.argv[3]) # Get tolerance to be used in comparisons
    res = compare_all_timeseries(sys.argv[1], sys.argv[2], eps=eps)
    
    #try:
    #    res = compare_all_timeseries(sys.argv[1], sys.argv[2], eps=eps)
    #except:
    #    print 'Failed'
    #    sys.exit(-1)
    #else:
    #    print 'Good', res
    #    sys.exit(0) 
