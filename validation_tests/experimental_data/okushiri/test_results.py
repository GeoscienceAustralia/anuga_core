"""Verify that simulation produced by ANUGA compares to published
validation timeseries ch5, ch7 and ch9 as well as the boundary timeseries.

RMS norm is printed and plots are produced as png files.
No plots are shown on screen.
"""

import sys

import numpy as num
import anuga
from anuga.file.netcdf import NetCDFFile

import project
from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.numerical_tools import cov
#from anuga.utilities.numerical_tools import get_machine_precision
from anuga.shallow_water.sww_interrogate import get_maximum_inundation_elevation
from anuga.shallow_water.sww_interrogate import get_maximum_inundation_location


args = anuga.get_args()
verbose = args.verbose


testing = True

#-------------------------
# Basic data
#-------------------------

finaltime = 22.5
timestep = 0.05

gauge_locations = [[0.000, 1.696]] # Boundary gauge
gauge_locations += [[4.521, 1.196],  [4.521, 1.696],  [4.521, 2.196]] #Ch 5-7-9
gauge_names = ['Boundary', 'ch5', 'ch7', 'ch9']

validation_data = {}
for key in gauge_names:
    validation_data[key] = []


# Expected values
expected_covariance = {'Boundary': 5.269569575007607815e-05, 
                       'ch5': 1.166277999581819919e-04,
                       'ch7': 1.127136457890861503e-04,
                       'ch9': 1.250659477418482129e-04}

expected_difference = {'Boundary': 8.350712673810733924e-04,
                       'ch5': 3.405426180525532483e-03,
                       'ch7': 2.852870417368218517e-03,
                       'ch9': 3.248778982037564891e-03}

expected_maximum = {'Boundary': 1.611749508386188523e-02,  
                    'ch5': 3.551308418158714147e-02,
                    'ch7': 3.858418457126511908e-02,
                    'ch9': 4.317962986578308127e-02}

expected_minimum = {'Boundary': -1.164547474575844919e-02, 
                    'ch5': -8.664439185502026408e-03,
                    'ch7': -2.726335488279797541e-03,
                    'ch9': -5.977581218447349659e-03}

expected_argmax = {'Boundary': 1.255000000000000071e+01, 
                   'ch5': 1.839999999999999858e+01,
                   'ch7': 1.700000000000000000e+01,
                   'ch9': 1.685000000000000142e+01}

expected_argmin = {'Boundary': 2.064999999999999858e+01, 
                   'ch5': 1.459999999999999964e+01,
                   'ch7': 1.230000000000000071e+01,
                   'ch9': 1.315000000000000036e+01}

#-------------------------
# Read validation data
#-------------------------

if verbose: print('Reading', project.boundary_filename)

fid = NetCDFFile(project.boundary_filename, 'r')
input_time = fid.variables['time'][:]
validation_data['Boundary'] = fid.variables['stage'][:]

reference_time = []
fid = open(project.validation_filename)
lines = fid.readlines()
fid.close()

for i, line in enumerate(lines[1:]):
    if i == len(input_time): break
    
    fields = line.split()

    reference_time.append(float(fields[0]))    # Record reference time
    for j, key in enumerate(gauge_names[1:]):  # Omit boundary gauge
        value = float(fields[1:][j])           # Omit time 
        validation_data[key].append(value/100) # Convert cm2m


# Checks
assert reference_time[0] == 0.0
assert reference_time[-1] == finaltime
assert num.allclose(reference_time, input_time)

for key in gauge_names:
    validation_data[key] = ensure_numeric(validation_data[key])

#--------------------------------------------------
# Read and interpolate model output
#--------------------------------------------------

#if len(sys.argv) > 1:
#    sww_filename = sys.argv[1]
#else:
sww_filename = project.output_filename
    
f = file_function(sww_filename,
                  quantities='stage',
                  interpolation_points=gauge_locations,
                  use_cache=False,
                  verbose=verbose)


def report_difference(name, computed_value, reference_value, rtol, atol):

    if abs(reference_value) > 0:
        msg = '%s (expected, computed):\n  (%.18e, %.18e):\n  Relative error=%.18e'\
              %(name, reference_value, computed_value,
                abs(reference_value-computed_value)/reference_value)
        print(msg)
        

    msg = '  Absolute error=%.18e'\
          %(abs(reference_value-computed_value))        
    print(msg)

    
    #print 'Allclose:', allclose(reference_value, computed_value,
    #                            rtol=rtol, atol=atol)
    if testing is True:
        assert num.allclose(reference_value, computed_value,
                            rtol=rtol, atol=atol), msg
    


#--------------------------------------------------
# Compare model output to validation data
#--------------------------------------------------


#eps = get_machine_precision()

# Tolerances  for 20,000 triangles
rtol = 2.0e-2
atol = 2.0e-2

# Tolerances  for 60,000 triangles
#rtol = 1.0e-2
#atol = 1.0e-2

if verbose: print('Precisions used: rtol=%e, atol=%e' %(rtol, atol))


#print reference_time
for k, name in enumerate(gauge_names):

    sqsum = 0
    denom = 0
    model = []
    if verbose: 
        print() 
        print('Validating ' + name)
    observed_timeseries = validation_data[name]
    for i, t in enumerate(reference_time):
        model.append(f(t, point_id=k)[0])

    # Covariance measure    
    res = cov(observed_timeseries, model)
    if verbose:
        report_difference('Covariance', res, expected_covariance[name], rtol, atol)
     
    # Difference measures    
    res = sum(abs(observed_timeseries-model))/len(model)
    if verbose:
        report_difference('Accumulated difference', res,
                      expected_difference[name], rtol, atol)    

    # Extrema
    res = max(model)
    if verbose:
        report_difference('Maximum', res, expected_maximum[name], rtol, atol)
    
    res = min(model)
    if verbose:
        report_difference('Minimum', res, expected_minimum[name], rtol, atol)    

    # Locations of extrema
    #i0 = argmax(observed_timeseries)
    i1 = num.argmax(model)
    res = reference_time[i1]
    if verbose:
        report_difference('Location of maximum', res, expected_argmax[name], rtol, atol)    
    

    if not name in ['ch7', 'ch9']:
        # Minima of ch7 and ch9 are very flat and hard to pinpoint
        i1 = num.argmin(model)
        res = reference_time[i1]
        if verbose:
            report_difference('Location of minimum', res, expected_argmin[name],
                          rtol, atol)        


# Check max runup

q = get_maximum_inundation_elevation(sww_filename)
loc = get_maximum_inundation_location(sww_filename)

if verbose:
    print('Max runup elevation: ', q)
    print('Max runup elevation (scaled by 400): ', q*400)
    print('Max runup location:  ', loc)


