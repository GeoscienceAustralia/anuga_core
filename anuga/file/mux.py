"""
    Read a mux2 file.
"""


from builtins import range
from anuga.utilities.numerical_tools import ensure_numeric
import numpy as num     
        
################################################################################
# READ MUX2 FILES line of points
################################################################################

WAVEHEIGHT_MUX_LABEL = '-z-mux'
EAST_VELOCITY_LABEL =  '-e-mux'
NORTH_VELOCITY_LABEL =  '-n-mux'

WAVEHEIGHT_MUX2_LABEL = '-z-mux2'
EAST_VELOCITY_MUX2_LABEL = '-e-mux2'
NORTH_VELOCITY_MUX2_LABEL = '-n-mux2'

def read_mux2_py(filenames,
                 weights=None,
                 permutation=None,
                 verbose=False):
    """Access the mux files specified in the filenames list. Combine the
       data found therin as a weighted linear sum as specifed by the weights.
       If permutation is None or empty extract timeseries data for all gauges
       within the files.

       Input:
           filenames:   List of filenames specifiying the file containing the
                        timeseries data (mux2 format) for each source
           weights:     Weighs associated with each source
                        (defaults to 1 for each source)
           permutation: Specifies the gauge numbers that for which data is to be
                        extracted
    """

    from .urs_ext import read_mux2

    numSrc = len(filenames)

    file_params = -1 * num.ones(3, float)                    # [nsta,dt,nt]

    # Convert verbose to int C flag
    if verbose:
        verbose = 1
    else:
        verbose = 0

    if weights is None:
        weights = num.ones(numSrc)

    if permutation is None:
        permutation = ensure_numeric([], int)

    # Call underlying C implementation urs2sts_ext.c
    cast_filenames = []
    for filename in filenames:
        cast_filenames.append(str(filename).encode())
    data = read_mux2(numSrc, cast_filenames, weights, file_params,
                     permutation, verbose)

    msg = 'File parameter values were not read in correctly from c file'
    assert len(num.compress(file_params > 0, file_params)) != 0, msg

    msg = 'The number of stations specifed in the c array and in the file ' \
          'are inconsistent'
    assert file_params[0] >= len(permutation), msg

    msg = 'The number of stations returned is inconsistent with ' \
          'the requested number'
    assert len(permutation) == 0 or len(permutation) == data.shape[0], msg

    nsta = int(file_params[0])
    msg = 'Must have at least one station'
    assert nsta > 0, msg

    dt = file_params[1]
    msg = 'Must have a postive timestep'
    assert dt > 0, msg

    nt = int(file_params[2])
    msg = 'Must have at least one gauge value'
    assert nt > 0, msg

    OFFSET = 5 # Number of site parameters p passed back with data
               # p = [geolat,geolon,depth,start_tstep,finish_tstep]

    # FIXME (Ole): What is the relationship with params and data.shape ?
    # It looks as if the following asserts should pass but they don't always
    #
    #msg = 'nt = %d, data.shape[1] == %d' %(nt, data.shape[1])
    #assert nt == data.shape[1] - OFFSET, msg
    #
    #msg = 'nsta = %d, data.shape[0] == %d' %(nsta, data.shape[0])
    #assert nsta == data.shape[0], msg

    # Number of stations in ordering file
    number_of_selected_stations = data.shape[0]

    # Index where data ends and parameters begin
    parameters_index = data.shape[1] - OFFSET

    times = dt * num.arange(parameters_index)
    latitudes = num.zeros(number_of_selected_stations, float)
    longitudes = num.zeros(number_of_selected_stations, float)
    elevation = num.zeros(number_of_selected_stations, float)
    quantity = num.zeros((number_of_selected_stations, parameters_index), \
                                                    float)

    starttime = 1e16
    for i in range(number_of_selected_stations):
        quantity[i][:] = data[i][:parameters_index]
        latitudes[i] = data[i][parameters_index]
        longitudes[i] = data[i][parameters_index+1]
        elevation[i] = -data[i][parameters_index+2]
        first_time_step = data[i][parameters_index+3]
        starttime = min(dt*first_time_step, starttime)

    return times, latitudes, longitudes, elevation, quantity, starttime



