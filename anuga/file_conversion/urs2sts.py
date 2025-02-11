




import numpy as num

import anuga.utilities.log as log

# ANUGA modules
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.coordinate_transforms.geo_reference import Geo_reference, \
     write_NetCDF_georeference, ensure_geo_reference

# local modules
from anuga.file.mux import read_mux2_py
from anuga.file.mux import WAVEHEIGHT_MUX2_LABEL, EAST_VELOCITY_MUX2_LABEL, \
                NORTH_VELOCITY_MUX2_LABEL
from anuga.file.sts import Write_sts                
from functools import reduce
                

def urs2sts(basename_in, basename_out=None,
            weights=None,
            verbose=False,
            origin=None,
            zone=None,
            central_meridian=None,            
            mean_stage=0.0,
            zscale=1.0,
            ordering_filename=None):
    """Convert URS mux2 format for wave propagation to sts format

    Also convert latitude and longitude to UTM. All coordinates are
    assumed to be given in the GDA94 datum

    origin is a 3-tuple with geo referenced
    UTM coordinates (zone, easting, northing)

    inputs:

    basename_in: list of source file prefixes

        These are combined with the extensions:
        WAVEHEIGHT_MUX2_LABEL = '-z-mux2' for stage
        EAST_VELOCITY_MUX2_LABEL = '-e-mux2' xmomentum
        NORTH_VELOCITY_MUX2_LABEL = '-n-mux2' and ymomentum

        to create a 2D list of mux2 file. The rows are associated with each
        quantity and must have the above extensions
        the columns are the list of file prefixes.

    ordering: a .txt file name specifying which mux2 gauge points are
              to be stored. This is indicated by the index of the gauge
              in the ordering file.

              ordering file format:
              1st line:    'index,longitude,latitude\n'
              other lines: index,longitude,latitude

              If ordering is None or ordering file is empty then
               all points are taken in the order they
              appear in the mux2 file.


    output:
      basename_out: name of sts file in which mux2 data is stored.
      
      
      
    NOTE: South is positive in mux files so sign of y-component of velocity is reverted
    """

    import os
    from anuga.file.netcdf import NetCDFFile
    from operator import __and__

    if not isinstance(basename_in, list):
        if verbose: log.critical('Reading single source')
        basename_in = [basename_in]

    # This is the value used in the mux file format to indicate NAN data
    # FIXME (Ole): This should be changed everywhere to IEEE NAN when
    #              we upgrade to Numpy
    NODATA = 99

    # Check that basename is a list of strings
    if not reduce(__and__, [isinstance(z,str) for z in basename_in]):
        msg= 'basename_in must be a string or list of strings'
        raise IOError(msg)

    # Find the number of sources to be used
    numSrc = len(basename_in)

    # A weight must be specified for each source
    if weights is None:
        # Default is equal weighting
        weights = num.ones(numSrc, float)/numSrc
    else:
        weights = ensure_numeric(weights)
        msg = 'When combining multiple sources must specify a weight for ' \
              'mux2 source file'
        assert len(weights) == numSrc, msg

    if verbose: log.critical('Weights used in urs2sts: %s' % str(weights))

    # Check output filename
    if basename_out is None:
        msg = 'STS filename must be specified as basename_out ' \
              'in function urs2sts'
        raise Exception(msg)

    if basename_out.endswith('.sts'):
        stsname = basename_out
    else:
        stsname = basename_out + '.sts'

    # Create input filenames from basenames and check their existence
    files_in = [[], [], []]
    for files in basename_in:
        files_in[0].append(files + WAVEHEIGHT_MUX2_LABEL),
        files_in[1].append(files + EAST_VELOCITY_MUX2_LABEL)
        files_in[2].append(files + NORTH_VELOCITY_MUX2_LABEL)

    quantities = ['HA','UA','VA'] # Quantity names used in the MUX2 format
    for i in range(len(quantities)):
        for file_in in files_in[i]:
            if (os.access(file_in, os.R_OK) == 0):
                msg = 'File %s does not exist or is not accessible' % file_in
                raise Exception(msg)

    # Establish permutation array
    if ordering_filename is not None:
        if verbose is True: log.critical('Reading ordering file %s'
                                         % ordering_filename)

        # Read ordering file
        try:
            fid = open(ordering_filename, 'r')
            file_header = fid.readline().split(',')
            ordering_lines = fid.readlines()
            fid.close()
        except:
            msg = 'Cannot open %s' % ordering_filename
            raise Exception(msg)

        reference_header = 'index, longitude, latitude\n'
        reference_header_split = reference_header.split(',')
        for i in range(3):
            if not file_header[i].strip() == reference_header_split[i].strip():
                msg = 'File must contain header: ' + reference_header
                raise Exception(msg)

        if len(ordering_lines) < 2:
            msg = 'File must contain at least two points'
            raise Exception(msg)

        permutation = [int(line.split(',')[0]) for line in ordering_lines]
        permutation = ensure_numeric(permutation)
    else:
        permutation = None

    # Read MUX2 files
    if (verbose): log.critical('reading mux2 file')

    mux={}
    times_old = 0.0
    latitudes_old = 0.0
    longitudes_old = 0.0
    elevation_old = 0.0
    starttime_old = 0.0
    
    for i, quantity in enumerate(quantities):
        # For each quantity read the associated list of source mux2 file with
        # extention associated with that quantity

        times, latitudes, longitudes, elevation, mux[quantity], starttime \
            = read_mux2_py(files_in[i], weights, permutation, verbose=verbose)

        # Check that all quantities have consistent time and space information
        if quantity != quantities[0]:
            msg = '%s, %s and %s have inconsistent gauge data' \
                  % (files_in[0], files_in[1], files_in[2])
            assert num.allclose(times, times_old), msg
            assert num.allclose(latitudes, latitudes_old), msg
            assert num.allclose(longitudes, longitudes_old), msg
            assert num.allclose(elevation, elevation_old), msg
            assert num.allclose(starttime, starttime_old), msg
        times_old = times
        latitudes_old = latitudes
        longitudes_old = longitudes
        elevation_old = elevation
        starttime_old = starttime

        # Self check - can be removed to improve speed
        #ref_longitudes = [float(line.split(',')[1]) for line in ordering_lines]
        #ref_latitudes = [float(line.split(',')[2]) for line in ordering_lines]
        #
        #msg = 'Longitudes specified in ordering file do not match those ' \
        #      'found in mux files. ' \
        #      'I got %s instead of %s (only beginning shown)' \
        #      % (str(longitudes[:10]) + '...',
        #         str(ref_longitudes[:10]) + '...')
        #assert allclose(longitudes, ref_longitudes), msg
        #
        #msg = 'Latitudes specified in ordering file do not match those ' \
        #      'found in mux files. '
        #      'I got %s instead of %s (only beginning shown)' \
        #      % (str(latitudes[:10]) + '...',
        #         str(ref_latitudes[:10]) + '...')
        #assert allclose(latitudes, ref_latitudes), msg

    # Store timeseries in STS file
    msg = 'File is empty and or clipped region not in file region'
    assert len(latitudes > 0), msg

    number_of_points = latitudes.shape[0]      # Number of stations retrieved
    number_of_times = times.shape[0]           # Number of timesteps
    number_of_latitudes = latitudes.shape[0]   # Number latitudes
    number_of_longitudes = longitudes.shape[0] # Number longitudes

    # The permutation vector of contains original indices
    # as given in ordering file or None in which case points
    # are assigned the trivial indices enumerating them from
    # 0 to number_of_points-1
    if permutation is None:
        permutation = num.arange(number_of_points, dtype=int)

    # NetCDF file definition
    outfile = NetCDFFile(stsname, netcdf_mode_w)

    description = 'Converted from URS mux2 files: %s' % basename_in

    # Create new file
    sts = Write_sts()
    sts.store_header(outfile,
                     times+starttime,
                     number_of_points,
                     description=description,
                     verbose=verbose,
                     sts_precision=netcdf_float)

    # Store
    from anuga.coordinate_transforms.redfearn import redfearn

    x = num.zeros(number_of_points, float)  # Easting
    y = num.zeros(number_of_points, float)  # Northing

    # Check zone boundaries
    if zone is None:
        refzone, _, _ = redfearn(latitudes[0], longitudes[0],
                                 central_meridian=central_meridian)
    else:
        refzone = zone

    old_zone = refzone
    old_easting = 0.0
    old_northing = 0.0

    for i in range(number_of_points):
        computed_zone, easting, northing = redfearn(latitudes[i], longitudes[i],
                                                    zone=zone,
                                                    central_meridian=central_meridian)
        x[i] = easting
        y[i] = northing
        if computed_zone != refzone:
            msg = 'All sts gauges need to be in the same zone. \n'
            msg += 'offending gauge:Zone %d,%.4f, %4f\n' \
                   % (computed_zone, easting, northing)
            msg += 'previous gauge:Zone %d,%.4f, %4f' \
                   % (old_zone, old_easting, old_northing)
            raise Exception(msg)
        old_zone = computed_zone
        old_easting = easting
        old_northing = northing

    if origin is None:
        origin = Geo_reference(refzone, min(x), min(y))
    geo_ref = write_NetCDF_georeference(origin, outfile)

    elevation = num.resize(elevation, outfile.variables['elevation'][:].shape)
    outfile.variables['permutation'][:] = permutation.astype(num.int32) # Opteron 64
    outfile.variables['x'][:] = x - geo_ref.get_xllcorner()
    outfile.variables['y'][:] = y - geo_ref.get_yllcorner()
    outfile.variables['elevation'][:] = elevation

    stage = outfile.variables['stage']
    xmomentum = outfile.variables['xmomentum']
    ymomentum = outfile.variables['ymomentum']

    if verbose: log.critical('Converting quantities')

    for j in range(len(times)):
        for i in range(number_of_points):
            ha = mux['HA'][i,j]
            ua = mux['UA'][i,j]
            va = mux['VA'][i,j]
            if ha == NODATA:
                if verbose:
                    msg = 'Setting nodata value %d to 0 at time = %f, ' \
                          'point = %d' % (ha, times[j], i)
                    log.critical(msg)
                ha = 0.0
                ua = 0.0
                va = 0.0

            w = zscale*ha + mean_stage
            h = w - elevation[i]
            stage[j,i] = w

            xmomentum[j,i] = ua * h
            ymomentum[j,i] = -va * h # South is positive in mux files


    outfile.close()
    
    if verbose:
        log.critical('Wrote sts file ' + stsname)    
    

