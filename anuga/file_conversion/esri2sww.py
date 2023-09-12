

class DataMissingValuesError(Exception):
    pass



def esri2sww(bath_dir,
                  elevation_dir,
                  ucur_dir,
                  vcur_dir,
                  sww_file,
                  minlat=None, maxlat=None,
                  minlon=None, maxlon=None,
                  zscale=1,
                  mean_stage=0,
                  fail_on_NaN=True,
                  elevation_NaN_filler=0,
                  bath_prefix='ba',
                  elevation_prefix='el',
                  verbose=False):
    """
    Produce an sww boundary file, from esri ascii data from CSIRO.

    Also convert latitude and longitude to UTM. All coordinates are
    assumed to be given in the GDA94 datum.

    assume:
    All files are in esri ascii format

    4 types of information
    bathymetry
    elevation
    u velocity
    v velocity

    Assumptions
    The metadata of all the files is the same
    Each type is in a seperate directory
    One bath file with extention .000
    The time period is less than 24hrs and uniform.
    """

    from anuga.file.netcdf import NetCDFFile

    from anuga.coordinate_transforms.redfearn import redfearn

    if sww_file[-4:] != '.sww':
        raise IOError('Output file %s should be of type .sww.' % sww_file)

    # So if we want to change the precision it's done here
    precision = netcdf_float 

    # go in to the bath dir and load the only file,
    bath_files = os.listdir(bath_dir)
    bath_file = bath_files[0]
    bath_dir_file =  bath_dir + os.sep + bath_file
    bath_metadata, bath_grid =  _read_asc(bath_dir_file)

    #Use the date.time of the bath file as a basis for
    #the start time for other files
    base_start = bath_file[-12:]

    #go into the elevation dir and load the 000 file
    elevation_dir_file = elevation_dir  + os.sep + elevation_prefix \
                         + base_start

    elevation_files = os.listdir(elevation_dir)
    ucur_files = os.listdir(ucur_dir)
    vcur_files = os.listdir(vcur_dir)
    elevation_files.sort()

    # the first elevation file should be the
    # file with the same base name as the bath data
    assert elevation_files[0] == 'el' + base_start

    number_of_latitudes = bath_grid.shape[0]
    number_of_longitudes = bath_grid.shape[1]
    number_of_volumes = (number_of_latitudes-1) * (number_of_longitudes-1) * 2

    longitudes = [bath_metadata['xllcorner'] + x*bath_metadata['cellsize'] \
                  for x in range(number_of_longitudes)]
    latitudes = [bath_metadata['yllcorner'] + y*bath_metadata['cellsize'] \
                 for y in range(number_of_latitudes)]

     # reverse order of lat, so the first lat represents the first grid row
    latitudes.reverse()

    kmin, kmax, lmin, lmax = get_min_max_indices(latitudes[:],longitudes[:],
                                                  minlat=minlat, maxlat=maxlat,
                                                  minlon=minlon, maxlon=maxlon)

    bath_grid = bath_grid[kmin:kmax,lmin:lmax]
    latitudes = latitudes[kmin:kmax]
    longitudes = longitudes[lmin:lmax]
    number_of_latitudes = len(latitudes)
    number_of_longitudes = len(longitudes)
    number_of_times = len(os.listdir(elevation_dir))
    number_of_points = number_of_latitudes * number_of_longitudes
    number_of_volumes = (number_of_latitudes-1) * (number_of_longitudes-1) * 2

    #Work out the times
    if len(elevation_files) > 1:
        # Assume: The time period is less than 24hrs.
        time_period = (int(elevation_files[1][-3:]) \
                       - int(elevation_files[0][-3:])) * 60*60
        times = [x*time_period for x in range(len(elevation_files))]
    else:
        times = [0.0]

    if verbose:
        log.critical('------------------------------------------------')
        log.critical('Statistics:')
        log.critical('  Extent (lat/lon):')
        log.critical('    lat in [%f, %f], len(lat) == %d'
                     % (min(latitudes), max(latitudes), len(latitudes)))
        log.critical('    lon in [%f, %f], len(lon) == %d'
                     % (min(longitudes), max(longitudes), len(longitudes)))
        log.critical('    t in [%f, %f], len(t) == %d'
                     % (min(times), max(times), len(times)))

    ######### WRITE THE SWW FILE #############

    # NetCDF file definition
    outfile = NetCDFFile(sww_file, netcdf_mode_w)

    #Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'Converted from XXX'

    #For sww compatibility
    outfile.smoothing = 'Yes'
    outfile.order = 1

    #Start time in seconds since the epoch (midnight 1/1/1970)
    outfile.starttime = starttime = times[0]

    # dimension definitions
    outfile.createDimension('number_of_volumes', number_of_volumes)
    outfile.createDimension('number_of_vertices', 3)
    outfile.createDimension('number_of_points', number_of_points)
    outfile.createDimension('number_of_timesteps', number_of_times)

    # variable definitions
    outfile.createVariable('x', precision, ('number_of_points',))
    outfile.createVariable('y', precision, ('number_of_points',))
    outfile.createVariable('elevation', precision, ('number_of_points',))

    #FIXME: Backwards compatibility
    #outfile.createVariable('z', precision, ('number_of_points',))
    #################################

    outfile.createVariable('volumes', netcdf_int, ('number_of_volumes',
                                                   'number_of_vertices'))

    outfile.createVariable('time', precision, ('number_of_timesteps',))

    outfile.createVariable('stage', precision, ('number_of_timesteps',
                                                'number_of_points'))

    outfile.createVariable('xmomentum', precision, ('number_of_timesteps',
                                                    'number_of_points'))

    outfile.createVariable('ymomentum', precision, ('number_of_timesteps',
                                                    'number_of_points'))

    #Store
    from anuga.coordinate_transforms.redfearn import redfearn

    x = num.zeros(number_of_points, float)  #Easting
    y = num.zeros(number_of_points, float)  #Northing

    if verbose: log.critical('Making triangular grid')

    #Get zone of 1st point.
    refzone, _, _ = redfearn(latitudes[0], longitudes[0])

    vertices = {}
    i = 0
    for k, lat in enumerate(latitudes):
        for l, lon in enumerate(longitudes):
            vertices[l,k] = i

            zone, easting, northing = redfearn(lat, lon)

            #msg = 'Zone boundary crossed at longitude =', lon
            #assert zone == refzone, msg
            #print '%7.2f %7.2f %8.2f %8.2f' %(lon, lat, easting, northing)
            x[i] = easting
            y[i] = northing
            i += 1

    #Construct 2 triangles per 'rectangular' element
    volumes = []
    for l in range(number_of_longitudes-1):    #X direction
        for k in range(number_of_latitudes-1): #Y direction
            v1 = vertices[l,k+1]
            v2 = vertices[l,k]
            v3 = vertices[l+1,k+1]
            v4 = vertices[l+1,k]

            #Note, this is different to the ferrit2sww code
            #since the order of the lats is reversed.
            volumes.append([v1,v3,v2]) #Upper element
            volumes.append([v4,v2,v3]) #Lower element

    volumes = num.array(volumes, int)      #array default#

    geo_ref = Geo_reference(refzone, min(x), min(y))
    geo_ref.write_NetCDF(outfile)

    # This will put the geo ref in the middle
    #geo_ref = Geo_reference(refzone, (max(x)+min(x))/2., (max(x)+min(y))/2.)

    if verbose:
        log.critical('------------------------------------------------')
        log.critical('More Statistics:')
        log.critical('  Extent (/lon):')
        log.critical('    x in [%f, %f], len(lat) == %d'
                     % (min(x), max(x), len(x)))
        log.critical('    y in [%f, %f], len(lon) == %d'
                     % (min(y), max(y), len(y)))
        log.critical('geo_ref: ', geo_ref)

    z = num.resize(bath_grid,outfile.variables['elevation'][:].shape)
    outfile.variables['x'][:] = x - geo_ref.get_xllcorner()
    outfile.variables['y'][:] = y - geo_ref.get_yllcorner()
    # FIXME (Ole): Remove once viewer has been recompiled and changed
    #              to use elevation instead of z
    #outfile.variables['z'][:] = z
    outfile.variables['elevation'][:] = z
    outfile.variables['volumes'][:] = volumes.astype(num.int32) # On Opteron 64

    stage = outfile.variables['stage']
    xmomentum = outfile.variables['xmomentum']
    ymomentum = outfile.variables['ymomentum']

    outfile.variables['time'][:] = times   #Store time relative

    if verbose: log.critical('Converting quantities')

    n = number_of_times
    for j in range(number_of_times):
        # load in files
        elevation_meta, elevation_grid = \
            _read_asc(elevation_dir + os.sep + elevation_files[j])

        _, u_momentum_grid = _read_asc(ucur_dir + os.sep + ucur_files[j])
        _, v_momentum_grid = _read_asc(vcur_dir + os.sep + vcur_files[j])

        #cut matrix to desired size
        elevation_grid = elevation_grid[kmin:kmax,lmin:lmax]
        u_momentum_grid = u_momentum_grid[kmin:kmax,lmin:lmax]
        v_momentum_grid = v_momentum_grid[kmin:kmax,lmin:lmax]

        # handle missing values
        missing = (elevation_grid == elevation_meta['NODATA_value'])
        if num.any (missing):
            if fail_on_NaN:
                msg = 'File %s contains missing values' \
                      % (elevation_files[j])
                raise(DataMissingValuesError, msg)
            else:
                elevation_grid = elevation_grid*(missing==0) \
                                 + missing*elevation_NaN_filler

        if verbose and j % ((n+10)//10) == 0: log.critical('  Doing %d of %d'
                                                          % (j, n))

        i = 0
        for k in range(number_of_latitudes):      #Y direction
            for l in range(number_of_longitudes): #X direction
                w = zscale*elevation_grid[k,l] + mean_stage
                stage[j,i] = w
                h = w - z[i]
                xmomentum[j,i] = u_momentum_grid[k,l]*h
                ymomentum[j,i] = v_momentum_grid[k,l]*h
                i += 1

    outfile.close()



def _read_asc(filename, verbose=False):
    """Read esri file from the following ASCII format (.asc)

    Example:

    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

    filename Path to the file to read.
    verbose True if this function is to be verbose.
    """

    datafile = open(filename)

    if verbose: log.critical('Reading DEM from %s' % filename)

    lines = datafile.readlines()
    datafile.close()

    if verbose: log.critical('Got %d lines' % len(lines))

    ncols = int(lines.pop(0).split()[1].strip())
    nrows = int(lines.pop(0).split()[1].strip())
    xllcorner = float(lines.pop(0).split()[1].strip())
    yllcorner = float(lines.pop(0).split()[1].strip())
    cellsize = float(lines.pop(0).split()[1].strip())
    NODATA_value = float(lines.pop(0).split()[1].strip())

    assert len(lines) == nrows

    #Store data
    grid = []

    n = len(lines)
    for i, line in enumerate(lines):
        cells = line.split()
        assert len(cells) == ncols
        grid.append(num.array([float(x) for x in cells]))
    grid = num.array(grid)

    return {'xllcorner':xllcorner,
            'yllcorner':yllcorner,
            'cellsize':cellsize,
            'NODATA_value':NODATA_value}, grid
