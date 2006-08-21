






###############################################
#OBSOLETE STUFF
#Native checkpoint format.
#Information needed to recreate a state is preserved
#FIXME: Rethink and maybe use netcdf format
def cpt_variable_writer(filename, t, v0, v1, v2):
    """Store all conserved quantities to file
    """

    M, N = v0.shape

    FN = create_filename(filename, 'cpt', M, t)
    #print 'Writing to %s' %FN

    fid = open(FN, 'w')
    for i in range(M):
        for j in range(N):
            fid.write('%.16e ' %v0[i,j])
        for j in range(N):
            fid.write('%.16e ' %v1[i,j])
        for j in range(N):
            fid.write('%.16e ' %v2[i,j])

        fid.write('\n')
    fid.close()


def cpt_variable_reader(filename, t, v0, v1, v2):
    """Store all conserved quantities to file
    """

    M, N = v0.shape

    FN = create_filename(filename, 'cpt', M, t)
    #print 'Reading from %s' %FN

    fid = open(FN)


    for i in range(M):
        values = fid.readline().split() #Get one line

        for j in range(N):
            v0[i,j] = float(values[j])
            v1[i,j] = float(values[3+j])
            v2[i,j] = float(values[6+j])

    fid.close()

def cpt_constant_writer(filename, X0, X1, X2, v0, v1, v2):
    """Writes x,y,z,z,z coordinates of triangles constituting the bed
    elevation.
    FIXME: Not in use pt
    """

    M, N = v0.shape


    print X0
    import sys; sys.exit()
    FN = create_filename(filename, 'cpt', M)
    print 'Writing to %s' %FN

    fid = open(FN, 'w')
    for i in range(M):
        for j in range(2):
            fid.write('%.16e ' %X0[i,j])   #x, y
        for j in range(N):
            fid.write('%.16e ' %v0[i,j])       #z,z,z,

        for j in range(2):
            fid.write('%.16e ' %X1[i,j])   #x, y
        for j in range(N):
            fid.write('%.16e ' %v1[i,j])

        for j in range(2):
            fid.write('%.16e ' %X2[i,j])   #x, y
        for j in range(N):
            fid.write('%.16e ' %v2[i,j])

        fid.write('\n')
    fid.close()



#Function for storing out to e.g. visualisation
#FIXME: Do we want this?
#FIXME: Not done yet for this version
def dat_constant_writer(filename, X0, X1, X2, v0, v1, v2):
    """Writes x,y,z coordinates of triangles constituting the bed elevation.
    """

    M, N = v0.shape

    FN = create_filename(filename, 'dat', M)
    #print 'Writing to %s' %FN

    fid = open(FN, 'w')
    for i in range(M):
        for j in range(2):
            fid.write('%f ' %X0[i,j])   #x, y
        fid.write('%f ' %v0[i,0])       #z

        for j in range(2):
            fid.write('%f ' %X1[i,j])   #x, y
        fid.write('%f ' %v1[i,0])       #z

        for j in range(2):
            fid.write('%f ' %X2[i,j])   #x, y
        fid.write('%f ' %v2[i,0])       #z

        fid.write('\n')
    fid.close()



def dat_variable_writer(filename, t, v0, v1, v2):
    """Store water height to file
    """

    M, N = v0.shape

    FN = create_filename(filename, 'dat', M, t)
    #print 'Writing to %s' %FN

    fid = open(FN, 'w')
    for i in range(M):
        fid.write('%.4f ' %v0[i,0])
        fid.write('%.4f ' %v1[i,0])
        fid.write('%.4f ' %v2[i,0])

        fid.write('\n')
    fid.close()


def read_sww(filename):
    """Read sww Net CDF file containing Shallow Water Wave simulation

    The integer array volumes is of shape Nx3 where N is the number of
    triangles in the mesh.

    Each entry in volumes is an index into the x,y arrays (the location).

    Quantities stage, elevation, xmomentum and ymomentum are all in arrays of dimensions
    number_of_timesteps, number_of_points.

    The momentum is not always stored.

    """
    from Scientific.IO.NetCDF import NetCDFFile
    print 'Reading from ', filename
    fid = NetCDFFile(filename, 'r')    #Open existing file for read
#latitude, longitude
    # Get the variables as Numeric arrays
    x = fid.variables['x']             #x-coordinates of vertices
    y = fid.variables['y']             #y-coordinates of vertices
    z = fid.variables['elevation']     #Elevation
    time = fid.variables['time']       #Timesteps
    stage = fid.variables['stage']     #Water level
    #xmomentum = fid.variables['xmomentum']   #Momentum in the x-direction
    #ymomentum = fid.variables['ymomentum']   #Momentum in the y-direction

    volumes = fid.variables['volumes'] #Connectivity

    #FIXME (Ole): What is this?
    #             Why isn't anything returned?
    #             Where's the unit test?



def sww2asc_obsolete(basename_in, basename_out = None,
            quantity = None,
            timestep = None,
            reduction = None,
            cellsize = 10,
            verbose = False,
            origin = None):
    """Read SWW file and convert to Digitial Elevation model format (.asc)

    Example:

    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

    Also write accompanying file with same basename_in but extension .prj
    used to fix the UTM zone, datum, false northings and eastings.

    The prj format is assumed to be as

    Projection    UTM
    Zone          56
    Datum         WGS84
    Zunits        NO
    Units         METERS
    Spheroid      WGS84
    Xshift        0.0000000000
    Yshift        10000000.0000000000
    Parameters


    if quantity is given, out values from quantity otherwise default to
    elevation

    if timestep (an index) is given, output quantity at that timestep

    if reduction is given use that to reduce quantity over all timesteps.

    """
    from Numeric import array, Float, concatenate, NewAxis, zeros,\
         sometrue
    from utilities.polygon import inside_polygon

    #FIXME: Should be variable
    datum = 'WGS84'
    false_easting = 500000
    false_northing = 10000000

    if quantity is None:
        quantity = 'elevation'

    if reduction is None:
        reduction = max

    if basename_out is None:
        basename_out = basename_in + '_%s' %quantity

    swwfile = basename_in + '.sww'
    ascfile = basename_out + '.asc'
    prjfile = basename_out + '.prj'


    if verbose: print 'Reading from %s' %swwfile
    #Read sww file
    from Scientific.IO.NetCDF import NetCDFFile
    fid = NetCDFFile(swwfile)

    #Get extent and reference
    x = fid.variables['x'][:]
    y = fid.variables['y'][:]
    volumes = fid.variables['volumes'][:]

    ymin = min(y); ymax = max(y)
    xmin = min(x); xmax = max(x)

    number_of_timesteps = fid.dimensions['number_of_timesteps']
    number_of_points = fid.dimensions['number_of_points']
    if origin is None:

        #Get geo_reference
        #sww files don't have to have a geo_ref
        try:
            geo_reference = Geo_reference(NetCDFObject=fid)
        except AttributeError, e:
            geo_reference = Geo_reference() #Default georef object

        xllcorner = geo_reference.get_xllcorner()
        yllcorner = geo_reference.get_yllcorner()
        zone = geo_reference.get_zone()
    else:
        zone = origin[0]
        xllcorner = origin[1]
        yllcorner = origin[2]


    #Get quantity and reduce if applicable
    if verbose: print 'Reading quantity %s' %quantity

    if quantity.lower() == 'depth':
        q = fid.variables['stage'][:] - fid.variables['elevation'][:]
    else:
        q = fid.variables[quantity][:]


    if len(q.shape) == 2:
        if verbose: print 'Reducing quantity %s' %quantity
        q_reduced = zeros( number_of_points, Float )

        for k in range(number_of_points):
            q_reduced[k] = reduction( q[:,k] )

        q = q_reduced

    #Now q has dimension: number_of_points

    #Create grid and update xll/yll corner and x,y
    if verbose: print 'Creating grid'
    ncols = int((xmax-xmin)/cellsize)+1
    nrows = int((ymax-ymin)/cellsize)+1

    newxllcorner = xmin+xllcorner
    newyllcorner = ymin+yllcorner

    x = x+xllcorner-newxllcorner
    y = y+yllcorner-newyllcorner

    vertex_points = concatenate ((x[:, NewAxis] ,y[:, NewAxis]), axis = 1)
    assert len(vertex_points.shape) == 2


    from Numeric import zeros, Float
    grid_points = zeros ( (ncols*nrows, 2), Float )


    for i in xrange(nrows):
        yg = i*cellsize
        for j in xrange(ncols):
            xg = j*cellsize
            k = i*ncols + j

            grid_points[k,0] = xg
            grid_points[k,1] = yg

    #Interpolate
    from least_squares import Interpolation


    #FIXME: This should be done with precrop = True, otherwise it'll
    #take forever. With expand_search  set to False, some grid points might
    #miss out....

    interp = Interpolation(vertex_points, volumes, grid_points, alpha=0.0,
                           precrop = False, expand_search = True,
                           verbose = verbose)

    #Interpolate using quantity values
    if verbose: print 'Interpolating'
    grid_values = interp.interpolate(q).flat

    #Write
    #Write prj file
    if verbose: print 'Writing %s' %prjfile
    prjid = open(prjfile, 'w')
    prjid.write('Projection    %s\n' %'UTM')
    prjid.write('Zone          %d\n' %zone)
    prjid.write('Datum         %s\n' %datum)
    prjid.write('Zunits        NO\n')
    prjid.write('Units         METERS\n')
    prjid.write('Spheroid      %s\n' %datum)
    prjid.write('Xshift        %d\n' %false_easting)
    prjid.write('Yshift        %d\n' %false_northing)
    prjid.write('Parameters\n')
    prjid.close()



    if verbose: print 'Writing %s' %ascfile
    NODATA_value = -9999

    ascid = open(ascfile, 'w')

    ascid.write('ncols         %d\n' %ncols)
    ascid.write('nrows         %d\n' %nrows)
    ascid.write('xllcorner     %d\n' %newxllcorner)
    ascid.write('yllcorner     %d\n' %newyllcorner)
    ascid.write('cellsize      %f\n' %cellsize)
    ascid.write('NODATA_value  %d\n' %NODATA_value)


    #Get bounding polygon from mesh
    P = interp.mesh.get_boundary_polygon()
    inside_indices = inside_polygon(grid_points, P)

    for i in range(nrows):
        if verbose and i%((nrows+10)/10)==0:
            print 'Doing row %d of %d' %(i, nrows)

        for j in range(ncols):
            index = (nrows-i-1)*ncols+j

            if sometrue(inside_indices == index):
                ascid.write('%f ' %grid_values[index])
            else:
                ascid.write('%d ' %NODATA_value)

        ascid.write('\n')

    #Close
    ascid.close()
    fid.close()

#********************
#*** END OF OBSOLETE FUNCTIONS
#***************
