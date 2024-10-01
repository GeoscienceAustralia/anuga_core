""" Random utilities for reading sww file data and for plotting
(in ipython, or in scripts)

    Functionality of note:

    plot_utils.get_outputs -- read the data from a single sww file
    into a single object
    
    plot_utils.combine_outputs -- read the data from a list of sww
    files into a single object
    
    plot_utils.near_transect -- for finding the indices of points
                          'near' to a given line, and
                          assigning these points a
                          coordinate along that line.

    This is useful for plotting outputs which are 'almost' along a
    transect (e.g. a channel cross-section) -- see example below

    plot_utils.sort_sww_filenames -- match sww filenames by a wildcard, and order
                               them according to their 'time'. This means that
                               they can be stuck together using
                               'combine_outputs' correctly

    plot_utils.triangle_areas -- compute the areas of every triangle
                           in a get_outputs object [ must be vertex-based]

    plot_utils.water_volume -- compute the water volume at every
                         time step in an sww file (needs both
                         vertex and centroid value input). 

    plot_utils.Make_Geotif -- convert sww centroids to a georeferenced tiff
 
    Here is an example ipython session which uses some of these functions:

    > from anuga import plot_utils
    > from matplotlib import pyplot as pyplot
    > p=plot_utils.get_output('myfile.sww',minimum_allowed_height=0.01)
    > p2=plot_utils.get_centroids(p,velocity_extrapolation=True)
    > xxx=plot_utils.near_transect(p,[95., 85.], [120.,68.],tol=2.) # Could equally well use p2
    > pyplot.ion() # Interactive plotting
    > pyplot.scatter(xxx[1],p.vel[140,xxx[0]],color='red') # Plot along the transect

    FIXME: TODO -- Convert to a single function 'get_output', which can either take a
          single filename, a list of filenames, or a wildcard defining a number of
          filenames, and ensure that in each case, the output will be as desired.

"""

from anuga.file.netcdf import NetCDFFile
import numpy
import copy
import matplotlib.cm

class combine_outputs(object):
    """
    Read in a list of filenames, and combine all their outputs into a single object.
    e.g.:

    p = util.combine_outputs(['file1.sww', 'file1_time_10000.sww', 'file1_time_20000.sww'], 0.01)
    
    will make an object p which has components p.x,p.y,p.time,p.stage, .... etc,
    where the values of stage / momentum / velocity from the sww files are concatenated as appropriate.

    This is nice for interactive interrogation of model outputs, or for sticking together outputs in scripts
   
    WARNING: It is easy to use lots of memory, if the sww files are large.

    Note: If you want the centroid values, then you could subsequently use:

    p2 = util.get_centroids(p,velocity_extrapolation=False)

    which would make an object p2 that is like p, but holds information at centroids
    """
    def __init__(self, filename_list, minimum_allowed_height=1.0e-03, verbose=False):
        #
        # Go through the sww files in 'filename_list', and combine them into one object.
        #

        for i, filename in enumerate(filename_list):
            if verbose: print(i, filename)
            # Store output from filename
            p_tmp = get_output(filename, minimum_allowed_height,verbose=verbose)
            if(i==0):
                # Create self
                p1=p_tmp
            else:
                # Append extra data to self
                # Note that p1.x, p1.y, p1.vols, p1.elev should not change
                assert (p1.x == p_tmp.x).all()
                assert (p1.y == p_tmp.y).all()
                assert (p1.vols ==p_tmp.vols).all()
                p1.time = numpy.append(p1.time, p_tmp.time)
                p1.stage = numpy.append(p1.stage, p_tmp.stage, axis=0)
                p1.height = numpy.append(p1.height, p_tmp.height, axis=0)
                p1.xmom = numpy.append(p1.xmom, p_tmp.xmom, axis=0)
                p1.ymom = numpy.append(p1.ymom, p_tmp.ymom, axis=0)
                p1.xvel = numpy.append(p1.xvel, p_tmp.xvel, axis=0)
                p1.yvel = numpy.append(p1.yvel, p_tmp.yvel, axis=0)
                p1.vel = numpy.append(p1.vel, p_tmp.vel, axis=0)
        
        self.x, self.y, self.time, self.vols, self.stage, \
                self.height, self.elev, self.friction, self.xmom, self.ymom, \
                self.xvel, self.yvel, self.vel, self.minimum_allowed_height,\
                self.xllcorner, self.yllcorner, self.timeSlices =\
                p1.x, p1.y, p1.time, p1.vols, p1.stage, \
                p1.height, p1.elev, p1.friction, p1.xmom, p1.ymom, \
                p1.xvel, p1.yvel, p1.vel, p1.minimum_allowed_height,\
                p1.xllcorner, p1.yllcorner, p1.timeSlices 

        self.filename = p1.filename
        self.verbose = p1.verbose


####################

def sort_sww_filenames(sww_wildcard):
    # Function to take a 'wildcard' sww filename, 
    # and return a list of all filenames of this type,
    # sorted by their time.
    # This can then be used efficiently in 'combine_outputs'
    # if you have many filenames starting with the same pattern
    import glob
    filenames=glob.glob(sww_wildcard)
    
    # Extract time from filenames
    file_time=list(range(len(filenames))) # Predefine
     
    for i,filename in enumerate(filenames):
        filesplit=filename.rsplit('_time_')
        if(len(filesplit)>1):
            file_time[i]=int(filesplit[1].split('_0.sww')[0])
        else:
            file_time[i]=0         
    
    name_and_time=list(zip(file_time,filenames))
    name_and_time.sort() # Sort by file_time
    
    output_times, output_names = list(zip(*name_and_time))
    
    return list(output_names)

#####################################################################

# FIXME (Ole): We should move this to e.g. the module sww.py as it has nothing to do with plotting ;-)
class get_output(object):
    """Read in data from an .sww file in a convenient form
       e.g. 
        p = plot_utils.get_output('channel3.sww', minimum_allowed_height=0.01)
        
       p then contains most relevant information as e.g., p.stage, p.elev, p.xmom, etc 
    """
    def __init__(self, filename, minimum_allowed_height=1.0e-03, timeSlices='all', verbose=False):
                # FIXME: verbose is not used
        self.x, self.y, self.time, self.vols, self.stage, \
                self.height, self.elev, self.friction, self.xmom, self.ymom, \
                self.xvel, self.yvel, self.vel, self.minimum_allowed_height,\
                self.xllcorner, self.yllcorner, self.timeSlices, self.starttime = \
                _read_output(filename, minimum_allowed_height, copy.copy(timeSlices))
        self.filename = filename
        self.verbose = verbose

####################################################################
def getInds(varIn, timeSlices, absMax=False):
    """
     Convenience function to get the indices we want in an array.
     There are a number of special cases that make this worthwhile
     having in its own function
    
     INPUT: varIn -- numpy array, either 1D (variables in space) or 2D
            (variables in time+space)
            timeSlices -- times that we want the variable, see read_output or get_output
            absMax -- if TRUE and timeSlices is 'max', then get max-absolute-values
     OUTPUT:
           
    """
    #import pdb
    #pdb.set_trace()

    if (len(varIn.shape)==2):
        # There are multiple time-slices
        if timeSlices == 'max':
            # Extract the maxima over time, assuming there are multiple
            # time-slices, and ensure the var is still a 2D array
            if( not absMax):
                var = (varIn[:]).max(axis=0, keepdims=True)
            else:
                # For variables xmom,ymom,xvel,yvel we want the 'maximum-absolute-value'
                varInds = abs(varIn[:]).argmax(axis=0)
                varNew = varInds * 0.
                for i in range(len(varInds)):
                    varNew[i] = varIn[varInds[i], i]
                var = varNew
                var=var.reshape((1, len(var)))
        else:
            var = numpy.zeros((len(timeSlices), varIn.shape[1]), dtype='float32')
            for i in range(len(timeSlices)):
                var[i,:]=varIn[timeSlices[i]]
            var.reshape((len(timeSlices), varIn.shape[1]))
    else:
        # There is 1 time slice only
        var = varIn[:]
    
    return var

############################################################################

def _read_output(filename, minimum_allowed_height, timeSlices):
    """
     Purpose: To read the sww file, and output a number of variables as arrays that 
              we can then e.g. plot, interrogate 

              See get_output for the typical interface, and get_centroids for
                working with centroids directly
    
     Input: filename -- The name of an .sww file to read data from,
                        e.g. read_sww('channel3.sww')
            minimum_allowed_height -- zero velocity when height < this
            timeSlices -- List of time indices to read (e.g. [100] or [0, 10, 21]), or 'all' or 'last' or 'max'
                          If 'max', the time-max of each variable will be computed. For xmom/ymom/xvel/yvel, the
                           one with maximum magnitude is reported
    
    
     Output: x, y, time, stage, height, elev, xmom, ymom, xvel, yvel, vel
             x,y are only stored at one time
             elevation may be stored at one or multiple times
             everything else is stored every time step for vertices
    """

    # Open ncdf connection
    fid=NetCDFFile(filename)
    
    time=fid.variables['time'][:]

    # Treat specification of timeSlices
    if(timeSlices == 'all'):
        inds = list(range(len(time)))
    elif(timeSlices == 'last'):
        inds = [len(time)-1]
    elif(timeSlices == 'max'):
        inds = 'max' #
    else:
        try:
            inds = list(timeSlices)
        except:
            inds = [timeSlices]
    
    if(inds != 'max'):
        time = time[inds]
    else:
        # We can't really assign a time to 'max', but I guess max(time) is
        # technically the right thing -- if not misleading
        time = time.max()

    
    # Get lower-left
    xllcorner = fid.xllcorner
    yllcorner = fid.yllcorner
    starttime = fid.starttime

    # Read variables
    x = fid.variables['x'][:]
    y = fid.variables['y'][:]

    stage = getInds(fid.variables['stage'], timeSlices=inds)
    elev = getInds(fid.variables['elevation'], timeSlices=inds)

    # Simple approach for volumes
    vols = fid.variables['volumes'][:]

    # Friction if it exists
    if('friction' in fid.variables):
        friction=getInds(fid.variables['friction'], timeSlices=inds) 
    else:
        # Set friction to nan if it is not stored
        friction = elev * 0. + numpy.nan

    # Trick to treat the case where inds == 'max'
    inds2 = copy.copy(inds)
    if inds == 'max':
        inds2 = list(range(len(fid.variables['time'])))
    
    # Get height
    if('height' in fid.variables):
        height = fid.variables['height'][inds2]
    else:
        # Back calculate height if it is not stored
        #height = fid.variables['stage'][inds2]+0.
        height = numpy.zeros((len(inds2), stage.shape[1]), dtype='float32')
        for i in range(len(inds2)):
            height[i,:] = fid.variables['stage'][inds2[i]]

        if(len(elev.shape)==2):
            height = height-elev
        else:
            for i in range(height.shape[0]):
                height[i,:] = height[i,:]-elev
    height = height*(height>0.)

    # Get xmom
    #xmom = fid.variables['xmomentum'][inds2]
    #ymom = fid.variables['ymomentum'][inds2]
    xmom = numpy.zeros((len(inds2), stage.shape[1]), dtype='float32')
    ymom = numpy.zeros((len(inds2), stage.shape[1]), dtype='float32')
    for i in range(len(inds2)):
        xmom[i,:] = fid.variables['xmomentum'][inds2[i]]
        ymom[i,:] = fid.variables['ymomentum'][inds2[i]]
    
    # Get vel
    h_inv = 1.0/(height+1.0e-12)
    hWet = (height > minimum_allowed_height)
    xvel = xmom*h_inv*hWet
    yvel = ymom*h_inv*hWet
    vel = (xmom**2 + ymom**2)**0.5*h_inv*hWet

    if inds == 'max':
        height = height.max(axis=0, keepdims=True)
        vel = vel.max(axis=0, keepdims=True)
        xvel = getInds(xvel, timeSlices=inds,absMax=True)
        yvel = getInds(yvel, timeSlices=inds,absMax=True)
        xmom = getInds(xmom, timeSlices=inds,absMax=True)
        ymom = getInds(ymom, timeSlices=inds,absMax=True)

    fid.close()

    return x, y, time, vols, stage, height, elev, friction, xmom, ymom,\
           xvel, yvel, vel, minimum_allowed_height, xllcorner,yllcorner, inds, starttime

######################################################################################

class get_centroids(object):
    """
    Extract centroid values from the output of get_output, OR from a
        filename  
    See _read_output or _get_centroid_values for further explanation of
        arguments
    e.g.
        # Case 1 -- get vertex values first, then centroids
        p = plot_utils.get_output('my_sww.sww', minimum_allowed_height=0.01) 
        pc=util.get_centroids(p, velocity_extrapolation=True) 

        # Case 2 -- get centroids directly
        pc=plot_utils.get_centroids('my_sww.sww', velocity_extrapolation=True) 

    NOTE: elevation is only stored once in the output, even if it was
          stored every timestep.
          Lots of existing plotting code assumes elevation is a 1D
          array. 
          But as a hack for the time being the elevation from the file 
          is available via elev_orig
    """
    def __init__(self, p, velocity_extrapolation=False, verbose=False,
                 timeSlices=None, minimum_allowed_height=1.0e-03):
        
        self.time, self.x, self.y, self.stage, self.xmom,\
            self.ymom, self.height, self.elev, self.elev_orig, self.friction, self.xvel,\
            self.yvel, self.vel, self.xllcorner, self.yllcorner, self.timeSlices= \
                _get_centroid_values(p, velocity_extrapolation,\
                                     timeSlices=copy.copy(timeSlices),\
                                     minimum_allowed_height=minimum_allowed_height,\
                                     verbose=verbose)

def _getCentVar(fid, varkey_c, time_indices, absMax=False,  vols = None, space_indices=None):
    """
        Convenience function used to get centroid variables from netCDF
        file connection fid

    """

    if vols is not None:
        vols0 = vols[:,0]
        vols1 = vols[:,1]
        vols2 = vols[:,2]

    if((varkey_c in fid.variables)==False):
        # It looks like centroid values are not stored
        # In this case, compute centroid values from vertex values
        assert (vols is not None), "Must specify vols since centroid quantity is not stored"

        newkey=varkey_c.replace('_c','')
        if time_indices != 'max':
            # Relatively efficient treatment is possible
            var_cent = fid.variables[newkey]
            if (len(var_cent.shape)>1):
                # array contain time slices
                var_cent = numpy.zeros((len(time_indices), fid.variables[newkey].shape[1]), dtype='float32')
                for i in range(len(time_indices)):
                    var_cent[i,:] = fid.variables[newkey][time_indices[i]]
                var_cent = (var_cent[:,vols0]+var_cent[:,vols1]+var_cent[:,vols2])/3.0
            else:
                var_cent = fid.variables[newkey][:]
                var_cent = (var_cent[vols0]+var_cent[vols1]+var_cent[vols2])/3.0
        else:
            # Requires reading all the data
            tmp = fid.variables[newkey][:]
            try: # array contain time slices
                tmp=(tmp[:,vols0]+tmp[:,vols1]+tmp[:,vols2])/3.0
            except:
                tmp=(tmp[vols0]+tmp[vols1]+tmp[vols2])/3.0
            var_cent=getInds(tmp, timeSlices=time_indices, absMax=absMax)
    else:
        if time_indices != 'max':
            if(len(fid.variables[varkey_c].shape)>1):
                var_cent = numpy.zeros((len(time_indices), fid.variables[varkey_c].shape[1]), dtype='float32')
                for i in range(len(time_indices)):
                    var_cent[i,:] = fid.variables[varkey_c][time_indices[i]]
            else:
                var_cent = fid.variables[varkey_c][:]
        else:
            var_cent=getInds(fid.variables[varkey_c][:], timeSlices=time_indices, absMax=absMax)

    if space_indices is not None:
        # Maybe only return particular space indices. Could do this more
        # efficiently by only reading those indices initially, if that proves
        # important
        if (len(var_cent.shape)>1):
            var_cent = var_cent[:,space_indices]
        else:
            var_cent = var_cent[space_indices]

    return var_cent

                                 
def _get_centroid_values(p, velocity_extrapolation, verbose, timeSlices, 
                         minimum_allowed_height):
    """
    Function to get centroid information -- main interface is through 
        get_centroids. 
        See get_centroids for usage examples, and read_output or get_output for further relevant info
     Input: 
           p --  EITHER:
                  The result of e.g. p=util.get_output('mysww.sww'). 
                  See the get_output class defined above. 
                 OR:
                  Alternatively, the name of an sww file
    
           velocity_extrapolation -- If true, and centroid values are not
            in the file, then compute centroid velocities from vertex velocities, and
            centroid momenta from centroid velocities. If false, and centroid values
            are not in the file, then compute centroid momenta from vertex momenta,
            and centroid velocities from centroid momenta
    
           timeSlices = list of integer indices when we want output for, or
                        'all' or 'last' or 'max'. See _read_output
    
           minimum_allowed_height = height at which velocities are zeroed. See _read_output
    
     Output: Values of x, y, Stage, xmom, ymom, elev, xvel, yvel, vel etc at centroids
    """

    # Figure out if p is a string (filename) or the output of get_output
    pIsFile = isinstance(p, str)
    if(pIsFile): 
        fid = NetCDFFile(p) 
    else:
        fid = NetCDFFile(p.filename)

    # UPDATE: 15/06/2014 -- below, we now get all variables directly from the file
    #         This is more flexible, and allows to get 'max' as well
    #         However, potentially it could have performance penalities vs the old approach (?)

    # Make 3 arrays, each containing one index of a vertex of every triangle.
    vols = fid.variables['volumes'][:]
    vols0 = vols[:,0]
    vols1 = vols[:,1]
    vols2 = vols[:,2]
    
    # Get lower-left offset
    xllcorner = fid.xllcorner
    yllcorner = fid.yllcorner
   
    #@ Get timeSlices 
    # It will be either a list of integers, or 'max'
    l = len(vols)
    time = fid.variables['time'][:]
    nts = len(time) # number of time slices in the file 
    if(timeSlices is None):
        if(pIsFile):
            # Assume all timeSlices
            timeSlices=list(range(nts))
        else:
            timeSlices=copy.copy(p.timeSlices)
    else:
        # Treat word-based special cases
        if(timeSlices == 'all'):
            timeSlices=list(range(nts))
        if(timeSlices == 'last'):
            timeSlices=[nts-1]

    #@ Get minimum_allowed_height
    if(minimum_allowed_height is None):
        if(pIsFile):
            minimum_allowed_height=0.
        else:
            minimum_allowed_height=copy.copy(p.minimum_allowed_height)

    # Treat specification of timeSlices
    if(timeSlices == 'all'):
        inds = list(range(len(time)))
    elif(timeSlices=='last'):
        inds = [len(time)-1]
    elif(timeSlices=='max'):
        inds = 'max' #
    else:
        try:
            inds = list(timeSlices)
        except:
            inds = [timeSlices]
    
    if(inds != 'max'):
        time = time[inds]
    else:
        # We can't really assign a time to 'max', but I guess max(time) is
        # technically the right thing -- if not misleading
        time = time.max()

    # Get coordinates
    x = fid.variables['x'][:]
    y = fid.variables['y'][:]
    x_cent = (x[vols0] + x[vols1] + x[vols2]) / 3.0
    y_cent = (y[vols0] + y[vols1] + y[vols2]) / 3.0

    # Stage and height and elevation
    stage_cent = _getCentVar(fid, 'stage_c', time_indices=inds, vols=vols)
    elev_cent = _getCentVar(fid, 'elevation_c', time_indices=inds, vols=vols)

    # Hack to allow refernece to time varying elevation
    elev_cent_orig = elev_cent
    
    if(len(elev_cent.shape) == 2):
        # Coerce to 1D array, since lots of our code assumes it is
        elev_cent = elev_cent[0,:]

    # Friction might not be stored at all
    try:
        friction_cent = _getCentVar(fid, 'friction_c', time_indices=inds, vols=vols)
    except:
        friction_cent = elev_cent*0.+numpy.nan
    
    # Trick to treat the case where inds == 'max'
    inds2 = copy.copy(inds)
    if inds == 'max':
        inds2 = list(range(len(fid.variables['time'])))
   
    # height
    height_cent = stage_cent + 0.
    for i in range(stage_cent.shape[0]):
        height_cent[i,:] = stage_cent[i,:] - elev_cent

    if 'xmomentum_c' in fid.variables:
        # The following commented out lines seem to only work on
        # some numpy/netcdf versions. So we loop
        #xmom_cent = fid.variables['xmomentum_c'][inds2]
        #ymom_cent = fid.variables['ymomentum_c'][inds2]
        xmom_cent = numpy.zeros((len(inds2), fid.variables['xmomentum_c'].shape[1]), dtype='float32')
        ymom_cent = numpy.zeros((len(inds2), fid.variables['ymomentum_c'].shape[1]), dtype='float32')
        height_c_tmp = numpy.zeros((len(inds2), fid.variables['stage_c'].shape[1]), dtype='float32')
        for i in range(len(inds2)):
            xmom_cent[i,:] = fid.variables['xmomentum_c'][inds2[i]]
            ymom_cent[i,:] = fid.variables['ymomentum_c'][inds2[i]]
            if 'height_c' in fid.variables:
                height_c_tmp[i,:] = fid.variables['height_c'][inds2[i]]
            else:
                height_c_tmp[i,:] = fid.variables['stage_c'][inds2[i]] - elev_cent

        # Vel
        hInv = 1.0/(height_c_tmp + 1.0e-12)
        hWet = (height_c_tmp > minimum_allowed_height)
        xvel_cent = xmom_cent*hInv*hWet
        yvel_cent = ymom_cent*hInv*hWet

    else:
        # Get important vertex variables
        xmom_v = numpy.zeros((len(inds2), fid.variables['xmomentum'].shape[1]), dtype='float32')
        ymom_v = numpy.zeros((len(inds2), fid.variables['ymomentum'].shape[1]), dtype='float32')
        stage_v = numpy.zeros((len(inds2), fid.variables['stage'].shape[1]), dtype='float32')
        for i in range(len(inds2)):
            xmom_v[i,:] = fid.variables['xmomentum'][inds2[i]]
            ymom_v[i,:] = fid.variables['ymomentum'][inds2[i]]
            stage_v[i,:] = fid.variables['stage'][inds2[i]]

        elev_v = fid.variables['elevation']
        # Fix elevation + get height at vertices
        if (len(elev_v.shape)>1):
            elev_v = numpy.zeros(elev_v.shape, dtype='float32')
            for i in range(elev_v.shape[0]):
                elev_v[i,:] = fid.variables['elevation'][inds2[i]]
            height_v = stage_v - elev_v
        else:
            elev_v = elev_v[:]
            height_v = stage_v + 0.
            for i in range(stage_v.shape[0]):
                height_v[i,:] = stage_v[i,:] - elev_v

        # Height at centroids        
        height_c_tmp = (height_v[:, vols0] + height_v[:,vols1] + height_v[:,vols2])/3.0
       
        # Compute xmom/xvel/ymom/yvel
        if velocity_extrapolation:

            xvel_v = xmom_v*0.
            yvel_v = ymom_v*0.

            hInv = 1.0/(height_v+1.0e-12)
            hWet = (height_v > minimum_allowed_height)

            xvel_v = xmom_v*hInv*hWet
            yvel_v = ymom_v*hInv*hWet

            # Final xmom/ymom centroid values
            xvel_cent = (xvel_v[:, vols0] + xvel_v[:,vols1] + xvel_v[:,vols2])/3.0
            xmom_cent = xvel_cent*height_c_tmp
            yvel_cent = (yvel_v[:, vols0] + yvel_v[:,vols1] + yvel_v[:,vols2])/3.0
            ymom_cent = yvel_cent*height_c_tmp

        else:
            hInv = 1.0/(height_c_tmp + 1.0e-12)
            hWet = (height_c_tmp > minimum_allowed_height)

            xmom_v = numpy.zeros((len(inds2), fid.variables['xmomentum'].shape[1]), dtype='float32')
            ymom_v = numpy.zeros((len(inds2), fid.variables['ymomentum'].shape[1]), dtype='float32')
            for i in range(len(inds2)):
                xmom_v[i,:] = fid.variables['xmomentum'][inds2[i]]
                ymom_v[i,:] = fid.variables['ymomentum'][inds2[i]]

            xmom_cent = (xmom_v[:,vols0] + xmom_v[:,vols1] + xmom_v[:,vols2])/3.0
            xvel_cent = xmom_cent*hInv*hWet
            ymom_cent = (ymom_v[:,vols0] + ymom_v[:,vols1] + ymom_v[:,vols2])/3.0
            yvel_cent = ymom_cent*hInv*hWet

    # Velocity
    vel_cent = (xvel_cent**2 + yvel_cent**2)**0.5

    if inds == 'max':
        vel_cent = vel_cent.max(axis=0, keepdims=True)
        #vel_cent = getInds(vel_cent, timeSlices=inds)
        xmom_cent = getInds(xmom_cent, timeSlices=inds, absMax=True)
        ymom_cent = getInds(ymom_cent, timeSlices=inds, absMax=True)
        xvel_cent = getInds(xvel_cent, timeSlices=inds, absMax=True)
        yvel_cent = getInds(yvel_cent, timeSlices=inds, absMax=True)

    fid.close()
    
    return time, x_cent, y_cent, stage_cent, xmom_cent,\
             ymom_cent, height_cent, elev_cent, elev_cent_orig, friction_cent,\
             xvel_cent, yvel_cent, vel_cent, xllcorner, yllcorner, inds


def animate_1D(time, var, x, ylab=' '): 
    """Animate a 2d array with a sequence of 1d plots

     Input: time = one-dimensional time vector;
            var =  array with first dimension = len(time) ;
            x = (optional) vector width dimension equal to var.shape[1];
            ylab = ylabel for plot
    """
    
    import pylab
    import numpy
   
    

    pylab.close()
    pylab.ion()

    # Initial plot
    vmin = var.min()
    vmax = var.max()
    line, = pylab.plot( (x.min(), x.max()), (vmin, vmax), 'o')

    # Lots of plots
    for i in range(len(time)):
        line.set_xdata(x)
        line.set_ydata(var[i,:])
        pylab.draw()
        pylab.xlabel('x')
        pylab.ylabel(ylab)
        pylab.title('time = ' + str(time[i]))
    
    return

def near_transect(p, point1, point2, tol=1.):
    # Function to get the indices of points in p less than 'tol' from the line
    # joining (x1,y1), and (x2,y2)
    # p comes from util.get_output('mysww.sww')
    #
    # e.g.
    # import util
    # from matplotlib import pyplot
    # p=util.get_output('merewether_1m.sww',0.01)
    # p2=util.get_centroids(p,velocity_extrapolation=True)
    # #xxx=transect_interpolate.near_transect(p,[95., 85.], [120.,68.],tol=2.)
    # xxx=util.near_transect(p,[95., 85.], [120.,68.],tol=2.)
    # pyplot.scatter(xxx[1],p.vel[140,xxx[0]],color='red')
    
    x1=point1[0]
    y1=point1[1]
    
    x2=point2[0]
    y2=point2[1]
    
    # Find line equation a*x + b*y + c = 0
    # based on y=gradient*x +intercept
    if x1!=x2:
        gradient= ((y2-y1) / (x2-x1))
        intercept = y1 - gradient*x1
        #
        a = -gradient
        b = 1.
        c = -intercept
    else:
        a=1.
        b=0.
        c=-x2 
    
    # Distance formula
    inv_denom = 1./(a**2 + b**2)**0.5
    distp = abs(p.x*a + p.y*b + c)*inv_denom
    
    near_points = (distp<tol).nonzero()[0]
    
    # Now find a 'local' coordinate for the point, projected onto the line
    # g1 = unit vector parallel to the line
    # g2 = vector joining (x1,y1) and (p.x,p.y)
    g1x = x2-x1 
    g1y = y2-y1
    g1_norm = (g1x**2 + g1y**2)**0.5
    g1x = g1x / g1_norm
    g1y = g1y / g1_norm
    
    g2x = p.x[near_points] - x1
    g2y = p.y[near_points] - y1
    
    # Dot product = projected distance == a local coordinate
    local_coord = g1x*g2x + g1y*g2y
    
    # only keep coordinates between zero and the distance along the line
    dl=((x1-x2)**2+(y1-y2)**2)**0.5
    keepers=(local_coord<=dl)*(local_coord>=0.)
    keepers=keepers.nonzero()
    
    return near_points[keepers], local_coord[keepers]


def triangle_areas(p, subset=None):
    # Compute areas of triangles in p -- assumes p contains vertex information
    # subset = vector of centroid indices to include in the computation. 

    if(subset is None):
        subset=list(range(len(p.vols[:,0])))
    
    x0=p.x[p.vols[subset,0]]
    x1=p.x[p.vols[subset,1]]
    x2=p.x[p.vols[subset,2]]
    
    y0=p.y[p.vols[subset,0]]
    y1=p.y[p.vols[subset,1]]
    y2=p.y[p.vols[subset,2]]
    
    # Vectors for cross-product
    v1_x=x0-x1
    v1_y=y0-y1
    #
    v2_x=x2-x1
    v2_y=y2-y1
    # Area
    area=(v1_x*v2_y-v1_y*v2_x)*0.5
    area=abs(area)
    return area


def water_volume(p, p2, per_unit_area=False, subset=None):
    # Compute the water volume from p(vertex values) and p2(centroid values)

    if(subset is None):
        subset=list(range(len(p2.x)))

    l=len(p2.time)
    area=triangle_areas(p, subset=subset)
    
    total_area=area.sum()
    volume=p2.time*0.
   
    # This accounts for how volume is measured in ANUGA 
    # Compute in 2 steps to reduce precision error from limited SWW precision
    # FIXME: Is this really needed?
    for i in range(l):
        #volume[i]=((p2.stage[i,subset]-p2.elev[subset])*(p2.stage[i,subset]>p2.elev[subset])*area).sum()
        volume[i]=((p2.stage[i,subset])*(p2.stage[i,subset]>p2.elev[subset])*area).sum()
        volume[i]=volume[i]+((-p2.elev[subset])*(p2.stage[i,subset]>p2.elev[subset])*area).sum()
    
    if(per_unit_area):
        volume = volume / total_area
    
    return volume


def get_triangle_containing_point(p, point, search_order=None):
    """
    Function to get the index of a triangle containing a point. 
    It loops over all points in the mesh until it finds on that contains the point.
    The search order (i.e. order in which triangles defined by p.vols are searched) can
    be provided. If it is not, it is estimated by computing the distance
    from the point to the first vertex of every triangle, and searching from smallest to largest.

    @param p Object containing mesh vertex information (e.g. from plot_utils.get_output)
    @param point A single point
    @param search_order An optional integer array giving the order in which to search the mesh triangles

    @return The index such that the triangle defined by p.vols[index,:] contains the point
    """

    V = p.vols

    x = p.x
    y = p.y

    from anuga.geometry.polygon import is_inside_polygon

    if search_order is None:
        # Estimate a good search order by finding the distance to the first
        # vertex of every triangle, and doing the search ordered by that
        # distance.
        point_distance2 = (x[V[:,0]] - point[0])**2 + (y[V[:,0]]-point[1])**2
        point_distance_order = point_distance2.argsort().tolist()
    else:
        point_distance_order = search_order

    for i in point_distance_order:
        i0 = V[i,0]
        i1 = V[i,1]
        i2 = V[i,2]
        poly = [ [x[i0], y[i0]], [x[i1], y[i1]], [x[i2], y[i2]] ]

        if is_inside_polygon(point, poly, closed=True):
            return i

    msg = 'Point %s not found within a triangle' %str(point)
    raise Exception(msg)

def get_triangle_near_point(p, point, tolerance=1.0e20):
    """
    Function to get the index of a triangle nearest to a point (as measured by distance to 
    centroid of the triangle).

    @param p Object containing mesh vertex information (e.g. from plot_utils.get_output)
    @param point A single point (absolute units)
    @param tolerance Raise an exception if "nearest" point is further that tolerance from the domain
    
    @return The index of the triangle "nearest" to point.
    """

    import numpy

    pc = get_centroids(p)
    
    xll_corner = p.xllcorner
    yll_corner = p.yllcorner

    X = pc.x[:] + xll_corner
    Y = pc.y[:] + yll_corner
    
    distance2 = (X - point[0])**2 + (Y - point[1])**2
    
    tid = numpy.argmin(distance2)

    if distance2[tid] > tolerance**2:
        msg = 'Point %s further than %g (m) from domain' % (str(point), tolerance)
        raise Exception(msg)
    else:
        return tid



def get_triangle_lookup_function(pv):
    """Return a function F(x,y) which gives the row index in pv.vols
    corresponding to  the triangle containing x,y. This function
    should be more efficient than get_triangle_containing_point
    if many points need to be looked-up

    @param pv object containing vertex information (e.g. from plot_utils.get_output)
    @return function F(x,y) which gives the index (or indices) in pv.vols
    corresponding to the triangle(s) containing x,y, where x,y can be numpy.arrays

    """
    import matplotlib.tri as tri

    # Get unique vertices for triangle hunting
    complex_verts, unique_inds = numpy.unique(pv.x + 1j*pv.y, return_inverse=True)

    reduced_x = numpy.real(complex_verts)
    reduced_y = numpy.imag(complex_verts)

    # Here the ordering is the same as pv.vols
    reduced_triangles = unique_inds[pv.vols]

    new_triangulation = tri.Triangulation(reduced_x, reduced_y, reduced_triangles)
    tri_lookup = new_triangulation.get_trifinder()

    return(tri_lookup)


def get_extent(p):

    import numpy

    x_min = numpy.min(p.x)
    x_max = numpy.max(p.x)

    y_min = numpy.min(p.y)
    y_max = numpy.max(p.y)

    return x_min, x_max, y_min, y_max



def make_grid(data, lats, lons, fileName, EPSG_CODE=None, proj4string=None, 
               creation_options=[]):
    """
        Convert data,lats,lons to a georeferenced raster tif
        INPUT: data -- array with desired raster cell values
               lats -- 1d array with 'latitude' or 'y' range
               lons -- 1D array with 'longitude' or 'x' range
               fileName -- name of file to write to
               EPSG_CODE -- Integer code with projection information in EPSG format 
               proj4string -- proj4string with projection information
               creation_options -- list of tif creation options for gdal (e.g. ["COMPRESS=DEFLATE"])

        NOTE: proj4string is used in preference to EPSG_CODE if available
    """

    try:
        import osgeo.gdal as gdal
        import osgeo.osr as osr
    except ImportError as e:
        msg='Failed to import gdal/ogr modules --'\
        + 'perhaps gdal python interface is not installed.'
        raise ImportError(msg)
    


    xres = lons[1] - lons[0]
    yres = lats[1] - lats[0]

    ysize = len(lats)
    xsize = len(lons)

    # Assume data/lats/longs refer to cell centres, and compute upper left coordinate
    ulx = lons[0] - (xres / 2.)
    uly = lats[lats.shape[0]-1] + (yres / 2.)

    # GDAL magic to make the tif
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(fileName, xsize, ysize, 1, gdal.GDT_Float32, 
                       creation_options)

    srs = osr.SpatialReference()
    if(proj4string is not None):
        srs.ImportFromProj4(proj4string)
    elif(EPSG_CODE is not None):
        srs.ImportFromEPSG(EPSG_CODE)
    else:
        raise Exception('No spatial reference information given')


    ds.SetProjection(srs.ExportToWkt())

    gt = [ulx, xres, 0, uly, 0, -yres ]
    #gt = [llx, xres, 0, lly, yres,0 ]
    ds.SetGeoTransform(gt)

    #import pdb
    #pdb.set_trace()
    import scipy

    outband = ds.GetRasterBand(1)
    outband.SetNoDataValue(numpy.nan)
    outband.WriteArray(data)

    ds = None
    return

##################################################################################

def Make_Geotif(swwFile=None, 
             output_quantities=['depth'],
             myTimeStep=0, CellSize=100.0, 
             lower_left=None, upper_right=None,
             EPSG_CODE=None, 
             proj4string=None,
             velocity_extrapolation=True,
             min_allowed_height=1.0e-05,
             output_dir='TIFS',
             bounding_polygon=None,
             internal_holes=None,
             verbose=False,
             k_nearest_neighbours=3,
             creation_options=[]):
    """
        Make a georeferenced tif by nearest-neighbour interpolation of sww file outputs (or a 3-column array with xyz Points)

        You must supply projection information as either a proj4string or an integer EPSG_CODE (but not both!)

        INPUTS: swwFile -- name of sww file, OR a 3-column array with x/y/z
                    points. In the latter case x and y are assumed to be in georeferenced
                    coordinates.  The output raster will contain 'z', and will have a name-tag
                    based on the name in 'output_quantities'.
                output_quantities -- list of quantitiies to plot, e.g.
                                ['depth', 'velocity', 'stage','elevation','depthIntegratedVelocity','friction']
                myTimeStep -- list containing time-index of swwFile to plot (e.g. [0, 10, 32] ) or 'last', or 'max', or 'all'
                CellSize -- approximate pixel size for output raster [adapted to fit lower_left / upper_right]
                lower_left -- [x0,y0] of lower left corner. If None, use extent of swwFile.
                upper_right -- [x1,y1] of upper right corner. If None, use extent of swwFile.
                EPSG_CODE -- Projection information as an integer EPSG code (e.g. 3123 for PRS92 Zone 3, 32756 for UTM Zone 56 S, etc). 
                             Google for info on EPSG Codes
                proj4string -- Projection information as a proj4string (e.g. '+init=epsg:3123')
                             Google for info on proj4strings. 
                velocity_extrapolation -- Compute velocity assuming the code extrapolates with velocity (instead of momentum)?
                min_allowed_height -- Minimum allowed height from ANUGA
                output_dir -- Write outputs to this directory
                bounding_polygon -- polygon (e.g. from read_polygon) If present, only set values of raster cells inside the bounding_polygon
                internal_holes -- a list of polygons. If present, do not set values of raster cells inside these polygons.
                k_nearest_neighbours -- how many neighbours to use in interpolation. If k>1, inverse-distance-weighted interpolation is used
                creation_options -- list of tif creation options for gdal, e.g. ['COMPRESS=DEFLATE']
    """

    import scipy.io
    import scipy.interpolate
    import scipy.spatial
    import anuga
    import os
    
    try:
        import osgeo.gdal as gdal
        import osgeo.osr as osr
    except ImportError as e:
        msg = 'Failed to import gdal/ogr modules --'\
        + 'perhaps gdal python interface is not installed.'
        raise ImportError(msg)

    # Check whether swwFile is an array, and if so, redefine various inputs to
    # make the code work
    if(type(swwFile) == numpy.ndarray):
        import copy
        xyzPoints = copy.copy(swwFile)
        swwFile = None

    if(((EPSG_CODE is None) & (proj4string is None) )|
       ((EPSG_CODE is not None) & (proj4string is not None))):
        raise Exception('Must specify EITHER an integer EPSG_CODE describing the file projection, OR a proj4string')


    # Make output_dir
    try:
        os.mkdir(output_dir)
    except:
        pass

    if(swwFile is not None):
        # Read in ANUGA outputs
            
        if(verbose):
            print('Reading sww File ...')
        p2 = get_centroids(swwFile, velocity_extrapolation, timeSlices=myTimeStep,
            minimum_allowed_height=min_allowed_height)
        xllcorner = p2.xllcorner
        yllcorner = p2.yllcorner

        myTimeStep_Orig = myTimeStep
        # Now, myTimeStep just holds indices we want to plot in p2
        if(myTimeStep != 'max'):
            myTimeStep = list(range(len(p2.time)))

        # Ensure myTimeStep is a list
        if type(myTimeStep) != list:
            myTimeStep = [myTimeStep]

        if(verbose):
            print('Extracting required data ...')
        # Get ANUGA points
        swwX = p2.x + xllcorner
        swwY = p2.y + yllcorner
    else:
        # Get the point data from the 3-column array
        if(xyzPoints.shape[1] != 3):
            raise Exception('If an array is passed, it must have exactly 3 columns')
        if(len(output_quantities) != 1):
            raise Exception('Can only have 1 output quantity when passing an array')
        swwX = xyzPoints[:,0]
        swwY = xyzPoints[:,1]
        myTimeStep = ['pointData']

    # Grid for meshing
    if(verbose):
        print('Computing grid of output locations...')
    # Get points where we want raster cells
    if(lower_left is None):
        lower_left = [swwX.min(), swwY.min()]
    if(upper_right is None):
        upper_right = [swwX.max(), swwY.max()]
    nx = int(round((upper_right[0]-lower_left[0])*1.0/(1.0*CellSize)) + 1)
    xres = (upper_right[0]-lower_left[0])*1.0/(1.0*(nx-1))
    desiredX = numpy.linspace(lower_left[0], upper_right[0],nx )
    ny = int(round((upper_right[1]-lower_left[1])*1.0/(1.0*CellSize)) + 1)
    yres = (upper_right[1]-lower_left[1])*1.0/(1.0*(ny-1))
    desiredY = numpy.linspace(lower_left[1], upper_right[1], ny)

    gridX, gridY = numpy.meshgrid(desiredX, desiredY)

    if(verbose):
        print('Making interpolation functions...')
    swwXY = numpy.array([swwX[:],swwY[:]]).transpose()

    # Get function to interpolate quantity onto gridXY_array
    gridXY_array = numpy.array([numpy.concatenate(gridX),
        numpy.concatenate(gridY)]).transpose()
    gridXY_array = numpy.ascontiguousarray(gridXY_array)

    # Create Interpolation function
    #basic_nearest_neighbour=False
    if(k_nearest_neighbours == 1):
        index_qFun = scipy.interpolate.NearestNDInterpolator(
            swwXY,
            numpy.arange(len(swwX),dtype='int64').transpose())
        gridqInd = index_qFun(gridXY_array).astype(int)

        #from pprint import pprint
        #print(72*"=")
        #pprint(gridqInd)
        
        # Function to do the interpolation
        def myInterpFun(quantity):
            return quantity[gridqInd]
    else:
        # Combined nearest neighbours and inverse-distance interpolation
        index_qFun = scipy.spatial.cKDTree(swwXY)
        NNInfo = index_qFun.query(gridXY_array, k=k_nearest_neighbours)
        # Weights for interpolation
        nn_wts = 1./(NNInfo[0]+1.0e-100)
        nn_inds = NNInfo[1]
        def myInterpFun(quantity):
            denom = 0.
            num = 0.
            for i in range(k_nearest_neighbours):
                denom += nn_wts[:,i]
                num += quantity[nn_inds[:,i]]*nn_wts[:,i]
            return num / denom

    if bounding_polygon is not None:
        # Find points to exclude (i.e. outside the bounding polygon)
        from anuga.geometry.polygon import outside_polygon
        cut_points = outside_polygon(gridXY_array, bounding_polygon)
        
    hole_points_list = []
    if internal_holes is not None:
        # Find points to exclude (i.e. inside the internal_holes)
        from anuga.geometry.polygon import inside_polygon
        for hole in internal_holes:
            cut_holes = inside_polygon(gridXY_array, hole)
            hole_points_list.append(cut_holes)

    # Loop over all output quantities and produce the output
    for myTSindex, myTSi in enumerate(myTimeStep):
        if(verbose):
            print('Reduction = ', myTSi)
        for output_quantity in output_quantities:
            if (verbose): print(f'output_quantity {output_quantity}')

            if(myTSi != 'max'):
                myTS = myTSi
            else:
                # We have already extracted the max, and e.g.
                # p2.stage is an array of dimension (1, number_of_pointS).
                myTS = 0

            if(type(myTS) == int):
                if(output_quantity == 'stage'):
                    gridq = myInterpFun(p2.stage[myTS,:])
                if(output_quantity == 'depth'):
                    gridq = p2.height[myTS,:]*(p2.height[myTS,:]>0.)# Force positive depth (tsunami alg)
                    gridq = myInterpFun(gridq)
                if(output_quantity == 'velocity'):
                    gridq = myInterpFun(p2.vel[myTS,:])
                if(output_quantity == 'friction'):
                    gridq = myInterpFun(p2.friction)
                if(output_quantity == 'depthIntegratedVelocity'):
                    swwDIVel = (p2.xmom[myTS,:]**2+p2.ymom[myTS,:]**2)**0.5
                    gridq = myInterpFun(swwDIVel)
                if(output_quantity == 'elevation'):
                    gridq = myInterpFun(p2.elev)
    
                if(myTSi == 'max'):
                    timestepString = 'max'
                else:
                    timestepString = str(myTimeStep[myTSindex])+'_Time_'+str(round(p2.time[myTS]))
            elif(myTS == 'pointData'):
                gridq = myInterpFun(xyzPoints[:,2])

            if ( (bounding_polygon is not None) and (len(cut_points)>0)):
                # Cut the points outside the bounding polygon
                gridq[cut_points] = numpy.nan

            if (internal_holes is not None) and (len(hole_points_list[0]) > 0):
                # Cut the points inside the hole polygons
                for hole_points in hole_points_list:
                    gridq[hole_points] = numpy.nan

            # Make name for output file
            if(myTS != 'pointData'):
                output_name = output_dir + '/' +\
                    os.path.splitext(os.path.basename(swwFile))[0] + '_' +\
                    output_quantity + '_' + timestepString + '.tif'
                            #'_'+str(myTS)+'.tif'
            else:
                output_name = output_dir+'/'+'PointData_'+output_quantity+'.tif'

            if(verbose):
                print('Making raster ...')

            gridq.shape = (len(desiredY),len(desiredX))
            make_grid(numpy.flipud(gridq), desiredY, desiredX, output_name, EPSG_CODE=EPSG_CODE, 
                      proj4string=proj4string, creation_options=creation_options)

    return

def plot_triangles(p, adjustLowerLeft=False, values=None, values_cmap=matplotlib.cm.jet, edgecolors='k'):
    """ Add mesh triangles to a pyplot plot
        
       @param p = object holding sww vertex information (from util.get_output)
       @param adjustLowerLeft = if TRUE, use spatial coordinates, otherwise use ANUGA internal coordinates     
       @param values = list or array of length(p.vols), or None. All triangles are assigned this value (for face plotting colors).
       @param values_cmap = colormap for faces [e.g. values_cmap = matplotlib.cm.get_cmap('spectral')]
       @param edgecolors = edge color for polygons (using matplotlib.colors notation). Use 'none' for no color
    """
    import matplotlib
    from matplotlib import pyplot as pyplot
    from matplotlib.collections import PolyCollection

    x0=p.xllcorner
    y0=p.yllcorner 

    # Make vertices for PolyCollection Object
    vertices = []
    for i in range(len(p.vols)):
        k1=p.vols[i][0]
        k2=p.vols[i][1]
        k3=p.vols[i][2]

        tri_coords = numpy.array([ [p.x[k1], p.y[k1]], [p.x[k2], p.y[k2]], [p.x[k3], p.y[k3]] ])
        if adjustLowerLeft:
            tri_coords[:,0] = tri_coords[:,0] + x0
            tri_coords[:,1] = tri_coords[:,1] + y0

        vertices.append(tri_coords)
     
    # Make PolyCollection 
    if values is None: 
        all_poly = PolyCollection( vertices, array = numpy.zeros(len(vertices)), 
            edgecolors=edgecolors)
        all_poly.set_facecolor('none')
    else:
        try:
            lv = len(values)
        except:
            values = numpy.array(len(p.vols)*[values])
            lv = len(values)

        msg = 'len(values) must be the same as len(p.vols) (or values can be a constant)'
        assert lv==len(p.vols), msg
        all_poly = PolyCollection( vertices, array = values, cmap = values_cmap, 
            edgecolors=edgecolors)

    # Add to plot
    # FIXME: To see the triangles, this might require that the user does
    # something else to the plot?
    pyplot.gca().add_collection(all_poly)

def find_neighbours(p,ind):
    """ 
        Find the triangles neighbouring triangle 'ind'
        p is an object from get_output containing mesh vertices
    """
    ind_nei=p.vols[ind]
    
    shared_nei0=p.vols[:,1]*0.0
    shared_nei1=p.vols[:,1]*0.0
    shared_nei2=p.vols[:,1]*0.0
    # Compute indices that match one of the vertices of triangle ind
    # Note: Each triangle can only match a vertex, at most, once
    for i in range(3):
        shared_nei0+=1*(p.x[p.vols[:,i]]==p.x[ind_nei[0]])*\
            1*(p.y[p.vols[:,i]]==p.y[ind_nei[0]])
        
        shared_nei1+=1*(p.x[p.vols[:,i]]==p.x[ind_nei[1]])*\
            1*(p.y[p.vols[:,i]]==p.y[ind_nei[1]])
        
        shared_nei2+=1*(p.x[p.vols[:,i]]==p.x[ind_nei[2]])*\
            1*(p.y[p.vols[:,i]]==p.y[ind_nei[2]])
    
    out=(shared_nei2 + shared_nei1 + shared_nei0)
    return((out==2).nonzero())

def calc_edge_elevations(p):
    """
        Compute the triangle edge elevations on p
        Return x,y,elev for edges
    """
    pe_x=p.x*0.
    pe_y=p.y*0.
    pe_el=p.elev*0.

   
    # Compute coordinates + elevations 
    pe_x[p.vols[:,0]] = 0.5*(p.x[p.vols[:,1]] + p.x[p.vols[:,2]])
    pe_y[p.vols[:,0]] = 0.5*(p.y[p.vols[:,1]] + p.y[p.vols[:,2]])
    pe_el[p.vols[:,0]] = 0.5*(p.elev[p.vols[:,1]] + p.elev[p.vols[:,2]])
    
    pe_x[p.vols[:,1]] = 0.5*(p.x[p.vols[:,0]] + p.x[p.vols[:,2]])
    pe_y[p.vols[:,1]] = 0.5*(p.y[p.vols[:,0]] + p.y[p.vols[:,2]])
    pe_el[p.vols[:,1]] = 0.5*(p.elev[p.vols[:,0]] + p.elev[p.vols[:,2]])

    pe_x[p.vols[:,2]] = 0.5*(p.x[p.vols[:,0]] + p.x[p.vols[:,1]])
    pe_y[p.vols[:,2]] = 0.5*(p.y[p.vols[:,0]] + p.y[p.vols[:,1]])
    pe_el[p.vols[:,2]] = 0.5*(p.elev[p.vols[:,0]] + p.elev[p.vols[:,1]])

    return [pe_x, pe_y, pe_el]


