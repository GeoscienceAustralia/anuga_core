'''
    Operations to extract information from an SWW file.
'''



from builtins import range
import os
import anuga.utilities.log as log
import numpy as num

from anuga.utilities.file_utils import get_all_swwfiles
from anuga.coordinate_transforms.geo_reference import Geo_reference 
from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.geometry.polygon import is_inside_polygon
from anuga.file.sww import get_mesh_and_quantities_from_file
from anuga.abstract_2d_finite_volumes.neighbour_mesh import segment_midpoints


def get_interpolated_quantities_at_polyline_midpoints(filename,
                                                      quantity_names=None,
                                                      polyline=None,
                                                      verbose=False):
    """Get values for quantities interpolated to polyline midpoints from SWW.

    filename        path to file to read
    quantity_names  quantity names to get
    polyline        representation of desired cross-section
                    may contain multiple sections allowing complex shapes
                    assume UTM coordinates
    verbose         True if this function is to be verbose

    Returns (segments, i_func)
    where segments is a list of Triangle_intersection instances
      and i_func is an instance of Interpolation_function.

    Note: For 'polyline' assume absolute UTM coordinates.

    This function is used by get_flow_through_cross_section and
    get_energy_through_cross_section.
    """

    from anuga.fit_interpolate.interpolate import Interpolation_function

    # Get mesh and quantities from sww file
    X = get_mesh_and_quantities_from_file(filename,
                                          quantities=quantity_names,
                                          verbose=verbose)
    mesh, quantities, time = X

    # Find all intersections and associated triangles.
    segments = mesh.get_intersecting_segments(polyline, verbose=verbose)

    # Get midpoints
    interpolation_points = segment_midpoints(segments)

    # Interpolate
    if verbose:
        log.critical('Interpolating - total number of interpolation points = %d'
                     % len(interpolation_points))

    I = Interpolation_function(time,
                               quantities,
                               quantity_names=quantity_names,
                               vertex_coordinates=mesh.nodes,
                               triangles=mesh.triangles,
                               interpolation_points=interpolation_points,
                               verbose=verbose)

    return segments, I


def get_flow_through_cross_section(filename, polyline, verbose=False):
    """Obtain flow (m^3/s) perpendicular to specified cross section.

    filename  path to SWW file to read
    polyline  representation of desired cross-section - it may contain
              multiple sections allowing for complex shapes. Assume
              absolute UTM coordinates.
              Format [[x0, y0], [x1, y1], ...]
    verbose   True if this function is to be verbose

    Return (time, Q)
    where time is a list of all stored times in SWW file
      and Q is a hydrograph of total flow across given segments for all
            stored times.

    The normal flow is computed for each triangle intersected by the polyline
    and added up.  Multiple segments at different angles are specified the
    normal flows may partially cancel each other.

    The typical usage of this function would be to get flow through a channel,
    and the polyline would then be a cross section perpendicular to the flow.
    """

    quantity_names =['elevation',
                     'stage',
                     'xmomentum',
                     'ymomentum']

    # Get values for quantities at each midpoint of poly line from sww file
    X = get_interpolated_quantities_at_polyline_midpoints(filename,
                                                          quantity_names=\
                                                              quantity_names,
                                                          polyline=polyline,
                                                          verbose=verbose)
    segments, interpolation_function = X

    # Get vectors for time and interpolation_points
    time = interpolation_function.time
    interpolation_points = interpolation_function.interpolation_points

    if verbose: log.critical('Computing hydrograph')

    # Compute hydrograph
    Q = []
    for t in time:
        total_flow = 0
        for i in range(len(interpolation_points)):
            elevation, stage, uh, vh = interpolation_function(t, point_id=i)
            normal = segments[i].normal

            # Inner product of momentum vector with segment normal [m^2/s]
            normal_momentum = uh*normal[0] + vh*normal[1]

            # Flow across this segment [m^3/s]
            segment_flow = normal_momentum * segments[i].length

            # Accumulate
            total_flow += segment_flow

        # Store flow at this timestep
        Q.append(total_flow)


    return time, Q

def get_flow_through_multiple_cross_sections(filename, polylines, verbose=False):
    """Obtain flow (m^3/s) perpendicular to specified cross sections.

    filename  path to SWW file to read
    polylines  representation of desired cross-sections - each cross section may contain
              multiple sections allowing for complex shapes. Assume
              absolute UTM coordinates.
              polyline format [[x0, y0], [x1, y1], ...]
              polylines format [polyline_0, polyline_1, ...,polyline_n-1]
    verbose   True if this function is to be verbose

    Return (time, Q)
    where time is a list of all stored times in SWW file
      and Q is a list of n hydrographs of total flow across given polylines for all
            stored times.

    The normal flow is computed for each triangle intersected by the polyline
    and added up.  Multiple segments at different angles are specified the
    normal flows may partially cancel each other.

    The typical usage of this function would be to get flow through a channel,
    and the polyline would then be a cross section perpendicular to the flow.
    """

    quantity_names =['elevation',
                     'stage',
                     'xmomentum',
                     'ymomentum']

    # Get values for quantities at each midpoint of poly line from sww file
    X = get_interpolated_quantities_at_multiple_polyline_midpoints(filename,
                                                          quantity_names=\
                                                              quantity_names,
                                                          polylines=polylines,
                                                          verbose=verbose)
    mult_segments, interpolation_function = X

    # Get vectors for time and interpolation_points
    time = interpolation_function.time
    interpolation_points = interpolation_function.interpolation_points

    if verbose: log.critical('Computing hydrographs')

    # Compute hydrograph
    mult_Q = []
    base_id = 0
    if verbose: print('', end=' ')
    for segments in mult_segments:
        if verbose: print('\b.', end=' ')
        Q = []
        for t in time:
            total_flow = 0.0
            i= base_id
            for segment in segments:
                elevation, stage, uh, vh = interpolation_function(t, point_id=i)
                normal = segment.normal

                # Inner product of momentum vector with segment normal [m^2/s]
                normal_momentum = uh*normal[0] + vh*normal[1]

                # Flow across this segment [m^3/s]
                segment_flow = normal_momentum * segment.length

                # Accumulate
                total_flow += segment_flow
                
                i=i+1

            # Store flow at this timestep
            Q.append(total_flow)
            
        base_id = base_id + len(segments)
        mult_Q.append(Q)

    if verbose: print()
    
    return time, mult_Q

def get_interpolated_quantities_at_multiple_polyline_midpoints(filename,
                                                      quantity_names=None,
                                                      polylines=None,
                                                      verbose=False):
    """Get values for quantities interpolated to multiple polyline midpoints from SWW.

    filename        path to file to read
    quantity_names  quantity names to get
    polylines       representation of desired cross-sections
                    may contain multiple cross sections with multiple 
                    sections allowing complex shapes
                    assume UTM coordinates
    verbose         True if this function is to be verbose

    Returns (segments, i_func)
    where a list of segments where segments are a list of Triangle_intersection instances
      and i_func is an instance of Interpolation_function.

    Note: For 'polylines' assume absolute UTM coordinates.

    This function is used by get_flow_through_multiple_cross_sections and
    get_energy_through_multipe_cross_sections.
    """

    from anuga.fit_interpolate.interpolate import Interpolation_function

    # Get mesh and quantities from sww file
    if verbose: print('Reading mesh and quantities from sww file')
    X = get_mesh_and_quantities_from_file(filename,
                                          quantities=quantity_names,
                                          verbose=verbose)
    mesh, quantities, time = X

    # Find all intersections and associated triangles.
    mult_segments = []
    if verbose: print('Intersecting segments with triangles')
    for polyline in polylines:
        mult_segments.append(mesh.get_intersecting_segments(polyline, verbose=verbose))

    # Get midpoints
    interpolation_points = []
    for segments in mult_segments:
        mid_points = segment_midpoints(segments)
        #print mid_points
        interpolation_points = interpolation_points + mid_points
     

    if verbose: 
        print('len interpolating points ', len(interpolation_points))

    # Interpolate
    if verbose:
        log.critical('Interpolating - total number of interpolation points = %d'
                     % len(interpolation_points))

    I = Interpolation_function(time,
                               quantities,
                               quantity_names=quantity_names,
                               vertex_coordinates=mesh.nodes,
                               triangles=mesh.triangles,
                               interpolation_points=interpolation_points,
                               verbose=verbose)

    #if verbose:
    #    print mult_segments
        
    return mult_segments, I



def get_energy_through_cross_section(filename,
                                     polyline,
                                     kind='total',
                                     verbose=False):
    """Obtain average energy head [m] across specified cross section.

    filename  path to file of interest
    polyline  representation of desired cross-section - it may contain multiple
              sections allowing for complex shapes. Assume absolute UTM
              coordinates.  Format [[x0, y0], [x1, y1], ...]
    kind      select energy to compute: 'specific' or 'total'
    verbose   True if this function is to be verbose

    Returns (time, E)
    where time is a list of timestep
    and E is a Average energy [m] across given segments for all stored times.

    The average velocity is computed for each triangle intersected by the
    polyline and averaged weighted by segment lengths.

    The typical usage of this function would be to get average energy of
    flow in a channel, and the polyline would then be a cross section
    perpendicular to the flow.

    #FIXME (Ole) - need name for this energy reflecting that its dimension
    is [m].
    """

    from anuga.config import g, epsilon, velocity_protection as h0

    quantity_names =['elevation',
                     'stage',
                     'xmomentum',
                     'ymomentum']

    # Get values for quantities at each midpoint of poly line from sww file
    X = get_interpolated_quantities_at_polyline_midpoints(filename,
                                                          quantity_names=\
                                                              quantity_names,
                                                          polyline=polyline,
                                                          verbose=verbose)
    segments, interpolation_function = X

    # Get vectors for time and interpolation_points
    time = interpolation_function.time
    interpolation_points = interpolation_function.interpolation_points

    if verbose: log.critical('Computing %s energy' % kind)

    # Compute total length of polyline for use with weighted averages
    total_line_length = 0.0
    for segment in segments:
        total_line_length += segment.length

    # Compute energy
    E = []
    for t in time:
        average_energy = 0.0
        for i, p in enumerate(interpolation_points):
            elevation, stage, uh, vh = interpolation_function(t, point_id=i)

            # Depth
            h = depth = stage-elevation

            # Average velocity across this segment
            if h > epsilon:
                # Use protection against degenerate velocities
                u = uh/ (h + h0/h)
                v = vh/ (h + h0/h)
            else:
                u = v = 0.0

            speed_squared = u*u + v*v
            kinetic_energy = 0.5 * speed_squared/ g

            if kind == 'specific':
                segment_energy = depth + kinetic_energy
            elif kind == 'total':
                segment_energy = stage + kinetic_energy
            else:
                msg = 'Energy kind must be either "specific" or "total". '
                msg += 'I got %s' % kind

            # Add to weighted average
            weigth = segments[i].length/ total_line_length
            average_energy += segment_energy * weigth

        # Store energy at this timestep
        E.append(average_energy)

    return time, E


def get_maximum_inundation_elevation(filename,
                                     polygon=None,
                                     time_interval=None,
                                     verbose=False):
    """Return highest elevation where depth > 0

    filename       path to SWW file containing ANUGA model output
    polygon        if specified restrict to points inside this polygon
    time_interval  if specified restrict to within the period specified
    verbose        True if this function is  to be verbose

    If no inundation is found within polygon and time_interval the return value
    is None signifying "No Runup" or "Everything is dry".

    Usage:
    max_runup = get_maximum_inundation_elevation(filename,
                                                 polygon=None,
                                                 time_interval=None,
                                                 verbose=False)

    See general function get_maximum_inundation_data for details.
    """

    runup, locatoin = get_maximum_inundation_data(filename,
                                           polygon=polygon,
                                           time_interval=time_interval,
                                           verbose=verbose)
    return runup


def get_maximum_inundation_location(filename,
                                    polygon=None,
                                    time_interval=None,
                                    verbose=False):
    """Return location of highest elevation where h > 0

    filename       path to SWW file containing ANUGA model output
    polygon        if specified restrict to points inside this polygon
    time_interval  if specified restrict to within the period specified
    verbose        True if this function is  to be verbose

    If no inundation is found within polygon and time_interval the return value
    is None signifying "No Runup" or "Everything is dry".

    Usage:
    max_runup_location = get_maximum_inundation_location(filename,
                                                         polygon=None,
                                                         time_interval=None,
                                                         verbose=False)

    See general function get_maximum_inundation_data for details.
    """

    runup, max_loc = get_maximum_inundation_data(filename,
                                             polygon=polygon,
                                             time_interval=time_interval,
                                             verbose=verbose)
    return max_loc


def get_maximum_inundation_data(filename, polygon=None, time_interval=None,
                                use_centroid_values=True,
                                return_time= False,
                                verbose=False):
    """Compute maximum run up height from sww file.

    filename             path to SWW file to read
    polygon              if specified resrict to points inside this polygon
                         assumed absolute coordinates and in same zone as
                         domain
    time_interval        if specified resrict to within the period specified
    use_centroid_values 
    verbose              True if this function is to be verbose

    Returns (maximal_runup, maximal_runup_location).

    Usage:
    runup, location = get_maximum_inundation_data(filename,
                                                  polygon=None,
                                                  time_interval=None,
                                                  verbose=False)

    Algorithm is as in get_maximum_inundation_elevation from
    shallow_water_domain except that this function works with the SWW file and
    computes the maximal runup height over multiple timesteps.

    If no inundation is found within polygon and time_interval the return value
    is None signifying "No Runup" or "Everything is dry".
    """

    # We are using nodal values here as that is what is stored in sww files.

    # Water depth below which it is considered to be 0 in the model
    # FIXME (Ole): Allow this to be specified as a keyword argument as well

    from anuga.geometry.polygon import inside_polygon
    from anuga.config import minimum_allowed_height
    from anuga.file.netcdf import NetCDFFile

    # Just find max inundation over one file
    dir, base = os.path.split(filename)
    #iterate_over = get_all_swwfiles(dir, base)
    iterate_over = [ filename[:-4] ]
    if verbose:
        print(iterate_over)
        
    # Read sww file
    if verbose: log.critical('Reading from %s' % filename)
    # FIXME: Use general swwstats (when done)

    maximal_runup = None
    maximal_runup_location = None
    maximal_time = None


    for _, swwfile in enumerate (iterate_over):
        # Read sww file
        filename = os.path.join(dir, swwfile+'.sww')

        if verbose: log.critical('Reading from %s' % filename)
        # FIXME: Use general swwstats (when done)

        fid = NetCDFFile(filename)

        # Get geo_reference
        # sww files don't have to have a geo_ref
        try:
            geo_reference = Geo_reference(NetCDFObject=fid)
        except AttributeError:
            geo_reference = Geo_reference() # Default georef object

        xllcorner = geo_reference.get_xllcorner()
        yllcorner = geo_reference.get_yllcorner()

        # Get extent
        volumes = fid.variables['volumes'][:]
        x = fid.variables['x'][:] + xllcorner
        y = fid.variables['y'][:] + yllcorner

        # Get the relevant quantities (Convert from single precison)
        try:
            elevation = num.array(fid.variables['elevation_c'][:], float)
            stage     = num.array(fid.variables['stage_c'][:], float)
            found_c_values = True
        except:
            elevation = num.array(fid.variables['elevation'][:], float)
            stage     = num.array(fid.variables['stage'][:], float)
            found_c_values = False

        if verbose:
            print('found c values ', found_c_values)
            print('stage.shape ',stage.shape)
            print('elevation.shape ',elevation.shape)
            
        # Here's where one could convert nodal information to centroid
        # information but is probably something we need to write in C.
        # Here's a Python thought which is NOT finished!!!
        if use_centroid_values is True:
            vols0=volumes[:,0]
            vols1=volumes[:,1]
            vols2=volumes[:,2]
            # Then use these to compute centroid location 
            x=(x[vols0]+x[vols1]+x[vols2])/3.0
            y=(y[vols0]+y[vols1]+y[vols2])/3.0

            if found_c_values:
                # don't have to do anything as found in sww file
                pass
            else:
                elevation=(elevation[vols0]+elevation[vols1]+elevation[vols2])/3.0
                stage=(stage[:,vols0]+stage[:,vols1]+stage[:,vols2])/3.0

        # Spatial restriction
        if polygon is not None:
            msg = 'polygon must be a sequence of points.'
            assert len(polygon[0]) == 2, msg
            # FIXME (Ole): Make a generic polygon input check in polygon.py
            # and call it here
            points = num.ascontiguousarray(num.concatenate((x[:, num.newaxis],
                                                            y[:, num.newaxis]),
                                                            axis=1))
            point_indices = inside_polygon(points, polygon)

            # Restrict quantities to polygon
            elevation = num.take(elevation, point_indices, axis=0)
            stage = num.take(stage, point_indices, axis=1)

            # Get info for location of maximal runup
            points_in_polygon = num.take(points, point_indices, axis=0)
            x = points_in_polygon[:,0]
            y = points_in_polygon[:,1]
        else:
            # Take all points
            point_indices = num.arange(len(x))

        # Temporal restriction
        time = fid.variables['time'][:]
        if verbose:
            print(time)
        all_timeindices = num.arange(len(time))
        
        if time_interval is not None:
            msg = 'time_interval must be a sequence of length 2.'
            assert len(time_interval) == 2, msg
            msg = 'time_interval %s must not be decreasing.' % time_interval
            assert time_interval[1] >= time_interval[0], msg
            msg = 'Specified time interval [%.8f:%.8f] ' % tuple(time_interval)
            msg += 'must does not match model time interval: [%.8f, %.8f]\n' \
                   % (time[0], time[-1])
            if time_interval[1] < time[0]:
                fid.close()
                raise ValueError(msg)
            if time_interval[0] > time[-1]:
                fid.close()
                raise ValueError(msg)

            # Take time indices corresponding to interval (& is bitwise AND)
            timesteps = num.compress((time_interval[0] <= time) \
                                     & (time <= time_interval[1]),
                                     all_timeindices)

            msg = 'time_interval %s did not include any model timesteps.' \
                  % time_interval
            assert not num.all(timesteps == 0), msg
        else:
            # Take them all
            timesteps = all_timeindices
        
        #print timesteps

        fid.close()

        # Compute maximal runup for each timestep
        #maximal_runup = None
        #maximal_runup_location = None
        #maximal_runups = [None]
        #maximal_runup_locations = [None]


        for i in timesteps:
            ## if use_centroid_values is True:
            ##     stage_i  = stage[i,:]
            ## else:
            ##     stage_i = stage[i,:]

            stage_i = stage[i,:]
            depth = stage_i - elevation

            if verbose:
                print('++++++++')
            # Get wet nodes i.e. nodes with depth>0 within given region
            # and timesteps
            wet_nodes = num.where(depth > 0.0)[0]


            if verbose:
                print(stage_i.shape)
                print(num.max(stage_i))
                #print max(wet_elevation)


            if num.all(wet_nodes == 0):
                runup = None
            else:
                # Find maximum elevation among wet nodes
                wet_elevation = num.take(elevation, wet_nodes, axis=0)


                if verbose:
                    pass
                    #print wet_elevation
                    
                runup_index = num.argmax(wet_elevation)
                runup = max(wet_elevation)
                if verbose:
                    print('max(wet_elevation) ',max(wet_elevation))
                assert wet_elevation[runup_index] == runup       # Must be True


            # Python3.8 no longer supports things like 3.333 > None or None > None
            # so this elegant conditional: if runup > maximal_runup
            # had to be unpacked thusly
            if maximal_runup is None and runup is None:
                # First time around
                pass
            elif maximal_runup is None or runup > maximal_runup:
                maximal_runup = runup      # works even if maximal_runup is None
                maximal_time = time[i]

                # Record location
                wet_x = num.take(x, wet_nodes, axis=0)
                wet_y = num.take(y, wet_nodes, axis=0)
                maximal_runup_location = [wet_x[runup_index],
                                          wet_y[runup_index]]
            else:
                pass
            
            if verbose:
                print(i, runup)

    if return_time:
        return maximal_runup, maximal_runup_location, maximal_time
    else:
        return maximal_runup, maximal_runup_location
        

