"""
Post-processing code to compute flux through cross-sections. Implementation
currently uses an approximate method

The user must hard-code the polylines directionay (defining the cross-sections) below.

See the help for 'get_approximate_discharge_timeseries' for info on how the discharge
is computed through the polylines (including how the sign is determined).

The user can either hard-code the sww_filename, knn, and desired_ds (below), or pass them
as command line arguments to the code

This calling approach will use the hard coded values
> python flow_through_cross_sections.py 

This calling approach will use the command line arguments [in order the sww
filename, knn, and desired_ds]
> python flow_through_cross_sections.py  MODEL_OUTPUTS/RUN_XXXX/mysww.sww 1 500.0

Gareth Davies, Geoscience Australia 2014+
"""
import sys, os
import pickle
import anuga
from anuga import plot_utils as util
import scipy.spatial
import numpy
from anuga.utilities import spatialInputUtil as su
import matplotlib

# NCI hack, since interactive plotting fails there
try:
    from matplotlib import pyplot as pyplot
except:
    matplotlib.use('Agg')
    from matplotlib import pyplot as pyplot


## USER INPUT ####

# polylines can have multiple segment
polylines = {
    'Offshore': [ [666779.0, 8253357.], [673906., 8096367.], [684597., 7983715.]], 
    'CoastalInlet': [ [393277., 8104579.], [395059., 8095373.]]
            }
    
# Hard-coded values, possibly overwritten by command line arguments
sww_filename = 'MODEL_OUTPUTS/RUN_20150625_110925_cairns_excel/cairns_excel.sww'
knn = 1 # How many neighbours to use for interpolation
desired_ds = 500.0 # Spacing of points along integration lines

## END USER INPUT ###



def get_approximate_discharge_timeseries(sww_filename, 
                                         polylines,
                                         desired_ds=0.5, 
                                         k_nearest_neighbours=1,
                                         search_mesh=True, 
                                         verbose=True):
    """Given an sww_filename and a dictionary of 1D polylines, estimate the
    discharge timeseries through each polyline by interpolating the centroid
    uh/vh onto evenly spaced points on the polyline (with spacing ~ desired_ds),
    computing the flux normal to the line, and using the trapezoidal rule to
    integrate it.

    The interpolation of centroid uh/vh onto the polyline points can be either
    based on 'k-nearest-neighbours', or a direct-search of the mesh triangles.
    The former can be fast and allow for smoothing, while the latter is often
    still fast enough, and might be more accurate.

    The positive/negative discharge direction is determined from the polyline.
    Consider a river channel. If the polyline begins on the left-bank and ends
    on the right bank (left/right defined when facing downstream) then
    discharge in the downstream direction is positive.

    WARNING: The result is approximate only because ANUGA's internal edge
    fluxes are derived differently (with the reimann solver), and because the
    interpolation does not follow ANUGA's, and because your transect might not
    be exactly perpendicular to the flow. None of the methods give an exact
    result at present. 

    Errors can be significant where the solution is changing rapidly.  It may
    be worth comparing multiple cross-sections in the vicinity of the site of
    interest [covering different mesh triangles, with slightly different
    orientations].

    @param sww_filename name of sww file
    @param polylines dictionary of polylines, e.g.
            polylines = {
                'Xsection1': [ [495., 1613.], [495., 1614.], [496., 1615.] ],
                'Xsection2': [ [496., 1614.], [4968., 1615.] ]
                        }
    @param desired_ds point spacing used for trapozoidal integration on
           polylines
    @param k_nearest_neighbours number of nearest neighbours used for
           interpolation of uh/vh onto polylines
    @param search_mesh If True AND k_nearest_neighbours=1, we search the
           mesh vertices to find the triangle containing our point. Otherwise
           do nearest-neighbours on the triangle centroids to estimate the
           'nearest' triangle
    @param verbose

    @return a list of length 2 with the output_times as a numpy array, and a
            dictionary with the flow timeseries

    """

    if (search_mesh) & (k_nearest_neighbours > 1):
        msg = 'k_nearest_neighbours must be 1 when search_mesh is true'
        raise Exception(msg)

    # 2 ways to associate transect points with triangle values
    # 1) knn on centroids, or
    # 2) directly search for mesh triangles containing transect points
    # 1 can be faster + allows for smoothing, but 2 might be usually better
    use_knn = (search_mesh == False) | (k_nearest_neighbours != 1)


    if use_knn:
        # Centroids are used for knn
        p = util.get_centroids(sww_filename, timeSlices=0)
        sww_xy = numpy.vstack([p.x+p.xllcorner, p.y+p.yllcorner]).transpose()
        point_index_kdtree = scipy.spatial.cKDTree(sww_xy)
    else:
        # Vertices are used for mesh search
        p = util.get_output(sww_filename, timeSlices=0)

    # To conserve memory read from netcdf directly
    from anuga.file.netcdf import NetCDFFile
    sww_nc = NetCDFFile(sww_filename)
    ud = sww_nc.variables['xmomentum_c']
    vd = sww_nc.variables['ymomentum_c']
    output_times = sww_nc.variables['time'][:]

    discharge_series = {}

    for pk in polylines.keys():

        if verbose: print pk

        pl_full = polylines[pk]

        for segment_num in range(len(pl_full)-1):

            pl = [ pl_full[segment_num], pl_full[segment_num + 1] ]

            segment_length = ( (pl[0][0] - pl[1][0])**2 +\
                               (pl[0][1] - pl[1][1])**2 )**0.5 

            # Normal vector
            n1 = (pl[0][1] - pl[1][1])/segment_length
            n2 = -(pl[0][0] - pl[1][0])/segment_length

            # Approximate segment as npts points
            npts = int(numpy.ceil( segment_length / (desired_ds) + 1.0))
            gridXY = numpy.vstack([scipy.linspace(pl[0][0], pl[1][0], num=npts), 
                                   scipy.linspace(pl[0][1], pl[1][1], num=npts)]
                                 ).transpose()

            # Actual distance between points
            ds = (numpy.diff(gridXY[:,0])**2 + numpy.diff(gridXY[:,1])**2)**0.5
            ds_trapz = numpy.hstack([ ds[0], (ds[0:-1] + ds[1:]), ds[-1]])*0.5

            if verbose: print 'Finding triangles containing point'

            if use_knn:
                point_distance, point_indices = point_index_kdtree.query(gridXY, 
                    k = k_nearest_neighbours)

            else:
                gridXY_offset = gridXY*0.
                gridXY_offset[:,0] = gridXY[:,0] - p.xllcorner
                gridXY_offset[:,1] = gridXY[:,1] - p.yllcorner
                point_indices = numpy.zeros( gridXY.shape[0]).astype(int)
                # Provide the order to search the points (might be faster?)
                v1 = p.vols[:,0]
                search_order_update_freq = 1 
                for i in range(gridXY.shape[0]):

                    # For efficiency, we don't recompute the order to search points
                    # everytime
                    if i%search_order_update_freq==0:
                        # Update the mesh triangle search order
                        first_vertex_d2 = (p.x[v1] - gridXY_offset[i,0])**2 +\
                                          (p.y[v1] - gridXY_offset[i,1])**2
                        search_order = first_vertex_d2.argsort().tolist()
                        # Estimate how often we should update the triangle ordering
                        # Use "distance of point to vertex" / "point spacing"
                        # Crude
                        search_order_update_freq = \
                            int(numpy.ceil((first_vertex_d2[search_order[0]]**0.5)/ds[0]))
                    point_indices[i] =\
                        util.get_triangle_containing_point(p, gridXY_offset[i,:],
                            search_order = search_order) 

            if verbose: print 'Computing the flux'

            if k_nearest_neighbours == 1:
                point_uh = ud[:][:, point_indices]    
                point_vh = vd[:][:, point_indices]
            else:
                point_uh = numpy.zeros( (len(output_times), 
                                         len(point_indices[:,0])))
                point_vh = numpy.zeros( (len(output_times), 
                                         len(point_indices[:,0])))
                # Compute the inverse distance weighted uh/vh 
                numerator = point_uh*0.
                denominator = point_uh*0.
                inv_dist = 1.0/(point_distance+1.0e-12) #Avoid zero division

                # uh
                for k in range(k_nearest_neighbours): 
                    ud_data = ud[:][:,point_indices[:,k]]
                    for ti in range(len(output_times)):
                        numerator[ti,:] += ud_data[ti,:]*inv_dist[:,k]
                        denominator[ti,:] += inv_dist[:,k]
                point_uh = numerator/denominator
                
                #vh
                numerator *= 0.
                denominator *= 0.
                for k in range(k_nearest_neighbours): 
                    vd_data = vd[:][:,point_indices[:,k]]
                    for ti in range(len(output_times)):
                        numerator[ti,:] += vd_data[ti,:]*inv_dist[:,k]
                        denominator[ti,:] += inv_dist[:,k]
                point_vh = numerator/denominator
                    
            Q = [ ((point_uh[i,:]*n1 + point_vh[i,:]*n2)*ds_trapz).sum() \
                    for i in range(len(output_times)) ]

            if segment_num == 0:
                discharge_series[pk] = numpy.array(Q)
            else:
                discharge_series[pk] += numpy.array(Q)

    return [output_times, discharge_series]



def plot_discharge_timeseries(discharge_series_in, output_times, subset=None):
    """Quick-and-dirty plot of the discharge timeseries

    """

    if subset is not None:
        discharge_series = discharge_series_subset(discharge_series_in, subset)
    else:
        discharge_series = discharge_series_in
    

    ## Plot all series
    site_order = discharge_series.keys()
    line_types = ['-', '-.', '--']
    site_order.sort()
    for i, pk in enumerate(site_order):
        pyplot.plot(output_times, discharge_series[pk],
                    line_types[i%3], label=pk)
    pyplot.legend(loc=3, fontsize='xx-small')
    pyplot.plot(output_times, output_times*0.,'--',color='black')

    return

def discharge_series_subset(discharge_series, river_name_pattern):
    """Make a new discharge_series dictionary from all sites which match a
       pattern

    """

    discharge_series_keys = discharge_series.keys()
    river_keys = [ discharge_series_keys[i] \
        for i in su.matchInds(river_name_pattern, discharge_series_keys) ]

    new_discharge_series = {}
    for rk in river_keys:
        new_discharge_series[rk] = discharge_series[rk]

    return new_discharge_series

###############################################################################

if __name__ == '__main__':

    # Parse command line arguments
    if len(sys.argv)>1:
        sww_filename = sys.argv[1]
        if len(sys.argv)>2:
            knn = int(sys.argv[2])
            if len(sys.argv) > 3:
                desired_ds = float(sys.argv[3])

    if knn==1:
        search_mesh = True
    else:
        search_mesh = False

    assert os.path.exists(sww_filename), 'sww_filename not found'
    print 'sww_filename: ' + sww_filename
    print 'knn: ' + str(knn)
    print 'desired_ds: ' + str(desired_ds)
    print ''

    output_times, discharge_series = get_approximate_discharge_timeseries(
        sww_filename, polylines, desired_ds=desired_ds, 
        k_nearest_neighbours=knn, search_mesh=search_mesh)
    
    # Pickle outputs
    output_pickle = os.path.join(os.path.dirname(sww_filename), 
        'discharge_series.pkl')
    pickle.dump([output_times, discharge_series], 
                open(output_pickle, 'w'))

    # Now write in text format
    for key in discharge_series.keys():
        temp_array = numpy.vstack([output_times, discharge_series[key]]).transpose()
        numpy.savetxt(
            (os.path.join(os.path.dirname(sww_filename), key + '.csv')), 
            temp_array,
            delimiter=',')

    # Lots of plots
    #try:
    #    pyplot.ion()
    #    pyplot.figure()
    #    plot_discharge_timeseries(discharge_series, output_times)

    #    rivers = ['Offshore', 'CoastalInlet']
    #    for river in rivers:
    #        pyplot.figure()
    #        plot_discharge_timeseries(discharge_series, output_times, river)
    #        pyplot.title(river)
    #except:
    #    print 'Interactive plotting failed (expected on NCI)'
        

