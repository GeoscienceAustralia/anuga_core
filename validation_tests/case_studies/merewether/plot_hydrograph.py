from anuga.utilities import plot_utils as util
from matplotlib import pyplot as pyplot
import numpy
from anuga.file.sww import get_mesh_and_quantities_from_file
from anuga.abstract_2d_finite_volumes.neighbour_mesh import segment_midpoints




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

    print 'In get_flow_through_cross_section'
    quantity_names =['elevation',
                     'stage',
                     'xmomentum',
                     'ymomentum']

    # Get values for quantities at each midpoint of poly line from sww file
#     X = get_interpolated_quantities_at_polyline_midpoints(filename,
#                                                           quantity_names=\
#                                                               quantity_names,
#                                                           polyline=polyline,
#                                                           verbose=verbose)
    from anuga.fit_interpolate.interpolate import Interpolation_function

    print 'In get_interpolated_quantities_at_polyline_midpoints'
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
        print 'Interpolating - total number of interpolation points = %d' % len(interpolation_points)

    interpolation_function = Interpolation_function(time,
                               quantities,
                               quantity_names=quantity_names,
                               vertex_coordinates=mesh.nodes,
                               triangles=mesh.triangles,
                               interpolation_points=interpolation_points,
                               verbose=verbose)

    # Get vectors for time and interpolation_points
    time = interpolation_function.time
    interpolation_points = interpolation_function.interpolation_points

    if verbose: print 'Computing hydrograph'

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


#===================================================================================

swwfile = 'merewether_1m.sww'

pyplot.clf()

xll = 382250.0
yll = 6354265.0
cross_section = [[103+xll, 100.+yll], [130.+xll,80.+yll]]

#from anuga import get_flow_through_cross_section
time, Q = get_flow_through_cross_section(swwfile,
                                         cross_section,
                                         verbose=True)

pyplot.plot(time,Q)

#pyplot.gca().set_aspect('equal')
pyplot.title("Flow through transect 1")

#pyplot.legend(loc='upper left')
pyplot.savefig('flowrate.png',dpi=100,bbox_inches='tight')


