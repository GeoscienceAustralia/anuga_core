"""Least squares interpolation.

   These functions and classes calculate a value at a particular point on
   the given mesh. It interpolates the values stored at the vertices of the
   mesh.

   For example, if you want to get the height of a terrain mesh at particular
   point, you pass the point to an Interpolate class. The point will intersect
   one of the triangles on the mesh, and the interpolated height will be an
   intermediate value between the three vertices of that triangle.
   This value is returned by the class.

   Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
   Geoscience Australia, 2004.

DESIGN ISSUES
* what variables should be global?
- if there are no global vars functions can be moved around alot easier

* The public interface to Interpolate
__init__
interpolate
interpolate_block

"""

import time
import os
import sys
from warnings import warn
from math import sqrt
from csv import writer, DictWriter

from anuga.caching.caching import cache
from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
from anuga.utilities.sparse import Sparse, Sparse_CSR
from anuga.utilities.cg_solve import conjugate_gradient, VectorShapeError
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.utilities.numerical_tools import ensure_numeric, NAN
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.geospatial_data.geospatial_data import ensure_absolute
from anuga.pmesh.mesh_quadtree import MeshQuadtree
from anuga.fit_interpolate.general_fit_interpolate import FitInterpolate
from anuga.abstract_2d_finite_volumes.file_function import file_function
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, epsilon
from anuga.geometry.polygon import interpolate_polyline, in_and_outside_polygon
import anuga.utilities.log as log


import numpy as num


# Interpolation specific exceptions

class Modeltime_too_late(BaseException): pass
class Modeltime_too_early(BaseException): pass


def interpolate(vertex_coordinates,
                triangles,
                vertex_values,
                interpolation_points,
                mesh_origin=None,
                start_blocking_len=500000,
                use_cache=False,
                verbose=False,
                output_centroids=False):
    """Interpolate vertex_values to interpolation points.

    Inputs (mandatory):


    vertex_coordinates: List of coordinate pairs [xi, eta] of
                        points constituting a mesh
                        (or an m x 2 numeric array or
                        a geospatial object)
                        Points may appear multiple times
                        (e.g. if vertices have discontinuities)

    triangles: List of 3-tuples (or a numeric array) of
               integers representing indices of all vertices
               in the mesh.

    vertex_values: Vector or array of data at the mesh vertices.
                   If array, interpolation will be done for each column as
                   per underlying matrix-matrix multiplication

    interpolation_points: Interpolate mesh data to these positions.
                          List of coordinate pairs [x, y] of
                          data points or an nx2 numeric array or a
                          Geospatial_data object

    Inputs (optional)

    mesh_origin: A geo_reference object or 3-tuples consisting of
                 UTM zone, easting and northing.
                 If specified vertex coordinates are assumed to be
                 relative to their respective origins.

                           Note: Don't supply a vertex coords as a geospatial
                           object and a mesh origin, since geospatial has its
                           own mesh origin.

    start_blocking_len: If the # of points is more or greater than this,
                        start blocking

    use_cache: True or False


    Output:

    Interpolated values at specified point_coordinates

    Note: This function is a simple shortcut for case where
    interpolation matrix is unnecessary
    Note: This function does not take blocking into account,
    but allows caching.

    """

    # FIXME(Ole): Probably obsolete since I is precomputed and
    #             interpolate_block caches

    from anuga.caching import cache

    # Create interpolation object with matrix
    args = (ensure_numeric(vertex_coordinates, float),
            ensure_numeric(triangles))
    kwargs = {'mesh_origin': mesh_origin,
              'verbose': verbose}

    if use_cache is True:
        I = cache(Interpolate, args, kwargs, verbose=verbose)
    else:
        I = Interpolate(*args, **kwargs)

    # Call interpolate method with interpolation points
    result = I.interpolate_block(vertex_values, interpolation_points,
                                 use_cache=use_cache,
                                 verbose=verbose,
                                 output_centroids=output_centroids)

    return result


class Interpolate (FitInterpolate):

    def __init__(self,
                 vertex_coordinates,
                 triangles,
                 mesh_origin=None,
                 verbose=False):

        """ Build interpolation matrix mapping from
        function values at vertices to function values at data points

        Inputs:
          vertex_coordinates: List of coordinate pairs [xi, eta] of
              points constituting a mesh (or an m x 2 numeric array or
              a geospatial object)
              Points may appear multiple times
              (e.g. if vertices have discontinuities)

          triangles: List of 3-tuples (or a numeric array) of
              integers representing indices of all vertices in the mesh.

          mesh_origin: A geo_reference object or 3-tuples consisting of
              UTM zone, easting and northing.
              If specified vertex coordinates are assumed to be
              relative to their respective origins.

          max_vertices_per_cell: Number of vertices in a quad tree cell
          at which the cell is split into 4.

          Note: Don't supply a vertex coords as a geospatial object and
              a mesh origin, since geospatial has its own mesh origin.
        """

        # FIXME (Ole): Need an input check

        FitInterpolate.__init__(self,
                                vertex_coordinates=vertex_coordinates,
                                triangles=triangles,
                                mesh_origin=mesh_origin,
                                verbose=verbose)

        # Initialise variables
        self._A_can_be_reused = False  # FIXME (Ole): Probably obsolete
        self._point_coordinates = None # FIXME (Ole): Probably obsolete
        self.interpolation_matrices = {} # Store precomputed matrices


    # FIXME: What is a good start_blocking_len value?
    def interpolate(self,
                    f,
                    point_coordinates=None,
                    start_blocking_len=500000,
                    NODATA_value = NAN,
                    verbose=False,
                    output_centroids=False):
        """Interpolate mesh data f to determine values, z, at points.

        f is the data on the mesh vertices.

        The mesh values representing a smooth surface are
        assumed to be specified in f.

        Inputs:
          f: Vector or array of data at the mesh vertices.
              If f is an array, interpolation will be done for each column as
              per underlying matrix-matrix multiplication

          point_coordinates: Interpolate mesh data to these positions.
              List of coordinate pairs [x, y] of
              data points or an nx2 numeric array or a Geospatial_data object

              If point_coordinates is absent, the points inputted last time
              this method was called are used, if possible.

          start_blocking_len: If the # of points is more or greater than this,
              start blocking

        Output:
          Interpolated values at inputted points (z).
        """

        # FIXME (Ole): Why is the interpolation matrix rebuilt everytime the
        # method is called even if interpolation points are unchanged.
        # This really should use some kind of caching in cases where
        # interpolation points are reused.
        #
        # This has now been addressed through an attempt in interpolate_block

        if verbose: log.critical('Build intepolation object')
        if isinstance(point_coordinates, Geospatial_data):
            point_coordinates = point_coordinates.get_data_points(absolute=True)

        # Can I interpolate, based on previous point_coordinates?
        if point_coordinates is None:
            if self._A_can_be_reused is True \
               and len(self._point_coordinates) < start_blocking_len:
                z = self._get_point_data_z(f, NODATA_value=NODATA_value, verbose=verbose)
            elif self._point_coordinates is not None:
                #     if verbose, give warning
                if verbose:
                    log.critical('WARNING: Recalculating A matrix, '
                                 'due to blocking.')
                point_coordinates = self._point_coordinates
            else:
                # There are no good point_coordinates. import sys; sys.exit()
                msg = 'ERROR (interpolate.py): No point_coordinates inputted'
                raise Exception(msg)

        if point_coordinates is not None:
            self._point_coordinates = point_coordinates
            if len(point_coordinates) < start_blocking_len \
               or start_blocking_len == 0:
                self._A_can_be_reused = True
                z = self.interpolate_block(f, point_coordinates, NODATA_value = NODATA_value,
                                           verbose=verbose, output_centroids=output_centroids)
            else:
                # Handle blocking
                self._A_can_be_reused = False
                start = 0
                # creating a dummy array to concatenate to.

                f = ensure_numeric(f, float)
                if len(f.shape) > 1:
                    z = num.zeros((0, f.shape[1]), int)     #array default#
                else:
                    z = num.zeros((0,), int)        #array default#

                for end in range(start_blocking_len,
                                 len(point_coordinates),
                                 start_blocking_len):
                    t = self.interpolate_block(f, point_coordinates[start:end], NODATA_value=NODATA_value,
                                               verbose=verbose, output_centroids=output_centroids)
                    z = num.concatenate((z, t), axis=0)    #??default#
                    start = end

                end = len(point_coordinates)
                t = self.interpolate_block(f, point_coordinates[start:end], NODATA_value=NODATA_value,
                                           verbose=verbose, output_centroids=output_centroids)
                z = num.concatenate((z, t), axis=0)    #??default#
        return z


    def interpolate_block(self, f, point_coordinates, NODATA_value=NAN,
                          use_cache=False, verbose=False, output_centroids=False):
        """
        Call this if you want to control the blocking or make sure blocking
        doesn't occur.

        Return the point data, z.

        See interpolate for doc info.
        """

        # FIXME (Ole): I reckon we should change the interface so that
        # the user can specify the interpolation matrix instead of the
        # interpolation points to save time.

        if isinstance(point_coordinates, Geospatial_data):
            point_coordinates = point_coordinates.get_data_points(absolute=True)

        # Convert lists to numeric arrays if necessary
        point_coordinates = ensure_numeric(point_coordinates, float)
        f = ensure_numeric(f, float)

        from anuga.caching import myhash
        import sys

        if use_cache is True:
            if sys.platform != 'win32':
                # FIXME (Ole): (Why doesn't this work on windoze?)
                # Still absolutely fails on Win 24 Oct 2008

                X = cache(self._build_interpolation_matrix_A,
                          args=(point_coordinates, output_centroids),
                          kwargs={'verbose': verbose},
                          verbose=verbose)
            else:
                # FIXME
                # Hash point_coordinates to memory location, reuse if possible
                # This will work on Linux as well if we want to use it there.
                key = myhash(point_coordinates)

                reuse_A = False

                if key in self.interpolation_matrices:
                    X, stored_points = self.interpolation_matrices[key]
                    if num.all(stored_points == point_coordinates):
                        reuse_A = True                # Reuse interpolation matrix

                if reuse_A is False:
                    X = self._build_interpolation_matrix_A(point_coordinates,
                                                           output_centroids,
                                                           verbose=verbose)
                    self.interpolation_matrices[key] = (X, point_coordinates)
        else:
            X = self._build_interpolation_matrix_A(point_coordinates, output_centroids,
                                                   verbose=verbose)

        # Unpack result
        self._A, self.inside_poly_indices, self.outside_poly_indices, self.centroids = X
        # Check that input dimensions are compatible
        msg = 'Two columns must be specified in point coordinates. ' \
              'I got shape=%s' % (str(point_coordinates.shape))
        assert point_coordinates.shape[1] == 2, msg

        msg = 'The number of rows in matrix A must be the same as the '
        msg += 'number of points supplied.'
        msg += ' I got %d points and %d matrix rows.' \
               % (point_coordinates.shape[0], self._A.shape[0])
        assert point_coordinates.shape[0] == self._A.shape[0], msg

        msg = 'The number of columns in matrix A must be the same as the '
        msg += 'number of mesh vertices.'
        msg += ' I got %d vertices and %d matrix columns.' \
               % (f.shape[0], self._A.shape[1])
        assert self._A.shape[1] == f.shape[0], msg

        # Compute Matrix vector product and return
        return self._get_point_data_z(f, NODATA_value=NODATA_value)


    def get_outside_poly_indices(self):
        """
        Return index of those data points outside (and in holes)
        the mesh

        Precondition: interpolation or interpolation_block has been called
        """
        return self.outside_poly_indices


    def _get_point_data_z(self, f, NODATA_value=NAN, verbose=False):
        """
        Return the point data, z.

        Precondition: The _A matrix has been created
        """

        z = self._A * f

        # Taking into account points outside the mesh.
        for i in self.outside_poly_indices:
            z[i] = NODATA_value
        return z


    def _build_interpolation_matrix_A(self,
                                      point_coordinates,
                                      output_centroids=False,
                                      verbose=False):
        """Build n x m interpolation matrix, where
        n is the number of data points and
        m is the number of basis functions phi_k (one per vertex)

        This algorithm uses a quad tree data structure for fast binning
        of data points
        origin is a 3-tuple consisting of UTM zone, easting and northing.
        If specified coordinates are assumed to be relative to this origin.

        This one will override any data_origin that may be specified in
        instance interpolation

        Preconditions:
            Point_coordindates and mesh vertices have the same origin.
        """

        if verbose: log.critical('Building interpolation matrix')

        # Convert point_coordinates to numeric arrays, in case it was a list.
        point_coordinates = ensure_numeric(point_coordinates, float)

        if verbose: log.critical('Getting indices inside mesh boundary')

        # Quick test against boundary, but will not deal with holes in the mesh,
        # that is done below
        inside_boundary_indices, outside_poly_indices = \
            in_and_outside_polygon(point_coordinates,
                                   self.mesh.get_boundary_polygon(),
                                   closed=True, verbose=verbose)

        # Build n x m interpolation matrix
        if verbose and len(outside_poly_indices) > 0:
            log.critical('WARNING: Points outside mesh boundary.')

        # Since you can block, throw a warning, not an error.
        if verbose and 0 == len(inside_boundary_indices):
            log.critical('WARNING: No points within the mesh!')

        m = self.mesh.number_of_nodes  # Nbr of basis functions (1/vertex)
        n = point_coordinates.shape[0] # Nbr of data points

        if verbose: log.critical('Number of datapoints: %d' % n)
        if verbose: log.critical('Number of basis functions: %d' % m)

        A = Sparse(n,m)

        n = len(inside_boundary_indices)

        centroids = []
        inside_poly_indices = []

        # Compute matrix elements for points inside the mesh
        if verbose: log.critical('Building interpolation matrix from %d points'
                                 % n)

        for d, i in enumerate(inside_boundary_indices):
            # For each data_coordinate point
            if verbose and d % ((n+10) / 10) == 0: log.critical('Doing %d of %d'
                                                                %(d, n))

            x = point_coordinates[i]
            element_found, sigma0, sigma1, sigma2, k = self.root.search_fast(x)
            # Update interpolation matrix A if necessary
            if element_found is True:

                #if verbose:
                #    print 'Point is within mesh:', d, i

                inside_poly_indices.append(i)

                # Assign values to matrix A
                j0 = self.mesh.triangles[k,0] # Global vertex id for sigma0
                j1 = self.mesh.triangles[k,1] # Global vertex id for sigma1
                j2 = self.mesh.triangles[k,2] # Global vertex id for sigma2
                js = [j0, j1, j2]

                if output_centroids is False:
                    # Weight each vertex according to its distance from x
                    sigmas = {j0:sigma0, j1:sigma1, j2:sigma2}
                    for j in js:
                        A[i, j] = sigmas[j]
                else:
                    # If centroids are needed, weight all 3 vertices equally
                    for j in js:
                        A[i, j] = 1.0/3.0
                    centroids.append(self.mesh.centroid_coordinates[k])
            else:
                if verbose:
                    log.critical('Mesh has a hole - moving this point to outside list')

                # This is a numpy arrays, so we need to do a slow transfer
                outside_poly_indices = num.append(outside_poly_indices, [i], axis=0)

        return A, inside_poly_indices, outside_poly_indices, centroids






def benchmark_interpolate(vertices,
                          vertex_attributes,
                          triangles, points,
                          max_points_per_cell=None,
                          start_blocking_len=500000,
                          mesh_origin=None):
    """
    points: Interpolate mesh data to these positions.
            List of coordinate pairs [x, y] of
            data points or an nx2 numeric array or a Geospatial_data object

    No test for this yet.
    Note, this has no time the input data has no time dimension.  Which is
    different from most of the data we interpolate, eg sww info.

    Output:
        Interpolated values at inputted points.
    """

    interp = Interpolate(vertices,
                         triangles,
                         max_vertices_per_cell=max_points_per_cell,
                         mesh_origin=mesh_origin)

    calc = interp.interpolate(vertex_attributes,
                              points,
                              start_blocking_len=start_blocking_len)




class CSV_files:
    """Wrapper class to handle the five separate CSV files used by function interpolate_sww2csv
       Notable, it'll make sure they are opened and closed properly.
    """

    def __init__(self,
                 depth_file,
                 velocity_x_file,
                 velocity_y_file,
                 stage_file=None,
                 froude_file=None):

        # Open files and assign to file handles
        self.depth_FH = open(depth_file, 'w', newline="")
        self.velocity_x_FH = open(velocity_x_file, 'w', newline="")
        self.velocity_y_FH = open(velocity_y_file, 'w', newline="")

        if stage_file is not None:
            self.stage_FH = open(stage_file, 'w', newline="")
        else:
            self.stage_FH = None

        if froude_file is not None:
            self.froude_FH = open(froude_file, 'w', newline="")
        else:
            self.froude_FH = None

        # Create CSV writers
        self.depth_writer = writer(self.depth_FH)
        self.velocity_x_writer = writer(self.velocity_x_FH)
        self.velocity_y_writer = writer(self.velocity_y_FH)
        if stage_file is not None:
            stage_writer = writer(self.stage_FH)
        if froude_file is not None:
            froude_writer = writer(self.froude_FH)

    def close_all(self):
        """Close all CSV file handles
        """

        self.depth_FH.close()
        self.velocity_x_FH.close()
        self.velocity_y_FH.close()

        if self.stage_FH is not None:
            self.stage_FH.close()
        if self.froude_FH is not None:
            self.froude_FH.close()

    def write_headings(self, heading):
        """Write heading to all CSV files"""

        self.depth_writer.writerow(heading)
        self.velocity_x_writer.writerow(heading)
        self.velocity_y_writer.writerow(heading)
        if self.stage_FH is not None:
            self.stage_writer.writerow(heading)
        if self.froude_FH is not None:
            self.froude_writer.writerow(heading)

    def write_row(self, depths, velocity_xs, velocity_ys, stages, froudes):
        """Write one row for each CSV file"""

        self.depth_writer.writerow(depths)
        self.velocity_x_writer.writerow(velocity_xs)
        self.velocity_y_writer.writerow(velocity_ys)
        if self.stage_FH is not None:
            self.stage_writer.writerow(stages)
        if self.froude_FH is not None:
            self.froude_writer.writerow(froudes)




def interpolate_sww2csv(sww_file,
                        points,
                        depth_file,
                        velocity_x_file,
                        velocity_y_file,
                        stage_file=None,
                        froude_file=None,
                        time_thinning=1,
                        g = 9.80665,
                        verbose=True,
                        use_cache = True):
    """
    Interpolate the quantities at a given set of locations, given
    an sww file.
    The results are written to csv files.

    sww_file is the input sww file.
    points is a list of the 'gauges' x,y location.
    depth_file is the name of the output depth file
    velocity_x_file is the name of the output x velocity file.
    velocity_y_file is the name of the output y velocity file.
    stage_file is the name of the output stage file.

    In the csv files columns represents the gauges and each row is a
    time slice.

    Time_thinning_number controls how many timesteps to use. Only
    timesteps with index%time_thinning_number == 0 will used, or
    in other words a value of 3, say, will cause the algorithm to
    use every third time step.

    In the future let points be a points file.
    And let the user choose the quantities.

    This is currently quite specific.
    If it is need to be more general, change things.
    """

    quantities =  ['stage', 'elevation', 'xmomentum', 'ymomentum']
    points = ensure_absolute(points)
    point_count = len(points)
    callable_sww = file_function(sww_file,
                                 quantities=quantities,
                                 interpolation_points=points,
                                 verbose=verbose,
                                 time_thinning=time_thinning,
                                 use_cache=use_cache)


    csv_files = CSV_files(depth_file, velocity_x_file, velocity_y_file,
                          stage_file=stage_file,
                          froude_file=froude_file)

    # Write heading
    heading = [str(x[0])+ ':' + str(x[1]) for x in points]
    heading.insert(0, 'time')

    csv_files.write_headings(heading)
    for time in callable_sww.get_time():
        depths = [time]
        velocity_xs = [time]
        velocity_ys = [time]

        stages = [time]   # May not be used if stage file is None, but makes code below simpler
        froudes = [time]  # May not be used if stage file is None, but makes code below simpler

        for point_i, point in enumerate(points):
            quantities = callable_sww(time,point_i)

            w = quantities[0]
            z = quantities[1]
            momentum_x = quantities[2]
            momentum_y = quantities[3]
            depth = w - z

            if w == NAN or z == NAN or momentum_x == NAN:
                velocity_x = NAN
            else:
                if depth > epsilon:
                    velocity_x = momentum_x / depth  # Absolute velocity
                else:
                    velocity_x = 0

            if w == NAN or z == NAN or momentum_y == NAN:
                velocity_y = NAN
            else:
                if depth > epsilon:
                    velocity_y = momentum_y / depth  # Absolute velocity
                else:
                    velocity_y = 0

            if depth < epsilon:
                froude = NAN
            else:
                froude = sqrt(velocity_x*velocity_x + velocity_y*velocity_y) / sqrt(depth * g) # gravity m/s/s

            # Append values to lists
            depths.append(depth)
            velocity_xs.append(velocity_x)
            velocity_ys.append(velocity_y)
            stages.append(w)
            froudes.append(froude)

        csv_files.write_row(depths, velocity_xs, velocity_ys, stages, froudes)

    # Clean up (force the file handles inside the writers to close
    csv_files.close_all()


class Interpolation_function(object):
    """Interpolation_interface - creates callable object f(t, id) or f(t, x, y)
    which is interpolated from time series defined at vertices of
    triangular mesh (such as those stored in sww files)

    Let m be the number of vertices, n the number of triangles
    and p the number of timesteps.
    Also, let N be the number of interpolation points.

    Mandatory input
        time:                 px1 array of monotonously increasing times (float)
        quantities:           Dictionary of arrays or 1 array (float)
                              The arrays must either have dimensions pxm or mx1.
                              The resulting function will be time dependent in
                              the former case while it will be constant with
                              respect to time in the latter case.

    Optional input:
        quantity_names:       List of keys into the quantities dictionary for
                              imposing a particular order on the output vector.
        vertex_coordinates:   mx2 array of coordinates (float)
        triangles:            nx3 array of indices into vertex_coordinates (int)
        interpolation_points: Nx2 array of coordinates to be interpolated to
        verbose:              Level of reporting

    The quantities returned by the callable object are specified by
    the list quantities which must contain the names of the
    quantities to be returned and also reflect the order, e.g. for
    the shallow water wave equation, on would have
    quantities = ['stage', 'xmomentum', 'ymomentum']

    The parameter interpolation_points decides at which points interpolated
    quantities are to be computed whenever object is called.
    If None, return average value

    FIXME (Ole): Need to allow vertex coordinates and interpolation points to
                 be geospatial data objects

    (FIXME (Ole): This comment should be removed)
    Time assumed to be relative to starttime
    All coordinates assume origin of (0,0) - e.g. georeferencing must be
    taken care of outside this function
    """

    def __init__(self,
                 time,
                 quantities,
                 quantity_names=None,
                 vertex_coordinates=None,
                 triangles=None,
                 interpolation_points=None,
                 time_thinning=1,
                 verbose=False,
                 gauge_neighbour_id=None,
                 output_centroids=False):
        """Initialise object and build spatial interpolation if required

        Time_thinning_number controls how many timesteps to use. Only timesteps
        with index%time_thinning_number == 0 will used, or in other words a
        value of 3, say, will cause the algorithm to use every third time step.
        """

        from anuga.config import time_format

        if verbose is True:
            log.critical('Interpolation_function: input checks')

        # Check temporal info
        time = ensure_numeric(time)

        if not num.all(time[1:] - time[:-1] >= 0):
            # This message is time consuming to form due to the conversion of
            msg = 'Time must be a monotonuosly increasing sequence %s' % time
            raise Exception(msg)

        # Check if quantities is a single array only
        if not isinstance(quantities, dict):
            quantities = ensure_numeric(quantities)
            quantity_names = ['Attribute']

            # Make it a dictionary
            quantities = {quantity_names[0]: quantities}

        # Use keys if no names are specified
        if quantity_names is None:
            quantity_names = list(quantities.keys())

        # Check spatial info
        if vertex_coordinates is None:
            self.spatial = False
        else:
            # FIXME (Ole): Try ensure_numeric here -
            #              this function knows nothing about georefering.
            vertex_coordinates = ensure_absolute(vertex_coordinates)

            if triangles is not None:
                triangles = ensure_numeric(triangles)
            self.spatial = True

        if verbose is True:
            log.critical('Interpolation_function: thinning by %d'
                         % time_thinning)


        # Thin timesteps if needed
        # Note array() is used to make the thinned arrays contiguous in memory
        self.time = num.array(time[::time_thinning])
        for name in quantity_names:
            if len(quantities[name].shape) == 2:
                quantities[name] = num.array(quantities[name][::time_thinning,:])

        if verbose is True:
            log.critical('Interpolation_function: precomputing')

        # Save for use with statistics
        self.quantities_range = {}
        for name in quantity_names:
            q = quantities[name][:].flatten()
            self.quantities_range[name] = [min(q), max(q)]

        self.quantity_names = quantity_names
        self.vertex_coordinates = vertex_coordinates
        self.interpolation_points = interpolation_points

        self.index = 0    # Initial time index
        self.precomputed_values = {}
        self.centroids = []

        # Precomputed spatial interpolation if requested
        if interpolation_points is not None:
            #no longer true. sts files have spatial = True but
            #if self.spatial is False:
            #    raise Exception('Triangles and vertex_coordinates must be specified')
            #
            try:
                self.interpolation_points = \
                    interpolation_points = ensure_numeric(interpolation_points)
            except:
                msg = 'Interpolation points must be an N x 2 numeric array ' \
                      'or a list of points\n'
                msg += 'Got: %s.' %(str(self.interpolation_points)[:60] + '...')
                raise Exception(msg)

            # Ensure 'mesh_boundary_polygon' is defined
            mesh_boundary_polygon = None

            if triangles is not None and vertex_coordinates is not None:
                # Check that all interpolation points fall within
                # mesh boundary as defined by triangles and vertex_coordinates.
                from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
                from anuga.geometry.polygon import outside_polygon

                # Create temporary mesh object from mesh info passed
                # into this function.
                mesh = Mesh(vertex_coordinates, triangles)
                mesh_boundary_polygon = mesh.get_boundary_polygon()

                indices = outside_polygon(interpolation_points,
                                          mesh_boundary_polygon)

                # Record result
                #self.mesh_boundary_polygon = mesh_boundary_polygon
                self.indices_outside_mesh = indices

                # Report
                if len(indices) > 0:
                    msg = 'Interpolation points in Interpolation function fall '
                    msg += 'outside specified mesh. Offending points:\n'
                    out_interp_pts = []
                    for i in indices:
                        msg += '%d: %s\n' % (i, interpolation_points[i])
                        out_interp_pts.append(
                                    ensure_numeric(interpolation_points[i]))

                    if verbose is True:
                        import sys
                        from anuga.geometry.polygon import plot_polygons
                        title = ('Interpolation points fall '
                                 'outside specified mesh')
                        plot_polygons([mesh_boundary_polygon,
                                       interpolation_points,
                                       out_interp_pts],
                                      ['line', 'point', 'outside'],
                                      figname='points_boundary_out',
                                      label=title)

                    # Joaquim Luis suggested this as an Exception, so
                    # that the user can now what the problem is rather than
                    # looking for NaN's. However, NANs are handy as they can
                    # be ignored leaving good points for continued processing.
                    if verbose:
                        log.critical(msg)
                    #raise Exception(msg)

            elif triangles is None and vertex_coordinates is not None:    #jj
                #Dealing with sts file
                pass
            else:
                raise Exception('Sww file function requires both triangles and '
                                'vertex_coordinates. sts file file function '
                                'requires the latter.')

            # Plot boundary and interpolation points,
            # but only if if 'mesh_boundary_polygon' has data.
            if verbose is True and mesh_boundary_polygon is not None:
                import sys
                if sys.platform == 'win32':
                    from anuga.geometry.polygon import plot_polygons
                    title = ('Interpolation function: '
                             'Polygon and interpolation points')
                    plot_polygons([mesh_boundary_polygon,
                                   interpolation_points],
                                  ['line', 'point'],
                                  figname='points_boundary',
                                  label=title)

            m = len(self.interpolation_points)
            p = len(self.time)

            for name in quantity_names:
                self.precomputed_values[name] = num.zeros((p, m), float)

            if verbose is True:
                log.critical('Build interpolator')


            # Build interpolator
            if triangles is not None and vertex_coordinates is not None:
                if verbose:
                    msg = 'Building interpolation matrix from source mesh '
                    msg += '(%d vertices, %d triangles)' \
                           % (vertex_coordinates.shape[0],
                              triangles.shape[0])
                    log.critical(msg)

                # This one is no longer needed for STS files
                interpol = Interpolate(vertex_coordinates,
                                       triangles,
                                       verbose=verbose)

            elif triangles is None and vertex_coordinates is not None:
                if verbose:
                    log.critical('Interpolation from STS file')



            if verbose:
                log.critical('Interpolating (%d interpolation points, %d timesteps).'
                             % (self.interpolation_points.shape[0], self.time.shape[0]))

                if time_thinning > 1:
                    log.critical('Timesteps were thinned by a factor of %d'
                                 % time_thinning)
                else:
                    log.critical()

            for i, t in enumerate(self.time):
                # Interpolate quantities at this timestep
                #if verbose and i%((p+10)/10) == 0:
                if verbose:
                    log.critical('  time step %d of %d' % (i, p))

                for name in quantity_names:
                    if len(quantities[name].shape) == 2:
                        Q = quantities[name][i,:] # Quantities at timestep i
                    else:
                        Q = quantities[name][:]   # No time dependency

                    #if verbose and i%((p+10)/10) == 0:
                    if verbose:
                        log.critical('    quantity %s, size=%d' % (name, len(Q)))

                    # Interpolate
                    if triangles is not None and vertex_coordinates is not None:
                        result = interpol.interpolate(Q,
                                                      point_coordinates=\
                                                      self.interpolation_points,
                                                      verbose=False,
                                                      output_centroids=output_centroids)
                        self.centroids = interpol.centroids
                    elif triangles is None and vertex_coordinates is not None:
                        result = interpolate_polyline(Q,
                                                      vertex_coordinates,
                                                      gauge_neighbour_id,
                                                      interpolation_points=\
                                                          self.interpolation_points)

                    #assert len(result), len(interpolation_points)
                    self.precomputed_values[name][i, :] = result

            # Report
            if verbose:
                log.critical(self.statistics())
        else:
            # Store quantitites as is
            for name in quantity_names:
                self.precomputed_values[name] = quantities[name]

#     def __repr__(self):
#         # return 'Interpolation function (spatio-temporal)'
#         return self.statistics()

    def __call__(self, t, point_id=None, x=None, y=None):
        """Evaluate f(t) or f(t, point_id)

        Inputs:
          t:        time - Model time. Must lie within existing timesteps
          point_id: index of one of the preprocessed points.

          If spatial info is present and all of point_id
          are None an exception is raised

          If no spatial info is present, point_id arguments are ignored
          making f a function of time only.

          FIXME: f(t, x, y) x, y could overrided location, point_id ignored
          FIXME: point_id could also be a slice
          FIXME: What if x and y are vectors?
          FIXME: What about f(x,y) without t?
        """

        from math import pi, cos, sin, sqrt

        if self.spatial is True:
            if point_id is None:
                if x is None or y is None:
                    msg = 'Either point_id or x and y must be specified'
                    raise Exception(msg)
            else:
                if self.interpolation_points is None:
                    msg = 'Interpolation_function must be instantiated ' + \
                          'with a list of interpolation points before ' + \
                          'parameter point_id can be used'
                    raise Exception(msg)

        msg = 'Model time %.16f' % t
        msg += ' is not contained in function domain [%.16f:%.16f].\n' % (self.time[0], self.time[-1])
        if t < self.time[0]: raise Modeltime_too_early(msg)
        if t > self.time[-1]: raise Modeltime_too_late(msg)

        oldindex = self.index #Time index

        # Find current time slot
        while t > self.time[self.index]: self.index += 1
        while t < self.time[self.index]: self.index -= 1

        if t == self.time[self.index]:
            # Protect against case where t == T[-1] (last time)
            #  - also works in general when t == T[i]
            ratio = 0
        else:
            # t is now between index and index+1
            ratio = (t - self.time[self.index]) / (self.time[self.index+1] - self.time[self.index])

        # Compute interpolated values
        q = num.zeros(len(self.quantity_names), float)
        for i, name in enumerate(self.quantity_names):
            Q = self.precomputed_values[name]

            if self.spatial is False:
                # If there is no spatial info
                assert len(Q.shape) == 1

                Q0 = Q[self.index]
                if ratio > 0: Q1 = Q[self.index+1]
            else:
                if x is not None and y is not None:
                    # Interpolate to x, y
                    raise Exception('x,y interpolation not yet implemented')
                else:
                    # Use precomputed point
                    Q0 = Q[self.index, point_id]
                    if ratio > 0:
                        Q1 = Q[self.index+1, point_id]

            # Linear temporal interpolation
            if ratio > 0:
                if Q0 == NAN and Q1 == NAN:
                    q[i] = Q0
                else:
                    q[i] = Q0 + ratio*(Q1 - Q0)
            else:
                q[i] = Q0

        # Return vector of interpolated values
        # FIXME:
        if self.spatial is True:
            return q
        else:
            # Replicate q according to x and y
            # This is e.g used for Wind_stress
            if x is None or y is None:
                return q
            else:
                try:
                    N = len(x)
                except:
                    return q
                else:
                    # x is a vector - Create one constant column for each value
                    N = len(x)
                    assert len(y) == N, 'x and y must have same length'
                    res = []
                    for col in q:
                        res.append(col*num.ones(N, float))

                return res

    def get_time(self):
        """Return model time as a vector of timesteps
        """
        return self.time

    def statistics(self):
        """Output statistics about interpolation_function
        """

        vertex_coordinates = self.vertex_coordinates
        interpolation_points = self.interpolation_points
        quantity_names = self.quantity_names
        #quantities = self.quantities
        precomputed_values = self.precomputed_values

        msg =  '------------------------------------------------\n'
        msg += 'Interpolation_function (spatio-temporal) statistics:\n'
        msg += '  Extent:\n'
        if vertex_coordinates is not None:
            x = vertex_coordinates[:,0]
            y = vertex_coordinates[:,1]

            msg += '    x in [%f, %f], len(x) == %d\n'\
                   %(min(x), max(x), len(x))
            msg += '    y in [%f, %f], len(y) == %d\n'\
                   %(min(y), max(y), len(y))



        msg += '    t in [%f, %f], len(t) == %d\n'\
               %(min(self.time), max(self.time), len(self.time))
        msg += '  Quantities:\n'
        for name in quantity_names:
            minq, maxq = self.quantities_range[name]
            msg += '    %s in [%f, %f]\n' %(name, minq, maxq)
            #q = quantities[name][:].flatten()
            #str += '    %s in [%f, %f]\n' %(name, min(q), max(q))

        if interpolation_points is not None:
            msg += '  Interpolation points (xi, eta):'\
                   ' number of points == %d\n' %interpolation_points.shape[0]
            msg += '    xi in [%f, %f]\n' %(min(interpolation_points[:,0]),
                                            max(interpolation_points[:,0]))
            msg += '    eta in [%f, %f]\n' %(min(interpolation_points[:,1]),
                                             max(interpolation_points[:,1]))
            msg += '  Interpolated quantities (over all timesteps):\n'

            for name in quantity_names:
                q = precomputed_values[name][:].flatten()
                msg += '    %s at interpolation points in [%f, %f]\n'\
                       %(name, min(q), max(q))
        msg += '------------------------------------------------\n'

        return msg


def interpolate_sww(sww_file, time, interpolation_points,
                    quantity_names=None, verbose=False):
    """
    obsolete.
    use file_function in utils
    """

    #open sww file
    x, y, volumes, time, quantities = read_sww(sww_file)
    log.critical("x=%s" % str(x))
    log.critical("y=%s" % str(y))

    log.critical("time=%s" % str(time))
    log.critical("quantities=%s" % str(quantities))

    #Add the x and y together
    vertex_coordinates = num.concatenate((x[:,num.newaxis], y[:,num.newaxis]),
                                         axis=1)

    #Will return the quantity values at the specified times and locations
    interp = Interpolation_interface(time,
                                     quantities,
                                     quantity_names=quantity_names,
                                     vertex_coordinates=vertex_coordinates,
                                     triangles=volumes,
                                     interpolation_points=interpolation_points,
                                     verbose=verbose)
