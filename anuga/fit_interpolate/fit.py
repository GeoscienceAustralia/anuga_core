"""Least squares fitting.

   Implements a penalised least-squares fit.
   putting point data onto the mesh.

   The penalty term (or smoothing term) is controlled by the smoothing
   parameter alpha.
   With a value of alpha=0, the fit function will attempt
   to interpolate as closely as possible in the least-squares sense.
   With values alpha > 0, a certain amount of smoothing will be applied.
   A positive alpha is essential in cases where there are too few
   data points.
   A negative alpha is not allowed.
   A typical value of alpha is 1.0e-6

   Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
   Geoscience Australia, 2004.
"""
import numpy as num
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from . import fitsmooth_ext as fitsmooth
import sys


from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
from anuga.caching import cache
from anuga.geospatial_data.geospatial_data import Geospatial_data, \
    ensure_absolute
from anuga.fit_interpolate.general_fit_interpolate import FitInterpolate

from anuga.utilities.sparse import Sparse_CSR
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.cg_solve import conjugate_gradient
from anuga.config import default_smoothing_parameter as DEFAULT_ALPHA
import anuga.utilities.log as log


class TooFewPointsError(Exception):
    pass


class VertsWithNoTrianglesError(Exception):
    pass


# ----------------------------------------------
# C code to build interpolation matrices
# ----------------------------------------------


class Fit(FitInterpolate):

    def __init__(self,
                 vertex_coordinates=None,
                 triangles=None,
                 mesh=None,
                 mesh_origin=None,
                 alpha=None,
                 verbose=False,
                 cg_precon='Jacobi',
                 use_c_cg=True):
        """
        Padarn Note 05/12/12: This documentation should probably
        be updated to account for the fact that the fitting is now
        done in C. I wasn't sure what details were necessary though.

        Fit data at points to the vertices of a mesh.

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

          Note: Don't supply a vertex coords as a geospatial object and
              a mesh origin, since geospatial has its own mesh origin.


        Usage,
        To use this in a blocking way, call  build_fit_subset, with z info,
        and then fit, with no point coord, z info.
        """
        # Initialise variabels
        if alpha is None:
            self.alpha = DEFAULT_ALPHA
        else:
            self.alpha = alpha  # may be a float or the string 'auto'

        FitInterpolate.__init__(self,
                                vertex_coordinates,
                                triangles,
                                mesh,
                                mesh_origin=mesh_origin,
                                verbose=verbose)

        self.AtA = None
        self.Atz = None
        self.D = None
        self.point_count = 0
        self.z_sq = None  # accumulated ||z||² for L-curve alpha selection

        # The smoothing matrix D is always built, even when alpha=0, because
        # the underlying C solver (fit_ext) expects the D array to be present.
        # Skipping construction when alpha=0 would require changes to the C
        # interface and is deferred.
        if verbose:
            log.info('Building smoothing matrix')
        self.D = self._build_smoothing_matrix_D()

        bd_poly = self.mesh.get_boundary_polygon()
        self.mesh_boundary_polygon = ensure_numeric(bd_poly)

        self.cg_precon = cg_precon
        self.use_c_cg = use_c_cg

    def _build_coefficient_matrix_B(self,
                                    verbose=False):
        """
        Build final coefficient matrix from AtA and D
        """

        msize = self.mesh.number_of_nodes

        self.B = fitsmooth.build_matrix_B(self.D,
                                          self.AtA, self.alpha)

        # Convert self.B matrix to CSR format
        self.B = Sparse_CSR(data=num.array(self.B[0]),
                            Colind=num.array(self.B[1]),
                            rowptr=num.array(self.B[2]),
                            m=msize, n=msize)
        # NOTE PADARN: The above step could potentially be removed
        # and the sparse matrix worked with directly in C. Not sure
        # if this would be worthwhile.

    def _build_smoothing_matrix_D(self):
        r"""Build m x m smoothing matrix, where
        m is the number of basis functions phi_k (one per vertex)

        The smoothing matrix is defined as

        D = D1 + D2

        where

        [D1]_{k,l} = \int_\Omega
           \frac{\partial \phi_k}{\partial x}
           \frac{\partial \phi_l}{\partial x}\,
           dx dy

        [D2]_{k,l} = \int_\Omega
           \frac{\partial \phi_k}{\partial y}
           \frac{\partial \phi_l}{\partial y}\,
           dx dy


        The derivatives \frac{\partial \phi_k}{\partial x},
        \frac{\partial \phi_k}{\partial x} for a particular triangle
        are obtained by computing the gradient a_k, b_k for basis function k

        NOTE PADARN: All of this is now done in an external C function, and the
        result is stored in a Capsule object, meaning the entries cannot be directly
        accessed.
        """

        # NOTE PADARN: Should the input arguments here be checked - making
        # sure that they are floats? Not sure if this is done elsewhere.
        # NOTE PADARN: Should global coordinates be used for the smoothing
        # matrix, or is thids not important?
        return fitsmooth.build_smoothing_matrix(self.mesh.triangles,
                                                self.mesh.areas, self.mesh.vertex_coordinates)

    # NOTE PADARN: This function was added to emulate behavior of the original
    # class not using external C functions. This method is dangerous as D could
    # be very large - it was added as it is used in a unit test.

    def get_D(self):
        return fitsmooth.return_full_D(self.D, self.mesh.number_of_nodes)

    # NOTE PADARN: This function was added to emulate behavior of the original
    # class so as to pass a unit test. It is completely unneeded.
    def build_fit_subset(self, point_coordinates, z=None, attribute_name=None,
                         verbose=False, output='dot'):
        self._build_matrix_AtA_Atz(
            point_coordinates, z, attribute_name, verbose, output)

    def _build_matrix_AtA_Atz(self, point_coordinates, z=None, attribute_name=None,
                              verbose=False, output='dot'):
        """Build:
        AtA  m x m  interpolation matrix, and,
        Atz  m x a  interpolation matrix where,
        m is the number of basis functions phi_k (one per vertex)
        a is the number of data attributes

        This algorithm uses a quad tree data structure for fast binning of
        data points.

        If Ata is None, the matrices AtA and Atz are created.

        This function can be called again and again, with sub-sets of
        the point coordinates.  Call fit to get the results.

        Preconditions
        z and points are numeric
        Point_coordindates and mesh vertices have the same origin.

        The number of attributes of the data points does not change
        """

        if isinstance(point_coordinates, Geospatial_data):
            point_coordinates = point_coordinates.get_data_points(
                absolute=True)

        # Convert input to numeric arrays
        if z is not None:
            z = ensure_numeric(z, float)
        else:
            msg = 'z not specified'
            assert isinstance(point_coordinates, Geospatial_data), msg
            z = point_coordinates.get_attributes(attribute_name)

        point_coordinates = ensure_numeric(point_coordinates, float)

        npts = len(z)
        z = num.array(z)
        # NOTE PADARN : This copy might be needed to
        # make sure memory is contig - would be better to read in C..
        z = z.copy()

        self.point_count += z.shape[0]

        # Accumulate ||z||² per attribute for L-curve alpha selection
        z_sq_block = num.sum(z**2, axis=0) if z.ndim > 1 else float(z @ z)
        self.z_sq = z_sq_block if self.z_sq is None else self.z_sq + z_sq_block

        zdim = 1
        if len(z.shape) != 1:
            zdim = z.shape[1]

        [AtA, Atz] = fitsmooth.build_matrix_AtA_Atz_points(self.root.root,
                                                           self.mesh.number_of_nodes,
                                                           self.mesh.triangles,
                                                           num.array(point_coordinates), z, zdim, npts)

        if verbose and output == 'dot':
            print('\b.', end=' ')
            sys.stdout.flush()
        if zdim == 1:
            Atz = num.array(Atz[0])
        else:
            Atz = num.array(Atz).transpose()

        if self.AtA is None and self.Atz is None:
            self.AtA = AtA
            self.Atz = Atz
        else:
            fitsmooth.combine_partial_AtA_Atz(self.AtA, AtA,
                                              self.Atz, Atz, zdim, self.mesh.number_of_nodes)

    def select_alpha(self, return_curve=False):
        """Select smoothing parameter alpha using the L-curve criterion.

        Tries alpha on a log scale (1e-6 … 1) and finds the corner of the
        L-curve (log residual norm vs log smoothing norm).  The corner
        corresponds to the best balance between data fit and smoothness.

        Must be called after _build_matrix_AtA_Atz (i.e. after data has been
        accumulated).

        Parameters
        ----------
        return_curve : bool
            If True, also return ``(log10_rss, log10_s, curvature)`` arrays.

        Returns
        -------
        float
            Optimal alpha value.
        tuple (optional)
            ``(log10_rss, log10_s, curvature)`` if return_curve is True.
        """
        assert self.AtA is not None, 'Call _build_matrix_AtA_Atz before select_alpha'
        assert self.z_sq is not None, 'z_sq not accumulated; re-run with current code'

        m = self.mesh.number_of_nodes
        # 20 candidates spanning 1e-6 … 100: wide enough to capture the L-corner
        # for both underdetermined (few data) and overdetermined (many data) cases.
        alphas = num.logspace(-6, 2, 20)

        # Extract AtA and D as scipy sparse matrices (non-destructive: dok_to_csr
        # only sorts the hash table, it does not modify values).
        #
        # Two row_ptr repairs are needed:
        # (a) sparse_dok.num_rows tracks only the highest row index seen, so
        #     row_ptr may be shorter than m+1 if high-indexed nodes have no data.
        #     Extend it with the sentinel (num_entries).
        # (b) convert_to_csr_ptr uses malloc (not calloc), so rows before the
        #     first populated row are uninitialized garbage.  Since entries are
        #     sorted, the first populated row always has rp[i_first]=0.  A
        #     backward min-scan propagates 0 back over any garbage entries.
        def _to_scipy_csr(cap):
            raw = fitsmooth.dok_to_csr(cap)
            data_   = num.array(raw[0], dtype=float)
            colind_ = num.array(raw[1], dtype=num.int64)
            rp      = num.array(raw[2], dtype=num.int64)
            n_ent   = len(data_)
            # (a) extend to m+1 if needed
            if len(rp) < m + 1:
                rp = num.concatenate(
                    [rp, num.full(m + 1 - len(rp), n_ent, dtype=num.int64)])
            # (b) repair garbage in leading uninitialized entries:
            #     clip to [0, n_ent], then backward scan to enforce monotonicity
            rp = num.clip(rp, 0, n_ent)
            for k in range(len(rp) - 2, -1, -1):
                if rp[k] > rp[k + 1]:
                    rp[k] = rp[k + 1]
            return sp.csr_matrix((data_, colind_, rp), shape=(m, m))

        AtA_sp = _to_scipy_csr(self.AtA)
        D_sp   = _to_scipy_csr(self.D)

        # Use first attribute only for alpha selection
        Atz = self.Atz
        if Atz.ndim > 1:
            Atz_1d = num.ascontiguousarray(Atz[:, 0])
            z_sq = float(self.z_sq[0])
        else:
            Atz_1d = Atz
            z_sq = float(self.z_sq)

        log_rss = num.zeros(len(alphas))
        log_s   = num.zeros(len(alphas))

        for i, alpha in enumerate(alphas):
            B = AtA_sp + float(alpha) * D_sp
            f = spla.spsolve(B, Atz_1d)

            # Numerically stable RSS: since (AtA + alpha*D)*f = Atz we have
            # f^T AtA f = f^T Atz - alpha*f_D_f, so
            # RSS = z_sq - 2*f·Atz + f^T AtA f = z_sq - f·Atz - alpha*f_D_f
            # This avoids catastrophic cancellation at small alpha.
            Atz_dot_f = float(Atz_1d @ f)
            f_D_f     = float(f @ (D_sp @ f))

            rss = z_sq - Atz_dot_f - float(alpha) * f_D_f
            rss = max(rss, 1e-300)
            s   = max(f_D_f, 1e-300)

            log_rss[i] = num.log10(rss)
            log_s[i]   = num.log10(s)

        # L-curve corner: maximum curvature of the log-log curve
        x, y = log_rss, log_s
        dx  = num.gradient(x)
        ddx = num.gradient(dx)
        dy  = num.gradient(y)
        ddy = num.gradient(dy)
        denom  = num.where((dx**2 + dy**2)**1.5 < 1e-300, 1e-300,
                           (dx**2 + dy**2)**1.5)
        kappa  = (ddx * dy - dx * ddy) / denom
        best   = int(num.argmax(kappa))

        # Degenerate cases: no clear interior corner.
        # (a) All-negative kappa: nearly vertical L-curve (overdetermined system)
        # (b) Corner lands on last candidate: L-curve hasn't turned within this range
        # In both cases fall back to DEFAULT_ALPHA which is a safe, conservative value.
        if kappa[best] <= 0 or best == len(alphas) - 1:
            opt_alpha = DEFAULT_ALPHA
        else:
            opt_alpha = float(alphas[best])

        if return_curve:
            return opt_alpha, (log_rss, log_s, kappa)
        return opt_alpha

    def fit(self, point_coordinates_or_filename=None, z=None,
            verbose=False,
            point_origin=None,
            attribute_name=None,
            max_read_lines=1e7):
        """Fit a smooth surface to given 1d array of data points z.

        The smooth surface is computed at each vertex in the underlying
        mesh using the formula given in the module doc string.

        Inputs:
        point_coordinates_or_filename: The co-ordinates of the data points.
              A filename of a .pts file or a
              List of coordinate pairs [x, y] of
              data points or an nx2 numeric array or a Geospatial_data object
              or points file filename
          z: Single 1d vector or array of data at the point_coordinates.

        """
        if isinstance(point_coordinates_or_filename, str):
            if point_coordinates_or_filename[-4:] != ".pts":
                use_blocking_option2 = False

        # NOTE PADARN 29/03/13: File reading from C has been removed. Now
        # the input is either a set of points, or a filename which is then
        # handled by the Geospatial_data object

        if verbose:
            print('Fit.fit: Initializing')

        # Use blocking to load in the point info
        if isinstance(point_coordinates_or_filename, str):
            msg = "Don't set a point origin when reading from a file"
            assert point_origin is None, msg
            filename = point_coordinates_or_filename

            G_data = Geospatial_data(filename,
                                     max_read_lines=max_read_lines,
                                     load_file_now=False,
                                     verbose=verbose)

            for i, geo_block in enumerate(G_data):

               # Build the array
                points = geo_block.get_data_points(absolute=True)
                z = geo_block.get_attributes(attribute_name=attribute_name)

                self._build_matrix_AtA_Atz(points, z, attribute_name, verbose)

            point_coordinates = None

            if verbose:
                print('')
        else:
            point_coordinates = point_coordinates_or_filename

        # This condition either means a filename was read or the function
        # recieved a None as input
        if point_coordinates is None:
            if verbose:
                log.warning('Fit.fit: Warning: no data points in fit')
            msg = 'No interpolation matrix.'
            assert self.AtA is not None, msg
            assert self.Atz is not None

        else:
            point_coordinates = ensure_absolute(point_coordinates,
                                                geo_reference=point_origin)
            # if isinstance(point_coordinates,Geospatial_data) and z is None:
            # z will come from the geo-ref

            self._build_matrix_AtA_Atz(
                point_coordinates, z, verbose=verbose, output='counter')

        # Check sanity
        m = self.mesh.number_of_nodes  # Nbr of basis functions (1/vertex)
        n = self.point_count
        if n < m and self.alpha == 0.0:
            msg = 'ERROR (least_squares): Too few data points\n'
            msg += 'There are only %d data points and alpha == 0. ' % n
            msg += 'Need at least %d\n' % m
            msg += 'Alternatively, set smoothing parameter alpha to a small '
            msg += 'positive value,\ne.g. 1.0e-3.'
            raise TooFewPointsError(msg)

        if self.alpha == 'auto':
            self.alpha = self.select_alpha()
            if verbose:
                log.info('L-curve selected alpha = %.2e' % self.alpha)

        self._build_coefficient_matrix_B(verbose)
        loners = self.mesh.get_lone_vertices()
        if len(loners) > 0:
            msg = 'WARNING: (least_squares): \nVertices with no triangles\n'
            msg += 'All vertices should be part of a triangle.\n'
            msg += 'In the future this will be inforced.\n'
            msg += 'The following vertices are not part of a triangle;\n'
            msg += str(loners)
            log.info(msg)

            #raise VertsWithNoTrianglesError(msg)
        return conjugate_gradient(self.B, self.Atz, self.Atz,
                                  imax=2 * len(self.Atz)+1000, use_c_cg=self.use_c_cg,
                                  precon=self.cg_precon)


# poin_coordiantes can also be a points file name

def fit_to_mesh(point_coordinates,
                vertex_coordinates=None,
                triangles=None,
                mesh=None,
                point_attributes=None,
                alpha=DEFAULT_ALPHA,
                verbose=False,
                mesh_origin=None,
                data_origin=None,
                max_read_lines=None,
                attribute_name=None,
                use_cache=False,
                cg_precon='Jacobi',
                use_c_cg=True):
    """Wrapper around internal function _fit_to_mesh for use with caching.
    """

    args = (point_coordinates, )
    kwargs = {'vertex_coordinates': vertex_coordinates,
              'triangles': triangles,
              'mesh': mesh,
              'point_attributes': point_attributes,
              'alpha': alpha,
              'verbose': verbose,
              'mesh_origin': mesh_origin,
              'data_origin': data_origin,
              'max_read_lines': max_read_lines,
              'attribute_name': attribute_name,
              'cg_precon': cg_precon,
              'use_c_cg': use_c_cg
              }

    if use_cache is True:
        if isinstance(point_coordinates, str):
            # We assume that point_coordinates is the name of a .csv/.txt
            # file which must be passed onto caching as a dependency
            # (in case it has changed on disk)
            dep = [point_coordinates]
        else:
            dep = None

        return cache(_fit_to_mesh,
                     args, kwargs,
                     verbose=verbose,
                     compression=False,
                     dependencies=dep)
    else:
        return _fit_to_mesh(*args, **kwargs)


# point_coordinates can also be a points file name

def _fit_to_mesh(point_coordinates,
                 vertex_coordinates=None,
                 triangles=None,
                 mesh=None,
                 point_attributes=None,
                 alpha=DEFAULT_ALPHA,
                 verbose=False,
                 mesh_origin=None,
                 data_origin=None,
                 max_read_lines=None,
                 attribute_name=None,
                 cg_precon='Jacobi',
                 use_c_cg=True):
    """
    Fit a smooth surface to a triangulation,
    given data points with attributes.


        Inputs:
        vertex_coordinates: List of coordinate pairs [xi, eta] of
        points constituting a mesh (or an m x 2 numeric array or
              a geospatial object)
              Points may appear multiple times
              (e.g. if vertices have discontinuities)

          triangles: List of 3-tuples (or a numeric array) of
          integers representing indices of all vertices in the mesh.

          point_coordinates: List of coordinate pairs [x, y] of data points
          (or an nx2 numeric array). This can also be a .csv/.txt/.pts
          file name.

          alpha: Smoothing parameter.

          mesh_origin: A geo_reference object or 3-tuples consisting of
              UTM zone, easting and northing.
              If specified vertex coordinates are assumed to be
              relative to their respective origins.

          point_attributes: Vector or array of data at the
                            point_coordinates.

    """

    if mesh is None:
        # Convert input to numeric arrays
        triangles = ensure_numeric(triangles, int)
        vertex_coordinates = ensure_absolute(vertex_coordinates,
                                             geo_reference=mesh_origin)

        if verbose:
            log.info('_fit_to_mesh: Building mesh')
        mesh = Mesh(vertex_coordinates, triangles)

        # Don't need this as we have just created the mesh
        # mesh.check_integrity()

    interp = Fit(mesh=mesh,
                 verbose=verbose,
                 alpha=alpha,
                 cg_precon=cg_precon,
                 use_c_cg=use_c_cg)

    vertex_attributes = interp.fit(point_coordinates,
                                   point_attributes,
                                   point_origin=data_origin,
                                   max_read_lines=max_read_lines,
                                   attribute_name=attribute_name,
                                   verbose=verbose)

    # Add the value checking stuff that's in least squares.
    # Maybe this stuff should get pushed down into Fit.
    # at least be a method of Fit.
    # Or intigrate it into the fit method, saving teh max and min's
    # as att's.

    return vertex_attributes


def fit_to_mesh_file(mesh_file, point_file, mesh_output_file,
                     alpha=DEFAULT_ALPHA, verbose=False,
                     expand_search=False,
                     precrop=False,
                     display_errors=True):
    """
    Given a mesh file (tsh) and a point attribute file, fit
    point attributes to the mesh and write a mesh file with the
    results.

    Note: the points file needs titles.  If you want anuga to use the tsh file,
    make sure the title is elevation.

    NOTE: Throws IOErrors, for a variety of file problems.

    """

    from anuga.load_mesh.loadASCII import import_mesh_file, \
        export_mesh_file, concatinate_attributelist

    try:
        mesh_dict = import_mesh_file(mesh_file)
    except OSError as e:
        if display_errors:
            log.warning("Could not load bad file: %s" % str(e))
        raise OSError  # Could not load bad mesh file.

    vertex_coordinates = mesh_dict['vertices']
    triangles = mesh_dict['triangles']
    if isinstance(mesh_dict['vertex_attributes'], num.ndarray):
        old_point_attributes = mesh_dict['vertex_attributes'].tolist()
    else:
        old_point_attributes = mesh_dict['vertex_attributes']

    if isinstance(mesh_dict['vertex_attribute_titles'], num.ndarray):
        old_title_list = mesh_dict['vertex_attribute_titles'].tolist()
    else:
        old_title_list = mesh_dict['vertex_attribute_titles']

    if verbose:
        log.info('tsh file %s loaded' % mesh_file)

    # load in the points file
    try:
        geo = Geospatial_data(point_file, verbose=verbose)
    except OSError as e:
        if display_errors:
            log.warning("Could not load bad file: %s" % str(e))
        raise OSError  # Re-raise exception

    point_coordinates = geo.get_data_points(absolute=True)
    title_list, point_attributes = concatinate_attributelist(
        geo.get_all_attributes())

    if 'geo_reference' in mesh_dict and \
            mesh_dict['geo_reference'] is not None:
        mesh_origin = mesh_dict['geo_reference'].get_origin()
    else:
        mesh_origin = None

    if verbose:
        log.info("points file loaded")
    if verbose:
        log.info("fitting to mesh")
    f = fit_to_mesh(point_coordinates,
                    vertex_coordinates,
                    triangles,
                    None,
                    point_attributes,
                    alpha=alpha,
                    verbose=verbose,
                    data_origin=None,
                    mesh_origin=mesh_origin)
    if verbose:
        log.info("finished fitting to mesh")

    # Put the newer attributes last
    if old_title_list != []:
        old_title_list.extend(title_list)
        combined = num.hstack([num.array(old_point_attributes), f])
        mesh_dict['vertex_attributes'] = combined.tolist()
        mesh_dict['vertex_attribute_titles'] = old_title_list
    else:
        mesh_dict['vertex_attributes'] = f.tolist()
        mesh_dict['vertex_attribute_titles'] = title_list

    if verbose:
        log.info("exporting to file %s" % mesh_output_file)

    try:
        export_mesh_file(mesh_output_file, mesh_dict)
    except OSError as e:
        if display_errors:
            log.warning("Could not write file %s" % str(e))
        raise OSError
