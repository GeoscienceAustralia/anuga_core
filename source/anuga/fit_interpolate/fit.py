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

   TO DO
   * test geo_ref, geo_spatial

   IDEAS
   * (DSG-) Change the interface of fit, so a domain object can
      be passed in. (I don't know if this is feasible). If could
      save time/memory.
"""
import types

from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
from anuga.caching import cache            
from anuga.geospatial_data.geospatial_data import Geospatial_data, \
     ensure_absolute
from anuga.fit_interpolate.general_fit_interpolate import FitInterpolate
from anuga.utilities.sparse import Sparse, Sparse_CSR
from anuga.utilities.polygon import inside_polygon, is_inside_polygon
from anuga.fit_interpolate.search_functions import search_tree_of_vertices

from anuga.utilities.cg_solve import conjugate_gradient
from anuga.utilities.numerical_tools import ensure_numeric, gradient
from anuga.config import default_smoothing_parameter as DEFAULT_ALPHA

import exceptions
class TooFewPointsError(exceptions.Exception): pass
class VertsWithNoTrianglesError(exceptions.Exception): pass

import Numeric as num


class Fit(FitInterpolate):
    
    def __init__(self,
                 vertex_coordinates=None,
                 triangles=None,
                 mesh=None,
                 mesh_origin=None,
                 alpha = None,
                 verbose=False,
                 max_vertices_per_cell=None):


        """
        Fit data at points to the vertices of a mesh.

        Inputs:

          vertex_coordinates: List of coordinate pairs [xi, eta] of
	      points constituting a mesh (or an m x 2 Numeric array or
              a geospatial object)
              Points may appear multiple times
              (e.g. if vertices have discontinuities)

          triangles: List of 3-tuples (or a Numeric array) of
              integers representing indices of all vertices in the mesh.

          mesh_origin: A geo_reference object or 3-tuples consisting of
              UTM zone, easting and northing.
              If specified vertex coordinates are assumed to be
              relative to their respective origins.

          max_vertices_per_cell: Number of vertices in a quad tree cell
          at which the cell is split into 4.

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
            self.alpha = alpha
            
        FitInterpolate.__init__(self,
                 vertex_coordinates,
                 triangles,
                 mesh,
                 mesh_origin,
                 verbose,
                 max_vertices_per_cell)
        
        m = self.mesh.number_of_nodes # Nbr of basis functions (vertices)
        
        self.AtA = None
        self.Atz = None

        self.point_count = 0
        if self.alpha <> 0:
            if verbose: print 'Building smoothing matrix'
            self._build_smoothing_matrix_D()
            
        bd_poly = self.mesh.get_boundary_polygon()    
        self.mesh_boundary_polygon = ensure_numeric(bd_poly)
            
    def _build_coefficient_matrix_B(self,
                                  verbose = False):
        """
        Build final coefficient matrix

        Precon
        If alpha is not zero, matrix D has been built
        Matrix Ata has been built
        """

        if self.alpha <> 0:
            #if verbose: print 'Building smoothing matrix'
            #self._build_smoothing_matrix_D()
            self.B = self.AtA + self.alpha*self.D
        else:
            self.B = self.AtA

        # Convert self.B matrix to CSR format for faster matrix vector
        self.B = Sparse_CSR(self.B)

    def _build_smoothing_matrix_D(self):
        """Build m x m smoothing matrix, where
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
        """
        
        # FIXME: algorithm might be optimised by computing local 9x9
        # "element stiffness matrices:

        m = self.mesh.number_of_nodes # Nbr of basis functions (1/vertex)

        self.D = Sparse(m,m)

        # For each triangle compute contributions to D = D1+D2
        for i in range(len(self.mesh)):

            # Get area
            area = self.mesh.areas[i]

            # Get global vertex indices
            v0 = self.mesh.triangles[i,0]
            v1 = self.mesh.triangles[i,1]
            v2 = self.mesh.triangles[i,2]

            # Get the three vertex_points
            xi0 = self.mesh.get_vertex_coordinate(i, 0)
            xi1 = self.mesh.get_vertex_coordinate(i, 1)
            xi2 = self.mesh.get_vertex_coordinate(i, 2)

            # Compute gradients for each vertex
            a0, b0 = gradient(xi0[0], xi0[1], xi1[0], xi1[1], xi2[0], xi2[1],
                              1, 0, 0)

            a1, b1 = gradient(xi0[0], xi0[1], xi1[0], xi1[1], xi2[0], xi2[1],
                              0, 1, 0)

            a2, b2 = gradient(xi0[0], xi0[1], xi1[0], xi1[1], xi2[0], xi2[1],
                              0, 0, 1)

            # Compute diagonal contributions
            self.D[v0,v0] += (a0*a0 + b0*b0)*area
            self.D[v1,v1] += (a1*a1 + b1*b1)*area
            self.D[v2,v2] += (a2*a2 + b2*b2)*area

            # Compute contributions for basis functions sharing edges
            e01 = (a0*a1 + b0*b1)*area
            self.D[v0,v1] += e01
            self.D[v1,v0] += e01

            e12 = (a1*a2 + b1*b2)*area
            self.D[v1,v2] += e12
            self.D[v2,v1] += e12

            e20 = (a2*a0 + b2*b0)*area
            self.D[v2,v0] += e20
            self.D[v0,v2] += e20

    def get_D(self):
        return self.D.todense()



    def _build_matrix_AtA_Atz(self,
                              point_coordinates,
                              z,
                              verbose = False):
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
        
        # Build n x m interpolation matrix
        if self.AtA == None:
            # AtA and Atz need to be initialised.
            m = self.mesh.number_of_nodes
            if len(z.shape) > 1:
                att_num = z.shape[1]
                self.Atz = num.zeros((m,att_num), num.Float)
            else:
                att_num = 1
                self.Atz = num.zeros((m,), num.Float)
            assert z.shape[0] == point_coordinates.shape[0] 

            AtA = Sparse(m,m)
            # The memory damage has been done by now.
        else:
             AtA = self.AtA #Did this for speed, did ~nothing
        self.point_count += point_coordinates.shape[0]


        inside_indices = inside_polygon(point_coordinates,
                                        self.mesh_boundary_polygon,
                                        closed=True,
                                        verbose=False) # Suppress output
        
        n = len(inside_indices)

        # Compute matrix elements for points inside the mesh
        triangles = self.mesh.triangles # Shorthand
        for d, i in enumerate(inside_indices):
            # For each data_coordinate point
            # if verbose and d%((n+10)/10)==0: print 'Doing %d of %d' %(d, n)
            x = point_coordinates[i]
            
            element_found, sigma0, sigma1, sigma2, k = \
                           search_tree_of_vertices(self.root, self.mesh, x)
            
            if element_found is True:
                j0 = triangles[k,0] #Global vertex id for sigma0
                j1 = triangles[k,1] #Global vertex id for sigma1
                j2 = triangles[k,2] #Global vertex id for sigma2

                sigmas = {j0:sigma0, j1:sigma1, j2:sigma2}
                js     = [j0,j1,j2]

                for j in js:
                    self.Atz[j] +=  sigmas[j]*z[i]
                    #print "self.Atz building", self.Atz
                    #print "self.Atz[j]", self.Atz[j]
                    #print " sigmas[j]", sigmas[j]
                    #print "z[i]",z[i]
                    #print "result", sigmas[j]*z[i]
                    
                    for k in js:
                        AtA[j,k] += sigmas[j]*sigmas[k]
            else:
                flag = is_inside_polygon(x,
                                         self.mesh_boundary_polygon,
                                         closed=True,
                                         verbose=False) # Suppress output
                msg = 'Point (%f, %f) is not inside mesh boundary' % tuple(x)
                assert flag is True, msg                
                
                # FIXME(Ole): This is the message referred to in ticket:314
                minx = min(self.mesh_boundary_polygon[:,0])
                maxx = max(self.mesh_boundary_polygon[:,0])                
                miny = min(self.mesh_boundary_polygon[:,1])
                maxy = max(self.mesh_boundary_polygon[:,1])
                msg = 'Could not find triangle for point %s. ' % str(x) 
                msg += 'Mesh boundary extent is (%.f, %.f), (%.f, %.f)'\
                    % (minx, maxx, miny, maxy)
                raise RuntimeError, msg
                
        self.AtA = AtA

        
    def fit(self, point_coordinates_or_filename=None, z=None,
            verbose=False,
            point_origin=None,
            attribute_name=None,
            max_read_lines=500):
        """Fit a smooth surface to given 1d array of data points z.

        The smooth surface is computed at each vertex in the underlying
        mesh using the formula given in the module doc string.

        Inputs:
        point_coordinates: The co-ordinates of the data points.
              List of coordinate pairs [x, y] of
	      data points or an nx2 Numeric array or a Geospatial_data object
              or points file filename 
          z: Single 1d vector or array of data at the point_coordinates.
          
        """
        # Use blocking to load in the point info
        if type(point_coordinates_or_filename) == types.StringType:
            msg = "Don't set a point origin when reading from a file"
            assert point_origin is None, msg
            filename = point_coordinates_or_filename
            G_data = Geospatial_data(filename,
                                     max_read_lines=max_read_lines,
                                     load_file_now=False,
                                     verbose=verbose)

            for i, geo_block in enumerate(G_data):
                if verbose is True and 0 == i%200: 
                    # The time this will take
                    # is dependant on the # of Triangles
                        
                    print 'Processing Block %d' %i
                    # FIXME (Ole): It would be good to say how many blocks
                    # there are here. But this is no longer necessary
                    # for pts files as they are reported in geospatial_data
                    # I suggest deleting this verbose output and make
                    # Geospatial_data more informative for txt files.
                    #
                    # I still think so (12/12/7, Ole).
            

                    
                # Build the array

                points = geo_block.get_data_points(absolute=True)
                z = geo_block.get_attributes(attribute_name=attribute_name)
                self.build_fit_subset(points, z, verbose=verbose)

                # FIXME(Ole): I thought this test would make sense here
                # See test_fitting_example_that_crashed_2 in test_shallow_water_domain.py
                # Committed 11 March 2009
                #msg = 'Matrix AtA was not built'
                #assert self.AtA is not None, msg
                
                #print 'Matrix was built OK'

                
            point_coordinates = None
        else:
            point_coordinates =  point_coordinates_or_filename
            
        if point_coordinates is None:
            if verbose: print 'Warning: no data points in fit'
            msg = 'No interpolation matrix.'
            assert self.AtA is not None, msg
            assert self.Atz is not None
            
            # FIXME (DSG) - do  a message
        else:
            point_coordinates = ensure_absolute(point_coordinates,
                                                geo_reference=point_origin)
            # if isinstance(point_coordinates,Geospatial_data) and z is None:
            # z will come from the geo-ref
            self.build_fit_subset(point_coordinates, z, verbose)

        # Check sanity
        m = self.mesh.number_of_nodes # Nbr of basis functions (1/vertex)
        n = self.point_count
        if n<m and self.alpha == 0.0:
            msg = 'ERROR (least_squares): Too few data points\n'
            msg += 'There are only %d data points and alpha == 0. ' %n
	    msg += 'Need at least %d\n' %m
            msg += 'Alternatively, set smoothing parameter alpha to a small '
	    msg += 'positive value,\ne.g. 1.0e-3.'
            raise TooFewPointsError(msg)

        self._build_coefficient_matrix_B(verbose)
        loners = self.mesh.get_lone_vertices()
        # FIXME  - make this as error message.
        # test with
        # Not_yet_test_smooth_att_to_mesh_with_excess_verts.
        if len(loners)>0:
            msg = 'WARNING: (least_squares): \nVertices with no triangles\n'
            msg += 'All vertices should be part of a triangle.\n'
            msg += 'In the future this will be inforced.\n'
	    msg += 'The following vertices are not part of a triangle;\n'
            msg += str(loners)
            print msg
            #raise VertsWithNoTrianglesError(msg)
        
        
        return conjugate_gradient(self.B, self.Atz, self.Atz,
                                  imax=2*len(self.Atz) )

        
    def build_fit_subset(self, point_coordinates, z=None, attribute_name=None,
                              verbose=False):
        """Fit a smooth surface to given 1d array of data points z.

        The smooth surface is computed at each vertex in the underlying
        mesh using the formula given in the module doc string.

        Inputs:
        point_coordinates: The co-ordinates of the data points.
              List of coordinate pairs [x, y] of
	      data points or an nx2 Numeric array or a Geospatial_data object
        z: Single 1d vector or array of data at the point_coordinates.
        attribute_name: Used to get the z values from the
              geospatial object if no attribute_name is specified,
              it's a bit of a lucky dip as to what attributes you get.
              If there is only one attribute it will be that one.

        """

        # FIXME(DSG-DSG): Check that the vert and point coords
        # have the same zone.
        if isinstance(point_coordinates,Geospatial_data):
            point_coordinates = point_coordinates.get_data_points( \
                absolute = True)
        
        # Convert input to Numeric arrays
        if z is not None:
            z = ensure_numeric(z, num.Float)
        else:
            msg = 'z not specified'
            assert isinstance(point_coordinates,Geospatial_data), msg
            z = point_coordinates.get_attributes(attribute_name)

        point_coordinates = ensure_numeric(point_coordinates, num.Float)
        self._build_matrix_AtA_Atz(point_coordinates, z, verbose)
        

############################################################################

def fit_to_mesh(point_coordinates, # this can also be a points file name
                vertex_coordinates=None,
                triangles=None,
                mesh=None,
                point_attributes=None,
                alpha=DEFAULT_ALPHA,
                verbose=False,
                acceptable_overshoot=1.01, 
                # FIXME: Move to config - this value is assumed in caching test
                # FIXME(Ole): Just realised that this was never implemented (29 Jan 2009). I suggest removing it altogether.
                mesh_origin=None,
                data_origin=None,
                max_read_lines=None,
                attribute_name=None,
                use_cache=False):
    """Wrapper around internal function _fit_to_mesh for use with caching.
    
    """
    
    args = (point_coordinates, )
    kwargs = {'vertex_coordinates': vertex_coordinates,
              'triangles': triangles,
              'mesh': mesh,
              'point_attributes': point_attributes,
              'alpha': alpha,
              'verbose': verbose,
              'acceptable_overshoot': acceptable_overshoot,
              'mesh_origin': mesh_origin,
              'data_origin': data_origin,
              'max_read_lines': max_read_lines,
              'attribute_name': attribute_name 
              }

    if use_cache is True:
        if isinstance(point_coordinates, basestring):
            # We assume that point_coordinates is the name of a .csv/.txt
            # file which must be passed onto caching as a dependency 
            # (in case it has changed on disk)
            dep = [point_coordinates]
        else:
            dep = None

            
        #from caching import myhash
        #import copy
        #print args
        #print kwargs
        #print 'hashing:'
        #print 'args', myhash( (args, kwargs) )
        #print 'again', myhash( copy.deepcopy( (args, kwargs)) )        
        
        #print 'mesh hash', myhash( kwargs['mesh'] )        
        
        #print '-------------------------'
        #print 'vertices hash', myhash( kwargs['mesh'].nodes )
        #print 'triangles hash', myhash( kwargs['mesh'].triangles )
        #print '-------------------------'        
        
        #for key in mesh.__dict__:
        #    print key, myhash(mesh.__dict__[key])
        
        #for key in mesh.quantities.keys():
        #    print key, myhash(mesh.quantities[key])
        
        #import sys; sys.exit() 
            
        return cache(_fit_to_mesh,
                     args, kwargs,
                     verbose=verbose,
                     compression=False,
                     dependencies=dep)
    else:
        return apply(_fit_to_mesh,
                     args, kwargs)

def _fit_to_mesh(point_coordinates, # this can also be a points file name
                 vertex_coordinates=None,
                 triangles=None,
                 mesh=None,
                 point_attributes=None,
                 alpha=DEFAULT_ALPHA,
                 verbose=False,
                 acceptable_overshoot=1.01,
                 mesh_origin=None,
                 data_origin=None,
                 max_read_lines=None,
                 attribute_name=None):
    """
    Fit a smooth surface to a triangulation,
    given data points with attributes.


        Inputs:
        vertex_coordinates: List of coordinate pairs [xi, eta] of
	      points constituting a mesh (or an m x 2 Numeric array or
              a geospatial object)
              Points may appear multiple times
              (e.g. if vertices have discontinuities)

          triangles: List of 3-tuples (or a Numeric array) of
          integers representing indices of all vertices in the mesh.

          point_coordinates: List of coordinate pairs [x, y] of data points
          (or an nx2 Numeric array). This can also be a .csv/.txt/.pts
          file name.

          alpha: Smoothing parameter.

          acceptable overshoot: NOT IMPLEMENTED
          controls the allowed factor by which
          fitted values
          may exceed the value of input data. The lower limit is defined
          as min(z) - acceptable_overshoot*delta z and upper limit
          as max(z) + acceptable_overshoot*delta z
          

          mesh_origin: A geo_reference object or 3-tuples consisting of
              UTM zone, easting and northing.
              If specified vertex coordinates are assumed to be
              relative to their respective origins.
          

          point_attributes: Vector or array of data at the
                            point_coordinates.

    """

    if mesh is None:
        # FIXME(DSG): Throw errors if triangles or vertex_coordinates
        # are None
            
        #Convert input to Numeric arrays
        triangles = ensure_numeric(triangles, num.Int)
        vertex_coordinates = ensure_absolute(vertex_coordinates,
                                             geo_reference = mesh_origin)

        if verbose: print 'FitInterpolate: Building mesh'        
        mesh = Mesh(vertex_coordinates, triangles)
        mesh.check_integrity()
    
    
    interp = Fit(mesh=mesh,
                 verbose=verbose,
                 alpha=alpha)

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


#def _fit(*args, **kwargs):
#    """Private function for use with caching. Reason is that classes
#    may change their byte code between runs which is annoying.
#    """
#    
#    return Fit(*args, **kwargs)


def fit_to_mesh_file(mesh_file, point_file, mesh_output_file,
                     alpha=DEFAULT_ALPHA, verbose= False,
                     expand_search = False,
                     precrop = False,
                     display_errors = True):
    """
    Given a mesh file (tsh) and a point attribute file, fit
    point attributes to the mesh and write a mesh file with the
    results.

    Note: the points file needs titles.  If you want anuga to use the tsh file,
    make sure the title is elevation.

    NOTE: Throws IOErrors, for a variety of file problems.
    
    """

    from load_mesh.loadASCII import import_mesh_file, \
         export_mesh_file, concatinate_attributelist


    try:
        mesh_dict = import_mesh_file(mesh_file)
    except IOError,e:
        if display_errors:
            print "Could not load bad file. ", e
        raise IOError  #Could not load bad mesh file.
    
    vertex_coordinates = mesh_dict['vertices']
    triangles = mesh_dict['triangles']
    if type(mesh_dict['vertex_attributes']) == num.ArrayType:
        old_point_attributes = mesh_dict['vertex_attributes'].tolist()
    else:
        old_point_attributes = mesh_dict['vertex_attributes']

    if type(mesh_dict['vertex_attribute_titles']) == num.ArrayType:
        old_title_list = mesh_dict['vertex_attribute_titles'].tolist()
    else:
        old_title_list = mesh_dict['vertex_attribute_titles']

    if verbose: print 'tsh file %s loaded' %mesh_file

    # load in the points file
    try:
        geo = Geospatial_data(point_file, verbose=verbose)
    except IOError,e:
        if display_errors:
            print "Could not load bad file. ", e
        raise IOError  #Re-raise exception  

    point_coordinates = geo.get_data_points(absolute=True)
    title_list,point_attributes = concatinate_attributelist( \
        geo.get_all_attributes())

    if mesh_dict.has_key('geo_reference') and \
           not mesh_dict['geo_reference'] is None:
        mesh_origin = mesh_dict['geo_reference'].get_origin()
    else:
        mesh_origin = None

    if verbose: print "points file loaded"
    if verbose: print "fitting to mesh"
    f = fit_to_mesh(point_coordinates,
                    vertex_coordinates,
                    triangles,
                    None,
                    point_attributes,
                    alpha = alpha,
                    verbose = verbose,
                    data_origin = None,
                    mesh_origin = mesh_origin)
    if verbose: print "finished fitting to mesh"

    # convert array to list of lists
    new_point_attributes = f.tolist()
    #FIXME have this overwrite attributes with the same title - DSG
    #Put the newer attributes last
    if old_title_list <> []:
        old_title_list.extend(title_list)
        #FIXME can this be done a faster way? - DSG
        for i in range(len(old_point_attributes)):
            old_point_attributes[i].extend(new_point_attributes[i])
        mesh_dict['vertex_attributes'] = old_point_attributes
        mesh_dict['vertex_attribute_titles'] = old_title_list
    else:
        mesh_dict['vertex_attributes'] = new_point_attributes
        mesh_dict['vertex_attribute_titles'] = title_list

    if verbose: print "exporting to file ", mesh_output_file

    try:
        export_mesh_file(mesh_output_file, mesh_dict)
    except IOError,e:
        if display_errors:
            print "Could not write file. ", e
        raise IOError
