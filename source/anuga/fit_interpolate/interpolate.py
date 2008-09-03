"""Least squares interpolation.

   Implements a least-squares interpolation.
   Putting mesh data onto points.

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
from warnings import warn
from math import sqrt
from csv import writer, DictWriter

from Numeric import zeros, array, Float, Int, dot, transpose, concatenate, \
     ArrayType, allclose, take, NewAxis, arange

from anuga.caching.caching import cache
from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
from anuga.utilities.sparse import Sparse, Sparse_CSR
from anuga.utilities.cg_solve import conjugate_gradient, VectorShapeError
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.utilities.numerical_tools import ensure_numeric, mean, NAN
from anuga.utilities.polygon import in_and_outside_polygon
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.geospatial_data.geospatial_data import ensure_absolute
from anuga.fit_interpolate.search_functions import search_tree_of_vertices
from anuga.fit_interpolate.general_fit_interpolate import FitInterpolate
from anuga.abstract_2d_finite_volumes.util import file_function

# Interpolation specific exceptions

class Modeltime_too_late(Exception): pass
class Modeltime_too_early(Exception): pass



class Interpolate (FitInterpolate):
        
    def __init__(self,
                 vertex_coordinates,
                 triangles,
                 mesh_origin=None,
                 verbose=False,
                 max_vertices_per_cell=None):


        """ Build interpolation matrix mapping from
        function values at vertices to function values at data points

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
        """

        # FIXME (Ole): Need an input check
        
        # Initialise variabels
        self._A_can_be_reused = False
        self._point_coordinates = None
        
        FitInterpolate.__init__(self,
                                vertex_coordinates=vertex_coordinates,
                                triangles=triangles,
                                mesh_origin=mesh_origin,
                                verbose=verbose,
                                max_vertices_per_cell=max_vertices_per_cell)

    def interpolate_polyline(self,
                             f,
                             vertex_coordinates,
                             gauge_neighbour_id,
                             point_coordinates=None,
                             verbose=False):
        """Interpolate linearly between values f on nodes (vertex coordinates) of a polyline to midpoints of triangles
        of boundary.

        f is the data on the polyline nodes.

        The mesh values representing a smooth surface are
        assumed to be specified in f.

        Inputs:
          f: Vector or array of data at the polyline nodes.
              If f is an array, interpolation will be done for each column as
              per underlying matrix-matrix multiplication
          point_coordinates: Interpolate polyline data to these positions.
              List of coordinate pairs [x, y] of
	      data points or an nx2 Numeric array or a Geospatial_data object
              
	Output:
	  Interpolated values at inputted points (z).
        """
        

        # FIXME: There is an option of passing a tolerance into this
        
        if isinstance(point_coordinates, Geospatial_data):
            point_coordinates = point_coordinates.get_data_points( \
                absolute = True)
 
        from utilities.polygon import point_on_line
        from Numeric import ones
        z=ones(len(point_coordinates),Float)

        msg='point coordinates are not given (interpolate.py)'
        assert point_coordinates is not None, msg
        msg='function value must be specified at every interpolation node'
        assert f.shape[0]==vertex_coordinates.shape[0], msg
        msg='Must define function value at one or more nodes'
        assert f.shape[0]>0, msg

        n=f.shape[0]
        if n==1:
            z=f*z
            msg = 'Polyline contained only one point. I need more. ', str(f)
            raise Exception, msg
            
        # FIXME (John): add unit test for only One vertex point. Exception should be thrown.

        
        elif n>1:
            for i in range(len(point_coordinates)):
                found = False
                for j in range(n):
                    if gauge_neighbour_id[j]>=0:
                        if point_on_line(point_coordinates[i],
                                         [vertex_coordinates[j], vertex_coordinates[gauge_neighbour_id[j]]],
                                         rtol=1.0e-6):
                            found=True
                            x0=vertex_coordinates[j][0]
                            y0=vertex_coordinates[j][1]
                            x1=vertex_coordinates[gauge_neighbour_id[j]][0]
                            y1=vertex_coordinates[gauge_neighbour_id[j]][1]
                            x2=point_coordinates[i][0]
                            y2=point_coordinates[i][1]
                            
                            segment_len=sqrt((x1-x0)**2+(y1-y0)**2)
                            dist=sqrt((x2-x0)**2+(y2-y0)**2)
                            z[i]=(f[gauge_neighbour_id[j]]-f[j])/segment_len*dist+f[j]
                            #print 'element found on segment'
                            break
                                  
                if not found:
                    z[i]=0.0
                    #print 'point not on urs boundary'
        return z

    # FIXME: What is a good start_blocking_len value?
    def interpolate(self,
                    f,
                    point_coordinates=None,
                    start_blocking_len=500000,
                    verbose=False):
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
	      data points or an nx2 Numeric array or a Geospatial_data object
              
	      If point_coordinates is absent, the points inputted last time
              this method was called are used, if possible.
          start_blocking_len: If the # of points is more or greater than this,
              start blocking 

	Output:
	  Interpolated values at inputted points (z).
        """

        # FIXME (Ole): Why is the interpolation matrix rebuilt everytime the
        # method is called even if interpolation points are unchanged.

        #print "point_coordinates interpolate.interpolate", point_coordinates
        if verbose: print 'Build intepolation object' 
        if isinstance(point_coordinates, Geospatial_data):
            point_coordinates = point_coordinates.get_data_points( \
                absolute = True)

        # Can I interpolate, based on previous point_coordinates?
        if point_coordinates is None:
            if self._A_can_be_reused is True and \
                   len(self._point_coordinates) < start_blocking_len:
                z = self._get_point_data_z(f,
                                           verbose=verbose)
            elif self._point_coordinates is not None:
                #     if verbose, give warning
                if verbose:
                    print 'WARNING: Recalculating A matrix, due to blocking.'
                point_coordinates = self._point_coordinates
            else:
                #There are no good point_coordinates. import sys; sys.exit()
                msg = 'ERROR (interpolate.py): No point_coordinates inputted'
                raise Exception(msg)
            
        if point_coordinates is not None:
            self._point_coordinates = point_coordinates
            if len(point_coordinates) < start_blocking_len or \
                   start_blocking_len == 0:
                self._A_can_be_reused = True
                z = self.interpolate_block(f, point_coordinates,
                                           verbose=verbose)
            else:
                #print 'BLOCKING'
                #Handle blocking
                self._A_can_be_reused = False
                start = 0
                # creating a dummy array to concatenate to.
                
                f = ensure_numeric(f, Float)
                #print "f.shape",f.shape 
                if len(f.shape) > 1:
                    z = zeros((0,f.shape[1]))
                else:
                    z = zeros((0,))
                    
                for end in range(start_blocking_len,
                                 len(point_coordinates),
                                 start_blocking_len):
                    
                    t = self.interpolate_block(f, point_coordinates[start:end],
                                               verbose=verbose)
                    #print "t", t
                    #print "z", z 
                    z = concatenate((z,t))
                    start = end
                    
                end = len(point_coordinates)
                t = self.interpolate_block(f, point_coordinates[start:end],
                                           verbose=verbose)
                z = concatenate((z,t))
        return z
    

    def interpolate_block(self, f, point_coordinates, verbose=False):
        """
        Call this if you want to control the blocking or make sure blocking
        doesn't occur.

        Return the point data, z.
        
        See interpolate for doc info.
        """
        if isinstance(point_coordinates,Geospatial_data):
            point_coordinates = point_coordinates.get_data_points(\
                absolute=True)

        # Convert lists to Numeric arrays if necessary
        point_coordinates = ensure_numeric(point_coordinates, Float)
        f = ensure_numeric(f, Float)        
            
        self._A = self._build_interpolation_matrix_A(point_coordinates,
                                                     verbose=verbose)


        # Check that input dimensions are compatible
        msg = 'Two colums must be specified in point coordinates. I got shape=%s'\
              %(str(point_coordinates.shape))
        assert point_coordinates.shape[1] == 2, msg

        msg = 'The number of rows in matrix A must be the same as the number of points supplied.'
        msg += ' I got %d points and %d matrix rows.'\
               %(point_coordinates.shape[0], self._A.shape[0])
        assert point_coordinates.shape[0] == self._A.shape[0], msg        

        msg = 'The number of columns in matrix A must be the same as the number of mesh vertices.'
        msg += ' I got %d vertices and %d matrix columns.'\
               %(f.shape[0], self._A.shape[1])        
        assert self._A.shape[1] == f.shape[0], msg

        # Compute Matrix vector product and return
        return self._get_point_data_z(f)
    

    def _get_point_data_z(self, f, verbose=False):
        """
        Return the point data, z.
        
        Precondition,
        The _A matrix has been created
        """

        z = self._A * f
        # Taking into account points outside the mesh.
        #print "self.outside_poly_indices", self.outside_poly_indices
        #print "self.inside_poly_indices", self.inside_poly_indices
        #print "z", z
        for i in self.outside_poly_indices: 
            z[i] = NAN
        return z

    def _build_interpolation_matrix_A(self,
                                      point_coordinates,
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

        Preconditions
        Point_coordindates and mesh vertices have the same origin.
        """

        if verbose: print 'Building interpolation matrix'

        # Convert point_coordinates to Numeric arrays, in case it was a list.
        point_coordinates = ensure_numeric(point_coordinates, Float)
        
        
        if verbose: print 'Getting indices inside mesh boundary'
        self.inside_poly_indices, self.outside_poly_indices  = \
                     in_and_outside_polygon(point_coordinates,
                                            self.mesh.get_boundary_polygon(),
                                            closed = True, verbose = verbose)
        
        #Build n x m interpolation matrix
        if verbose and len(self.outside_poly_indices) > 0:
            print '\n WARNING: Points outside mesh boundary. \n'
        # Since you can block, throw a warning, not an error.
        if verbose and 0 == len(self.inside_poly_indices):
            print '\n WARNING: No points within the mesh! \n'
            
        m = self.mesh.number_of_nodes  # Nbr of basis functions (1/vertex)
        n = point_coordinates.shape[0] # Nbr of data points

        if verbose: print 'Number of datapoints: %d' %n
        if verbose: print 'Number of basis functions: %d' %m

        A = Sparse(n,m)

        n = len(self.inside_poly_indices)
        #Compute matrix elements for points inside the mesh
        if verbose: print 'Building interpolation matrix from %d points' %n
        for d, i in enumerate(self.inside_poly_indices):
            # For each data_coordinate point
            if verbose and d%((n+10)/10)==0: print 'Doing %d of %d' %(d, n)
            x = point_coordinates[i]
            element_found, sigma0, sigma1, sigma2, k = \
                           search_tree_of_vertices(self.root, self.mesh, x)
            
	    # Update interpolation matrix A if necessary
            if element_found is True:
                # Assign values to matrix A

                j0 = self.mesh.triangles[k,0] # Global vertex id for sigma0
                j1 = self.mesh.triangles[k,1] # Global vertex id for sigma1
                j2 = self.mesh.triangles[k,2] # Global vertex id for sigma2

                sigmas = {j0:sigma0, j1:sigma1, j2:sigma2}
                js     = [j0,j1,j2]

                for j in js:
                    A[i,j] = sigmas[j]
            else:
                msg = 'Could not find triangle for point', x 
                raise Exception(msg)
        return A

def benchmark_interpolate(vertices,
                          vertex_attributes,
                          triangles, points,
                          max_points_per_cell=None,
                          start_blocking_len=500000,
                          mesh_origin=None):
    """
    points: Interpolate mesh data to these positions.
              List of coordinate pairs [x, y] of
	      data points or an nx2 Numeric array or a Geospatial_data object
              
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
            
    calc = interp.interpolate(vertex_attributes
                              ,points
                              ,start_blocking_len=start_blocking_len)
    #print "calc", calc
    
def interpolate_sww2csv(sww_file,
                        points,
                        depth_file,
                        velocity_x_file,
                        velocity_y_file,
                        stage_file=None,
                        froude_file=None,
                        #quantities = ['depth', 'velocity'],
                        time_thinning=1,
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
    #print "points",points 
    points = ensure_absolute(points)
    point_count = len(points)
    callable_sww = file_function(sww_file,
                                 quantities=quantities,
                                 interpolation_points=points,
                                 verbose=verbose,
                                 time_thinning=time_thinning,
                                 use_cache=use_cache)
    
    depth_writer = writer(file(depth_file, "wb"))
    velocity_x_writer = writer(file(velocity_x_file, "wb"))
    velocity_y_writer = writer(file(velocity_y_file, "wb"))
    if stage_file is not None:
        stage_writer = writer(file(stage_file, "wb"))
    if froude_file is not None:
        froude_writer = writer(file(froude_file, "wb"))
    # Write heading
    heading = [str(x[0])+ ':' + str(x[1]) for x in points]
    heading.insert(0, "time")
    depth_writer.writerow(heading)
    velocity_x_writer.writerow(heading)
    velocity_y_writer.writerow(heading)
    if stage_file is not None:
        stage_writer.writerow(heading) 
    if froude_file is not None:
        froude_writer.writerow(heading)     
    
    for time in callable_sww.get_time():
        depths = [time]
        velocity_xs = [time]
        velocity_ys = [time]
        if stage_file is not None:  
            stages = [time]  
        if froude_file is not None:  
            froudes = [time]  
        for point_i, point in enumerate(points):
            quantities = callable_sww(time,point_i)
            #print "quantities", quantities
            
            w = quantities[0]
            z = quantities[1]
            momentum_x = quantities[2]
            momentum_y = quantities[3]
            depth = w - z 
              
            if w == NAN or z == NAN or momentum_x == NAN:
                velocity_x = NAN
            else:
                if depth > 1.e-30: # use epsilon
                    velocity_x = momentum_x / depth  #Absolute velocity
                else:
                    velocity_x = 0
            if w == NAN or z == NAN or momentum_y == NAN:
                velocity_y = NAN
            else:
                if depth > 1.e-30: # use epsilon
                    velocity_y = momentum_y / depth  #Absolute velocity
                else:
                    velocity_y = 0
            if depth < 1.e-30: # use epsilon
                froude = NAN
            else:
                froude = sqrt(velocity_x*velocity_x + velocity_y*velocity_y)/ \
                         sqrt(depth * 9.8066) # gravity m/s/s
            depths.append(depth)
            velocity_xs.append(velocity_x)
            velocity_ys.append(velocity_y)
            if stage_file is not None:
                stages.append(w)
            if froude_file is not None:
                froudes.append(froude)
        depth_writer.writerow(depths)
        velocity_x_writer.writerow(velocity_xs)
        velocity_y_writer.writerow(velocity_ys)
        if stage_file is not None:
            stage_writer.writerow(stages)   
        if froude_file is not None:
            froude_writer.writerow(froudes)            


class Interpolation_function:
    """Interpolation_interface - creates callable object f(t, id) or f(t,x,y)
    which is interpolated from time series defined at vertices of
    triangular mesh (such as those stored in sww files)

    Let m be the number of vertices, n the number of triangles
    and p the number of timesteps.
    Also, let N be the number of interpolation points.

    Mandatory input
        time:               px1 array of monotonously increasing times (Float)
        quantities:         Dictionary of arrays or 1 array (Float) 
                            The arrays must either have dimensions pxm or mx1.
                            The resulting function will be time dependent in
                            the former case while it will be constant with
                            respect to time in the latter case.
        
    Optional input:
        quantity_names:     List of keys into the quantities dictionary for
                            imposing a particular order on the output vector.
        vertex_coordinates: mx2 array of coordinates (Float)
        triangles:          nx3 array of indices into vertex_coordinates (Int)
        interpolation_points: Nx2 array of coordinates to be interpolated to 
        verbose:            Level of reporting
    
    
    The quantities returned by the callable object are specified by
    the list quantities which must contain the names of the
    quantities to be returned and also reflect the order, e.g. for
    the shallow water wave equation, on would have
    quantities = ['stage', 'xmomentum', 'ymomentum']

    The parameter interpolation_points decides at which points interpolated
    quantities are to be computed whenever object is called.
    If None, return average value

    FIXME (Ole): Need to allow vertex coordinates and interpolation points to be
    geospatial data objects

    Time assumed to be relative to starttime (FIXME (Ole): This comment should be removed)
    All coordinates assume origin of (0,0) - e.g. georeferencing must be taken care of
    outside this function
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
                 gauge_neighbour_id=None):
        """Initialise object and build spatial interpolation if required

        Time_thinning_number controls how many timesteps to use. Only timesteps with
        index%time_thinning_number == 0 will used, or in other words a value of 3, say,
        will cause the algorithm to use every third time step.
        """

        from Numeric import array, zeros, Float, alltrue, concatenate,\
             reshape, ArrayType


        from anuga.config import time_format
        import types


	# Check temporal info
        time = ensure_numeric(time)        
        msg = 'Time must be a monotonuosly '
        msg += 'increasing sequence %s' %time
        assert alltrue(time[1:] - time[:-1] >= 0 ), msg


        # Check if quantities is a single array only
        if type(quantities) != types.DictType:
            quantities = ensure_numeric(quantities)
            quantity_names = ['Attribute']

            # Make it a dictionary
            quantities = {quantity_names[0]: quantities}


        # Use keys if no names are specified
        if quantity_names is None:
            quantity_names = quantities.keys()


        # Check spatial info
        if vertex_coordinates is None:
            self.spatial = False
        else:
            # FIXME (Ole): Try ensure_numeric here -
            #this function knows nothing about georefering.
            vertex_coordinates = ensure_absolute(vertex_coordinates)

            if triangles is not None:
                triangles = ensure_numeric(triangles)
            self.spatial = True          

        # Thin timesteps if needed
        # Note array() is used to make the thinned arrays contiguous in memory
        self.time = array(time[::time_thinning])          
        for name in quantity_names:
            if len(quantities[name].shape) == 2:
                quantities[name] = array(quantities[name][::time_thinning,:])
             
        # Save for use with statistics
        self.quantities_range = {}
        for name in quantity_names:
            q = quantities[name][:].flat
            self.quantities_range[name] = [min(q), max(q)]
        
        self.quantity_names = quantity_names        
        self.vertex_coordinates = vertex_coordinates 
        self.interpolation_points = interpolation_points
        

        self.index = 0    # Initial time index
        self.precomputed_values = {}
        
            
        # Precomputed spatial interpolation if requested
        if interpolation_points is not None:
            #no longer true. sts files have spatial = True but
            #if self.spatial is False:
            #    raise 'Triangles and vertex_coordinates must be specified'
            # 
            try:
	        self.interpolation_points = interpolation_points = ensure_numeric(interpolation_points)
            except:
	        msg = 'Interpolation points must be an N x 2 Numeric array '+\
                      'or a list of points\n'
		msg += 'I got: %s.' %(str(self.interpolation_points)[:60] +\
                                      '...')
                raise msg

            if triangles is not None and vertex_coordinates is not None:
                # Check that all interpolation points fall within
                # mesh boundary as defined by triangles and vertex_coordinates.
                from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
                from anuga.utilities.polygon import outside_polygon            

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
                    msg += 'outside specified mesh. '
                    msg += 'Offending points:\n'
                    out_interp_pts = []
                    for i in indices:
                        msg += '%d: %s\n' %(i, interpolation_points[i])
                        out_interp_pts.append(ensure_numeric(interpolation_points[i]))

                    if verbose is True:
                        import sys
                        if sys.platform == 'win32': # FIXME (Ole): Why only Windoze?
                            from anuga.utilities.polygon import plot_polygons
                            #out_interp_pts = take(interpolation_points,[indices])
                            title = 'Interpolation points fall outside specified mesh'
                            plot_polygons([mesh_boundary_polygon,
                                           interpolation_points,
                                           out_interp_pts],
                                          ['line','point','outside'],
                                          figname='points_boundary_out',
                                          label=title,
                                          verbose=verbose)

                    # Joaquim Luis suggested this as an Exception, so
                    # that the user can now what the problem is rather than
                    # looking for NaN's. However, NANs are handy as they can
                    # be ignored leaving good points for continued processing.
                    if verbose:
                        print msg
                    #raise Exception(msg)
                    
            elif triangles is None and vertex_coordinates is not None:#jj
                #Dealing with sts file
                pass
            else:
                msg = 'Sww file function requires both triangles and vertex_coordinates. sts file file function requires the later.'
                raise Exception(msg)

            # Plot boundary and interpolation points
            if verbose is True:
                import sys
                if sys.platform == 'win32':
                    from anuga.utilities.polygon import plot_polygons
                    title = 'Interpolation function: Polygon and interpolation points'
                    plot_polygons([mesh_boundary_polygon,
                                   interpolation_points],
                                  ['line','point'],
                                  figname='points_boundary',
                                  label=title,
                                  verbose=verbose)

            m = len(self.interpolation_points)
            p = len(self.time)
            
	    for name in quantity_names:
                self.precomputed_values[name] = zeros((p, m), Float)

            # Build interpolator
            if verbose:
                if triangles is not None and vertex_coordinates is not None:
                    msg = 'Building interpolation matrix from source mesh '
                    msg += '(%d vertices, %d triangles)' %(vertex_coordinates.shape[0],
                                                           triangles.shape[0])
                elif triangles is None and vertex_coordinates is not None:
                    msg = 'Building interpolation matrix from source points'
                
                print msg

                
            interpol = Interpolate(vertex_coordinates,
                                   triangles,
                                   verbose=verbose)

            if verbose:
                print 'Interpolating (%d interpolation points, %d timesteps).'\
                      %(self.interpolation_points.shape[0], self.time.shape[0]),
            
                if time_thinning > 1:
                    print 'Timesteps were thinned by a factor of %d' %time_thinning
                else:
                    print

	    for i, t in enumerate(self.time):
                # Interpolate quantities at this timestep
                if verbose and i%((p+10)/10)==0:
                    print '  time step %d of %d' %(i, p)
                    
                for name in quantity_names:
                    if len(quantities[name].shape) == 2:
                        Q = quantities[name][i,:] # Quantities at timestep i
                    else:
                        Q = quantities[name][:]   # No time dependency

                    if verbose and i%((p+10)/10)==0:
                        print '    quantity %s, size=%d' %(name, len(Q))
                        
                    # Interpolate 
                    if triangles is not None and vertex_coordinates is not None:   
                        result = interpol.interpolate(Q,
                                                      point_coordinates=\
                                                      self.interpolation_points,
                                                      verbose=False) # Don't clutter
                    elif triangles is None and vertex_coordinates is not None:
                        result=interpol.interpolate_polyline(Q,vertex_coordinates,gauge_neighbour_id,point_coordinates=self.interpolation_points)

                    #assert len(result), len(interpolation_points)
                    self.precomputed_values[name][i, :] = result

                   
            # Report
            if verbose:
                print self.statistics()
                #self.print_statistics()
            
        else:
            # Store quantitites as is
            for name in quantity_names:
                self.precomputed_values[name] = quantities[name]

    def __repr__(self):
        # return 'Interpolation function (spatio-temporal)'
        return self.statistics()
 
    def __call__(self, t, point_id=None, x=None, y=None):
        """Evaluate f(t) or f(t, point_id)
        
	Inputs:
	  t: time - Model time. Must lie within existing timesteps
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
        from Numeric import zeros, Float
        from anuga.abstract_2d_finite_volumes.util import mean        

        if self.spatial is True:
            if point_id is None:
                if x is None or y is None:
                    msg = 'Either point_id or x and y must be specified'
                    raise Exception(msg)
            else:
                if self.interpolation_points is None:
                    msg = 'Interpolation_function must be instantiated ' +\
                          'with a list of interpolation points before parameter ' +\
                          'point_id can be used'
                    raise Exception(msg)

        msg = 'Time interval [%.16f:%.16f]' %(self.time[0], self.time[-1])
        msg += ' does not match model time: %.16f\n' %t
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
            ratio = (t - self.time[self.index])/\
                    (self.time[self.index+1] - self.time[self.index])

        # Compute interpolated values
        q = zeros(len(self.quantity_names), Float)
        # print "self.precomputed_values", self.precomputed_values
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
                    
                    raise 'x,y interpolation not yet implemented'
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
        # if len(q) == 1:
        #     return q[0]
        # else:
        #     return q


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
                    from Numeric import ones, Float
                    # x is a vector - Create one constant column for each value
                    N = len(x)
                    assert len(y) == N, 'x and y must have same length'
                    res = []
                    for col in q:
                        res.append(col*ones(N, Float))
                        
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
                
        x = vertex_coordinates[:,0]
        y = vertex_coordinates[:,1]                

        str =  '------------------------------------------------\n'
        str += 'Interpolation_function (spatio-temporal) statistics:\n'
        str += '  Extent:\n'
        str += '    x in [%f, %f], len(x) == %d\n'\
               %(min(x), max(x), len(x))
        str += '    y in [%f, %f], len(y) == %d\n'\
               %(min(y), max(y), len(y))
        str += '    t in [%f, %f], len(t) == %d\n'\
               %(min(self.time), max(self.time), len(self.time))
        str += '  Quantities:\n'
        for name in quantity_names:
            minq, maxq = self.quantities_range[name]
            str += '    %s in [%f, %f]\n' %(name, minq, maxq)            
            #q = quantities[name][:].flat
            #str += '    %s in [%f, %f]\n' %(name, min(q), max(q))

        if interpolation_points is not None:    
            str += '  Interpolation points (xi, eta):'\
                   ' number of points == %d\n' %interpolation_points.shape[0]
            str += '    xi in [%f, %f]\n' %(min(interpolation_points[:,0]),
                                            max(interpolation_points[:,0]))
            str += '    eta in [%f, %f]\n' %(min(interpolation_points[:,1]),
                                             max(interpolation_points[:,1]))
            str += '  Interpolated quantities (over all timesteps):\n'
        
            for name in quantity_names:
                q = precomputed_values[name][:].flat
                str += '    %s at interpolation points in [%f, %f]\n'\
                       %(name, min(q), max(q))
        str += '------------------------------------------------\n'

        return str


def interpolate_sww(sww_file, time, interpolation_points,
                    quantity_names=None, verbose=False):
    """
    obsolete.
    use file_function in utils
    """
    #open sww file
    x, y, volumes, time, quantities = read_sww(sww_file)
    print "x",x
    print "y",y
    
    print "time", time
    print "quantities", quantities

    #Add the x and y together
    vertex_coordinates = concatenate((x[:,NewAxis], y[:,NewAxis]),axis=1)

    #Will return the quantity values at the specified times and locations
    interp = Interpolation_interface(time,
                                     quantities,
                                     quantity_names=quantity_names,  
                                     vertex_coordinates=vertex_coordinates,
                                     triangles=volumes,
                                     interpolation_points=interpolation_points,
                                     verbose=verbose)


def read_sww(file_name):
    """
    obsolete - Nothing should be calling this
    
    Read in an sww file.
    
    Input;
    file_name - the sww file
    
    Output;
    x - Vector of x values
    y - Vector of y values
    z - Vector of bed elevation
    volumes - Array.  Each row has 3 values, representing
    the vertices that define the volume
    time - Vector of the times where there is stage information
    stage - array with respect to time and vertices (x,y)
    """

    msg = 'Function read_sww in interpolat.py is obsolete'
    raise Exception, msg
    
    #FIXME Have this reader as part of data_manager?
    
    from Scientific.IO.NetCDF import NetCDFFile     
    import tempfile
    import sys
    import os
    
    #Check contents
    #Get NetCDF
    
    # see if the file is there.  Throw a QUIET IO error if it isn't
    # I don't think I could implement the above
    
    #throws prints to screen if file not present
    junk = tempfile.mktemp(".txt")
    fd = open(junk,'w')
    stdout = sys.stdout
    sys.stdout = fd
    fid = NetCDFFile(file_name, 'r') 
    sys.stdout = stdout
    fd.close()
    #clean up
    os.remove(junk)    	
      
    # Get the variables
    x = fid.variables['x'][:]
    y = fid.variables['y'][:]
    volumes = fid.variables['volumes'][:] 
    time = fid.variables['time'][:]

    keys = fid.variables.keys()
    keys.remove("x")
    keys.remove("y")
    keys.remove("volumes")
    keys.remove("time")
     #Turn NetCDF objects into Numeric arrays
    quantities = {}
    for name in keys:
        quantities[name] = fid.variables[name][:]
            
    fid.close()
    return x, y, volumes, time, quantities


#-------------------------------------------------------------
if __name__ == "__main__":
    names = ["x","y"]
    someiterable = [[1,2],[3,4]]
    csvwriter = writer(file("some.csv", "wb"))
    csvwriter.writerow(names)
    for row in someiterable:
        csvwriter.writerow(row)
