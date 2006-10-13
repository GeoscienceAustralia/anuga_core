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
from anuga.utilities.quad import build_quadtree
from anuga.utilities.numerical_tools import ensure_numeric, mean, NAN
from anuga.utilities.polygon import in_and_outside_polygon
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.geospatial_data.geospatial_data import ensure_absolute
from anuga.fit_interpolate.search_functions import search_tree_of_vertices
from anuga.fit_interpolate.general_fit_interpolate import FitInterpolate
from anuga.abstract_2d_finite_volumes.util import file_function


class Interpolate (FitInterpolate):
        
    def __init__(self,
                 vertex_coordinates,
                 triangles,
                 mesh_origin=None,
                 verbose=False,
                 max_vertices_per_cell=30):


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
                                vertex_coordinates,
                                triangles,
                                mesh_origin,
                                verbose,
                                max_vertices_per_cell)

    # FIXME: What is a good start_blocking_len value?
    def interpolate(self, f, point_coordinates = None,
                    start_blocking_len = 500000, verbose=False):
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

        
        # FIXME (Ole): Need an input check that dimensions are compatible

        # FIXME (Ole): Why is the interpolation matrix rebuilt everytime the
        # method is called even if interpolation points are unchanged.
        
        #print "point_coordinates interpolate.interpolate", point_coordinates 
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

    def interpolate_block(self, f, point_coordinates = None, verbose=False):
        """
        Call this if you want to control the blocking or make sure blocking
        doesn't occur.

        Return the point data, z.
        
        See interpolate for doc info.
        """
        if isinstance(point_coordinates,Geospatial_data):
            point_coordinates = point_coordinates.get_data_points( \
                absolute = True)
        if point_coordinates is not None:
            self._A = self._build_interpolation_matrix_A(point_coordinates,
                                                         verbose=verbose)
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
                                      verbose = False):
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

        #print 'Building interpolation matrix'
        
        #Convert point_coordinates to Numeric arrays, in case it was a list.
        point_coordinates = ensure_numeric(point_coordinates, Float)
	
        if verbose: print 'Getting indices inside mesh boundary'
        self.inside_poly_indices, self.outside_poly_indices  = \
                     in_and_outside_polygon(point_coordinates,
                                            self.mesh.get_boundary_polygon(),
                                            closed = True, verbose = verbose)
        #print "self.inside_poly_indices",self.inside_poly_indices
        #print "self.outside_poly_indices",self.outside_poly_indices 
        #Build n x m interpolation matrix
        if verbose and len(self.outside_poly_indices) >0:
            print '\n WARNING: Points outside mesh boundary. \n'
        # Since you can block, throw a warning, not an error.
        if verbose and 0 == len(self.inside_poly_indices):
            print '\n WARNING: No points within the mesh! \n'
            
        m = self.mesh.coordinates.shape[0] #Nbr of basis functions (1/vertex)
        n = point_coordinates.shape[0]     #Nbr of data points

        if verbose: print 'Number of datapoints: %d' %n
        if verbose: print 'Number of basis functions: %d' %m

        A = Sparse(n,m)

        n = len(self.inside_poly_indices)
        #Compute matrix elements for points inside the mesh
        if verbose: print 'Building interpolation matrix fram %d points' %n
        for k, i in enumerate(self.inside_poly_indices):
            #For each data_coordinate point
            if verbose and k%((n+10)/10)==0: print 'Doing %d of %d' %(k, n)
            x = point_coordinates[i]
            element_found, sigma0, sigma1, sigma2, k = \
                           search_tree_of_vertices(self.root, self.mesh, x)
	    #Update interpolation matrix A if necessary
            if element_found is True:
                #Assign values to matrix A

                j0 = self.mesh.triangles[k,0] #Global vertex id for sigma0
                j1 = self.mesh.triangles[k,1] #Global vertex id for sigma1
                j2 = self.mesh.triangles[k,2] #Global vertex id for sigma2

                sigmas = {j0:sigma0, j1:sigma1, j2:sigma2}
                js     = [j0,j1,j2]

                for j in js:
                    A[i,j] = sigmas[j]
            else:
                msg = 'Could not find triangle for point', x 
                raise Exception(msg)
        return A

def interpolate_sww2csv(sww_file,
                        points,
                        depth_file,
                        velocity_x_file,
                        velocity_y_file,
                        #quantities = ['depth', 'velocity'],
                        verbose=True,
                        use_cache = True):
    """
    Interpolate the quantities at a given set of locations, given
    an sww file.
    The results are written to a csv file.

    In the future let points be a csv or xya file.
    And the user choose the quantities.

    This is currently quite specific.
    If it need to be more general, chagne things.

    This is really returning speed, not velocity.
    """
    quantities =  ['stage', 'elevation', 'xmomentum', 'ymomentum']
    #print "points",points 
    points = ensure_absolute(points)
    point_count = len(points)
    callable_sww = file_function(sww_file,
                                 quantities=quantities,
                                 interpolation_points=points,
                                 verbose=verbose,
                                 use_cache=use_cache)
    
    depth_writer = writer(file(depth_file, "wb"))
    velocity_x_writer = writer(file(velocity_x_file, "wb"))
    velocity_y_writer = writer(file(velocity_y_file, "wb"))
    # Write heading
    heading = [str(x[0])+ ':' + str(x[1]) for x in points]
    heading.insert(0, "time")
    depth_writer.writerow(heading)
    velocity_x_writer.writerow(heading)
    velocity_y_writer.writerow(heading)
    
    for time in callable_sww.get_time():
        depths = [time]
        velocity_xs = [time]
        velocity_ys = [time]
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
            depths.append(depth)
            velocity_xs.append(velocity_x)
            velocity_ys.append(velocity_y)
        depth_writer.writerow(depths)
        velocity_x_writer.writerow(velocity_xs)
        velocity_y_writer.writerow(velocity_ys)


class Interpolation_function:
    """Interpolation_interface - creates callable object f(t, id) or f(t,x,y)
    which is interpolated from time series defined at vertices of
    triangular mesh (such as those stored in sww files)

    Let m be the number of vertices, n the number of triangles
    and p the number of timesteps. 

    Mandatory input
        time:               px1 array of monotonously increasing times (Float)
        quantities:         Dictionary of arrays or 1 array (Float) 
                            The arrays must either have dimensions pxm or mx1.
                            The resulting function will be time dependent in
                            the former case while it will be constan with
                            respect to time in the latter case.
        
    Optional input:
        quantity_names:     List of keys into the quantities dictionary 
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
    """
  
    
    def __init__(self,
                 time,
                 quantities,
                 quantity_names = None,  
                 vertex_coordinates = None,
                 triangles = None,
                 interpolation_points = None,
                 verbose = False):
        """Initialise object and build spatial interpolation if required
        """

        from Numeric import array, zeros, Float, alltrue, concatenate,\
             reshape, ArrayType


        #from anuga.abstract_2d_finite_volumes.util import mean, ensure_numeric
        from anuga.config import time_format
        import types


	#Check temporal info
        time = ensure_numeric(time)        
        msg = 'Time must be a monotonuosly '
        msg += 'increasing sequence %s' %time
        assert alltrue(time[1:] - time[:-1] >= 0 ), msg


        #Check if quantities is a single array only
        if type(quantities) != types.DictType:
            quantities = ensure_numeric(quantities)
            quantity_names = ['Attribute']

            #Make it a dictionary
            quantities = {quantity_names[0]: quantities}


        #Use keys if no names are specified
        if quantity_names is None:
            quantity_names = quantities.keys()


        #Check spatial info
        if vertex_coordinates is None:
            self.spatial = False
        else:    
            vertex_coordinates = ensure_numeric(vertex_coordinates)

            assert triangles is not None, 'Triangles array must be specified'
            triangles = ensure_numeric(triangles)
            self.spatial = True            

             
        #Save for use with statistics

        self.quantities_range = {}
        for name in quantity_names:
            q = quantities[name][:].flat
            self.quantities_range[name] = [min(q), max(q)]
        
        self.quantity_names = quantity_names        
        #self.quantities = quantities  #Takes too much memory      
        self.vertex_coordinates = vertex_coordinates 
        self.interpolation_points = interpolation_points
        self.time = time[:]  # Time assumed to be relative to starttime
        self.index = 0    # Initial time index
        self.precomputed_values = {}
        
            
        #Precomputed spatial interpolation if requested
        if interpolation_points is not None:
            if self.spatial is False:
                raise 'Triangles and vertex_coordinates must be specified'
            
            try:
	        self.interpolation_points = ensure_numeric(interpolation_points)
            except:
	        msg = 'Interpolation points must be an N x 2 Numeric array '+\
                      'or a list of points\n'
		msg += 'I got: %s.' %(str(self.interpolation_points)[:60] +\
                                      '...')
                raise msg


            m = len(self.interpolation_points)
            p = len(self.time)
            
	    for name in quantity_names:
                self.precomputed_values[name] = zeros((p, m), Float)

            #Build interpolator
            interpol = Interpolate(vertex_coordinates,
                                   triangles,
                                   #point_coordinates = \
                                   #self.interpolation_points,
                                   #alpha = 0,
                                   verbose = verbose)

            if verbose: print 'Interpolate'
	    for i, t in enumerate(self.time):
                # Interpolate quantities at this timestep
                if verbose and i%((p+10)/10)==0:
                    print ' time step %d of %d' %(i, p)
                    
                for name in quantity_names:
                    if len(quantities[name].shape) == 2:
                        Q = quantities[name][i,:] # Quantities at timestep i
                    else:
                        Q = quantities[name][:]   # No time dependency
                       
                    # Interpolate    
                    result = interpol.interpolate(Q,
                                                  point_coordinates=\
                                                  self.interpolation_points)
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
        """Evaluate f(t), f(t, point_id) or f(t, x, y)

	Inputs:
	  t: time - Model time. Must lie within existing timesteps
	  point_id: index of one of the preprocessed points.
          x, y:     Overrides location, point_id ignored
          
	  If spatial info is present and all of x,y,point_id
          are None an exception is raised 
                    
          If no spatial info is present, point_id and x,y arguments are ignored
          making f a function of time only.

          
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
        if t < self.time[0]: raise Exception(msg)
        if t > self.time[-1]: raise Exception(msg)

        oldindex = self.index #Time index

        #Find current time slot
        while t > self.time[self.index]: self.index += 1
        while t < self.time[self.index]: self.index -= 1

        if t == self.time[self.index]:
            #Protect against case where t == T[-1] (last time)
            # - also works in general when t == T[i]
            ratio = 0
        else:
            #t is now between index and index+1
            ratio = (t - self.time[self.index])/\
                    (self.time[self.index+1] - self.time[self.index])

        #Compute interpolated values
        q = zeros(len(self.quantity_names), Float)
        #print "self.precomputed_values", self.precomputed_values
	for i, name in enumerate(self.quantity_names):
            Q = self.precomputed_values[name]

            if self.spatial is False:
                #If there is no spatial info                
                assert len(Q.shape) == 1

                Q0 = Q[self.index]
                if ratio > 0: Q1 = Q[self.index+1]

            else:
                if x is not None and y is not None:
                    #Interpolate to x, y
                    
                    raise 'x,y interpolation not yet implemented'
                else:
                    #Use precomputed point
                    Q0 = Q[self.index, point_id]
                    if ratio > 0:
                        Q1 = Q[self.index+1, point_id]
            
            #Linear temporal interpolation    
            if ratio > 0:
                if Q0 == NAN and Q1 == NAN:
                    q[i]  = Q0
                else:
                    q[i] = Q0 + ratio*(Q1 - Q0)
            else:
                q[i] = Q0


        #Return vector of interpolated values
        #if len(q) == 1:
        #    return q[0]
        #else:
        #    return q


        #Return vector of interpolated values
        #FIXME:
        if self.spatial is True:
            return q
        else:
            #Replicate q according to x and y
            #This is e.g used for Wind_stress
            if x is None or y is None: 
                return q
            else:
                try:
                    N = len(x)
                except:
                    return q
                else:
                    from Numeric import ones, Float
                    #x is a vector - Create one constant column for each value
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
                    quantity_names = None, verbose = False):
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
