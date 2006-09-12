"""Least squares smooting and interpolation.

   Implements a penalised least-squares fit and associated interpolations.

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

import exceptions
class ShapeError(exceptions.Exception): pass
class FittingError(exceptions.Exception): pass


#from general_mesh import General_mesh
from Numeric import zeros, array, Float, Int, transpose, concatenate, ArrayType, NewAxis
from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh

from Numeric import dot, zeros, take, compress, array, Float, Int, transpose, concatenate, ArrayType
from anuga.utilities.sparse import Sparse, Sparse_CSR
from anuga.utilities.cg_solve import conjugate_gradient, VectorShapeError
from anuga.utilities.numerical_tools import ensure_numeric, mean, gradient


from anuga.coordinate_transforms.geo_reference import Geo_reference

import time

DEFAULT_ALPHA = 0.001

def fit_to_mesh_file(mesh_file, point_file, mesh_output_file,
                     alpha=DEFAULT_ALPHA, verbose= False,
                     expand_search = False,
                     data_origin = None,
                     mesh_origin = None,
                     precrop = False,
                     display_errors = True):
    """
    Given a mesh file (tsh) and a point attribute file (xya), fit
    point attributes to the mesh and write a mesh file with the
    results.


    If data_origin is not None it is assumed to be
    a 3-tuple with geo referenced
    UTM coordinates (zone, easting, northing)

    NOTE: Throws IOErrors, for a variety of file problems.
    
    mesh_origin is the same but refers to the input tsh file.
    FIXME: When the tsh format contains it own origin, these parameters can go.
    FIXME: And both origins should be obtained from the specified files.
    """

    from load_mesh.loadASCII import import_mesh_file, \
                 import_points_file, export_mesh_file, \
                 concatinate_attributelist


    try:
        mesh_dict = import_mesh_file(mesh_file)
    except IOError,e:
        if display_errors:
            print "Could not load bad file. ", e
        raise IOError  #Re-raise exception
        
    vertex_coordinates = mesh_dict['vertices']
    triangles = mesh_dict['triangles']
    if type(mesh_dict['vertex_attributes']) == ArrayType:
        old_point_attributes = mesh_dict['vertex_attributes'].tolist()
    else:
        old_point_attributes = mesh_dict['vertex_attributes']

    if type(mesh_dict['vertex_attribute_titles']) == ArrayType:
        old_title_list = mesh_dict['vertex_attribute_titles'].tolist()
    else:
        old_title_list = mesh_dict['vertex_attribute_titles']

    if verbose: print 'tsh file %s loaded' %mesh_file

    # load in the .pts file
    try:
        point_dict = import_points_file(point_file, verbose=verbose)
    except IOError,e:
        if display_errors:
            print "Could not load bad file. ", e
        raise IOError  #Re-raise exception  

    point_coordinates = point_dict['pointlist']
    title_list,point_attributes = concatinate_attributelist(point_dict['attributelist'])

    if point_dict.has_key('geo_reference') and not point_dict['geo_reference'] is None:
        data_origin = point_dict['geo_reference'].get_origin()
    else:
        data_origin = (56, 0, 0) #FIXME(DSG-DSG)

    if mesh_dict.has_key('geo_reference') and not mesh_dict['geo_reference'] is None:
        mesh_origin = mesh_dict['geo_reference'].get_origin()
    else:
        mesh_origin = (56, 0, 0) #FIXME(DSG-DSG)

    if verbose: print "points file loaded"
    if verbose: print "fitting to mesh"
    f = fit_to_mesh(vertex_coordinates,
                    triangles,
                    point_coordinates,
                    point_attributes,
                    alpha = alpha,
                    verbose = verbose,
                    expand_search = expand_search,
                    data_origin = data_origin,
                    mesh_origin = mesh_origin,
                    precrop = precrop)
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

    #FIXME (Ole): Remember to output mesh_origin as well
    if verbose: print "exporting to file ", mesh_output_file

    try:
        export_mesh_file(mesh_output_file, mesh_dict)
    except IOError,e:
        if display_errors:
            print "Could not write file. ", e
        raise IOError

def fit_to_mesh(vertex_coordinates,
                triangles,
                point_coordinates,
                point_attributes,
                alpha = DEFAULT_ALPHA,
                verbose = False,
                acceptable_overshoot = 1.01,
                expand_search = False,
                data_origin = None,
                mesh_origin = None,
                precrop = False,
                use_cache = False):
    """
    Fit a smooth surface to a triangulation,
    given data points with attributes.


        Inputs:

          vertex_coordinates: List of coordinate pairs [xi, eta] of points
          constituting mesh (or a an m x 2 Numeric array)

          triangles: List of 3-tuples (or a Numeric array) of
          integers representing indices of all vertices in the mesh.

          point_coordinates: List of coordinate pairs [x, y] of data points
          (or an nx2 Numeric array)

          alpha: Smoothing parameter.

          acceptable overshoot: controls the allowed factor by which fitted values
          may exceed the value of input data. The lower limit is defined
          as min(z) - acceptable_overshoot*delta z and upper limit
          as max(z) + acceptable_overshoot*delta z
          

          point_attributes: Vector or array of data at the point_coordinates.

          data_origin and mesh_origin are 3-tuples consisting of
          UTM zone, easting and northing. If specified
          point coordinates and vertex coordinates are assumed to be
          relative to their respective origins.

    """
    if use_cache is True:
        from anuga.caching.caching import cache
        interp = cache(_interpolation,
                       (vertex_coordinates,
                        triangles,
                        point_coordinates),
                       {'alpha': alpha,
                        'verbose': verbose,
                        'expand_search': expand_search,
                        'data_origin': data_origin,
                        'mesh_origin': mesh_origin,
                        'precrop': precrop},
                       verbose = verbose)        
        
    else:
        interp = Interpolation(vertex_coordinates,
                               triangles,
                               point_coordinates,
                               alpha = alpha,
                               verbose = verbose,
                               expand_search = expand_search,
                               data_origin = data_origin,
                               mesh_origin = mesh_origin,
                               precrop = precrop)

    vertex_attributes = interp.fit_points(point_attributes, verbose = verbose)


    #Sanity check
    point_coordinates = ensure_numeric(point_coordinates)
    vertex_coordinates = ensure_numeric(vertex_coordinates)

    #Data points
    X = point_coordinates[:,0]
    Y = point_coordinates[:,1]	
    Z = ensure_numeric(point_attributes)
    if len(Z.shape) == 1:
        Z = Z[:, NewAxis]
        

    #Data points inside mesh boundary
    indices = interp.point_indices
    if indices is not None:    
        Xc = take(X, indices)
        Yc = take(Y, indices)	
        Zc = take(Z, indices)
    else:
        Xc = X
        Yc = Y	
        Zc = Z        
    
    #Vertex coordinates
    Xi = vertex_coordinates[:,0]
    Eta = vertex_coordinates[:,1]	
    Zeta = ensure_numeric(vertex_attributes)
    if len(Zeta.shape) == 1:
        Zeta = Zeta[:, NewAxis]    

    for i in range(Zeta.shape[1]): #For each attribute
        zeta = Zeta[:,i]
        z = Z[:,i]                
        zc = Zc[:,i]

        max_zc = max(zc)
        min_zc = min(zc)
        delta_zc = max_zc-min_zc
        upper_limit = max_zc + delta_zc*acceptable_overshoot
        lower_limit = min_zc - delta_zc*acceptable_overshoot        
        

        if max(zeta) > upper_limit or min(zeta) < lower_limit:
            msg = 'Least sqares produced values outside the allowed '
            msg += 'range [%f, %f].\n' %(lower_limit, upper_limit)
            msg += 'z in [%f, %f], zeta in [%f, %f].\n' %(min_zc, max_zc,
                                                          min(zeta), max(zeta))
            msg += 'If greater range is needed, increase the value of '
            msg += 'acceptable_fit_overshoot (currently %.2f).\n' %(acceptable_overshoot)


            offending_vertices = (zeta > upper_limit or zeta < lower_limit)
            Xi_c = compress(offending_vertices, Xi)
            Eta_c = compress(offending_vertices, Eta)
            offending_coordinates = concatenate((Xi_c[:, NewAxis],
                                                 Eta_c[:, NewAxis]),
                                                axis=1)

            msg += 'Offending locations:\n %s' %(offending_coordinates)
            
            raise FittingError, msg


    
        if verbose:
            print '+------------------------------------------------'
            print 'Least squares statistics'
            print '+------------------------------------------------'	
            print 'points: %d points' %(len(z))
            print '    x in [%f, %f]'%(min(X), max(X))
            print '    y in [%f, %f]'%(min(Y), max(Y))
            print '    z in [%f, %f]'%(min(z), max(z))
            print

            if indices is not None:
                print 'Cropped points: %d points' %(len(zc))
                print '    x in [%f, %f]'%(min(Xc), max(Xc))
                print '    y in [%f, %f]'%(min(Yc), max(Yc))
                print '    z in [%f, %f]'%(min(zc), max(zc))
                print
            

            print 'Mesh: %d vertices' %(len(zeta))
            print '    xi in [%f, %f]'%(min(Xi), max(Xi))
            print '    eta in [%f, %f]'%(min(Eta), max(Eta))
            print '    zeta in [%f, %f]'%(min(zeta), max(zeta))
            print '+------------------------------------------------'

    return vertex_attributes



def pts2rectangular(pts_name, M, N, alpha = DEFAULT_ALPHA,
                    verbose = False, reduction = 1):
    """Fits attributes from pts file to MxN rectangular mesh

    Read pts file and create rectangular mesh of resolution MxN such that
    it covers all points specified in pts file.

    FIXME: This may be a temporary function until we decide on
    netcdf formats etc

    FIXME: Uses elevation hardwired
    """

    import  mesh_factory
    from load_mesh.loadASCII import import_points_file
    
    if verbose: print 'Read pts'
    points_dict = import_points_file(pts_name)
    #points, attributes = util.read_xya(pts_name)

    #Reduce number of points a bit
    points = points_dict['pointlist'][::reduction]
    elevation = points_dict['attributelist']['elevation']  #Must be elevation
    elevation = elevation[::reduction]

    if verbose: print 'Got %d data points' %len(points)

    if verbose: print 'Create mesh'
    #Find extent
    max_x = min_x = points[0][0]
    max_y = min_y = points[0][1]
    for point in points[1:]:
        x = point[0]
        if x > max_x: max_x = x
        if x < min_x: min_x = x
        y = point[1]
        if y > max_y: max_y = y
        if y < min_y: min_y = y

    #Create appropriate mesh
    vertex_coordinates, triangles, boundary =\
    	 mesh_factory.rectangular(M, N, max_x-min_x, max_y-min_y,
	 			(min_x, min_y))

    #Fit attributes to mesh
    vertex_attributes = fit_to_mesh(vertex_coordinates,
                	triangles,
                	points,
                	elevation, alpha=alpha, verbose=verbose)



    return vertex_coordinates, triangles, boundary, vertex_attributes


def _interpolation(*args, **kwargs):
    """Private function for use with caching. Reason is that classes
    may change their byte code between runs which is annoying.
    """
    
    return Interpolation(*args, **kwargs)


class Interpolation:

    def __init__(self,
                 vertex_coordinates,
                 triangles,
                 point_coordinates = None,
                 alpha = None,
                 verbose = False,
                 expand_search = True,
                 interp_only = False,
                 max_points_per_cell = 30,
                 mesh_origin = None,
                 data_origin = None,
                 precrop = False):


        """ Build interpolation matrix mapping from
        function values at vertices to function values at data points

        Inputs:

          vertex_coordinates: List of coordinate pairs [xi, eta] of
	  points constituting mesh (or a an m x 2 Numeric array)
          Points may appear multiple times
          (e.g. if vertices have discontinuities)

          triangles: List of 3-tuples (or a Numeric array) of
          integers representing indices of all vertices in the mesh.

          point_coordinates: List of coordinate pairs [x, y] of
	  data points (or an nx2 Numeric array)
	  If point_coordinates is absent, only smoothing matrix will
	  be built

          alpha: Smoothing parameter

          data_origin and mesh_origin are 3-tuples consisting of
          UTM zone, easting and northing. If specified
          point coordinates and vertex coordinates are assumed to be
          relative to their respective origins.

        """
        #Convert input to Numeric arrays
        triangles = ensure_numeric(triangles, Int)
        vertex_coordinates = ensure_numeric(vertex_coordinates, Float)

        #Build underlying mesh
        if verbose: print 'Building mesh'
        #self.mesh = General_mesh(vertex_coordinates, triangles,
        #FIXME: Trying the normal mesh while testing precrop,
        #       The functionality of boundary_polygon is needed for that

        #FIXME - geo ref does not have to go into mesh.
        # Change the point co-ords to conform to the
        # mesh co-ords early in the code
        if mesh_origin is None:
            geo = None
        else:
            geo = Geo_reference(mesh_origin[0],mesh_origin[1],mesh_origin[2])

                
        self.mesh = Mesh(vertex_coordinates, triangles,
                         geo_reference = geo)

        if verbose: print 'Checking mesh integrity'
        self.mesh.check_integrity()

        if verbose: print 'Mesh integrity checked'
        
        self.data_origin = data_origin

        self.point_indices = None

        #Smoothing parameter
        if alpha is None:
            self.alpha = DEFAULT_ALPHA
        else:    
            self.alpha = alpha


        if point_coordinates is not None:
            if verbose: print 'Building interpolation matrix'
            self.build_interpolation_matrix_A(point_coordinates,
                                              verbose = verbose,
                                              expand_search = expand_search,
                                              interp_only = interp_only, 
                                              max_points_per_cell =\
                                              max_points_per_cell,
                                              data_origin = data_origin,
                                              precrop = precrop)
        #Build coefficient matrices
        if interp_only == False:
            self.build_coefficient_matrix_B(point_coordinates,
                                        verbose = verbose,
                                        expand_search = expand_search,
                                        max_points_per_cell =\
                                        max_points_per_cell,
                                        data_origin = data_origin,
                                        precrop = precrop)
        if verbose: print 'Finished interpolation'

    def set_point_coordinates(self, point_coordinates,
                              data_origin = None,
                              verbose = False,
                              precrop = True):
        """
        A public interface to setting the point co-ordinates.
        """
        if point_coordinates is not None:
            if verbose: print 'Building interpolation matrix'
            self.build_interpolation_matrix_A(point_coordinates,
                                              verbose = verbose,
                                              data_origin = data_origin,
                                              precrop = precrop)
        self.build_coefficient_matrix_B(point_coordinates, data_origin)

    def build_coefficient_matrix_B(self, point_coordinates=None,
                                   verbose = False, expand_search = True,
                                   max_points_per_cell=30,
                                   data_origin = None,
                                   precrop = False):
        """Build final coefficient matrix"""


        if self.alpha <> 0:
            if verbose: print 'Building smoothing matrix'
            self.build_smoothing_matrix_D()

        if point_coordinates is not None:
            if self.alpha <> 0:
                self.B = self.AtA + self.alpha*self.D
            else:
                self.B = self.AtA

            #Convert self.B matrix to CSR format for faster matrix vector
            self.B = Sparse_CSR(self.B)

    def build_interpolation_matrix_A(self, point_coordinates,
                                     verbose = False, expand_search = True,
                                     max_points_per_cell=30,
                                     data_origin = None,
                                     precrop = False,
                                     interp_only = False):
        """Build n x m interpolation matrix, where
        n is the number of data points and
        m is the number of basis functions phi_k (one per vertex)

        This algorithm uses a quad tree data structure for fast binning of data points
        origin is a 3-tuple consisting of UTM zone, easting and northing.
        If specified coordinates are assumed to be relative to this origin.

        This one will override any data_origin that may be specified in
        interpolation instance

        """



        #FIXME (Ole): Check that this function is memeory efficient.
        #6 million datapoints and 300000 basis functions
        #causes out-of-memory situation
        #First thing to check is whether there is room for self.A and self.AtA
        #
        #Maybe we need some sort of blocking

        from anuga.abstract_2d_finite_volumes.quad import build_quadtree
        from anuga.utilities.polygon import inside_polygon
	

        if data_origin is None:
            data_origin = self.data_origin #Use the one from
                                           #interpolation instance

        #Convert input to Numeric arrays just in case.
        point_coordinates = ensure_numeric(point_coordinates, Float)

        #Keep track of discarded points (if any).
        #This is only registered if precrop is True
        self.cropped_points = False

        #Shift data points to same origin as mesh (if specified)

        #FIXME this will shift if there was no geo_ref.
        #But all this should be removed anyhow.
        #change coords before this point
        mesh_origin = self.mesh.geo_reference.get_origin()
        if point_coordinates is not None:
            if data_origin is not None:
                if mesh_origin is not None:

                    #Transformation:
                    #
                    #Let x_0 be the reference point of the point coordinates
                    #and xi_0 the reference point of the mesh.
                    #
                    #A point coordinate (x + x_0) is then made relative
                    #to xi_0 by
                    #
                    # x_new = x + x_0 - xi_0
                    #
                    #and similarly for eta

                    x_offset = data_origin[1] - mesh_origin[1]
                    y_offset = data_origin[2] - mesh_origin[2]
                else: #Shift back to a zero origin
                    x_offset = data_origin[1]
                    y_offset = data_origin[2]

                point_coordinates[:,0] += x_offset
                point_coordinates[:,1] += y_offset
            else:
                if mesh_origin is not None:
                    #Use mesh origin for data points
                    point_coordinates[:,0] -= mesh_origin[1]
                    point_coordinates[:,1] -= mesh_origin[2]



        #Remove points falling outside mesh boundary
        #This reduced one example from 1356 seconds to 825 seconds

        
        if precrop is True:
            from Numeric import take

            if verbose: print 'Getting boundary polygon'
            P = self.mesh.get_boundary_polygon()

            if verbose: print 'Getting indices inside mesh boundary'
            indices = inside_polygon(point_coordinates, P, verbose = verbose)


            if len(indices) != point_coordinates.shape[0]:
                self.cropped_points = True
                if verbose:
                    print 'Done - %d points outside mesh have been cropped.'\
                          %(point_coordinates.shape[0] - len(indices))

            point_coordinates = take(point_coordinates, indices)
            self.point_indices = indices




        #Build n x m interpolation matrix
        m = self.mesh.coordinates.shape[0] #Nbr of basis functions (1/vertex)
        n = point_coordinates.shape[0]     #Nbr of data points

        if verbose: print 'Number of datapoints: %d' %n
        if verbose: print 'Number of basis functions: %d' %m

        #FIXME (Ole): We should use CSR here since mat-mat mult is now OK.
        #However, Sparse_CSR does not have the same methods as Sparse yet
        #The tests will reveal what needs to be done

        #
        #self.A = Sparse_CSR(Sparse(n,m))
        #self.AtA = Sparse_CSR(Sparse(m,m))
        self.A = Sparse(n,m)
        self.AtA = Sparse(m,m)

        #Build quad tree of vertices (FIXME: Is this the right spot for that?)
        root = build_quadtree(self.mesh,
                              max_points_per_cell = max_points_per_cell)
        #root.show()
        self.expanded_quad_searches = []
        #Compute matrix elements
        for i in range(n):
            #For each data_coordinate point

            if verbose and i%((n+10)/10)==0: print 'Doing %d of %d' %(i, n)
            x = point_coordinates[i]

            #Find vertices near x
            candidate_vertices = root.search(x[0], x[1])
            is_more_elements = True


            element_found, sigma0, sigma1, sigma2, k = \
                           self.search_triangles_of_vertices(candidate_vertices, x)
            first_expansion = True
            while not element_found and is_more_elements and expand_search:
                #if verbose: print 'Expanding search'
                if first_expansion == True:
                    self.expanded_quad_searches.append(1)
                    first_expansion = False
                else:
                    end = len(self.expanded_quad_searches) - 1
                    assert end >= 0
                    self.expanded_quad_searches[end] += 1
                candidate_vertices, branch = root.expand_search()
                if branch == []:
                    # Searching all the verts from the root cell that haven't
                    # been searched.  This is the last try
                    element_found, sigma0, sigma1, sigma2, k = \
                      self.search_triangles_of_vertices(candidate_vertices, x)
                    is_more_elements = False
                else:
                    element_found, sigma0, sigma1, sigma2, k = \
                      self.search_triangles_of_vertices(candidate_vertices, x)

                
	    #Update interpolation matrix A if necessary
            if element_found is True:
                #Assign values to matrix A

                j0 = self.mesh.triangles[k,0] #Global vertex id for sigma0
                j1 = self.mesh.triangles[k,1] #Global vertex id for sigma1
                j2 = self.mesh.triangles[k,2] #Global vertex id for sigma2

                sigmas = {j0:sigma0, j1:sigma1, j2:sigma2}
                js     = [j0,j1,j2]

                for j in js:
                    self.A[i,j] = sigmas[j]
                    for k in js:
                        if interp_only == False:
                            self.AtA[j,k] += sigmas[j]*sigmas[k]
            else:
                pass
                #Ok if there is no triangle for datapoint
		#(as in brute force version)
                #raise 'Could not find triangle for point', x



    def search_triangles_of_vertices(self, candidate_vertices, x):

            
            #Find triangle containing x:
            element_found = False

            # This will be returned if element_found = False
            sigma2 = -10.0
            sigma0 = -10.0
            sigma1 = -10.0
            k = -10.0
            #print "*$* candidate_vertices", candidate_vertices
	    #For all vertices in same cell as point x
            for v in candidate_vertices:
                #FIXME (DSG-DSG): this catches verts with no triangle.
                #Currently pmesh is producing these.
                #this should be stopped, 
                if self.mesh.vertexlist[v] is None:
                    continue
                #for each triangle id (k) which has v as a vertex
		for k, _ in self.mesh.vertexlist[v]:

                    #Get the three vertex_points of candidate triangle
                    xi0 = self.mesh.get_vertex_coordinate(k, 0)
                    xi1 = self.mesh.get_vertex_coordinate(k, 1)
                    xi2 = self.mesh.get_vertex_coordinate(k, 2)

                    #print "PDSG - k", k
                    #print "PDSG - xi0", xi0
                    #print "PDSG - xi1", xi1
                    #print "PDSG - xi2", xi2
                    #print "PDSG element %i verts((%f, %f),(%f, %f),(%f, %f))"\
                    #   % (k, xi0[0], xi0[1], xi1[0], xi1[1], xi2[0], xi2[1])

                    #Get the three normals
                    n0 = self.mesh.get_normal(k, 0)
                    n1 = self.mesh.get_normal(k, 1)
                    n2 = self.mesh.get_normal(k, 2)


                    #Compute interpolation
                    sigma2 = dot((x-xi0), n2)/dot((xi2-xi0), n2)
                    sigma0 = dot((x-xi1), n0)/dot((xi0-xi1), n0)
                    sigma1 = dot((x-xi2), n1)/dot((xi1-xi2), n1)

                    #print "PDSG - sigma0", sigma0
                    #print "PDSG - sigma1", sigma1
                    #print "PDSG - sigma2", sigma2

                    #FIXME: Maybe move out to test or something
                    epsilon = 1.0e-6
                    assert abs(sigma0 + sigma1 + sigma2 - 1.0) < epsilon

                    #Check that this triangle contains the data point

                    #Sigmas can get negative within
                    #machine precision on some machines (e.g nautilus)
                    #Hence the small eps
                    eps = 1.0e-15
                    if sigma0 >= -eps and sigma1 >= -eps and sigma2 >= -eps:
                        element_found = True
                        break

                if element_found is True:
                    #Don't look for any other triangle
                    break
            return element_found, sigma0, sigma1, sigma2, k



    def build_interpolation_matrix_A_brute(self, point_coordinates):
        """Build n x m interpolation matrix, where
        n is the number of data points and
        m is the number of basis functions phi_k (one per vertex)

        This is the brute force which is too slow for large problems,
	but could be used for testing
        """


        #Convert input to Numeric arrays
        point_coordinates = ensure_numeric(point_coordinates, Float)

        #Build n x m interpolation matrix
        m = self.mesh.coordinates.shape[0] #Nbr of basis functions (1/vertex)
        n = point_coordinates.shape[0]     #Nbr of data points

        self.A = Sparse(n,m)
        self.AtA = Sparse(m,m)

        #Compute matrix elements
        for i in range(n):
            #For each data_coordinate point

            x = point_coordinates[i]
            element_found = False
            k = 0
            while not element_found and k < len(self.mesh):
                #For each triangle (brute force)
                #FIXME: Real algorithm should only visit relevant triangles

                #Get the three vertex_points
                xi0 = self.mesh.get_vertex_coordinate(k, 0)
                xi1 = self.mesh.get_vertex_coordinate(k, 1)
                xi2 = self.mesh.get_vertex_coordinate(k, 2)

                #Get the three normals
                n0 = self.mesh.get_normal(k, 0)
                n1 = self.mesh.get_normal(k, 1)
                n2 = self.mesh.get_normal(k, 2)

                #Compute interpolation
                sigma2 = dot((x-xi0), n2)/dot((xi2-xi0), n2)
                sigma0 = dot((x-xi1), n0)/dot((xi0-xi1), n0)
                sigma1 = dot((x-xi2), n1)/dot((xi1-xi2), n1)

                #FIXME: Maybe move out to test or something
                epsilon = 1.0e-6
                assert abs(sigma0 + sigma1 + sigma2 - 1.0) < epsilon

                #Check that this triangle contains data point
                if sigma0 >= 0 and sigma1 >= 0 and sigma2 >= 0:
                    element_found = True
                    #Assign values to matrix A

                    j0 = self.mesh.triangles[k,0] #Global vertex id
                    #self.A[i, j0] = sigma0

                    j1 = self.mesh.triangles[k,1] #Global vertex id
                    #self.A[i, j1] = sigma1

                    j2 = self.mesh.triangles[k,2] #Global vertex id
                    #self.A[i, j2] = sigma2

                    sigmas = {j0:sigma0, j1:sigma1, j2:sigma2}
                    js     = [j0,j1,j2]

                    for j in js:
                        self.A[i,j] = sigmas[j]
                        for k in js:
                            self.AtA[j,k] += sigmas[j]*sigmas[k]
                k = k+1



    def get_A(self):
        return self.A.todense()

    def get_B(self):
        return self.B.todense()

    def get_D(self):
        return self.D.todense()

        #FIXME: Remember to re-introduce the 1/n factor in the
        #interpolation term

    def build_smoothing_matrix_D(self):
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

        #FIXME: algorithm might be optimised by computing local 9x9
        #"element stiffness matrices:

        m = self.mesh.coordinates.shape[0] #Nbr of basis functions (1/vertex)

        self.D = Sparse(m,m)

        #For each triangle compute contributions to D = D1+D2
        for i in range(len(self.mesh)):

            #Get area
            area = self.mesh.areas[i]

            #Get global vertex indices
            v0 = self.mesh.triangles[i,0]
            v1 = self.mesh.triangles[i,1]
            v2 = self.mesh.triangles[i,2]

            #Get the three vertex_points
            xi0 = self.mesh.get_vertex_coordinate(i, 0)
            xi1 = self.mesh.get_vertex_coordinate(i, 1)
            xi2 = self.mesh.get_vertex_coordinate(i, 2)

            #Compute gradients for each vertex
            a0, b0 = gradient(xi0[0], xi0[1], xi1[0], xi1[1], xi2[0], xi2[1],
                              1, 0, 0)

            a1, b1 = gradient(xi0[0], xi0[1], xi1[0], xi1[1], xi2[0], xi2[1],
                              0, 1, 0)

            a2, b2 = gradient(xi0[0], xi0[1], xi1[0], xi1[1], xi2[0], xi2[1],
                              0, 0, 1)

            #Compute diagonal contributions
            self.D[v0,v0] += (a0*a0 + b0*b0)*area
            self.D[v1,v1] += (a1*a1 + b1*b1)*area
            self.D[v2,v2] += (a2*a2 + b2*b2)*area

            #Compute contributions for basis functions sharing edges
            e01 = (a0*a1 + b0*b1)*area
            self.D[v0,v1] += e01
            self.D[v1,v0] += e01

            e12 = (a1*a2 + b1*b2)*area
            self.D[v1,v2] += e12
            self.D[v2,v1] += e12

            e20 = (a2*a0 + b2*b0)*area
            self.D[v2,v0] += e20
            self.D[v0,v2] += e20


    def fit(self, z):
        """Fit a smooth surface to given 1d array of data points z.

        The smooth surface is computed at each vertex in the underlying
        mesh using the formula given in the module doc string.

        Pre Condition:
          self.A, self.AtA and self.B have been initialised

        Inputs:
          z: Single 1d vector or array of data at the point_coordinates.
        """
        
        #Convert input to Numeric arrays
        z = ensure_numeric(z, Float)

        if len(z.shape) > 1 :
            raise VectorShapeError, 'Can only deal with 1d data vector'

        if self.point_indices is not None:
            #Remove values for any points that were outside mesh
            z = take(z, self.point_indices)

        #Compute right hand side based on data
        #FIXME (DSG-DsG): could Sparse_CSR be used here?  Use this format
        # after a matrix is built, before calcs.
        Atz = self.A.trans_mult(z)


        #Check sanity
        n, m = self.A.shape
        if n<m and self.alpha == 0.0:
            msg = 'ERROR (least_squares): Too few data points\n'
            msg += 'There are only %d data points and alpha == 0. ' %n
	    msg += 'Need at least %d\n' %m
            msg += 'Alternatively, set smoothing parameter alpha to a small '
	    msg += 'positive value,\ne.g. 1.0e-3.'
            raise Exception(msg)



        return conjugate_gradient(self.B, Atz, Atz, imax=2*len(Atz) )
        #FIXME: Should we store the result here for later use? (ON)


    def fit_points(self, z, verbose=False):
        """Like fit, but more robust when each point has two or more attributes
	FIXME (Ole): The name fit_points doesn't carry any meaning
	for me. How about something like fit_multiple or fit_columns?
        """

        try:
            if verbose: print 'Solving penalised least_squares problem'
            return self.fit(z)
        except VectorShapeError, e:
            # broadcasting is not supported.

            #Convert input to Numeric arrays
            z = ensure_numeric(z, Float)

            #Build n x m interpolation matrix
            m = self.mesh.coordinates.shape[0] #Number of vertices
            n = z.shape[1]                     #Number of data points

            f = zeros((m,n), Float) #Resulting columns

            for i in range(z.shape[1]):
                f[:,i] = self.fit(z[:,i])

            return f


    def interpolate(self, f):
        """Evaluate smooth surface f at data points implied in self.A.

        The mesh values representing a smooth surface are
        assumed to be specified in f. This argument could,
        for example have been obtained from the method self.fit()

        Pre Condition:
          self.A has been initialised

        Inputs:
          f: Vector or array of data at the mesh vertices.
          If f is an array, interpolation will be done for each column as
          per underlying matrix-matrix multiplication

	Output:
	  Interpolated values at data points implied in self.A

        """
        print "obsolete in least_squares, use fit_interpolate.interpolate"
        return self.A * f

    def cull_outsiders(self, f):
        pass




class Interpolation_function:
    """Interpolation_function - creates callable object f(t, id) or f(t,x,y)
    which is interpolated from time series defined at vertices of
    triangular mesh (such as those stored in sww files)

    Let m be the number of vertices, n the number of triangles
    and p the number of timesteps. 

    Mandatory input
        time:               px1 array of monotonously increasing times (Float)
        quantities:         Dictionary of arrays or 1 array (Float) 
                            The arrays must either have dimensions pxm or mx1.
                            The resulting function will be time dependent in
                            the former case while it will be constant with
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
        self.quantity_names = quantity_names        
        self.quantities = quantities        
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
            interpol = Interpolation(vertex_coordinates,
                                     triangles,
                                     point_coordinates = \
                                     self.interpolation_points,
                                     alpha = 0,
                                     precrop = False, 
                                     verbose = verbose)

            if verbose: print 'Interpolate'
	    for i, t in enumerate(self.time):
                #Interpolate quantities at this timestep
                if verbose and i%((p+10)/10)==0:
                    print ' time step %d of %d' %(i, p)
                    
                for name in quantity_names:
                    if len(quantities[name].shape) == 2:
                        result = interpol.interpolate(quantities[name][i,:])
                    else:
                       #Assume no time dependency 
                       result = interpol.interpolate(quantities[name][:])
                       
                    self.precomputed_values[name][i, :] = result
                   
                        

            #Report
            if verbose:
                print self.statistics()
                #self.print_statistics()
            
        else:
            #Store quantitites as is
	    for name in quantity_names:
                self.precomputed_values[name] = quantities[name]


        #else:
        #    #Return an average, making this a time series
	#    for name in quantity_names:
        #        self.values[name] = zeros(len(self.time), Float)
        #
        #    if verbose: print 'Compute mean values'
	#    for i, t in enumerate(self.time):
        #        if verbose: print ' time step %d of %d' %(i, len(self.time))
        #        for name in quantity_names:
	#	    self.values[name][i] = mean(quantities[name][i,:])




    def __repr__(self):
        #return 'Interpolation function (spatio-temporal)'
        return self.statistics()
    

    def __call__(self, t, point_id = None, x = None, y = None):
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
        from anuga.utilities.numerical_tools import mean        

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


        msg = 'Time interval [%s:%s]' %(self.time[0], self.time[-1])
        msg += ' does not match model time: %s\n' %t
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
                    if ratio > 0: Q1 = Q[self.index+1, point_id]

            #Linear temporal interpolation    
            if ratio > 0:
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
        quantities = self.quantities
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
            q = quantities[name][:].flat
            str += '    %s in [%f, %f]\n' %(name, min(q), max(q))

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

        #FIXME: Delete
        #print '------------------------------------------------'
        #print 'Interpolation_function statistics:'
        #print '  Extent:'
        #print '    x in [%f, %f], len(x) == %d'\
        #      %(min(x), max(x), len(x))
        #print '    y in [%f, %f], len(y) == %d'\
        #      %(min(y), max(y), len(y))
        #print '    t in [%f, %f], len(t) == %d'\
        #      %(min(self.time), max(self.time), len(self.time))
        #print '  Quantities:'
        #for name in quantity_names:
        #    q = quantities[name][:].flat
        #    print '    %s in [%f, %f]' %(name, min(q), max(q))
        #print '  Interpolation points (xi, eta):'\
        #      ' number of points == %d ' %interpolation_points.shape[0]
        #print '    xi in [%f, %f]' %(min(interpolation_points[:,0]),
        #                             max(interpolation_points[:,0]))
        #print '    eta in [%f, %f]' %(min(interpolation_points[:,1]),
        #                              max(interpolation_points[:,1]))
        #print '  Interpolated quantities (over all timesteps):'
        #
        #for name in quantity_names:
        #    q = precomputed_values[name][:].flat
        #    print '    %s at interpolation points in [%f, %f]'\
        #          %(name, min(q), max(q))
        #print '------------------------------------------------'


#-------------------------------------------------------------
if __name__ == "__main__":
    """
    Load in a mesh and data points with attributes.
    Fit the attributes to the mesh.
    Save a new mesh file.
    """
    import os, sys
    usage = "usage: %s mesh_input.tsh point.xya mesh_output.tsh [expand|no_expand][vervose|non_verbose] [alpha] [display_errors|no_display_errors]"\
            %os.path.basename(sys.argv[0])

    if len(sys.argv) < 4:
        print usage
    else:
        mesh_file = sys.argv[1]
        point_file = sys.argv[2]
        mesh_output_file = sys.argv[3]

        expand_search = False
        if len(sys.argv) > 4:
            if sys.argv[4][0] == "e" or sys.argv[4][0] == "E":
                expand_search = True
            else:
                expand_search = False

        verbose = False
        if len(sys.argv) > 5:
            if sys.argv[5][0] == "n" or sys.argv[5][0] == "N":
                verbose = False
            else:
                verbose = True

        if len(sys.argv) > 6:
            alpha = sys.argv[6]
        else:
            alpha = DEFAULT_ALPHA

        # This is used more for testing
        if len(sys.argv) > 7:
            if sys.argv[7][0] == "n" or sys.argv[5][0] == "N":
                display_errors = False
            else:
                display_errors = True
            
        t0 = time.time()
        try:
            fit_to_mesh_file(mesh_file,
                         point_file,
                         mesh_output_file,
                         alpha,
                         verbose= verbose,
                         expand_search = expand_search,
                         display_errors = display_errors)
        except IOError,e:
            import sys; sys.exit(1)

        print 'That took %.2f seconds' %(time.time()-t0)

