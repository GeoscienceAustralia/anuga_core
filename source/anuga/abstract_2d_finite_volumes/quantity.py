"""Class Quantity - Implements values at each triangular element

To create:

   Quantity(domain, vertex_values)

   domain: Associated domain structure. Required.

   vertex_values: N x 3 array of values at each vertex for each element.
                  Default None

   If vertex_values are None Create array of zeros compatible with domain.
   Otherwise check that it is compatible with dimenions of domain.
   Otherwise raise an exception
"""

from Numeric import array, zeros, Float, less, concatenate, NewAxis, argmax, allclose
from anuga.utilities.numerical_tools import ensure_numeric

class Quantity:

    def __init__(self, domain, vertex_values=None):

        from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh

        msg = 'First argument in Quantity.__init__ '
        msg += 'must be of class Mesh (or a subclass thereof)'
        assert isinstance(domain, Mesh), msg

        if vertex_values is None:
            N = len(domain) # number_of_elements
            self.vertex_values = zeros((N, 3), Float)
        else:
            self.vertex_values = array(vertex_values).astype(Float)

            N, V = self.vertex_values.shape
            assert V == 3,\
                   'Three vertex values per element must be specified'


            msg = 'Number of vertex values (%d) must be consistent with'\
                  %N
            msg += 'number of elements in specified domain (%d).'\
                   %len(domain)

            assert N == len(domain), msg

        self.domain = domain

        #Allocate space for other quantities
        self.centroid_values = zeros(N, Float)
        self.edge_values = zeros((N, 3), Float)

        #Intialise centroid and edge_values
        self.interpolate()



    #Methods for operator overloading
    def __len__(self):
        return self.centroid_values.shape[0]


    def __neg__(self):
        """Negate all values in this quantity giving meaning to the
        expression -Q where Q is an instance of class Quantity
        """

        Q = Quantity(self.domain)
        Q.set_values(-self.vertex_values)
        return Q


    def __add__(self, other):
        """Add to self anything that could populate a quantity

        E.g other can be a constant, an array, a function, another quantity
        (except for a filename or points, attributes (for now))
        - see set_values for details
        """

        Q = Quantity(self.domain)
        Q.set_values(other)

        result = Quantity(self.domain)
        result.set_values(self.vertex_values + Q.vertex_values)
        return result

    def __radd__(self, other):
        """Handle cases like 7+Q, where Q is an instance of class Quantity
        """
        return self + other


    def __sub__(self, other):
        return self + -other  #Invoke __neg__

    def __mul__(self, other):
        """Multiply self with anything that could populate a quantity

        E.g other can be a constant, an array, a function, another quantity
        (except for a filename or points, attributes (for now))
        - see set_values for details

        Note that if two quantitites q1 and q2 are multiplied,
        vertex values are multiplied entry by entry
        while centroid and edge values are re-interpolated.
        Hence they won't be the product of centroid or edge values
        from q1 and q2.
        """

        Q = Quantity(self.domain)
        Q.set_values(other)

        result = Quantity(self.domain)
        result.set_values(self.vertex_values * Q.vertex_values)
        return result

    def __rmul__(self, other):
        """Handle cases like 3*Q, where Q is an instance of class Quantity
        """
        return self * other

    def __pow__(self, other):
        """Raise quantity to (numerical) power

        As with __mul__ vertex values are processed entry by entry
        while centroid and edge values are re-interpolated.

        Example using __pow__:
          Q = (Q1**2 + Q2**2)**0.5

        """

        result = Quantity(self.domain)
        result.set_values(self.vertex_values**other)
        return result



    def interpolate(self):
        """Compute interpolated values at edges and centroid
        Pre-condition: vertex_values have been set
        """

        N = self.vertex_values.shape[0]
        for i in range(N):
            v0 = self.vertex_values[i, 0]
            v1 = self.vertex_values[i, 1]
            v2 = self.vertex_values[i, 2]

            self.centroid_values[i] = (v0 + v1 + v2)/3

        self.interpolate_from_vertices_to_edges()


    def interpolate_from_vertices_to_edges(self):
        #Call correct module function
        #(either from this module or C-extension)
        interpolate_from_vertices_to_edges(self)




    #New leaner interface to setting values
    def set_values(self,
                   numeric = None,    # List, numeric array or constant
                   quantity = None,   # Another quantity
                   function = None,   # Callable object: f(x,y)
                   geospatial_data = None, #Arbitrary dataset
                   points = None, values = None, data_georef = None, #Input
                   # for fit (obsoleted by use of geo_spatial object)
                   filename = None, attribute_name = None, #Input from file
                   alpha = None,
                   location = 'vertices',
                   indices = None,
                   verbose = False,
                   use_cache = False):

        """Set values for quantity based on different sources.

        numeric:
          Compatible list, Numeric array (see below) or constant.
          If callable it will treated as a function (see below)
          If instance of another Quantity it will be treated as such.
          If geo_spatial object it will be treated as such

        quantity:
          Another quantity (compatible quantity, e.g. obtained as a
          linear combination of quantities)

        function:
          Any callable object that takes two 1d arrays x and y
          each of length N and returns an array also of length N.
          The function will be evaluated at points determined by
          location and indices in the underlying mesh.

        geospatial_data:
          Arbitrary geo spatial dataset in the form of the class
          Geospatial_data. Mesh points are populated using
          fit_interpolate.fit fitting

        points:
          Nx2 array of data points for use with fit_interpolate.fit
          If points are present, an N array of attribute
          values corresponding to
          each data point must be present.

        values:
          If points is specified, values is an array of length N containing
          attribute values for each point.

        data_georef:
          If points is specified, geo_reference applies to each point.

        filename:
          Name of a .pts file containing data points and attributes for
          use with fit_interpolate.fit.

        attribute_name:
          If specified, any array matching that name
          will be used. from file or geospatial_data.
          Otherwise a default will be used.

        alpha:
          Smoothing parameter to be used with fit_interpolate.fit.
          See module fit_interpolate.fit for further details about alpha.
          Alpha will only be used with points, values or filename.
          Otherwise it will be ignored.


        location: Where values are to be stored.
                  Permissible options are: vertices, edges, centroids
                  Default is 'vertices'

                  In case of location == 'centroids' the dimension values must
                  be a list of a Numerical array of length N,
                  N being the number of elements.
                  Otherwise it must be of dimension Nx3


                  The values will be stored in elements following their
                  internal ordering.

                  If location is not 'unique vertices' Indices is the
                  set of element ids that the operation applies to.
                  If location is 'unique vertices' Indices is the set
                  of vertex ids that the operation applies to.

                  If selected location is vertices, values for
                  centroid and edges will be assigned interpolated
                  values.  In any other case, only values for the
                  specified locations will be assigned and the others
                  will be left undefined.

        verbose: True means that output to stdout is generated

        use_cache: True means that caching of intermediate results is
                   attempted for fit_interpolate.fit.




        Exactly one of the arguments
          numeric, quantity, function, points, filename
        must be present.
        """

        from anuga.geospatial_data.geospatial_data import Geospatial_data
        from types import FloatType, IntType, LongType, ListType, NoneType
        from Numeric import ArrayType

        #General input checks
        L = [numeric, quantity, function, geospatial_data, points, filename]
        msg = 'Exactly one of the arguments '+\
              'numeric, quantity, function, geospatial_data, points, '+\
              'or filename must be present.'
        assert L.count(None) == len(L)-1, msg


        if location not in ['vertices', 'centroids', 'edges',
                            'unique vertices']:
            msg = 'Invalid location: %s' %location
            raise msg


        msg = 'Indices must be a list or None'
        assert type(indices) in [ListType, NoneType, ArrayType], msg


        if not(points is None and values is None and data_georef is None):
            from warnings import warn

            msg = 'Using points, values or data_georef with set_quantity '
            msg += 'is obsolete. Please use a Geospatial_data object instead.'
            warn(msg, DeprecationWarning)



        #Determine which 'set_values_from_...' to use

        if numeric is not None:
            if type(numeric) in [FloatType, IntType, LongType]:
                self.set_values_from_constant(numeric,
                                              location, indices, verbose)
            elif type(numeric) in [ArrayType, ListType]:
                self.set_values_from_array(numeric,
                                           location, indices, verbose)
            elif callable(numeric):
                self.set_values_from_function(numeric,
                                              location, indices, verbose)
            elif isinstance(numeric, Quantity):
                self.set_values_from_quantity(numeric,
                                              location, indices, verbose)
            elif isinstance(numeric, Geospatial_data):
                self.set_values_from_geospatial_data(numeric,
                                                     alpha,
                                                     location, indices,
                                                     verbose = verbose,
                                                     use_cache = use_cache)
            else:
                msg = 'Illegal type for argument numeric: %s' %str(numeric)
                raise msg

        elif quantity is not None:
            self.set_values_from_quantity(quantity,
                                          location, indices, verbose)
        elif function is not None:
            msg = 'Argument function must be callable'
            assert callable(function), msg
            self.set_values_from_function(function,
                                          location, indices, verbose)
        elif geospatial_data is not None:
                self.set_values_from_geospatial_data(geospatial_data,
                                                     alpha,
                                                     location, indices,
                                                     verbose = verbose,
                                                     use_cache = use_cache)
        elif points is not None:
            print 'The usage of points in set_values will be deprecated.' +\
                  'Please use the geospatial_data object.'

            msg = 'When points are specified, associated values must also be.'
            assert values is not None, msg
            self.set_values_from_points(points, values, alpha,
                                        location, indices,
                                        data_georef = data_georef,
                                        verbose = verbose,
                                        use_cache = use_cache)
        elif filename is not None:
            self.set_values_from_file(filename, attribute_name, alpha,
                                      location, indices,
                                      verbose = verbose,
                                      use_cache = use_cache)
        else:
            raise 'This can\'t happen :-)'


        #Update all locations in triangles
        if location == 'vertices' or location == 'unique vertices':
            #Intialise centroid and edge_values
            self.interpolate()

        if location == 'centroids':
            #Extrapolate 1st order - to capture notion of area being specified
            self.extrapolate_first_order()



    #Specific functions for setting values
    def set_values_from_constant(self, X,
                                 location, indices, verbose):
        """Set quantity values from specified constant X
        """


        if location == 'centroids':
            if (indices ==  None):
                self.centroid_values[:] = X
            else:
                #Brute force
                for i in indices:
                    self.centroid_values[i,:] = X

        elif location == 'edges':
            if (indices ==  None):
                self.edge_values[:] = X
            else:
                #Brute force
                for i in indices:
                    self.edge_values[i,:] = X

        elif location == 'unique vertices':
            if (indices ==  None):
                self.edge_values[:] = X
            else:

                #Go through list of unique vertices
                for unique_vert_id in indices:
                    triangles = self.domain.vertexlist[unique_vert_id]

                    #In case there are unused points
                    if triangles is None: continue

                    #Go through all triangle, vertex pairs
                    #and set corresponding vertex value
                    for triangle_id, vertex_id in triangles:
                        self.vertex_values[triangle_id, vertex_id] = X

                    #Intialise centroid and edge_values
                    self.interpolate()
        else:
            if (indices ==  None):
                self.vertex_values[:] = X
            else:
                #Brute force
                for i_vertex in indices:
                    self.vertex_values[i_vertex,:] = X






    def set_values_from_array(self, values,
                              location, indices, verbose):
        """Set values for quantity

        values: Numeric array
        location: Where values are to be stored.
        Permissible options are: vertices, edges, centroid, unique vertices
        Default is 'vertices'

        indices - if this action is carried out on a subset of
        elements or unique vertices
        The element/unique vertex indices are specified here.

        In case of location == 'centroid' the dimension values must
        be a list of a Numerical array of length N, N being the number
        of elements.

        Otherwise it must be of dimension Nx3

        The values will be stored in elements following their
        internal ordering.

        If selected location is vertices, values for centroid and edges
        will be assigned interpolated values.
        In any other case, only values for the specified locations
        will be assigned and the others will be left undefined.
        """

        from Numeric import array, Float, Int, allclose

        values = array(values).astype(Float)

        if indices is not None:
            indices = array(indices).astype(Int)
            msg = 'Number of values must match number of indices'
            assert values.shape[0] == indices.shape[0], msg

        N = self.centroid_values.shape[0]

        if location == 'centroids':
            assert len(values.shape) == 1, 'Values array must be 1d'

            if indices is None:
                msg = 'Number of values must match number of elements'
                assert values.shape[0] == N, msg

                self.centroid_values = values
            else:
                msg = 'Number of values must match number of indices'
                assert values.shape[0] == indices.shape[0], msg

                #Brute force
                for i in range(len(indices)):
                    self.centroid_values[indices[i]] = values[i]

        elif location == 'edges':
            assert len(values.shape) == 2, 'Values array must be 2d'

            msg = 'Number of values must match number of elements'
            assert values.shape[0] == N, msg

            msg = 'Array must be N x 3'
            assert values.shape[1] == 3, msg

            self.edge_values = values

        elif location == 'unique vertices':
            assert len(values.shape) == 1 or allclose(values.shape[1:], 1),\
                   'Values array must be 1d'

            self.set_vertex_values(values.flat, indices=indices)
        else:
            if len(values.shape) == 1:
                self.set_vertex_values(values, indices=indices)

            elif len(values.shape) == 2:
                #Vertex values are given as a triplet for each triangle

                msg = 'Array must be N x 3'
                assert values.shape[1] == 3, msg

                if indices == None:
                    self.vertex_values = values
                else:
                    for element_index, value in map(None, indices, values):
                        self.vertex_values[element_index] = value
            else:
                msg = 'Values array must be 1d or 2d'
                raise msg

    def set_values_from_quantity(self, q,
                                 location, indices, verbose):
        """Set quantity values from specified quantity instance q

        Location is ignored
        """


        A = q.vertex_values

        from Numeric import allclose
        msg = 'Quantities are defined on different meshes. '+\
              'This might be a case for implementing interpolation '+\
              'between different meshes.'
        assert allclose(A.shape, self.vertex_values.shape), msg

        self.set_values(A, location='vertices',
                        indices=indices,
                        verbose=verbose)


    def set_values_from_function(self, f,
                                 location, indices, verbose):
        """Set values for quantity using specified function

        f: x, y -> z Function where x, y and z are arrays
        location: Where values are to be stored.
                  Permissible options are: vertices, centroid, edges,
                  unique vertices
                  Default is "vertices"
        """

        #FIXME: Should check that function returns something sensible and
        #raise a meaningfull exception if it returns None for example

        #FIXME: Should supply absolute coordinates

        from Numeric import take

        if (indices is None):
            indices = range(len(self))
            is_subset = False
        else:
            is_subset = True

        if location == 'centroids':
            P = take(self.domain.centroid_coordinates, indices)
            if is_subset:
                self.set_values(f(P[:,0], P[:,1]),
                                location = location,
                                indices = indices)
            else:
                self.set_values(f(P[:,0], P[:,1]), location = location)
        elif location == 'vertices':
            P = self.domain.vertex_coordinates
            if is_subset:
                #Brute force
                for e in indices:
                    for i in range(3):
                        self.vertex_values[e,i] = f(P[e,2*i], P[e,2*i+1])
            else:
                for i in range(3):
                    self.vertex_values[:,i] = f(P[:,2*i], P[:,2*i+1])
        else:
            raise 'Not implemented: %s' %location



    def set_values_from_geospatial_data(self, geospatial_data, alpha,
                                        location, indices,
                                        verbose = False,
                                        use_cache = False):

        #FIXME: Use this function for the time being. Later move code in here

        points = geospatial_data.get_data_points(absolute = False)
        values = geospatial_data.get_attributes()
        data_georef = geospatial_data.get_geo_reference()



        self.set_values_from_points(points, values, alpha,
                                    location, indices,
                                    data_georef = data_georef,
                                    verbose = verbose,
                                    use_cache = use_cache)



    def set_values_from_points(self, points, values, alpha,
                               location, indices,
                               data_georef = None,
                               verbose = False,
                               use_cache = False):
        """
        Set quantity values from arbitray data points using
        fit_interpolate.fit
        """


        from anuga.fit_interpolate.fit import fit_to_mesh
        from anuga.coordinate_transforms.geo_reference import Geo_reference


        points = ensure_numeric(points, Float)
        values = ensure_numeric(values, Float)

        if location != 'vertices':
            msg = 'set_values_from_points is only defined for '+\
                  'location=\'vertices\''
            raise ms

        coordinates = self.domain.coordinates
        triangles = self.domain.triangles


        #Take care of georeferencing
        if data_georef is None:
            data_georef = Geo_reference()


        mesh_georef = self.domain.geo_reference

        #print mesh_georef
        #print data_georef
        #print points


        #Call fit_interpolate.fit function
        args = (coordinates, triangles, points, values)
        kwargs = {'data_origin': data_georef.get_origin(),
                  'mesh_origin': mesh_georef.get_origin(),
                  'alpha': alpha,
                  'verbose': verbose}

        #print kwargs

        if use_cache is True:
            try:
                from caching import cache
            except:
                msg = 'Caching was requested, but caching module'+\
                      'could not be imported'
                raise msg

            vertex_attributes = cache(fit_to_mesh,
                                      args, kwargs,
                                      verbose=verbose,
                                      compression=False)
        else:

            vertex_attributes = apply(fit_to_mesh,
                                      args, kwargs)

        #Call underlying method using array values
        self.set_values_from_array(vertex_attributes,
                                   location, indices, verbose)



    def set_values_from_file(self, filename, attribute_name, alpha,
                             location, indices,
                             verbose = False,
                             use_cache = False):
        """Set quantity based on arbitrary points in .pts file
        using attribute_name selects name of attribute
        present in file.
        If not specified try to use whatever is available in file.
        """

        from load_mesh.loadASCII import import_points_file
        from anuga.geospatial_data.geospatial_data import\
             points_dictionary2geospatial_data

        from types import StringType
        msg = 'Filename must be a text string'
        assert type(filename) == StringType, msg


        # Read from (NetCDF) file
        # FIXME (Ole): This function should really return a
        # Geospatial_data object.
        points_dict = import_points_file(filename)
        points = points_dict['pointlist']
        attributes = points_dict['attributelist']

        if attribute_name is None:
            names = attributes.keys()
            attribute_name = names[0]

        msg = 'Attribute_name must be a text string'
        assert type(attribute_name) == StringType, msg


        if verbose:
            print 'Using attribute %s from file %s' %(attribute_name, filename)
            print 'Available attributes: %s' %(names)

        #try:
        #    z = attributes[attribute_name]
        #except:
        #    msg = 'Could not extract attribute %s from file %s'\
        #          %(attribute_name, filename)
        #    raise msg


        #Take care of georeferencing
        if points_dict.has_key('geo_reference') and \
               points_dict['geo_reference'] is not None:
            data_georef = points_dict['geo_reference']
        else:
            data_georef = None



        #Call underlying method for geospatial data
        geospatial_data = points_dictionary2geospatial_data(points_dict)
        geospatial_data.set_default_attribute_name(attribute_name)

        self.set_values_from_geospatial_data(geospatial_data,
                                             alpha,
                                             location, indices,
                                             verbose = verbose,
                                             use_cache = use_cache)

   
    def get_maximum_index(self, indices=None):
        """Return index for maximum value of quantity (on centroids)

        Optional argument:
            indices is the set of element ids that the operation applies to.

        Usage:
            i = get_maximum_index()

        Notes:
            We do not seek the maximum at vertices as each vertex can
            have multiple values - one for each triangle sharing it.

            If there are multiple cells with same maximum value, the
            first cell encountered in the triangle array is returned.
        """

        V = self.get_values(location='centroids', indices=indices)

        # Always return absolute indices
        i = argmax(V)

        if indices is None:
            return i
        else:
            return indices[i]

        
    def get_maximum_value(self, indices=None):
        """Return maximum value of quantity (on centroids)

        Optional argument:
            indices is the set of element ids that the operation applies to.

        Usage:
            v = get_maximum_value()

        Note, we do not seek the maximum at vertices as each vertex can
        have multiple values - one for each triangle sharing it            
        """


        i = self.get_maximum_index(indices)
        V = self.get_values(location='centroids') #, indices=indices)
        
        return V[i]
        

    def get_maximum_location(self, indices=None):
        """Return location of maximum value of quantity (on centroids)

        Optional argument:
            indices is the set of element ids that the operation applies to.

        Usage:
            x, y = get_maximum_location()


        Notes:
            We do not seek the maximum at vertices as each vertex can
            have multiple values - one for each triangle sharing it.

            If there are multiple cells with same maximum value, the
            first cell encountered in the triangle array is returned.       
        """

        i = self.get_maximum_index(indices)
        x, y = self.domain.get_centroid_coordinates()[i]

        return x, y




    def get_interpolated_values(self, interpolation_points):

        # Interpolation object based on internal (discontinuous triangles)
        x, y, vertex_values, triangles = self.get_vertex_values(xy=True,
                                                                smooth=False)
        # FIXME: This concat should roll into get_vertex_values
        vertex_coordinates = concatenate((x[:, NewAxis], y[:, NewAxis]),
                                         axis=1)

        can_reuse = False
        if hasattr(self, 'interpolation_object'):
            # Reuse to save time
            I = self.interpolation_object

            if allclose(interpolation_points, I._point_coordinates):
                can_reuse = True
                

        if can_reuse is True:
            # Use absence of points to indicate reuse in I.interpolate
            result = I.interpolate(vertex_values) 
        else:    
            from anuga.fit_interpolate.interpolate import Interpolate

            # Create interpolation object with matrix
            I = Interpolate(vertex_coordinates, triangles)
            self.interpolation_object = I

            # Call interpolate with points the first time
            interpolation_points = ensure_numeric(interpolation_points, Float)
            result = I.interpolate(vertex_values, interpolation_points)     

        return result


    def get_values(self, interpolation_points=None,
                   location='vertices',
                   indices = None):
        """get values for quantity

        return X, Compatible list, Numeric array (see below)
        interpolation_points: List of x, y coordinates where value is
        sought (using interpolation). If points are given, values of
        location and indices are ignored
        
        location: Where values are to be stored.
                  Permissible options are: vertices, edges, centroid
                  and unique vertices. Default is 'vertices'


        The returned values with be a list the length of indices
        (N if indices = None).

        In case of location == 'centroids' the dimension of returned
        values will be a list or a Numerical array of length N, N being
        the number of elements.
        
        In case of location == 'vertices' or 'edges' the dimension of
        returned values will be of dimension Nx3

        In case of location == 'unique vertices' the average value at
        each vertex will be returned and the dimension of returned values
        will be a 1d array of length "number of vertices" 
        
        Indices is the set of element ids that the operation applies to.

        The values will be stored in elements following their
        internal ordering.

        """
        from Numeric import take

        if interpolation_points is not None:
            return self.get_interpolated_values(interpolation_points)
        
        

        if location not in ['vertices', 'centroids', 'edges',
                            'unique vertices']:
            msg = 'Invalid location: %s' %location
            raise msg

        import types, Numeric
        assert type(indices) in [types.ListType, types.NoneType,
                                 Numeric.ArrayType],\
                                 'Indices must be a list or None'

        if location == 'centroids':
            if (indices ==  None):
                indices = range(len(self))
            return take(self.centroid_values,indices)
        elif location == 'edges':
            if (indices ==  None):
                indices = range(len(self))
            return take(self.edge_values,indices)
        elif location == 'unique vertices':
            if (indices ==  None):
                indices=range(self.domain.coordinates.shape[0])
            vert_values = []
            #Go through list of unique vertices
            for unique_vert_id in indices:
                triangles = self.domain.vertexlist[unique_vert_id]

                #In case there are unused points
                if triangles is None:
                    msg = 'Unique vertex not associated with triangles'
                    raise msg

                # Go through all triangle, vertex pairs
                # Average the values
                
                # FIXME (Ole): Should we merge this with get_vertex_values
                # and use the concept of a reduction operator?
                sum = 0
                for triangle_id, vertex_id in triangles:
                    sum += self.vertex_values[triangle_id, vertex_id]
                vert_values.append(sum/len(triangles))
            return Numeric.array(vert_values)
        else:
            if (indices ==  None):
                indices = range(len(self))
            return take(self.vertex_values,indices)



    def set_vertex_values(self, A, indices = None):
        """Set vertex values for all unique vertices based on input array A
        which has one entry per unique vertex, i.e.
        one value for each row in array self.domain.coordinates or
        one value for each row in vertexlist.

        indices is the list of vertex_id's that will be set.

        This function is used by set_values_from_array
        """

        from Numeric import array, Float

        #Assert that A can be converted to a Numeric array of appropriate dim
        A = array(A, Float)

        #print 'SHAPE A', A.shape
        assert len(A.shape) == 1

        if indices == None:
            assert A.shape[0] == self.domain.coordinates.shape[0]
            vertex_list = range(A.shape[0])
        else:
            assert A.shape[0] == len(indices)
            vertex_list = indices

        #Go through list of unique vertices
        for i_index, unique_vert_id in enumerate(vertex_list):
            triangles = self.domain.vertexlist[unique_vert_id]

            if triangles is None: continue #In case there are unused points

            #Go through all triangle, vertex pairs
            #touching vertex unique_vert_id and set corresponding vertex value
            for triangle_id, vertex_id in triangles:
                self.vertex_values[triangle_id, vertex_id] = A[i_index]

        #Intialise centroid and edge_values
        self.interpolate()


    def smooth_vertex_values(self, value_array='field_values',
                             precision = None):
        """ Smooths field_values or conserved_quantities data.
        TODO: be able to smooth individual fields
        NOTE:  This function does not have a test.
        FIXME: NOT DONE - do we need it?
        FIXME: this function isn't called by anything.
               Maybe it should be removed..-DSG
        """

        from Numeric import concatenate, zeros, Float, Int, array, reshape


        A,V = self.get_vertex_values(xy=False,
                                     value_array=value_array,
                                     smooth = True,
                                     precision = precision)

        #Set some field values
        for volume in self:
            for i,v in enumerate(volume.vertices):
                if value_array == 'field_values':
                    volume.set_field_values('vertex', i, A[v,:])
                elif value_array == 'conserved_quantities':
                    volume.set_conserved_quantities('vertex', i, A[v,:])

        if value_array == 'field_values':
            self.precompute()
        elif value_array == 'conserved_quantities':
            Volume.interpolate_conserved_quantities()


    #Method for outputting model results
    #FIXME: Split up into geometric and numeric stuff.
    #FIXME: Geometric (X,Y,V) should live in mesh.py
    #FIXME: STill remember to move XY to mesh
    def get_vertex_values(self,
                          xy=True,
                          smooth = None,
                          precision = None,
                          reduction = None):
        """Return vertex values like an OBJ format

        The vertex values are returned as one sequence in the 1D float array A.
        If requested the coordinates will be returned in 1D arrays X and Y.

        The connectivity is represented as an integer array, V, of dimension
        M x 3, where M is the number of volumes. Each row has three indices
        into the X, Y, A arrays defining the triangle.

        if smooth is True, vertex values corresponding to one common
        coordinate set will be smoothed according to the given
        reduction operator. In this case vertex coordinates will be
        de-duplicated.

        If no smoothings is required, vertex coordinates and values will
        be aggregated as a concatenation of values at
        vertices 0, vertices 1 and vertices 2


        Calling convention
        if xy is True:
           X,Y,A,V = get_vertex_values
        else:
           A,V = get_vertex_values

        """

        from Numeric import concatenate, zeros, Float, Int, array, reshape


        if smooth is None:
            smooth = self.domain.smooth

        if precision is None:
            precision = Float

        #Create connectivity

        if smooth == True:
            
            if reduction is None:
                reduction = self.domain.reduction

            V = self.domain.get_vertices()
            N = len(self.domain.vertexlist)
            A = zeros(N, precision)

            #Smoothing loop
            for k in range(N):
                L = self.domain.vertexlist[k]

                #Go through all triangle, vertex pairs
                #contributing to vertex k and register vertex value

                if L is None: continue #In case there are unused points

                contributions = []
                for volume_id, vertex_id in L:
                    v = self.vertex_values[volume_id, vertex_id]
                    contributions.append(v)

                A[k] = reduction(contributions)


            if xy is True:
                X = self.domain.coordinates[:,0].astype(precision)
                Y = self.domain.coordinates[:,1].astype(precision)

                return X, Y, A, V
            else:
                return A, V
        else:
            #Don't smooth
            #obj machinery moved to general_mesh

            # Create a V like [[0 1 2], [3 4 5]....[3*m-2 3*m-1 3*m]]
            # These vert_id's will relate to the verts created below
            #m = len(self.domain)  #Number of volumes
            #M = 3*m        #Total number of unique vertices
            #V = reshape(array(range(M)).astype(Int), (m,3))

            V = self.domain.get_triangles(obj=True)
            #FIXME use get_vertices, when ready

            A = self.vertex_values.flat

            #Do vertex coordinates
            if xy is True:
                C = self.domain.get_vertex_coordinates()

                X = C[:,0:6:2].copy()
                Y = C[:,1:6:2].copy()

                return X.flat, Y.flat, A, V
            else:
                return A, V


    def extrapolate_first_order(self):
        """Extrapolate conserved quantities from centroid to
        vertices for each volume using
        first order scheme.
        """

        qc = self.centroid_values
        qv = self.vertex_values

        for i in range(3):
            qv[:,i] = qc


    def get_integral(self):
        """Compute the integral of quantity across entire domain
        """
        integral = 0
        for k in range(len(self.domain)):
            area = self.domain.areas[k]
            qc = self.centroid_values[k]
            integral += qc*area

        return integral




class Conserved_quantity(Quantity):
    """Class conserved quantity adds to Quantity:

    boundary values, storage and method for updating, and
    methods for (second order) extrapolation from centroid to vertices inluding
    gradients and limiters
    """

    def __init__(self, domain, vertex_values=None):
        Quantity.__init__(self, domain, vertex_values)

        from Numeric import zeros, Float

        #Allocate space for boundary values
        L = len(domain.boundary)
        self.boundary_values = zeros(L, Float)

        #Allocate space for updates of conserved quantities by
        #flux calculations and forcing functions

        N = len(domain) # number_of_triangles
        self.explicit_update = zeros(N, Float )
        self.semi_implicit_update = zeros(N, Float )


    def update(self, timestep):
        #Call correct module function
        #(either from this module or C-extension)
        return update(self, timestep)


    def compute_gradients(self):
        #Call correct module function
        #(either from this module or C-extension)
        return compute_gradients(self)


    def limit(self):
        #Call correct module function
        #(either from this module or C-extension)
        limit(self)


    def extrapolate_second_order(self):
        #Call correct module function
        #(either from this module or C-extension)
        extrapolate_second_order(self)


def update(quantity, timestep):
    """Update centroid values based on values stored in
    explicit_update and semi_implicit_update as well as given timestep

    Function implementing forcing terms must take on argument
    which is the domain and they must update either explicit
    or implicit updates, e,g,:

    def gravity(domain):
        ....
        domain.quantities['xmomentum'].explicit_update = ...
        domain.quantities['ymomentum'].explicit_update = ...



    Explicit terms must have the form

        G(q, t)

    and explicit scheme is

       q^{(n+1}) = q^{(n)} + delta_t G(q^{n}, n delta_t)


    Semi implicit forcing terms are assumed to have the form

       G(q, t) = H(q, t) q

    and the semi implicit scheme will then be

      q^{(n+1}) = q^{(n)} + delta_t H(q^{n}, n delta_t) q^{(n+1})


    """

    from Numeric import sum, equal, ones, exp, Float

    N = quantity.centroid_values.shape[0]


    #Divide H by conserved quantity to obtain G (see docstring above)


    for k in range(N):
        x = quantity.centroid_values[k]
        if x == 0.0:
            #FIXME: Is this right
            quantity.semi_implicit_update[k] = 0.0
        else:
            quantity.semi_implicit_update[k] /= x


    #Semi implicit updates
    denominator = ones(N, Float)-timestep*quantity.semi_implicit_update

    if sum(less(denominator, 1.0)) > 0.0:
        msg = 'denominator < 1.0 in semi implicit update. Call Stephen :-)'
        raise msg

    if sum(equal(denominator, 0.0)) > 0.0:
        msg = 'Zero division in semi implicit update. Call Stephen :-)'
        raise msg
    else:
        #Update conserved_quantities from semi implicit updates
        quantity.centroid_values /= denominator

#    quantity.centroid_values = exp(timestep*quantity.semi_implicit_update)*quantity.centroid_values

    #Explicit updates
    quantity.centroid_values += timestep*quantity.explicit_update

def interpolate_from_vertices_to_edges(quantity):
    """Compute edge values from vertex values using linear interpolation
    """

    for k in range(quantity.vertex_values.shape[0]):
        q0 = quantity.vertex_values[k, 0]
        q1 = quantity.vertex_values[k, 1]
        q2 = quantity.vertex_values[k, 2]

        quantity.edge_values[k, 0] = 0.5*(q1+q2)
        quantity.edge_values[k, 1] = 0.5*(q0+q2)
        quantity.edge_values[k, 2] = 0.5*(q0+q1)



def extrapolate_second_order(quantity):
    """Extrapolate conserved quantities from centroid to
    vertices for each volume using
    second order scheme.
    """

    a, b = quantity.compute_gradients()

    X = quantity.domain.get_vertex_coordinates()
    qc = quantity.centroid_values
    qv = quantity.vertex_values

    #Check each triangle
    for k in range(len(quantity.domain)):
        #Centroid coordinates
        x, y = quantity.domain.centroid_coordinates[k]

        #vertex coordinates
        x0, y0, x1, y1, x2, y2 = X[k,:]

        #Extrapolate
        qv[k,0] = qc[k] + a[k]*(x0-x) + b[k]*(y0-y)
        qv[k,1] = qc[k] + a[k]*(x1-x) + b[k]*(y1-y)
        qv[k,2] = qc[k] + a[k]*(x2-x) + b[k]*(y2-y)


def compute_gradients(quantity):
    """Compute gradients of triangle surfaces defined by centroids of
    neighbouring volumes.
    If one edge is on the boundary, use own centroid as neighbour centroid.
    If two or more are on the boundary, fall back to first order scheme.
    """

    from Numeric import zeros, Float
    from utilitites.numerical_tools import gradient

    centroid_coordinates = quantity.domain.centroid_coordinates
    surrogate_neighbours = quantity.domain.surrogate_neighbours
    centroid_values = quantity.centroid_values
    number_of_boundaries = quantity.domain.number_of_boundaries

    N = centroid_values.shape[0]

    a = zeros(N, Float)
    b = zeros(N, Float)

    for k in range(N):
        if number_of_boundaries[k] < 2:
            #Two or three true neighbours

            #Get indices of neighbours (or self when used as surrogate)
            k0, k1, k2 = surrogate_neighbours[k,:]

            #Get data
            q0 = centroid_values[k0]
            q1 = centroid_values[k1]
            q2 = centroid_values[k2]

            x0, y0 = centroid_coordinates[k0] #V0 centroid
            x1, y1 = centroid_coordinates[k1] #V1 centroid
            x2, y2 = centroid_coordinates[k2] #V2 centroid

            #Gradient
            a[k], b[k] = gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2)

        elif number_of_boundaries[k] == 2:
            #One true neighbour

            #Get index of the one neighbour
            for k0 in surrogate_neighbours[k,:]:
                if k0 != k: break
            assert k0 != k

            k1 = k  #self

            #Get data
            q0 = centroid_values[k0]
            q1 = centroid_values[k1]

            x0, y0 = centroid_coordinates[k0] #V0 centroid
            x1, y1 = centroid_coordinates[k1] #V1 centroid

            #Gradient
            a[k], b[k] = gradient2(x0, y0, x1, y1, q0, q1)
        else:
            #No true neighbours -
            #Fall back to first order scheme
            pass


    return a, b



def limit(quantity):
    """Limit slopes for each volume to eliminate artificial variance
    introduced by e.g. second order extrapolator

    This is an unsophisticated limiter as it does not take into
    account dependencies among quantities.

    precondition:
    vertex values are estimated from gradient
    postcondition:
    vertex values are updated
    """

    from Numeric import zeros, Float

    N = len(quantity.domain)

    beta_w = quantity.domain.beta_w

    qc = quantity.centroid_values
    qv = quantity.vertex_values

    #Find min and max of this and neighbour's centroid values
    qmax = zeros(qc.shape, Float)
    qmin = zeros(qc.shape, Float)

    for k in range(N):
        qmax[k] = qc[k]
        qmin[k] = qc[k]
        for i in range(3):
            n = quantity.domain.neighbours[k,i]
            if n >= 0:
                qn = qc[n] #Neighbour's centroid value

                qmin[k] = min(qmin[k], qn)
                qmax[k] = max(qmax[k], qn)
        qmax[k] = min(qmax[k], 2.0*qc[k])
        qmin[k] = max(qmin[k], 0.5*qc[k])


    #Diffences between centroids and maxima/minima
    dqmax = qmax - qc
    dqmin = qmin - qc

    #Deltas between vertex and centroid values
    dq = zeros(qv.shape, Float)
    for i in range(3):
        dq[:,i] = qv[:,i] - qc

    #Phi limiter
    for k in range(N):

        #Find the gradient limiter (phi) across vertices
        phi = 1.0
        for i in range(3):
            r = 1.0
            if (dq[k,i] > 0): r = dqmax[k]/dq[k,i]
            if (dq[k,i] < 0): r = dqmin[k]/dq[k,i]

            phi = min( min(r*beta_w, 1), phi )

        #Then update using phi limiter
        for i in range(3):
            qv[k,i] = qc[k] + phi*dq[k,i]



from anuga.utilities import compile
if compile.can_use_C_extension('quantity_ext.c'):
    #Replace python version with c implementations

    from quantity_ext import compute_gradients, limit,\
    extrapolate_second_order, interpolate_from_vertices_to_edges, update
