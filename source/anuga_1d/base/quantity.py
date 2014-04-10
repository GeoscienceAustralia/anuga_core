"""
Class Quantity - Implements values at each 1d element

To create:

   Quantity(domain, vertex_values)

   domain: Associated domain structure. Required.

   vertex_values: N x 2 array of values at each vertex for each element.
                  Default None

   If vertex_values are None Create array of zeros compatible with domain.
   Otherwise check that it is compatible with dimenions of domain.
   Otherwise raise an exception
"""

import numpy



class Quantity:

    
    def __init__(self, domain, vertex_values=None):
        #Initialise Quantity using optional vertex values.
        
        from anuga_1d.base.generic_domain import Generic_domain

        msg = 'First argument in Quantity.__init__ '
        msg += 'must be of class Generic_domain (or a subclass thereof)'
        assert isinstance(domain, Generic_domain), msg

        if vertex_values is None:
            N = domain.number_of_elements
            self.vertex_values = numpy.zeros((N, 2), numpy.float)
        else:
            self.vertex_values = numpy.array(vertex_values, numpy.float)

            N, V = self.vertex_values.shape
            assert V == 2,\
                   'Two vertex values per element must be specified'


            msg = 'Number of vertex values (%d) must be consistent with'\
                  %N
            msg += 'number of elements in specified domain (%d).'\
                   %domain.number_of_elements

            assert N == domain.number_of_elements, msg

        self.domain = domain

        #Allocate space for other quantities
        self.centroid_values = numpy.zeros(N, numpy.float)
        self.centroid_backup_values = numpy.zeros(N, numpy.float)
        #self.edge_values = numpy.zeros((N, 2), numpy.float)
        #edge values are values of the ends of each interval
      
        #Intialise centroid values
        self.interpolate()


        
        #Allocate space for boundary values
        #L = len(domain.boundary)
        self.boundary_values = numpy.zeros(2, numpy.float) #assumes no parrellism

        #Allocate space for updates of conserved quantities by
        #flux calculations and forcing functions
        
        self.N = domain.number_of_elements
        N = self.N
        self.explicit_update = numpy.zeros(N, numpy.float )
        self.semi_implicit_update = numpy.zeros(N, numpy.float )

        self.gradients = numpy.zeros(N, numpy.float)
        self.qmax = numpy.zeros(self.centroid_values.shape, numpy.float)
        self.qmin = numpy.zeros(self.centroid_values.shape, numpy.float)

        #These are taken from domain but can be set for each quantity
        # if so desired
        self.beta = 0.0
        self.limiter = 'vanleer'


        self.beta_p = numpy.zeros(N,numpy.float)
        self.beta_m = numpy.zeros(N,numpy.float)
        self.beta_x = numpy.zeros(N,numpy.float)


        self.dx  = numpy.zeros((N,2), numpy.float)
        self.phi = numpy.zeros(N, numpy.float)


    def __len__(self):
        """
        Returns number of intervals.
        """
        return self.centroid_values.shape[0]

    def __neg__(self):
        """Negate all values in this quantity giving meaning to the
        expression -Q where Q is an instance of class Quantity
        """

        Q = Quantity(self.domain)
        Q.set_values_from_numeric(-self.vertex_values)
        return Q



    def __add__(self, other):
        """Add to self anything that could populate a quantity

        E.g other can be a constant, an array, a function, another quantity
        (except for a filename or points, attributes (for now))
        - see set_values_from_numeric for details
        """

        Q = Quantity(self.domain)
        Q.set_values_from_numeric(other)

        result = Quantity(self.domain)
        result.set_values_from_numeric(self.vertex_values + Q.vertex_values)
        return result

    def __radd__(self, other):
        """Handle cases like 7+Q, where Q is an instance of class Quantity
        """

        return self + other

    def __sub__(self, other):
        return self + -other            # Invoke self.__neg__()

    def __mul__(self, other):
        """Multiply self with anything that could populate a quantity

        E.g other can be a constant, an array, a function, another quantity
        (except for a filename or points, attributes (for now))
        - see set_values_from_numeric for details
        """

        if isinstance(other, Quantity):
            Q = other
        else:
            Q = Quantity(self.domain)
            Q.set_values_from_numeric(other)

        result = Quantity(self.domain)

        # The product of vertex_values, edge_values and centroid_values
        # are calculated and assigned directly without using
        # set_values_from_numeric (which calls interpolate). Otherwise
        # centroid values wouldn't be products from q1 and q2
        result.vertex_values = self.vertex_values * Q.vertex_values
        result.centroid_values = self.centroid_values * Q.centroid_values

        return result

    def __rmul__(self, other):
        """Handle cases like 3*Q, where Q is an instance of class Quantity
        """

        return self * other

    def __div__(self, other):
        """Divide self with anything that could populate a quantity

        E.g other can be a constant, an array, a function, another quantity
        (except for a filename or points, attributes (for now))
        - see set_values_from_numeric for details

        Zero division is dealt with by adding an epsilon to the divisore
        FIXME (Ole): Replace this with native INF once we migrate to NumPy
        """

        if isinstance(other, Quantity):
            Q = other
        else:
            Q = Quantity(self.domain)
            Q.set_values_from_numeric(other)

        result = Quantity(self.domain)

        # The quotient of vertex_values, edge_values and centroid_values
        # are calculated and assigned directly without using
        # set_values_from_numeric (which calls interpolate). Otherwise
        # centroid values wouldn't be quotient of q1 and q2
        result.vertex_values = self.vertex_values/(Q.vertex_values + epsilon)
        result.centroid_values = self.centroid_values/(Q.centroid_values + epsilon)

        return result

    def __rdiv__(self, other):
        """Handle cases like 3/Q, where Q is an instance of class Quantity
        """

        return self / other

    def __pow__(self, other):
        """Raise quantity to (numerical) power

        As with __mul__ vertex values are processed entry by entry
        while centroid and edge values are re-interpolated.

        Example using __pow__:
          Q = (Q1**2 + Q2**2)**0.5
        """

        if isinstance(other, Quantity):
            Q = other
        else:
            Q = Quantity(self.domain)
            Q.set_values_from_numeric(other)

        result = Quantity(self.domain)

        # The power of vertex_values, edge_values and centroid_values
        # are calculated and assigned directly without using
        # set_values_from_numeric (which calls interpolate). Otherwise
        # centroid values wouldn't be correct
        result.vertex_values = self.vertex_values ** other
        result.centroid_values = self.centroid_values ** other

        return result

    def set_values_from_numeric(self, numeric):
        
        x = numpy.array([1.0,2.0])
        y = [1.0,2.0]

        if type(numeric) == type(y):
            self.set_values_from_array(numeric)
        elif type(numeric) == type(x):
            self.set_values_from_array(numeric)
        elif callable(numeric):
            self.set_values_from_function(numeric)
        elif isinstance(numeric, Quantity):
            self.set_values_from_quantity(numeric)
        else:   # see if it's coercible to a float (float, int or long, etc)
            try:
                numeric = float(numeric)
            except ValueError:
                msg = ("Illegal type for variable 'numeric': %s" % type(numeric))
                raise Exception(msg)
            self.set_values_from_constant(numeric)

    def set_values_from_constant(self,numeric):

        self.vertex_values[:,:]   = numeric
        self.centroid_values[:,] = numeric


    def set_values_from_array(self,numeric):

        self.vertex_values[:,:]   = numeric
        self.interpolate()


    def set_values_from_quantity(self,quantity):

        self.vertex_values[:,:]   = quantity.vertex_values
        self.centroid_values[:,] = quantity.centroid_values

    def set_values_from_function(self,function):

        self.vertex_values[:,:]   = map(function, self.domain.vertices)
        self.centroid_values[:,] = map(function, self.domain.centroids)


    def interpolate(self):
        """
        Compute interpolated values at centroid
        Pre-condition: vertex_values have been set
        """

        N = self.vertex_values.shape[0]
        for i in range(N):
            v0 = self.vertex_values[i, 0]
            v1 = self.vertex_values[i, 1]

            self.centroid_values[i] = (v0 + v1)/2.0




    def set_values(self, X, location='vertices'):
        """Set values for quantity

        X: Compatible list, Numeric array (see below), constant or function
        location: Where values are to be stored.
                  Permissible options are: vertices, centroid
                  Default is "vertices"

        In case of location == 'centroid' the dimension values must
        be a list of a Numerical array of length N, N being the number
        of elements in the mesh. Otherwise it must be of dimension Nx2

        The values will be stored in elements following their
        internal ordering.

        If values are described a function, it will be evaluated at specified points

        If selected location is vertices, values for centroid and edges
        will be assigned interpolated values.
        In any other case, only values for the specified locations
        will be assigned and the others will be left undefined.
        """

        if location not in ['vertices', 'centroids', 'unique_vertices']:
            msg = 'Invalid location: %s, (possible choices vertices, centroids, unique_vertices)' %location
            raise Exception(msg)

        if X is None:
            msg = 'Given values are None'
            raise Exception(msg)

        import types

        if callable(X):
            #Use function specific method
            self.set_function_values(X, location)
          
        elif type(X) in [types.FloatType, types.IntType, types.LongType]:
            
            self.centroid_values[:,] = float(X)
            self.vertex_values[:,:] = float(X)

        elif isinstance(X, Quantity):
            self.set_array_values(X.vertex_values, location)

        else:
            #Use array specific method
            self.set_array_values(X, location)

        if location == 'vertices' or location == 'unique_vertices':
            #Intialise centroid
            self.interpolate()

        if location == 'centroid':
            self.extrapolate_first_order()





    def set_function_values(self, f, location='vertices'):
        """Set values for quantity using specified function

        f: x -> z Function where x and z are arrays
        location: Where values are to be stored.
                  Permissible options are: vertices, centroid
                  Default is "vertices"
        """
        
        if location == 'centroids':
         
            P = self.domain.centroids
            self.set_values(f(P), location)
        else:
            #Vertices
            
            P = self.domain.get_vertices()
            
            for i in range(2):
               
                self.vertex_values[:,i] = f(P[:,i])
               
    def set_array_values(self, values, location='vertices'):
        """Set values for quantity

        values: Numeric array
        location: Where values are to be stored.
                  Permissible options are: vertices, centroid, unique_vertices
                  Default is "vertices"

        In case of location == 'centroid' the dimension values must
        be a list of a Numerical array of length N, N being the number
        of elements in the mesh. If location == 'unique_vertices' the
        dimension values must be a list or a Numeric array of length N+1.
        Otherwise it must be of dimension Nx2

        The values will be stored in elements following their
        internal ordering.

        If selected location is vertices, values for centroid
        will be assigned interpolated values.
        In any other case, only values for the specified locations
        will be assigned and the others will be left undefined.
        """


        values = numpy.array(values).astype(numpy.float)

        N = self.centroid_values.shape[0]


        if location == 'centroids':
            msg = 'Number of values must match number of elements'
            assert values.shape[0] == N, msg
            assert len(values.shape) == 1, 'Values array must be 1d'

            for i in range(N):
                self.centroid_values[i] = values[i]

            self.vertex_values[:,0] = values
            self.vertex_values[:,1] = values
 
        if location == 'vertices':
            msg = 'Number of values must match number of elements'
            assert values.shape[0] == N, msg
            assert len(values.shape) == 2, 'Values array must be 2d'
            msg = 'Array must be N x 2'
            assert values.shape[1] == 2, msg

            self.vertex_values[:,:] = values[:,:]

        if location == 'unique_vertices':
            msg = 'Number of values must match number of elements +1'
            assert values.shape[0] == N+1, msg
            assert len(values.shape) == 1, 'Values array must be 1d'

            self.vertex_values[:,0] = values[0:-1]
            self.vertex_values[:,1] = values[1:N+1]  


            


    def get_values(self, location='vertices', indices = None):
        """get values for quantity

        return X, Compatible list, Numeric array (see below)
        location: Where values are to be stored.
                  Permissible options are: vertices, centroid
                  and unique vertices. Default is 'vertices'

        In case of location == 'centroids' the dimension values must
        be a list of a Numerical array of length N, N being the number
        of elements. Otherwise it must be of dimension Nx3

        The returned values with be a list the length of indices
        (N if indices = None).  Each value will be a list of the three
        vertex values for this quantity.

        Indices is the set of element ids that the operation applies to.

        """


        if location not in ['vertices', 'centroids', 'unique vertices']:
            msg = 'Invalid location: %s' %location
            raise Exception(msg)

        import types, numpy
        assert type(indices) in [types.ListType, types.NoneType,
                                 numpy.array],\
                                 'Indices must be a list or None'

        if location == 'centroids':
            if (indices ==  None):
                indices = range(len(self))
            return numpy.take(self.centroid_values, indices)
        elif location == 'unique vertices':
            if (indices ==  None):
                indices=range(self.domain.coordinates.shape[0])
            vert_values = []
            #Go through list of unique vertices
            for unique_vert_id in indices:
                cells = self.domain.vertexlist[unique_vert_id]

                #In case there are unused points
                if cells is None:
                    msg = 'Unique vertex not associated with cells'
                    raise Exception(msg)

                # Go through all cells, vertex pairs
                # Average the values
                sum = 0
                for cell_id, vertex_id in cells:
                    sum += self.vertex_values[cell_id, vertex_id]
                vert_values.append(sum/len(cells))
            return numpy.array(vert_values)
        else:
            if (indices ==  None):
                indices = range(len(self))
            return numpy.take(self.vertex_values,indices)


    def get_vertex_values(self,
                          x=True,
                          smooth = None,
                          precision = None,
                          reduction = None):
        """Return vertex values like an OBJ format

        The vertex values are returned as one sequence in the 1D float array A.
        If requested the coordinates will be returned in 1D arrays X.

        The connectivity is represented as an integer array, V, of dimension
        M x 2, where M is the number of volumes. Each row has two indices
        into the X, A arrays defining the element.

        if smooth is True, vertex values corresponding to one common
        coordinate set will be smoothed according to the given
        reduction operator. In this case vertex coordinates will be
        de-duplicated.

        If no smoothings is required, vertex coordinates and values will
        be aggregated as a concatenation of values at
        vertices 0, vertices 1


        Calling convention
        if x is True:
           X,A,V = get_vertex_values
        else:
           A,V = get_vertex_values

        """




        if smooth is None:
            smooth = self.domain.smooth

        if precision is None:
            precision = numpy.float

        if reduction is None:
            reduction = self.domain.reduction

        #Create connectivity

        if smooth == True:

            V = self.domain.get_vertices()
            N = len(self.domain.vertexlist)
            #N = len(self.domain.vertices)
            A = numpy.zeros(N, precision)

            #Smoothing loop
            for k in range(N):
                L = self.domain.vertexlist[k]
                #L = self.domain.vertices[k]

                #Go through all triangle, vertex pairs
                #contributing to vertex k and register vertex value

                if L is None: continue #In case there are unused points

                contributions = []
                for volume_id, vertex_id in L:
                    v = self.vertex_values[volume_id, vertex_id]
                    contributions.append(v)

                A[k] = reduction(contributions)

            if x is True:
                 #X = self.domain.coordinates[:,0].astype(precision)
                 X = self.domain.coordinates[:].astype(precision)
                 #Y = self.domain.coordinates[:,1].astype(precision)

                 #return X, Y, A, V
                 return X, A, V
            
            #else:
            return A, V
        else:
            #Don't smooth
            #obj machinery moved to general_mesh

            # Create a V like [[0 1 2], [3 4 5]....[3*m-2 3*m-1 3*m]]
            # These vert_id's will relate to the verts created below
            #m = len(self.domain)  #Number of volumes
            #M = 3*m        #Total number of unique vertices
            #V = reshape(array(range(M)).astype(Int), (m,3))

            #V = self.domain.get_triangles(obj=True)
            V = self.domain.get_vertices
            #FIXME use get_vertices, when ready

            A = self.vertex_values.flat

            #Do vertex coordinates
            if x is True:
                X = self.domain.get_vertex_coordinates()

                #X = C[:,0:6:2].copy()
                #Y = C[:,1:6:2].copy()

                return X.flat, A, V
            else:
                return A, V

    def get_integral(self):
        """Compute the integral of quantity across entire domain
        """
        integral = 0
        for k in range(self.domain.number_of_elements):
            area = self.domain.areas[k]
            qc = self.centroid_values[k]
            integral += qc*area

        return integral

    def get_beta(self,beta):
        """Get limiting parameter
        """
        
        return self.beta

    def set_beta(self,beta):
        """Set limiting parameter
        """

        #Probably should test that it is not too large
        self.beta = beta


    def get_limiter(self):
        return self.limiter

    def set_limiter(self,limiter):

        possible_limiters = \
        ['pyvolution', 'minmod_steve', 'minmod', 'minmod_kurganov', 'superbee', 'vanleer', 'vanalbada']

        if limiter in possible_limiters:
            self.limiter = limiter
            return

        msg = '%s is an incorrect limiter type.\n'% limiter
        msg += 'Possible types are: '+ ", ".join(["%s" % el for el in possible_limiters])
        raise Exception, msg

    def update(self, timestep):
        """Update centroid values based on values stored in
        explicit_update and semi_implicit_update as well as given timestep
        """


        #Explicit updates
        self.centroid_values += timestep*self.explicit_update
        
        #Semi implicit updates
        denominator = 1.0-timestep*self.semi_implicit_update

#        if sum(numpy.equal(denominator, 0.0)) > 0.0:
#            msg = 'Zero division in semi implicit update. Call Stephen :-)'
#            raise Exception(msg)
#        else:
#            #Update conserved_quantities from semi implicit updates
#            self.centroid_values /= denominator
#            

        #Update conserved_quantities from semi implicit updates
        self.centroid_values /= denominator


    def compute_gradients(self):
        """Compute gradients of piecewise linear function defined by centroids of
        neighbouring volumes.
        """


        N = self.centroid_values.shape[0]


        G = self.gradients
        Q = self.centroid_values
        X = self.domain.centroids

    	# first element

        k = 0

        #Get data
        k0 = k
        k1 = k+1
        k2 = k+2

        q0 = Q[k0]
        q1 = Q[k1]
        q2 = Q[k2]

        x0 = X[k0] #V0 centroid
        x1 = X[k1] #V1 centroid
        x2 = X[k2]

        #Gradient
        #G[k] = (q1 - q0)/(x1 - x0)

        G[k] = (q1 - q0)*(x2 - x0)*(x2 - x0) - (q2 - q0)*(x1 - x0)*(x1 - x0)
        G[k] /= (x1 - x0)*(x2 - x0)*(x2 - x1)

        #last element
        k = N-1


        k0 = k
        k1 = k-1
        k2 = k-2

        q0 = Q[k0]
        q1 = Q[k1]
        q2 = Q[k2]

        x0 = X[k0] #V0 centroid
        x1 = X[k1] #V1 centroid
        x2 = X[k2]

        #Gradient
        #G[k] = (q1 - q0)/(x1 - x0)

        G[k] = (q1 - q0)*(x2 - x0)*(x2 - x0) - (q2 - q0)*(x1 - x0)*(x1 - x0)
        G[k] /= (x1 - x0)*(x2 - x0)*(x2 - x1)


        #Interior Volume (2 neighbours)


        q0 = Q[0:-2]
        q1 = Q[1:-1]
        q2 = Q[2:]

        x0 = X[0:-2] #V0 centroid
        x1 = X[1:-1]  #V1 centroid (Self)
        x2 = X[2:] #V2 centroid

        #Gradient
        #G[k] = (q2-q0)/(x2-x0)
        G[1:-1] = ((q0-q1)/(x0-x1)*(x2-x1) - (q2-q1)/(x2-x1)*(x0-x1))/(x2-x0)



    def compute_minmod_gradients(self):
        """Compute gradients of piecewise linear function defined by centroids of
        neighbouring volumes.
        """

        #print 'compute_minmod_gradients'
        from numpy import sign

        
        def xmin(a,b):
            from numpy import sign, minimum
            return 0.5*(sign(a)+sign(b))*minimum(abs(a),abs(b))

        def xmic(t,a,b):
            return xmin(t*xmin(a,b), 0.50*(a+b) )



        N = self.centroid_values.shape[0]


        G = self.gradients
        Q = self.centroid_values
        X = self.domain.centroids

        #-----------------
        #first element
        #-----------------
        k = 0

        k0 = k
        k1 = k+1
        k2 = k+2

        q0 = Q[k0]
        q1 = Q[k1]
        q2 = Q[k2]

        x0 = X[k0] #V0 centroid
        x1 = X[k1] #V1 centroid
        x2 = X[k2]

        #Gradient
        #G[k] = (q1 - q0)/(x1 - x0)

        G[k] = (q1 - q0)*(x2 - x0)*(x2 - x0) - (q2 - q0)*(x1 - x0)*(x1 - x0)
        G[k] /= (x1 - x0)*(x2 - x0)*(x2 - x1)

        #-------------------
        # Last element
        #-------------------
        k = N-1

        k0 = k
        k1 = k-1
        k2 = k-2

        q0 = Q[k0]
        q1 = Q[k1]
        q2 = Q[k2]

        x0 = X[k0] #V0 centroid
        x1 = X[k1] #V1 centroid
        x2 = X[k2]

        #Gradient
        #G[k] = (q1 - q0)/(x1 - x0)

        G[k] = (q1 - q0)*(x2 - x0)*(x2 - x0) - (q2 - q0)*(x1 - x0)*(x1 - x0)
        G[k] /= (x1 - x0)*(x2 - x0)*(x2 - x1)



        #------------------------------
        #Interior Volume (2 neighbours)
        #------------------------------

        q0 = Q[0:-2]
        q1 = Q[1:-1]
        q2 = Q[2:]

        x0 = X[0:-2] #V0 centroid
        x1 = X[1:-1]  #V1 centroid (Self)
        x2 = X[2:] #V2 centroid

        # assuming uniform grid
        d1 = (q1 - q0)/(x1-x0)
        d2 = (q2 - q1)/(x2-x1)

        #Gradient
        G[1:-1] = xmic( self.beta, d1, d2 )      

        

    def extrapolate_first_order(self):
        """Extrapolate conserved quantities from centroid to
        vertices for each volume using
        first order scheme.
        """

        qc = self.centroid_values
        qv = self.vertex_values

        for i in range(2):
            qv[:,i] = qc


    def extrapolate_second_order(self):
        """Extrapolate conserved quantities from centroid to
        vertices for each volume using
        second order scheme.
        """

        if self.limiter == "pyvolution":
            self.limit_pyvolution()

        elif self.limiter == "minmod_steve":
            self.limit_minmod()

        else:
            self.limit_range()
        




    def find_qmax_qmin(self):
        """ Find min and max of this and neighbour's centroid values"""

        from numpy import maximum, minimum

        qmax = self.qmax
        qmin = self.qmin

        qc = self.centroid_values
        
        qmax[:] = qc
        qmin[:] = qc
        
        # check left
        qmax[1:] = maximum(qmax[1:], qc[0:-1])
        qmin[1:] = minimum(qmin[1:], qc[0:-1])
        
        # check right
        qmax[0:-1] = maximum(qmax[0:-1], qc[1:])
        qmin[0:-1] = minimum(qmin[0:-1], qc[1:])



#        for k in range(N):
#            qmax[k] = qmin[k] = qc[k]
#            for i in range(2):
#                n = self.domain.neighbours[k,i]
#                if n >= 0:
#                    qn = qc[n] #Neighbour's centroid value
#
#                    qmin[k] = min(qmin[k], qn)
#                    qmax[k] = max(qmax[k], qn)



    def limit_minmod(self):
        #Z = self.gradients
        #print "gradients 1",Z
        self.compute_minmod_gradients()
        #print "gradients 2", Z

        G = self.gradients
        V = self.domain.vertices
        qc = self.centroid_values
        qv = self.vertex_values        

        x = self.domain.centroids

        x0 = V[:,0]
        x1 = V[:,1]

        #Extrapolate
        qv[:,0] = qc + G*(x0-x)
        qv[:,1] = qc + G*(x1-x)

#        #Check each triangle
#        for k in range(self.domain.number_of_elements):
#            #Centroid coordinates
#            x = self.domain.centroids[k]
#
#            #vertex coordinates
#            x0, x1 = V[k,:]
#
#            #Extrapolate
#            qv[k,0] = qc[k] + G[k]*(x0-x)
#            qv[k,1] = qc[k] + G[k]*(x1-x)

 
    def limit_pyvolution(self):
        """
        Limit slopes for each volume to eliminate artificial variance
        introduced by e.g. second order extrapolator

        This is an unsophisticated limiter as it does not take into
        account dependencies among quantities.

        precondition:
        vertex values are estimated from gradient
        postcondition:
        vertex values are updated
        """


        N = self.domain.number_of_elements
        beta = self.domain.beta
        #beta = 0.8

        self.compute_gradients()


        G = self.gradients
        V = self.domain.vertices
        C = self.domain.centroids
        qc = self.centroid_values
        qv = self.vertex_values

        V0 = V[:,0]
        V1 = V[:,1]

        # Extrapolate to vertices
        qv[:,0] = qc + G*(V0-C)
        qv[:,1] = qc + G*(V1-C)


        # Find max and min values
        self.find_qmax_qmin()

        qmax = self.qmax
        qmin = self.qmin

        #Diffences between centroids and maxima/minima
        dqmax = qmax - qc
        dqmin = qmin - qc

        #Deltas between vertex and centroid values
        dq = numpy.zeros(qv.shape, numpy.float)

        dq[:,0] = qv[:,0] - qc
        dq[:,1] = qv[:,1] - qc

        phi = numpy.ones(qc.shape, numpy.float)

        r0 = numpy.where(dq[:,0]>0.0,dqmax/dq[:,0],1.0)
        r0 = numpy.where(dq[:,0]<0.0,dqmin/dq[:,0],r0)

        r1 = numpy.where(dq[:,1]>0.0,dqmax/dq[:,1],1.0)
        r1 = numpy.where(dq[:,1]<0.0,dqmin/dq[:,1],r1)

        phi = numpy.min(r0*beta,numpy.min(r1*beta,1.0))

        qv[:,0] = qc + phi*dq[:,0]
        qv[:,1] = qc + phi*dq[:,1]

#        #Phi limiter
#        for k in range(N):
#
#            #Find the gradient limiter (phi) across vertices
#            phi = 1.0
#            for i in range(2):
#                r = 1.0
#                if (dq[k,i] > 0): r = dqmax[k]/dq[k,i]
#                if (dq[k,i] < 0): r = dqmin[k]/dq[k,i]
#
#                phi = min( min(r*beta, 1), phi )
#
#            #Then update using phi limiter
#            for i in range(2):
#                qv[k,i] = qc[k] + phi*dq[k,i]

    def limit_range(self):
        import sys

        from limiters_python import minmod, minmod_kurganov, minmod_kurganov_old, maxmod, vanleer, vanalbada

        limiter = self.get_limiter()
        #print limiter
        
        #print 'limit_range'
        N = self.N
        qc = self.centroid_values
        qv = self.vertex_values
        xc = self.domain.centroids
        x0 = self.domain.vertices[:,0]
        x1 = self.domain.vertices[:,1]

        beta_p = self.beta_p
        beta_m = self.beta_m
        beta_x = self.beta_x
        phi = self.phi
        dx  = self.dx


        beta_p[1:]  = (qc[1:]-qc[:-1])/(xc[1:]-xc[:-1])
        beta_m[:-1] = beta_p[1:]
        beta_x[1:-1] = (qc[2:]-qc[:-2])/(xc[2:]-xc[:-2])

        dx[:,0] = x0 - xc
        dx[:,1] = x1 - xc

        phi[0] = (qc[1] - qc[0])/(xc[1] - xc[0])
        phi[N-1] = (qc[N-1] - qc[N-2])/(xc[N-1] - xc[N-2])


        if limiter == "minmod":
            phi[1:-1] = minmod(beta_p[1:-1],beta_m[1:-1])

        elif limiter == "vanleer":
            phi[1:-1] = vanleer(beta_p[1:-1],beta_m[1:-1])
            
        elif limiter == "vanalbada":
            phi[1:-1] = vanalbada(beta_p[1:-1],beta_m[1:-1])

        elif limiter == "minmod_kurganov":
            theta = self.beta
            phi[1:-1] = minmod_kurganov(theta*beta_p[1:-1],theta*beta_m[1:-1], beta_x[1:-1])

        elif limiter == "superbee":
            slope1 = minmod(beta_m[1:-1],2.0*beta_p[1:-1])
            slope2 = minmod(2.0*beta_m[1:-1],beta_p[1:-1])
            phi[1:-1] = maxmod(slope1,slope2)

        else:
            msg = 'Unknown limiter'
            raise Exception, msg


        
        qv[:,0] = qc + phi*dx[:,0]
        qv[:,1] = qc + phi*dx[:,1]




    def limit_steve_slope(self):    

        import sys

        from util import minmod, minmod_kurganov, maxmod, vanleer

        N = self.domain.number_of_elements
        limiter = self.domain.limiter
        limiter_type = self.domain.limiter_type
            
        qc = self.centroid_values
        qv = self.vertex_values

        #Find min and max of this and neighbour's centroid values
        beta_p = numpy.zeros(N,numpy.float)
        beta_m = numpy.zeros(N,numpy.float)
        beta_x = numpy.zeros(N,numpy.float)
        C = self.domain.centroids
        X = self.domain.vertices

        for k in range(N):
        
            n0 = self.domain.neighbours[k,0]
            n1 = self.domain.neighbours[k,1]
            
            if (n0 >= 0) & (n1 >= 0):
                # Check denominator not zero
                if (qc[k+1]-qc[k]) == 0.0:
                    beta_p[k] = float(sys.maxint)
                    beta_m[k] = float(sys.maxint)
                else:
                    #STEVE LIMIT
                    beta_p[k] = (qc[k]-qc[k-1])/(qc[k+1]-qc[k])
                    beta_m[k] = (qc[k+2]-qc[k+1])/(qc[k+1]-qc[k])

        #Deltas between vertex and centroid values
        dq = numpy.zeros(qv.shape, numpy.float)
        for i in range(2):
            dq[:,i] =self.domain.vertices[:,i]-self.domain.centroids
            
        #Phi limiter
        for k in range(N):
                
            phi = 0.0
            if limiter == "flux_minmod":
                #FLUX MINMOD
                phi = minmod_kurganov(1.0,beta_m[k],beta_p[k])
            elif limiter == "flux_superbee":
                #FLUX SUPERBEE
                phi = max(0.0,min(1.0,2.0*beta_m[k]),min(2.0,beta_m[k]))+max(0.0,min(1.0,2.0*beta_p[k]),min(2.0,beta_p[k]))-1.0
            elif limiter == "flux_muscl":
                #FLUX MUSCL
                phi = max(0.0,min(2.0,2.0*beta_m[k],2.0*beta_p[k],0.5*(beta_m[k]+beta_p[k])))
            elif limiter == "flux_vanleer":
                #FLUX VAN LEER
                phi = (beta_m[k]+abs(beta_m[k]))/(1.0+abs(beta_m[k]))+(beta_p[k]+abs(beta_p[k]))/(1.0+abs(beta_p[k]))-1.0
                
                #Then update using phi limiter
                n = self.domain.neighbours[k,1]
                if n>=0:
                    #qv[k,0] = qc[k] - 0.5*phi*(qc[k+1]-qc[k])
                    #qv[k,1] = qc[k] + 0.5*phi*(qc[k+1]-qc[k])
                    qv[k,0] = qc[k] + 0.5*phi*(qv[k,0]-qc[k])
                    qv[k,1] = qc[k] + 0.5*phi*(qv[k,1]-qc[k])
                else:
                    qv[k,i] = qc[k]

    def backup_centroid_values(self):
        # Call correct module function
        # (either from this module or C-extension)
        #backup_centroid_values(self)

        self.centroid_backup_values[:,] = (self.centroid_values).astype('f')

    def saxpy_centroid_values(self,a,b):
        # Call correct module function
        # (either from this module or C-extension)
        self.centroid_values[:,] = (a*self.centroid_values + b*self.centroid_backup_values).astype('f')
        



    
if __name__ == "__main__":
    #from domain import Domain
    
    from anuga_1d.base.generic_domain import Generic_domain as Domain


    def newLinePlot(title='Simple Plot'):
        import pylab as g
        g.ion()
        g.hold(False)
        g.title(title)
        g.xlabel('x')
        g.ylabel('y')
        g.show()


    def linePlot(x,y):
        import pylab as g
        g.plot(x.flat,y.flat)
        g.show()


    def closePlots():
        import pylab as g
        g.close('all')

        
    import pylab as g

    points1 = [0.0, 1.0, 2.0, 3.0]
    vertex_values = [[1.0,2.0],[4.0,5.0],[-1.0,2.0]]

    D1 = Domain(points1)

    Q1 = Quantity(D1, vertex_values)

    print Q1.vertex_values
    print Q1.centroid_values

    new_vertex_values = [[2.0,1.0],[3.0,4.0],[-2.0,4.0]]

    Q1.set_values(new_vertex_values)

    print Q1.vertex_values
    print Q1.centroid_values

    new_centroid_values = [20,30,40]
    Q1.set_values(new_centroid_values,'centroids')

    print Q1.vertex_values
    print Q1.centroid_values

    class FunClass:
        def __init__(self,value):
            self.value = value

        def __call__(self,x):
            return self.value*(x**2)


    fun = FunClass(1.0)
    Q1.set_values(fun,'vertices')

    print Q1.vertex_values
    print Q1.centroid_values

    Xc = Q1.domain.vertices
    Qc = Q1.vertex_values
    print Xc
    print Qc

    Qc[1,0] = 3

    Q1.extrapolate_second_order()
    #Q1.limit_minmod()

    newLinePlot('plots')
    linePlot(Xc,Qc)
    raw_input('press return')

    points2 = numpy.arange(10)
    D2 = Domain(points2)

    Q2 = Quantity(D2)
    Q2.set_values(fun,'vertices')
    Xc = Q2.domain.vertices
    Qc = Q2.vertex_values
    linePlot(Xc,Qc)
    raw_input('press return')

    
    Q2.extrapolate_second_order()
    #Q2.limit_minmod()
    Xc = Q2.domain.vertices
    Qc = Q2.vertex_values
    print Q2.centroid_values
    print Qc
    linePlot(Xc,Qc)
    raw_input('press return')


    for i in range(10):
        import pylab as g
        g.hold(True)
        fun = FunClass(i/10.0)
        Q2.set_values(fun,'centroids')
        Q2.extrapolate_second_order()
        #Q2.limit_minmod()
        Qc = Q2.vertex_values
        linePlot(Xc,Qc)
        g.show()
        raw_input('press return')

    raw_input('press return to quit')

    closePlots()
