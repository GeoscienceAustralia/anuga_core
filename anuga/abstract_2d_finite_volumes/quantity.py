
"""Class Quantity - Implements values at each triangular element

To create:

   Quantity(domain, vertex_values)

   domain: Associated domain structure. Required.

   vertex_values: N x 3 array of values at each vertex for each element.
                  Default None

   If vertex_values are None Create array of zeros compatible with domain.
   Otherwise check that it is compatible with dimensions of domain.
   Otherwise raise an exception
   
   For Quantities that need to be saved during checkpointing, set register=True. Registered
   Quantities can be found in the dictionary domain.quantities (note, other Quantities can 
   exist). 
"""

import types
import os.path

from anuga.utilities.numerical_tools import ensure_numeric, is_scalar
from anuga.geometry.polygon import inside_polygon
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.fit_interpolate.fit import fit_to_mesh
from anuga.config import points_file_block_line_size as default_block_line_size
from anuga.config import epsilon
from anuga.caching import cache
import anuga.utilities.log as log


import numpy as num


class Quantity:


    counter = 0

    def __init__(self, domain, vertex_values=None, name=None, register=False):
        from anuga.abstract_2d_finite_volumes.generic_domain \
                            import Generic_Domain

        msg = ('First argument in Quantity.__init__() must be of class %s '
               '(or a subclass thereof). I got %s.' 
               % (str(Generic_Domain.__name__),str(domain.__class__)))
        assert isinstance(domain, Generic_Domain), msg

        if vertex_values is None:
            N = len(domain)             # number_of_elements
            self.vertex_values = num.zeros((N, 3), num.float)
        else:
            self.vertex_values = num.array(vertex_values, num.float)

            N, V = self.vertex_values.shape
            assert V == 3, 'Three vertex values per element must be specified'

            msg = 'Number of vertex values (%d) must be consistent with' % N
            msg += 'number of elements in specified domain (%d).' % len(domain)
            assert N == len(domain), msg

        self.domain = domain

        # Allocate space for other quantities
        self.centroid_values = num.zeros(N, num.float)
        self.edge_values = num.zeros((N, 3), num.float)

        # Allocate space for Gradient
        self.x_gradient = num.zeros(N, num.float)
        self.y_gradient = num.zeros(N, num.float)

        # Allocate space for Limiter Phi
        self.phi = num.zeros(N, num.float)

        # Intialise centroid and edge_values
        self.interpolate()

        # Allocate space for boundary values
        #self.boundary_length = domain.boundary_length
        self.boundary_length = L = self.domain.boundary_length
        self.boundary_values = num.zeros(L, num.float)

        # Allocate space for updates of conserved quantities by
        # flux calculations and forcing functions

        # Allocate space for update fields
        self.explicit_update = num.zeros(N, num.float )
        self.semi_implicit_update = num.zeros(N, num.float )
        self.centroid_backup_values = num.zeros(N, num.float)

        self.set_beta(1.0)

        # Keep count of quantities
        Quantity.counter += 1

        self.set_name(name)
        
        # Lets register the quantity with the quantity (useful for checkpointing)
        if register:
            self.domain.quantities[self.name] = self

    ############################################################################
    # Methods for operator overloading
    ############################################################################

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
        return self + -other            # Invoke self.__neg__()

    def __mul__(self, other):
        """Multiply self with anything that could populate a quantity

        E.g other can be a constant, an array, a function, another quantity
        (except for a filename or points, attributes (for now))
        - see set_values for details
        """

        if isinstance(other, Quantity):
            Q = other
        else:
            Q = Quantity(self.domain)
            Q.set_values(other)

        result = Quantity(self.domain)

        # The product of vertex_values, edge_values and centroid_values
        # are calculated and assigned directly without using
        # set_values (which calls interpolate). Otherwise
        # edge and centroid values wouldn't be products from q1 and q2
        result.vertex_values = self.vertex_values * Q.vertex_values
        result.edge_values = self.edge_values * Q.edge_values
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
        - see set_values for details

        Zero division is dealt with by adding an epsilon to the divisore
        FIXME (Ole): Replace this with native INF once we migrate to NumPy
        """

        if isinstance(other, Quantity):
            Q = other
        else:
            Q = Quantity(self.domain)
            Q.set_values(other)

        result = Quantity(self.domain)

        # The quotient of vertex_values, edge_values and centroid_values
        # are calculated and assigned directly without using
        # set_values (which calls interpolate). Otherwise
        # edge and centroid values wouldn't be quotient of q1 and q2
        result.vertex_values = self.vertex_values/(Q.vertex_values + epsilon)
        result.edge_values = self.edge_values/(Q.edge_values + epsilon)
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
            Q.set_values(other)

        result = Quantity(self.domain)

        # The power of vertex_values, edge_values and centroid_values
        # are calculated and assigned directly without using
        # set_values (which calls interpolate). Otherwise
        # edge and centroid values wouldn't be correct
        result.vertex_values = self.vertex_values ** other
        result.edge_values = self.edge_values ** other
        result.centroid_values = self.centroid_values ** other

        return result





    def maximum(self, other):
        """Max of self with anything that could populate a quantity

        E.g other can be a constant, an array, a function, another quantity
        (except for a filename or points, attributes (for now))
        - see set_values for details
        """

        if isinstance(other, Quantity):
            Q = other
        else:
            Q = Quantity(self.domain)
            Q.set_values(other)


        # The maximum of vertex_values, edge_values and centroid_values
        # are calculated and assigned directly without using
        # set_values (which calls interpolate). Otherwise
        # edge and centroid values wouldn't be max from q1 and q2
        self.vertex_values[:] = num.maximum(self.vertex_values, Q.vertex_values)
        self.edge_values[:] = num.maximum(self.edge_values, Q.edge_values)
        self.centroid_values[:] = num.maximum(self.centroid_values, Q.centroid_values)

        return self

        


    def minimum(self, other):
        """Max of self with anything that could populate a quantity

        E.g other can be a constant, an array, a function, another quantity
        (except for a filename or points, attributes (for now))
        - see set_values for details
        """

        if isinstance(other, Quantity):
            Q = other
        else:
            Q = Quantity(self.domain)
            Q.set_values(other)


        # The minimum of vertex_values, edge_values and centroid_values
        # are calculated and assigned directly without using
        # set_values (which calls interpolate). Otherwise
        self.vertex_values[:]= num.minimum(self.vertex_values, Q.vertex_values)
        self.edge_values[:] = num.minimum(self.edge_values, Q.edge_values)
        self.centroid_values[:] = num.minimum(self.centroid_values, Q.centroid_values)

        return self


    ############################################################################
    # Setters/Getters
    ############################################################################
    
    def save_centroid_data_to_csv(self,filename=None):


        #FIXME SR: Should add code to deal with parallel

        ids = self.domain.tri_full_flag == 1
        c_v = self.centroid_values[ids].reshape((-1,1))
        c_x = self.domain.centroid_coordinates[ids,0].reshape((-1,1))
        c_y = self.domain.centroid_coordinates[ids,1].reshape((-1,1))

        import numpy

        c_xyv = numpy.hstack((c_x, c_y, c_v))

        if filename is None:
            filename= self.name

        if self.domain.parallel:
            fullfilename = filename+'_centroid_data_P%g_%g.csv'% \
                    ( self.domain.numproc, self.domain.processor)
        else:
            fullfilename = filename+'_centroid_data.csv'

        numpy.savetxt(fullfilename, c_xyv, delimiter=',', fmt =  ['%.15e', '%.15e', '%.15e' ])


        if self.domain.parallel:
            import pypar
            pypar.barrier()

            # On processor 0 catenate the files
            if self.domain.processor == 0:
                import shutil
                import os
                destination = open(filename+'_centroid_data.csv','wb')
                np = self.domain.numproc
                files = [ filename+'_centroid_data'+"_P"+str(np)+"_"+str(v)+".csv" for v in range(np)]
                for file in files:
                    shutil.copyfileobj(open(file,'rb'), destination)
                destination.close()
                for file in files:
                    os.remove(file)


    def save_data_to_dem(self,filename=None):


        #FIXME SR: Should add code to deal with parallel

        ids = self.domain.tri_full_flag == 1
        c_v = self.centroid_values[ids].reshape((-1,1))
        c_x = self.domain.centroid_coordinates[ids,0].reshape((-1,1))
        c_y = self.domain.centroid_coordinates[ids,1].reshape((-1,1))

        import numpy

        c_xyv = numpy.hstack((c_x, c_y, c_v))

        if filename is None:
            filename= self.name

        if self.domain.parallel:
            fullfilename = filename+'_centroid_data_P%g_%g.csv'% \
                    ( self.domain.numproc, self.domain.processor)
        else:
            fullfilename = filename+'_centroid_data.csv'

        numpy.savetxt(fullfilename, c_xyv, delimiter=',', fmt =  ['%.15e', '%.15e', '%.15e' ])


        if self.domain.parallel:
            import pypar
            pypar.barrier()

            # On processor 0 catenate the files
            if self.domain.processor == 0:
                import shutil
                import os
                destination = open(filename+'_centroid_data.csv','wb')
                np = self.domain.numproc
                files = [ filename+'_centroid_data'+"_P"+str(np)+"_"+str(v)+".csv" for v in range(np)]
                for file in files:
                    shutil.copyfileobj(open(file,'rb'), destination)
                destination.close()
                for file in files:
                    os.remove(file)




    def plot_quantity(self, filename=None, show=True):

        X, Y, A, V = self.get_vertex_values(smooth=True)

        import matplotlib.pyplot as plt
        import numpy as np

        plt.clf()
        plt.tripcolor(X, Y, V, A)
        plt.colorbar()

        if filename is None:
            filename = self.name+'.png'

        plt.savefig(filename)
        
        if show:
            plt.show()
        

    def save_to_array(self,
                cellsize=None,
                NODATA_value=-9999.0,
                smooth=None,
                easting_min=None,
                easting_max=None,
                northing_min=None,
                northing_max=None,
                origin=None,
                verbose=False):
        
        """Interpolate quantity to an array
        """
    

        verbose = False
        
    
        #Get extent and reference
        
        domain = self.domain

        volumes = domain.triangles
        
        #smooth = True

        x,y,a,v= self.get_vertex_values(xy=True, smooth=smooth)
        
        
        false_easting = 500000
        false_northing = 10000000
    
    
    
        geo_ref = self.domain.geo_reference
        
        xllcorner = geo_ref.get_xllcorner()
        yllcorner = geo_ref.get_yllcorner()
        
        if verbose: 
            print 
            print xllcorner
            print yllcorner
            print x
            print y
            print node_coordinates[:,0]
            print node_coordinates[:,1]
            
        
        
        # Create grid and update xll/yll corner and x,y
        # Relative extent
        if easting_min is None:
            xmin = min(x)
        else:
            xmin = easting_min - xllcorner
    
        if easting_max is None:
            xmax = max(x)
        else:
            xmax = easting_max - xllcorner
    
        if northing_min is None:
            ymin = min(y)
        else:
            ymin = northing_min - yllcorner
    
        if northing_max is None:
            ymax = max(y)
        else:
            ymax = northing_max - yllcorner
    
        msg = 'Implementation of Quantity.save_to_array() is not completed'
        #raise Exception, msg
    
        
        msg = 'xmax must be greater than or equal to xmin.\n'
        msg += 'I got xmin = %f, xmax = %f' %(xmin, xmax)
        assert xmax >= xmin, msg
    
        msg = 'ymax must be greater than or equal to xmin.\n'
        msg += 'I got ymin = %f, ymax = %f' %(ymin, ymax)
        assert ymax >= ymin, msg
        
        
        if verbose: log.critical('Creating grid')
        
        xrange = xmax-xmin
        yrange = ymax-ymin
        
        if cellsize == None:
            cellsize = max(xrange,yrange)/10.0
    
 
        ncols = int(xrange/cellsize) + 1
        nrows = int(yrange/cellsize) + 1
    
        # New absolute reference and coordinates
        newxllcorner = xmin + xllcorner
        newyllcorner = ymin + yllcorner
    
        x = x + xllcorner - newxllcorner
        y = y + yllcorner - newyllcorner
    
    
        grid_values = num.zeros( (nrows*ncols, ), num.float)

    
        num_tri =  len(v)
        norms = num.zeros(6*num_tri, num.float)

    
        #Use fast method to calc grid values
        from anuga.file_conversion.calc_grid_values_ext import calc_grid_values
    
        calc_grid_values(nrows, ncols, cellsize, NODATA_value,
                         x,y, norms, v, a, grid_values)
    
    
        y_g = num.arange(nrows)*cellsize + yllcorner - newyllcorner
        x_g = num.arange(ncols)*cellsize + xllcorner - newxllcorner
        #print outside_indices
    
        if verbose:
            log.critical('Interpolated values are in [%f, %f]'
                         % (num.min(grid_values), num.max(grid_values)))
    
    
        return x_g,y_g, grid_values.reshape(nrows,ncols)#[::-1,:]



    def set_name(self, name=None):
        
        if name is not None:
            self.name = name
        else:
            self.name = 'quantity_%g' % Quantity.counter



    def get_name(self):

        return self.name

    
    def set_beta(self, beta):
        """Set default beta value for limiting """

        if beta < 0.0:
            log.critical('WARNING: setting beta < 0.0')
        if beta > 2.0:
            log.critical('WARNING: setting beta > 2.0')

        self.beta = beta

    ##
    # @brief Get the current beta value.
    # @return The current beta value.
    def get_beta(self):
        """Get default beta value for limiting"""

        return self.beta



    ##
    # @brief Set boundary values using a function or array or scalar
    # @param numeric: function or array or scalar
    def set_boundary_values(self, numeric = 0.0):
        """Set boundary values """

        if isinstance(numeric, (list, num.ndarray)):
            self._set_boundary_values_from_array(numeric)
        elif callable(numeric):
            self._set_boundary_values_from_function(numeric)
        else:   # see if it's coercible to a float (float, int or long, etc)
            try:
                numeric = float(numeric)
            except ValueError:
                msg = ("Illegal type for variable 'numeric': %s"
                           % type(numeric))
                raise Exception(msg)
            self._set_boundary_values_from_constant(numeric)

    def set_boundary_values_from_edges(self):
        """Set boundary values by simply extrapolating
        from the cells
        """

#        for j in range(self.boundary_length):
#            vol_id  = self.domain.boundary_cells[j]
#            edge_id = self.domain.boundary_edges[j]
#            self.boundary_values[j] = self.edge_values[vol_id,edge_id]

        vol_id  = self.domain.boundary_cells
        edge_id = self.domain.boundary_edges
        self.boundary_values = (self.edge_values.flat)[3*vol_id+edge_id]

    ##
    # @brief Set boundary values using a function
    # @param numeric: function
    def _set_boundary_values_from_function(self, function):
        """Set boundary values from function of x,y """

        for j in range(self.boundary_length):
            vol_id  = self.domain.boundary_cells[j]
            edge_id = self.domain.boundary_edges[j]
            [x,y] = self.domain.get_edge_midpoint_coordinates(vol_id)[edge_id]
            self.boundary_values[j] = function(x,y)

    ##
    # @brief Set boundary values using a scalar
    # @param numeric: scalar
    def _set_boundary_values_from_constant(self, scalar):
        """Set boundary values from scalar """

        #for (vol_id, edge_id) , j in self.domain.boundary_enumeration.items():
        #[x,y] = self.domain.get_edge_midpoint_coordinates(vol_id)[edge_id]
        self.boundary_values[:] = scalar

    ##
    # @brief Set boundary values using a scalar
    # @param numeric: scalar
    def _set_boundary_values_from_array(self, array):
        """Set boundary values from array """

        assert len(array)  == len(self.domain.boundary_enumeration)
        
        #for (vol_id, edge_id) , j in self.domain.boundary_enumeration.items():
        #[x,y] = self.domain.get_edge_midpoint_coordinates(vol_id)[edge_id]
        self.boundary_values[:] = array

    ##
    # @brief Compute interpolated values at edges and centroid.
    # @note vertex_values must be set before calling this.
    def interpolate_old(self):
        """Compute interpolated values at edges and centroid
        Pre-condition: vertex_values have been set
        """

        # FIXME (Ole): Maybe this function
        # should move to the C-interface?
        # However, it isn't called by validate_all.py, so it
        # may not be that important to optimise it?

#        N = self.vertex_values.shape[0]
#        for i in range(N):
#            v0 = self.vertex_values[i, 0]
#            v1 = self.vertex_values[i, 1]
#            v2 = self.vertex_values[i, 2]
#
#            self.centroid_values[i] = (v0 + v1 + v2)/3

        v0 = self.vertex_values[:, 0]
        v1 = self.vertex_values[:, 1]
        v2 = self.vertex_values[:, 2]

        self.centroid_values[:] = (v0 + v1 + v2)/3

        self.interpolate_from_vertices_to_edges()


    ##
    # @brief Compute interpolated values at edges and centroid.
    # @note vertex_values must be set before calling this.
    def interpolate(self):
        """Compute interpolated values at edges and centroid
        Pre-condition: vertex_values have been set
        """
        from quantity_ext import interpolate
        interpolate(self)


    def interpolate_from_vertices_to_edges(self):
        # Call correct module function (either from this module or C-extension)

        from quantity_ext import interpolate_from_vertices_to_edges
        interpolate_from_vertices_to_edges(self)

    def interpolate_from_edges_to_vertices(self):
        # Call correct module function (either from this module or C-extension)

        from quantity_ext import interpolate_from_edges_to_vertices
        interpolate_from_edges_to_vertices(self)

    #---------------------------------------------
    # Public interface for setting quantity values
    #---------------------------------------------

    def set_values(self, numeric=None,         # List, numeric array or constant
                         quantity=None,        # Another quantity
                         function=None,        # Callable object: f(x,y)
                         geospatial_data=None, # Arbitrary dataset
                         filename=None,
                         attribute_name=None,  # Input from file
                         alpha=None,
                         location='vertices',
                         polygon=None,
                         indices=None,
                         smooth=False,
                         verbose=False,
                         use_cache=False):
        """Set values for quantity based on different sources.

        numeric:
          Compatible list, numeric array (see below) or constant.
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

        filename:
          Name of a points file containing data points and attributes for
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
                  be a list of a numerical array of length N,
                  N being the number of elements.
                  Otherwise it must be of dimension Nx3


                  The values will be stored in elements following their
                  internal ordering.

                  If location is 'unique vertices' indices refers the set
                  of node ids that the operation applies to.
                  If location is not 'unique vertices' indices refers the
                  set of triangle ids that the operation applies to.


                  If selected location is vertices, values for
                  centroid and edges will be assigned interpolated
                  values.  In any other case, only values for the
                  specified locations will be assigned and the others
                  will be left undefined.


        polygon: Restrict update of quantity to locations that fall
                 inside polygon. Polygon works by selecting indices
                 and calling set_values recursively.
                 Polygon mode has only been implemented for
                 constant values so far.

        indices: Restrict update of quantity to locations that are
                 identified by indices (e.g. node ids if location
                 is 'unique vertices' or triangle ids otherwise).

        verbose: True means that output to stdout is generated

        use_cache: True means that caching of intermediate results is
                   attempted for fit_interpolate.fit.




        Exactly one of the arguments
          numeric, quantity, function, filename
        must be present.
        """

        from anuga.geospatial_data.geospatial_data import Geospatial_data

        # Treat special case: Polygon situation
        # Location will be ignored and set to 'centroids'
        # FIXME (Ole): This needs to be generalised and
        # perhaps the notion of location and indices simplified

        # FIXME (Ole): Need to compute indices based on polygon
        # (and location) and use existing code after that.

        # See ticket:275, ticket:250, ticeket:254 for refactoring plan

        if polygon is not None:
            if indices is not None:
                msg = 'Only one of polygon and indices can be specified'
                raise Exception(msg)

            msg = 'With polygon selected, set_quantity must provide '
            msg += 'the keyword numeric and it must (currently) be '
            msg += 'a constant.'
            if numeric is None:
                raise Exception(msg)
            else:
                # Check that numeric is as constant
                assert isinstance(numeric, (float, int, long)), msg

            location = 'centroids'

            points = self.domain.get_centroid_coordinates(absolute=True)
            indices = inside_polygon(points, polygon)

            self.set_values_from_constant(numeric, location, indices, verbose)

            self.extrapolate_first_order()

            if smooth:
                self.smooth_vertex_values(use_cache=use_cache,
                                          verbose=verbose)

            return

        # General input checks
        L = [numeric, quantity, function, geospatial_data, filename]
        msg = ('Exactly one of the arguments numeric, quantity, function, '
               'geospatial_data, or filename must be present.')

        count = 0
        for entry in L:
            if entry is None:
                count = count + 1
                
        assert count == len(L)-1, msg

        if location == 'edges':
            msg = 'edges has been deprecated as valid location'
            raise Exception(msg)

        if location not in ['vertices', 'centroids', 'unique vertices']:
            msg = 'Invalid location: %s' % location
            raise Exception(msg)

        msg = 'Indices must be a list, array or None'
        assert isinstance(indices, (types.NoneType, list, num.ndarray)), msg

        # Determine which 'set_values_from_...' to use
        if numeric is not None:
            if isinstance(numeric, (list, num.ndarray)):
                self.set_values_from_array(numeric, location, indices,
                                           use_cache=use_cache, verbose=verbose)
            elif callable(numeric):
                self.set_values_from_function(numeric, location, indices,
                                              use_cache=use_cache,
                                              verbose=verbose)
            elif isinstance(numeric, Quantity):
                self.set_values_from_quantity(numeric, location, indices,
                                              verbose=verbose)
            elif isinstance(numeric, Geospatial_data):
                self.set_values_from_geospatial_data(numeric, alpha, location,
                                                     indices, verbose=verbose,
                                                     use_cache=use_cache)
            else:   # see if it's coercible to a float (float, int or long, etc)
                try:
                    numeric = float(numeric)
                except ValueError:
                    msg = ("Illegal type for variable 'numeric': %s"
                           % type(numeric))
                    raise Exception(msg)
                self.set_values_from_constant(numeric, location,
                                              indices, verbose)
        elif quantity is not None:
            self.set_values_from_quantity(quantity, location, indices, verbose)
        elif function is not None:
            msg = 'Argument function must be callable'
            assert callable(function), msg
            self.set_values_from_function(function, location, indices,
                                          use_cache=use_cache, verbose=verbose)
        elif geospatial_data is not None:
                self.set_values_from_geospatial_data(geospatial_data, alpha,
                                                     location, indices,
                                                     verbose=verbose,
                                                     use_cache=use_cache)
        elif filename is not None:

            if hasattr(self.domain, 'points_file_block_line_size'):
                max_read_lines = self.domain.points_file_block_line_size
            else:
                max_read_lines = default_block_line_size
            self.set_values_from_file(filename, attribute_name, alpha, location,
                                      indices, verbose=verbose,
                                      max_read_lines=max_read_lines,
                                      use_cache=use_cache)
        else:
            raise Exception("This can't happen :-)")

        # Update all locations in triangles
        if location == 'vertices' or location == 'unique vertices':
            # Intialise centroid and edge_values
            self.interpolate()

        if location == 'centroids':
            # Extrapolate 1st order - to capture notion of area being specified
            self.extrapolate_first_order()

    ############################################################################
    # Specific internal functions for setting values based on type
    ############################################################################

    def set_values_from_constant(self, X, location, indices, verbose):
        """Set quantity values from specified constant X"""

        # FIXME (Ole): Somehow indices refer to centroids
        # rather than vertices as default. See unit test
        # test_set_vertex_values_using_general_interface_with_subset(self):

        if location == 'centroids':
            if indices is None:
                self.centroid_values[:] = X
            else:
                # Brute force
                for i in indices:
                    self.centroid_values[i] = X
        elif location == 'unique vertices':
            if indices is None:
                self.edge_values[:] = X  #FIXME (Ole): Shouldn't this be vertex_values?
            else:
                # Go through list of unique vertices
                for unique_vert_id in indices:
                    triangles = \
                        self.domain.get_triangles_and_vertices_per_node(node=unique_vert_id)

                    # In case there are unused points
                    if len(triangles) == 0:
                        continue

                    # Go through all triangle, vertex pairs
                    # and set corresponding vertex value
                    for triangle_id, vertex_id in triangles:
                        self.vertex_values[triangle_id, vertex_id] = X

                    # Intialise centroid and edge_values
                    self.interpolate()
        else:
            if indices is None:
                self.vertex_values[:] = X
            else:
                # Brute force
                for i_vertex in indices:
                    self.vertex_values[i_vertex] = X

    def set_values_from_array(self, values,
                                    location='vertices',
                                    indices=None,
                                    use_cache=False,
                                    verbose=False):
        """Set values for quantity

        values: numeric array
        location: Where values are to be stored.
        Permissible options are: vertices, centroid, unique vertices
        Default is 'vertices'

        indices - if this action is carried out on a subset of
        elements or unique vertices
        The element/unique vertex indices are specified here.

        In case of location == 'centroid' the dimension values must
        be a list of a numerical array of length N, N being the number
        of elements.

        Otherwise it must be of dimension Nx3

        The values will be stored in elements following their
        internal ordering.

        If selected location is vertices, values for centroid and edges
        will be assigned interpolated values.
        In any other case, only values for the specified locations
        will be assigned and the others will be left undefined.
        """

        values = num.array(values, num.float)

        if indices is not None:
            indices = num.array(indices, num.int)
            msg = ('Number of values must match number of indices: You '
                   'specified %d values and %d indices'
                   % (values.shape[0], indices.shape[0]))
            assert values.shape[0] == indices.shape[0], msg

        N = self.centroid_values.shape[0]

        if location == 'centroids':
            assert len(values.shape) == 1, 'Values array must be 1d'

            if indices is None:
                msg = 'Number of values must match number of elements'
                assert values.shape[0] == N, msg

                self.centroid_values[:] = values
            else:
                msg = 'Number of values must match number of indices'
                assert values.shape[0] == indices.shape[0], msg

                # Brute force
                self.centroid_values[indices] = values
                #for i in range(len(indices)):
                #    self.centroid_values[indices[i]] = values[i]
        elif location == 'unique vertices':
            msg = 'Values array must be 1d'
            assert len(values.shape) == 1 or num.allclose(values.shape[1:], 1), msg

            self.set_vertex_values(values.flat, indices=indices,
                                   use_cache=use_cache, verbose=verbose)
        else:
            # Location vertices
            if len(values.shape) == 1:
                # This is the common case arising from fitted
                # values (e.g. from pts file).
                self.set_vertex_values(values, indices=indices,
                                       use_cache=use_cache, verbose=verbose)
            elif len(values.shape) == 2:
                # Vertex values are given as a triplet for each triangle
                msg = 'Array must be N x 3'
                assert values.shape[1] == 3, msg

                if indices is None:
                    self.vertex_values[:] = values
                else:
                    self.vertex_values[indices] = values
                    #for element_index, value in map(None, indices, values):
                    #    self.vertex_values[element_index] = value
            else:
                msg = 'Values array must be 1d or 2d'
                raise Exception(msg)

    def set_values_from_quantity(self, q, location, indices, verbose):
        """Set quantity values from specified quantity instance q

        Location is ignored - vertices will always be used here.
        """


        A = q.vertex_values

        msg = 'Quantities are defined on different meshes. '+\
              'This might be a case for implementing interpolation '+\
              'between different meshes.'
        assert num.allclose(A.shape, self.vertex_values.shape), msg

        self.set_values(A, location='vertices',
                        indices=indices, verbose=verbose)

    def set_values_from_function(self,
                                 f,
                                 location='vertices',
                                 indices=None,
                                 use_cache=False,
                                 verbose=False):
        """Set values for quantity using specified function

        Input
        f: x, y -> z Function where x, y and z are arrays
        location: Where values are to be stored.
                  Permissible options are: vertices, centroid,
                  unique vertices
                  Default is "vertices"
        indices:
        """

        # FIXME: Should check that function returns something sensible and
        # raise a meaningful exception if it returns None for example

        # FIXME: Should supply absolute coordinates

        # Compute the function values and call set_values again
        if location == 'centroids':
            if indices is None:
                indices = range(len(self))

            V = num.take(self.domain.get_centroid_coordinates(), indices, axis=0)
            x = V[:,0]; y = V[:,1]
            if use_cache is True:
                res = cache(f, (x, y), verbose=verbose)
            else:
                res = f(x, y)

            self.set_values(res, location=location, indices=indices)
        elif location == 'vertices':
            # This is the default branch taken by set_quantity
            M = self.domain.number_of_triangles
            V = self.domain.get_vertex_coordinates()

            x = V[:,0];
            y = V[:,1]
            if use_cache is True:
                values = cache(f, (x, y), verbose=verbose)                
            else:
                if verbose is True:
                    log.critical('Evaluating function in set_values')
                values = f(x, y)

            # FIXME (Ole): This code should replace all the
            # rest of this function and it would work, except
            # one unit test in test_region fails.
            # If that could be resolved this one will be
            # more robust and simple.

            # This should be removed
            if is_scalar(values):
                # Function returned a constant value
                self.set_values_from_constant(values, location,
                                              indices, verbose)
                return

            # This should be removed
            if indices is None:
                for j in range(3):
                    self.vertex_values[:, j] = values[j::3]
            else:
                # Brute force
                for i in indices:
                    for j in range(3):
                        self.vertex_values[i, j] = values[3*i + j]
        else:
            raise Exception('Not implemented: %s' % location)

    def set_values_from_geospatial_data(self,
                                        geospatial_data,
                                        alpha,
                                        location,
                                        indices,
                                        verbose=False,
                                        use_cache=False):
        """Set values based on geo referenced geospatial data object."""

        from anuga.coordinate_transforms.geo_reference import Geo_reference

        points = geospatial_data.get_data_points(absolute=False)
        values = geospatial_data.get_attributes()
        data_georef = geospatial_data.get_geo_reference()

        from anuga.coordinate_transforms.geo_reference import Geo_reference

        points = ensure_numeric(points, num.float)
        values = ensure_numeric(values, num.float)

        if location != 'vertices':
            msg = ("set_values_from_points is only defined for "
                   "location='vertices'")
            raise Exception(msg)

        # Take care of georeferencing
        if data_georef is None:
            data_georef = Geo_reference()

        mesh_georef = self.domain.geo_reference

        # Call fit_interpolate.fit function
        args = (points, )
        kwargs = {'vertex_coordinates': None,
                  'triangles': None,
                  'mesh': self.domain.mesh,
                  'point_attributes': values,
                  'data_origin': data_georef.get_origin(),
                  'mesh_origin': mesh_georef.get_origin(),
                  'alpha': alpha,
                  'verbose': verbose}

        vertex_attributes = apply(fit_to_mesh, args, kwargs)

        # Call underlying method using array values
        self.set_values_from_array(vertex_attributes, location, indices,
                                   use_cache=use_cache, verbose=verbose)

    def set_values_from_points(self,
                               points,
                               values,
                               alpha,
                               location,
                               indices,
                               data_georef=None,
                               verbose=False,
                               use_cache=False):
        """Set quantity values from arbitray data points using fit_interpolate.fit"""

        raise Exception('set_values_from_points is obsolete, use geospatial data object instead')

    def set_values_from_file(self,
                             filename,
                             attribute_name,
                             alpha,
                             location,
                             indices,
                             verbose=False,
                             use_cache=False,
                             max_read_lines=None):
        """Set quantity based on arbitrary points in a points file using
        attribute_name selects name of attribute present in file.
        If attribute_name is not specified, use first available attribute
        as defined in geospatial_data.
        """


        msg = 'Filename must be a text string'
        assert isinstance(filename, basestring), msg
        
        msg = 'Extension should be .pts .dem, .csv, or txt'
        assert os.path.splitext(filename)[1] in ['.pts', '.dem', '.csv', '.txt'], msg
        

        if location != 'vertices':
            msg = "set_values_from_file is only defined for location='vertices'"
            raise Exception(msg)

            
        if True:
            # Use mesh as defined by domain
            # This used to cause problems for caching due to quantities
            # changing, but it now works using the appropriate Mesh object.
            # This addressed ticket:242 and was made to work when bug
            # in ticket:314 was fixed 18 March 2009.
            vertex_attributes = fit_to_mesh(filename,
                                            mesh=self.domain.mesh,
                                            alpha=alpha,
                                            attribute_name=attribute_name,
                                            use_cache=use_cache,
                                            verbose=verbose,
                                            max_read_lines=max_read_lines)
        else:
            # This variant will cause Mesh object to be recreated
            # in fit_to_mesh thus doubling up on the neighbour structure
            # FIXME(Ole): This is now obsolete 19 Jan 2009 except for bug 
            # (ticket:314) which was fixed 18 March 2009.
            nodes = self.domain.get_nodes(absolute=True)
            triangles = self.domain.get_triangles()
            vertex_attributes = fit_to_mesh(filename,
                                            nodes, triangles,
                                            mesh=None,
                                            alpha=alpha,
                                            attribute_name=attribute_name,
                                            use_cache=use_cache,
                                            verbose=verbose,
                                            max_read_lines=max_read_lines)

        # Call underlying method using array values
        if verbose:
            log.critical('Applying fitted data to domain')
        self.set_values_from_array(vertex_attributes, location,
                                   indices, use_cache=use_cache,
                                   verbose=verbose)
        
        
        
        
        
    def set_values_from_utm_grid_file(self,
                             filename,
                             location='vertices',
                             indices=None,
                             verbose=False):
        
        """Read Digital Elevation model from the following ASCII format (.asc or .grd)
    
        Example:
        ncols         3121
        nrows         1800
        xllcorner     722000
        yllcorner     5893000
        cellsize      25
        NODATA_value  -9999
        138.3698 137.4194 136.5062 135.5558 ..........
    
    
        An accompanying file with same basename but extension .prj must exist
        and is used to fix the UTM zone, datum, false northings and eastings.
    
        The prj format is assumed to be as
    
        Projection    UTM
        Zone          56
        Datum         WGS84
        Zunits        NO
        Units         METERS
        Spheroid      WGS84
        Xshift        0.0000000000
        Yshift        10000000.0000000000
        Parameters
        """
        
        
        
        from anuga.file_conversion.grd2array import grd2array
        
        x,y,Z = grd2array(filename)
           
        
        if location == 'centroids':
            points = self.domain.centroid_coordinates
        else:
            points = self.domain.vertex_coordinates
            
        from anuga.geospatial_data.geospatial_data import Geospatial_data,  ensure_absolute
        
        points = ensure_absolute(points, geo_reference=self.domain.geo_reference)        
        
        from  anuga.fit_interpolate.interpolate2d import interpolate2d 
        
        #print points 
        values = interpolate2d(x, y, Z, points, mode='linear', bounds_error=False)
        
        #print values

        # Call underlying method using array values
        if verbose:
            log.critical('Applying fitted data to quantity')
            
        
        if location == 'centroids':
            if indices is None:
                msg = 'Number of values must match number of elements'
                #assert values.shape[0] == N, msg

                self.centroid_values[:] = values
            else:
                msg = 'Number of values must match number of indices'
                assert values.shape[0] == indices.shape[0], msg

                # Brute force
                self.centroid_values[indices] = values
        else:
            if indices is None:
                msg = 'Number of values must match number of elements'
                #assert values.shape[0] == N, msg

                #print values.shape
                #print self.vertex_values.shape 
                self.vertex_values[:] = values.reshape((-1,3))
            else:
                msg = 'Number of values must match number of indices'
                assert values.shape[0] == indices.shape[0], msg

                # Brute force
                self.vertex_values[indices] = values.reshape((-1,3))
            # Cleanup centroid values
            self.interpolate()
            
            
            
            
    def set_values_from_lat_long_grid_file(self,
                             filename,
                             location='vertices',
                             indices=None,
                             verbose=False):
        
        """Read Digital Elevation model from the following ASCII format (.asc or .grd)
    
        Example:
        ncols         3121
        nrows         1800
        xllcorner     722000
        yllcorner     5893000
        cellsize      25
        NODATA_value  -9999
        138.3698 137.4194 136.5062 135.5558 ..........
    
    
        An accompanying file with same basename but extension .prj must exist
        and is used to fix the UTM zone, datum, false northings and eastings.
    
        The prj format is assumed to be as
    
        Projection    UTM
        Zone          56
        Datum         WGS84
        Zunits        NO
        Units         METERS
        Spheroid      WGS84
        Xshift        0.0000000000
        Yshift        10000000.0000000000
        Parameters
        """


        msg = "Function not implemented yet"
        raise Exception, msg

        msg = 'Filename must be a text string'
        assert isinstance(filename, basestring), msg
        
        msg = 'Extension should be .grd or asc'
        assert os.path.splitext(filename)[1] in ['.grd', '.asc'], msg
        

        msg = "set_values_from_grd_file is only defined for location='vertices' or 'centroids'"
        assert location in ['vertices', 'centroids'], msg

    
        root = filename[:-4]
    
    
        #Read DEM data
        datafile = open(filename)
    
        if verbose: log.critical('Reading data from %s' % (filename))
        
        lines = datafile.readlines()
        datafile.close()
    
    
        if verbose: log.critical('Got %d lines' % len(lines))
    
        # Parse the line data
        ncols = int(lines[0].split()[1].strip())
        nrows = int(lines[1].split()[1].strip())

    
        # Do cellsize (line 4) before line 2 and 3
        cellsize = float(lines[4].split()[1].strip())
    
        # Checks suggested by Joaquim Luis
        # Our internal representation of xllcorner
        # and yllcorner is non-standard.
        xref = lines[2].split()
        if xref[0].strip() == 'xllcorner':
            xllcorner = float(xref[1].strip()) # + 0.5*cellsize # Correct offset
        elif xref[0].strip() == 'xllcenter':
            xllcorner = float(xref[1].strip())
        else:
            msg = 'Unknown keyword: %s' % xref[0].strip()
            raise Exception, msg
    
        yref = lines[3].split()
        if yref[0].strip() == 'yllcorner':
            yllcorner = float(yref[1].strip()) # + 0.5*cellsize # Correct offset
        elif yref[0].strip() == 'yllcenter':
            yllcorner = float(yref[1].strip())
        else:
            msg = 'Unknown keyword: %s' % yref[0].strip()
            raise Exception, msg
    
        NODATA_value = int(float(lines[5].split()[1].strip()))
    
        assert len(lines) == nrows + 6
    
  
        #Store data
        import numpy
    
        datafile = open(filename)
        Z = numpy.loadtxt(datafile, skiprows=6)
        datafile.close()
        
        #print Z.shape
        #print Z
        
        # For raster data we need to a flip and transpose
        Z = numpy.flipud(Z)

        # Transpose z to have y coordinates along the first axis and x coordinates
        # along the second axis
        Z = Z.transpose()
    
        x = num.linspace(xllcorner, xllcorner+cellsize*(ncols-1), ncols)
        y = num.linspace(yllcorner, yllcorner+cellsize*(nrows-1), nrows)
        
        
        if location == 'centroids':
            points = self.domain.centroid_coordinates
        else:
            points = self.domain.vertex_coordinates
            
        from anuga.geospatial_data.geospatial_data import Geospatial_data,  ensure_absolute
        points = ensure_absolute(points, geo_reference=self.domain.geo_reference)
            
#         print numpy.max(points[:,0])
#         print numpy.min(points[:,0])
#         print numpy.max(points[:,1])
#         print numpy.min(points[:,1])
#         
#         print numpy.max(x)
#         print numpy.min(x)
#         print numpy.max(y)
#         print numpy.min(y)
        
        
        #print x.shape, x
        #print y.shape, y
        
        
        
        from  anuga.fit_interpolate.interpolate2d import interpolate2d 
        
        #print points 
        values = interpolate2d(x, y, Z, points, mode='linear', bounds_error=False)
        
        #print values

        # Call underlying method using array values
        if verbose:
            log.critical('Applying fitted data to quantity')
            
        
        if location == 'centroids':
            if indices is None:
                msg = 'Number of values must match number of elements'
                #assert values.shape[0] == N, msg

                self.centroid_values[:] = values
            else:
                msg = 'Number of values must match number of indices'
                assert values.shape[0] == indices.shape[0], msg

                # Brute force
                self.centroid_values[indices] = values
        else:
            if indices is None:
                msg = 'Number of values must match number of elements'
                #assert values.shape[0] == N, msg

                #print values.shape
                #print self.vertex_values.shape 
                self.vertex_values[:] = values.reshape((-1,3))
            else:
                msg = 'Number of values must match number of indices'
                assert values.shape[0] == indices.shape[0], msg

                # Brute force
                self.vertex_values[indices] = values.reshape((-1,3))
            # Cleanup centroid values
            self.interpolate()
                            
    def get_extremum_index(self, mode=None, indices=None):
        """Return index for maximum or minimum value of quantity (on centroids)

        Optional arguments:
            mode is either 'max'(default) or 'min'.
            indices is the set of element ids that the operation applies to.

        Usage:
            i = Q.get_extreme_index()

        Notes:
            We do not seek the extremum at vertices as each vertex can
            have multiple values - one for each triangle sharing it.

            If there are multiple cells with same maximum value, the
            first cell encountered in the triangle array is returned.
        """

        V = self.get_values(location='centroids', indices=indices)

        # Always return absolute indices
        if mode is None or mode == 'max':
            i = num.argmax(V)
        elif mode == 'min':
            i = num.argmin(V)
        else:
            raise ValueError('Bad mode value, got: %s' % str(mode))

        if indices is None:
            return i
        else:
            return indices[i]

    def get_maximum_index(self, indices=None):
        """See get extreme index for details"""

        return self.get_extremum_index(mode='max', indices=indices)

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
        V = self.get_values(location='centroids')      #, indices=indices)

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

    def get_minimum_index(self, indices=None):
        """See get extreme index for details"""

        return self.get_extremum_index(mode='min', indices=indices)

    def get_minimum_value(self, indices=None):
        """Return minimum value of quantity (on centroids)

        Optional argument:
            indices is the set of element ids that the operation applies to.

        Usage:
            v = get_minimum_value()

        See get_maximum_value for more details.
        """

        i = self.get_minimum_index(indices)
        V = self.get_values(location='centroids')

        return V[i]


    def get_minimum_location(self, indices=None):
        """Return location of minimum value of quantity (on centroids)

        Optional argument:
            indices is the set of element ids that the operation applies to.

        Usage:
            x, y = get_minimum_location()

        Notes:
            We do not seek the maximum at vertices as each vertex can
            have multiple values - one for each triangle sharing it.

            If there are multiple cells with same maximum value, the
            first cell encountered in the triangle array is returned.
        """

        i = self.get_minimum_index(indices)
        x, y = self.domain.get_centroid_coordinates()[i]

        return x, y

    def get_interpolated_values(self,
                                interpolation_points,
                                use_cache=False,
                                verbose=False):
        """Get values at interpolation points

        The argument interpolation points must be given as either a
        list of absolute UTM coordinates or a geospatial data object.
        """

        # FIXME (Ole): Points might be converted to coordinates relative to mesh origin
        # This could all be refactored using the
        # 'change_points_geo_ref' method of Class geo_reference.
        # The purpose is to make interpolation points relative
        # to the mesh origin.
        #
        # Speed is also a consideration here.

        # Ensure that interpolation points is either a list of
        # points, Nx2 array, or geospatial and convert to numeric array
        if isinstance(interpolation_points, Geospatial_data):
            # Ensure interpolation points are in absolute UTM coordinates
            interpolation_points = \
                interpolation_points.get_data_points(absolute=True)

        # Reconcile interpolation points with georeference of domain
        interpolation_points = \
            self.domain.geo_reference.get_relative(interpolation_points)
        interpolation_points = ensure_numeric(interpolation_points)


        # Get internal representation (disconnected) of vertex values
        vertex_values, triangles = self.get_vertex_values(xy=False,
                                                          smooth=False)

        # Get possibly precomputed interpolation object
        I = self.domain.get_interpolation_object()

        # Call interpolate method with interpolation points
        result = I.interpolate_block(vertex_values, interpolation_points,
                                     use_cache=use_cache, verbose=verbose)

        return result

    def get_values(self,
                   interpolation_points=None,
                   location='vertices',
                   indices=None,
                   use_cache=False,
                   verbose=False):
        """Get values for quantity

        Extract values for quantity as a numeric array.

        Inputs:
           interpolation_points: List of x, y coordinates where value is
                                 sought (using interpolation). If points
                                 are given, values of location and indices
                                 are ignored. Assume either absolute UTM
                                 coordinates or geospatial data object.

           location: Where values are to be stored.
                     Permissible options are: vertices, edges, centroids
                     and unique vertices. Default is 'vertices'


        The returned values will have the leading dimension equal to 
        length of the indices list or N (all values) if indices is None.

        In case of location == 'centroids' the dimension of returned
        values will be a list or a numerical array of length N, N being
        the number of elements.

        In case of location == 'vertices' or 'edges' the dimension of
        returned values will be of dimension Nx3

        In case of location == 'unique vertices' the average value at
        each vertex will be returned and the dimension of returned values
        will be a 1d array of length "number of vertices"

        Indices is the set of element ids that the operation applies to.
        If indices is None (or omitted) all elements are returned as
        a copy of the relevant array. If performance is critical, 
        use arrays domain.centroid_values, domain.vertex_values and 
        domain.edge_values directly.

        The values will be stored in elements following their
        internal ordering.
        """

        # FIXME (Ole): I reckon we should have the option of passing a
        #              polygon into get_values. The question becomes how
        #              resulting values should be ordered.

        if verbose is True:
            log.critical('Getting values from %s' % location)

        if interpolation_points is not None:
            return self.get_interpolated_values(interpolation_points,
                                                use_cache=use_cache,
                                                verbose=verbose)

        # FIXME (Ole): Consider deprecating 'edges' - but not if it is used
        # elsewhere in ANUGA.
        # Edges have already been deprecated in set_values, see changeset:5521,
        # but *might* be useful in get_values. Any thoughts anyone?
        # YES (Ole): Edge values are necessary for volumetric balance 
        # check and inflow boundary. Keep them!

        if location not in ['vertices', 'centroids',
                            'edges', 'unique vertices']:
            msg = 'Invalid location: %s' % location
            raise Exception(msg)


        msg = '\'indices\' must be a list, array or None'
        assert isinstance(indices, (types.NoneType, list, num.ndarray)), msg

        if location == 'centroids':
            if indices is None:
                return self.centroid_values.copy()
            else:
                return num.take(self.centroid_values, indices, axis=0)
        elif location == 'edges':
            if indices is  None:
                return self.edge_values.copy()
            else:    
                return num.take(self.edge_values, indices, axis=0)
        elif location == 'unique vertices':
            if indices is None:
                indices=range(self.domain.get_number_of_nodes())
            vert_values = []

            # Go through list of unique vertices
            for unique_vert_id in indices:
                triangles = self.domain.get_triangles_and_vertices_per_node(node=unique_vert_id)

                # In case there are unused points
                if len(triangles) == 0:
                    msg = 'Unique vertex not associated with triangles'
                    raise Exception(msg)

                # Go through all triangle, vertex pairs
                # Average the values
                # FIXME (Ole): Should we merge this with get_vertex_values
                sum = 0
                for triangle_id, vertex_id in triangles:
                    sum += self.vertex_values[triangle_id, vertex_id]
                vert_values.append(sum / len(triangles))
            return num.array(vert_values, num.float)
        else:
            if indices is None:
                return self.vertex_values.copy()
            else:    
                return num.take(self.vertex_values, indices, axis=0)

    def set_vertex_values(self,
                          A,
                          indices=None,
                          use_cache=False,
                          verbose=False):
        """Set vertex values for all unique vertices based on input array A
        which has one entry per unique vertex, i.e. one value for each row in
        array self.domain.nodes.

        indices is the list of vertex_id's that will be set.

        This function is used by set_values_from_array
        """

        # Check that A can be converted to array and is of appropriate dim
        A = ensure_numeric(A, num.float)
        assert len(A.shape) == 1

        if indices is None:
            assert A.shape[0] == self.domain.get_nodes().shape[0]
            vertex_list = range(A.shape[0])
        else:
            assert A.shape[0] == len(indices)
            vertex_list = indices

        #FIXME(Ole): This function ought to be faster.
        # We need to get the triangles_and_vertices list
        # from domain in one hit, then cache the computation of the
        # Nx3 array of vertex values that can then be assigned using
        # set_values_from_array.
        #
        # Alternatively, some C code would be handy
        #
        self._set_vertex_values(vertex_list, A)



    # Note Padarn 27/11/12:
    # This function has been changed and now uses an external c function
    # to set the 'vertex_values' instead of a python for loop. The function
    # 'get_triangles_and_vertices_per_node' has been removed and replaced by
    # 'build_inverted_triangle_structure. This now adds extra stored array to
    # the mesh object - this could be removed after the c function below uses
    #them.
    # Note, the naming of this function seems confusing - it seems to actually
    # update the 'node values' given a list of vertices.
    def _set_vertex_values(self, vertex_list, A):
        """Go through list of unique vertices
        This is the common case e.g. when values
        are obtained from a pts file through fitting
        """
        # If required, set up required arrays storing information about the triangle
        # vertex structure of the mesh.
        if not (hasattr(self.domain.mesh, 'number_of_triangles_per_node') and \
                hasattr(self.domain.mesh, 'vertex_value_indices') and \
                hasattr(self.domain.mesh, 'node_index')):
            self.build_inverted_triangle_structure()

        set_vertex_values_c(self, num.array(vertex_list), A)
        self.interpolate()

    def smooth_vertex_values(self, use_cache=False, verbose=False):
        """Smooths vertex values."""

        A, V = self.get_vertex_values(xy=False, smooth=True)
        self.set_vertex_values(A, use_cache=use_cache, verbose=verbose)



    ############################################################################
    # Methods for outputting model results
    ############################################################################

    def get_vertex_values(self, xy=True, smooth=None, precision=None):
        """Return vertex values like an OBJ format i.e. one value per node.

        The vertex values are returned as one sequence in the 1D float array A.
        If requested the coordinates will be returned in 1D arrays X and Y.

        The connectivity is represented as an integer array, V, of dimension
        Mx3, where M is the number of triangles. Each row has three indices
        defining the triangle and they correspond to elements in the arrays
        X, Y and A.

        If smooth is True, vertex values corresponding to one common coordinate
        set will be smoothed by taking the average of vertex values for each
        node.  In this case vertex coordinates will be de-duplicated
        corresponding to the original nodes as obtained from the method
        general_mesh.get_nodes()

        If no smoothings is required, vertex coordinates and values will be
        aggregated as a concatenation of values at vertices 0, vertices 1 and
        vertices 2.  This corresponds to the node coordinates obtained from the
        method general_mesh.get_vertex_coordinates()

        Calling convention
        if xy is True:
           X, Y, A, V = get_vertex_values
        else:
           A, V = get_vertex_values
        """

        if smooth is None:
            # Take default from domain
            try:
                smooth = self.domain.smooth
            except:
                smooth = False

        if precision is None:
            precision = num.float

        if smooth is True:
            # Ensure continuous vertex values by averaging values at each node
            V = self.domain.get_triangles()
            N = self.domain.number_of_full_nodes # Ignore ghost nodes if any
            A = num.zeros(N, num.float)
            points = self.domain.get_nodes()

            if 1:
                # Fast C version
                average_vertex_values(ensure_numeric(self.domain.vertex_value_indices),
                                      ensure_numeric(self.domain.number_of_triangles_per_node),
                                      ensure_numeric(self.vertex_values),
                                      A)
                A = A.astype(precision)
            else:
                # Slow Python version
                current_node = 0
                k = 0 # Track triangles touching on node
                total = 0.0
                for index in self.domain.vertex_value_indices:
                    if current_node == N:
                        msg = 'Current node exceeding number of nodes (%d) ' % N
                        raise Exception(msg)

                    k += 1

                    volume_id = index / 3
                    vertex_id = index % 3

                    v = self.vertex_values[volume_id, vertex_id]
                    total += v

                    if self.domain.number_of_triangles_per_node[current_node] == k:
                        A[current_node] = total/k

                        # Move on to next node
                        total = 0.0
                        k = 0
                        current_node += 1
        else:
            # Return disconnected internal vertex values
            V = self.domain.get_disconnected_triangles()
            points = self.domain.get_vertex_coordinates()
            A = self.vertex_values.flatten().astype(precision)

        # Return
        if xy is True:
            X = points[:,0].astype(precision)
            Y = points[:,1].astype(precision)

            return X, Y, A, V
        else:
            return A, V

    def extrapolate_first_order(self):
        """Extrapolate conserved quantities from centroid to vertices and edges
        for each volume using first order scheme.
        """

        qc = self.centroid_values
        qv = self.vertex_values
        qe = self.edge_values

        #for i in range(3):
        #    qe[:,i] = qc

        qe[:] = qc[:,num.newaxis]
        qv[:] = qe

        self.x_gradient[:] = 0.0
        self.y_gradient[:] = 0.0

    def get_integral(self, full_only=True, region=None, indices=None):

        """Compute the integral of quantity across entire domain,
        or over a region. Eg
        my_region = anuga.Region(polygon = user_polygon)

        Q.get_integral(region = my_region)"""

        assert region is None or indices is None

        
        areas = self.domain.get_areas()

        if region is None and indices is None:
            if full_only:
                indices = num.where(self.domain.tri_full_flag ==1)[0]
                return num.sum(areas[indices]*self.centroid_values[indices])
            else:
                return num.sum(areas*self.centroid_values)

        if indices is None:
            indices = region.get_indices(full_only)

        if indices == []:
            return 0.0
        elif indices is None:
            return num.sum(areas*self.centroid_values)
        else:
            return num.sum(areas[indices]*self.centroid_values[indices])


    def get_gradients(self):
        """Provide gradients. Use compute_gradients first."""

        return self.x_gradient, self.y_gradient

    def update(self, timestep):
        # Call correct module function
        # (either from this module or C-extension)
        return update(self, timestep)

    def compute_gradients(self):
        # Call correct module function
        # (either from this module or C-extension)
        return compute_gradients(self)


    def compute_local_gradients(self):
        # Call correct module function
        # (either from this module or C-extension)
        return compute_local_gradients(self)




    def limit(self):
        # Call correct module depending on whether
        # basing limit calculations on edges or vertices
        limit_old(self)

    def limit_vertices_by_all_neighbours(self):
        # Call correct module function
        # (either from this module or C-extension)
        limit_vertices_by_all_neighbours(self)

    def limit_edges_by_all_neighbours(self):
        # Call correct module function
        # (either from this module or C-extension)
        limit_edges_by_all_neighbours(self)

    def limit_edges_by_neighbour(self):
        # Call correct module function
        # (either from this module or C-extension)
        limit_edges_by_neighbour(self)

    def extrapolate_second_order(self):
        # Call correct module function
        # (either from this module or C-extension)
        compute_gradients(self)
        extrapolate_from_gradient(self)

    def extrapolate_second_order_and_limit_by_edge(self):
        # Call correct module function
        # (either from this module or C-extension)
        extrapolate_second_order_and_limit_by_edge(self)

    def extrapolate_second_order_and_limit_by_vertex(self):
        # Call correct module function
        # (either from this module or C-extension)
        extrapolate_second_order_and_limit_by_vertex(self)

    def bound_vertices_below_by_constant(self, bound):
        # Call correct module function
        # (either from this module or C-extension)
        bound_vertices_below_by_constant(self, bound)

    def bound_vertices_below_by_quantity(self, quantity):
        # Call correct module function
        # (either from this module or C-extension)

        # check consistency
        assert self.domain == quantity.domain
        bound_vertices_below_by_quantity(self, quantity)

    def backup_centroid_values(self):
        # Call correct module function
        # (either from this module or C-extension)
        backup_centroid_values(self)

    def saxpy_centroid_values(self, a, b):
        # Call correct module function
        # (either from this module or C-extension)
        saxpy_centroid_values(self, a, b)


class Conserved_quantity(Quantity):
    """Class conserved quantity being removed, use Quantity."""

    def __init__(self, domain, vertex_values=None):
        msg = 'ERROR: Use Quantity instead of Conserved_quantity'
        raise Exception(msg)


######
# Prepare the C extensions.
######
from quantity_ext import \
         average_vertex_values,\
         backup_centroid_values,\
         saxpy_centroid_values,\
         compute_gradients,\
         compute_local_gradients,\
         limit_old,\
         limit_vertices_by_all_neighbours,\
         limit_edges_by_all_neighbours,\
         limit_edges_by_neighbour,\
         limit_gradient_by_neighbour,\
         extrapolate_from_gradient,\
         extrapolate_second_order_and_limit_by_edge,\
         extrapolate_second_order_and_limit_by_vertex,\
         bound_vertices_below_by_constant,\
         bound_vertices_below_by_quantity,\
         interpolate_from_vertices_to_edges,\
         interpolate_from_edges_to_vertices,\
         set_vertex_values_c, \
         update
