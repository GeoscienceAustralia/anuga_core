
"""
boundary.py - Classes for implementing boundary conditions
"""


from builtins import str
from builtins import object
from warnings import warn

from anuga.utilities.numerical_tools import NAN    
from anuga.fit_interpolate.interpolate import Modeltime_too_late
from anuga.fit_interpolate.interpolate import Modeltime_too_early
import anuga.utilities.log as log

import numpy as num


class Boundary(object):
    """Base class for boundary conditions.
       Specific boundary conditions must provide values for
       the conserved_quantities

       A boundary object has one neighbour; the one it
       serves.
    """

    def __init__(self):
        self.verbose = True
        self.default_boundary = None
        self.default_boundary_invoked = False
        

    def evaluate(self, vol_id=None, edge_id=None):
        msg = 'Generic class Boundary must be subclassed'
        raise Exception(msg)


    def evaluate_segment(self, domain=None, segment_edges=None):
        """
        Evaluate boundary condition at edges of a domain in a list
        defined by segment_edges
        
        segment_edges are a sublist of the list of edges definded by the
        arrays domain.boundary_cells  and domain.boundary_edges

        Go through list of boundary objects and update boundary values
        for all conserved quantities on boundary.
        It is assumed that the ordering of conserved quantities is
        consistent between the domain and the boundary object, i.e.
        the jth element of vector q must correspond to the jth conserved
        quantity in domain.

        This implementation uses the evaluate proceudure, but for efficiency
        it is recommended to overload this procedure with a numpy or C based
        implementation.
        """

        if segment_edges is None:
            return
        if domain is None:
            return
        
        for i in segment_edges:
            vol_id  = domain.boundary_cells[i]
            edge_id = domain.boundary_edges[i]


            q_bdry = self.evaluate(vol_id, edge_id)

            if len(q_bdry) == len(domain.evolved_quantities):
                # conserved and evolved quantities are the same
                q_evol = q_bdry
            elif len(q_bdry) == len(domain.conserved_quantities):
                # boundary just returns conserved quantities
                # Need to calculate all the evolved quantities
                # Use default conversion

                q_evol = domain.get_evolved_quantities(vol_id, edge = edge_id)

                q_evol = domain.conserved_values_to_evolved_values \
                                                        (q_bdry, q_evol)
            else:
                msg = 'Boundary must return array of either conserved'
                msg += ' or evolved quantities'
                raise Exception(msg)

            for j, name in enumerate(domain.evolved_quantities):
                Q = domain.quantities[name]
                Q.boundary_values[i] = q_evol[j]


    def get_time(self):

        return self.domain.get_time()

        
    def get_boundary_values(self, t=None):

        if t is None:
            t = self.get_time()

        try:
            res = self.function(t)
        except Modeltime_too_early as e:
            raise Modeltime_too_early(e)
        except Modeltime_too_late as e:
            if self.default_boundary is None:
                raise Modeltime_too_late(e) # Reraise exception
            else:
                # Pass control to default boundary
                res = self.default_boundary

                # Ensure that result cannot be manipulated
                # This is a real danger in case the 
                # default_boundary is a Dirichlet type 
                # for instance.
                from copy import deepcopy 
                res = deepcopy(res) 

                if self.default_boundary_invoked is False:
                    if self.verbose:                
                        # Issue warning the first time
                        msg = '%s' %str(e)
                        msg += 'Instead I will use the default boundary value: %s\n'\
                            %str(self.default_boundary) 
                        msg += 'Note: Further warnings will be suppressed'
                        log.critical(msg)

                    # FIXME (Ole): Replace this crude flag with
                    # Python's ability to print warnings only once.
                    # See http://docs.python.org/lib/warning-filter.html
                    self.default_boundary_invoked = True

        return res



class Transmissive_boundary(Boundary):
    """Transmissive boundary returns same conserved quantities as
    those present in its neighbour volume.

    Underlying domain must be specified when boundary is instantiated
    """

    def __init__(self, domain = None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for transmissive boundary'
            raise Exception(msg)

        self.domain = domain

    def __repr__(self):
        return 'Transmissive_boundary(%s)' %self.domain

    def evaluate(self, vol_id, edge_id):
        """Transmissive boundaries return the edge values
        of the volume they serve.
        """

        if self.domain.centroid_transmissive_bc :
            q = self.domain.get_conserved_quantities(vol_id)
        else:
            q = self.domain.get_conserved_quantities(vol_id, edge = edge_id)


        return q



    def evaluate_segment(self, domain, segment_edges):

        if segment_edges is None:
            return
        if domain is None:
            return


        ids = segment_edges
        vol_ids  = domain.boundary_cells[ids]
        edge_ids = domain.boundary_edges[ids]


        if domain.centroid_transmissive_bc :
            for j, name in enumerate(domain.evolved_quantities):
                Q = domain.quantities[name]
                Q.boundary_values[ids] = Q.centroid_values[vol_ids]
        else:
            for j, name in enumerate(domain.evolved_quantities):
                Q = domain.quantities[name]
                Q.boundary_values[ids] = Q.edge_values[vol_ids,edge_ids]



class Dirichlet_boundary(Boundary):
    """Dirichlet boundary returns constant values for the
    conserved quantities
    """


    def __init__(self, dirichlet_values=None):
        Boundary.__init__(self)

        if dirichlet_values is None:
            msg = 'Must specify one value for each quantity'
            raise Exception(msg)

        self.dirichlet_values=num.array(dirichlet_values, float)



    def __repr__(self):
        return 'Dirichlet boundary (%s)' %self.dirichlet_values


    def evaluate(self, vol_id=None, edge_id=None):
        return self.dirichlet_values

    def evaluate_segment(self, domain, segment_edges):

        if segment_edges is None:
            return
        if domain is None:
            return


        ids = segment_edges


        vol_ids  = domain.boundary_cells[ids]
        edge_ids = domain.boundary_edges[ids]

        q_bdry = self.dirichlet_values

        conserved_quantities = True
        if len(q_bdry) == len(domain.evolved_quantities):
            # enough dirichlet values to set evolved quantities
            conserved_quantities = False

        #--------------------------------------------------
        # First populate all the boundary values with
        # interior edge values
        #--------------------------------------------------
        if  conserved_quantities:
            for j, name in enumerate(domain.evolved_quantities):
                Q = domain.quantities[name]
                Q.boundary_values[ids] = Q.edge_values[vol_ids,edge_ids]

        #--------------------------------------------------
        # Now over write with constant Dirichlet value
        #--------------------------------------------------
        if conserved_quantities:
            quantities = domain.conserved_quantities
        else:
            quantities = domain.evolved_quantities

        #-------------------------------------------------
        # Now update to Dirichlet values
        #-------------------------------------------------
        for j, name in enumerate(quantities):
            Q = domain.quantities[name]
            Q.boundary_values[ids] = q_bdry[j]

class Compute_fluxes_boundary(Boundary):
    """ Associate the Compute_fluxes_boundary BC
    to each boundary segment where a special calculation
    of the fluxes will occur.
    """

    def __init__(self):
        Boundary.__init__(self)
        """ Instantiate a Compute_fluxes_boundary.
            domain: underlying domain
            """


    def __repr__(self):
        """ Return a representation of this instance. """
        return 'Compute_fluxes_boundary(%s)' % self.domain

    def evaluate(self, vol_id, edge_id):
        """Do something basic"""

        return None

    def evaluate_segment(self, domain, segment_edges):
        """Don't do any calculation when applying this BC.
        The work will be done in compute_fluxes.
        """

        return



class Time_boundary(Boundary):
    """Time dependent boundary returns values for the
    conserved quantities as a function of time.
    Must specify domain to get access to model time and a function of t
    which must return conserved quantities as a function time.
    
    Example:
      B = Time_boundary(domain, 
                        function=lambda t: [(60<t<3660)*2, 0, 0])
      
      This will produce a boundary condition with is a 2m high square wave
      starting 60 seconds into the simulation and lasting one hour.
      Momentum applied will be 0 at all times.
                        
    """

    def __init__(self, domain=None,
                 #f=None, # Should be removed and replaced by function below
                 function=None,
                 default_boundary=None,
                 verbose=False):
        Boundary.__init__(self)

        self.default_boundary = default_boundary
        self.default_boundary_invoked = False    # Flag
        self.domain = domain
        self.verbose = verbose

        if domain is None:
            raise Exception('You must specify a domain to Time_boundary')
            
        if function is None:
            raise Exception('You must specify a function to Time_boundary')
            
        
        try:
            q = function(0.0)
        except Exception as e:
            msg = 'Function for time boundary could not be executed:\n%s' %e
            raise Exception(msg)


        try:
            q = num.array(q, float)
        except:
            msg = 'Return value from time boundary function could '
            msg += 'not be converted into a numeric array of floats.\n'
            msg += 'Specified function should return either list or array.\n'
            msg += 'I got %s' %str(q)
            raise Exception(msg)

        msg = 'ERROR: Time boundary function must return a 1d list or array '
        assert len(q.shape) == 1, msg

        d = len(domain.conserved_quantities)
        msg = 'Return value for function must be a list or an array of length %d' %d
        assert len(q) == d, msg

        self.function = function
        self.domain = domain

    def __repr__(self):
        return 'Time boundary'

    def get_time(self):

        return self.domain.get_time()

    def evaluate(self, vol_id=None, edge_id=None):

        return self.get_boundary_values()


    def evaluate_segment(self, domain, segment_edges):

        if segment_edges is None:
            return
        if domain is None:
            return

        ids = segment_edges

        vol_ids  = domain.boundary_cells[ids]
        edge_ids = domain.boundary_edges[ids]

        q_bdry = self.get_boundary_values()

        conserved_quantities = True
        if len(q_bdry) == len(domain.evolved_quantities):
            # enough dirichlet values to set evolved quantities
            conserved_quantities = False

        #--------------------------------------------------
        # First populate all the boundary values with
        # interior edge values
        #--------------------------------------------------
        if  conserved_quantities:
            for j, name in enumerate(domain.evolved_quantities):
                Q = domain.quantities[name]
                Q.boundary_values[ids] = Q.edge_values[vol_ids,edge_ids]

        #--------------------------------------------------
        # Now over write with constant Dirichlet value
        #--------------------------------------------------
        if conserved_quantities:
            quantities = domain.conserved_quantities
        else:
            quantities = domain.evolved_quantities

        #-------------------------------------------------
        # Now update to Dirichlet values
        #-------------------------------------------------
        for j, name in enumerate(quantities):
            Q = domain.quantities[name]
            Q.boundary_values[ids] = q_bdry[j]
            






class Time_space_boundary(Boundary):
    """Time and spatially dependent boundary returns values for the
    conserved quantities as a function of time and space.
    Must specify domain to get access to model time and a function of t,x,y
    which must return conserved quantities at specified time and location.
    
    Example:
      B = Time_space_boundary(domain, 
                              function=lambda t,x,y: [(60<t<3660)*2, 0, 0])
      
      This will produce a boundary condition with is a 2m high square wave
      starting 60 seconds into the simulation and lasting one hour.
      Momentum applied will be 0 at all times.
                        
    """


    def __init__(self, domain=None,
                 function=None,
                 default_boundary=None,
                 verbose=False):
        Boundary.__init__(self)

        self.default_boundary = default_boundary
        self.default_boundary_invoked = False    # Flag
        self.domain = domain
        self.verbose = verbose

        
        if function is None:
            raise Exception('You must specify a function to Time_space_boundary')
            
        try:
            q = function(0.0, 0.0, 0.0)
        except Exception as e:
            msg = 'Function for time_space_boundary could not be executed:\n%s' %e
            raise Exception(msg)


        try:
            q = num.array(q, float)
        except:
            msg = 'Return value from time_space_boundary function could '
            msg += 'not be converted into a numeric array of floats.\n'
            msg += 'Specified function should return either list or array.\n'
            msg += 'I got %s' %str(q)
            raise Exception(msg)

        msg = 'ERROR: Time_space_boundary function must return a 1d list or array '
        assert len(q.shape) == 1, msg

        d = len(domain.conserved_quantities)
        msg = 'Return value for function must be a list or an array of length %d' %d
        assert len(q) == d, msg

        self.function = function
        self.domain = domain

    def __repr__(self):
        return 'Time space boundary'

    def evaluate(self, vol_id=None, edge_id=None):
        # FIXME (Ole): I think this should be get_time(), see ticket:306

        [x,y] = self.domain.get_edge_midpoint_coordinate(vol_id, edge_id)

        try:
            res = self.function(self.domain.get_time(), x, y)
        except Modeltime_too_early as e:
            raise Modeltime_too_early(e)
        except Modeltime_too_late as e:
            if self.default_boundary is None:
                raise Exception(e) # Reraise exception
            else:
                # Pass control to default boundary
                res = self.default_boundary.evaluate(vol_id, edge_id)
                
                # Ensure that result cannot be manipulated
                # This is a real danger in case the 
                # default_boundary is a Dirichlet type 
                # for instance. 
                res = res.copy() 
                
                if self.default_boundary_invoked is False:
                    if self.verbose:                
                        # Issue warning the first time
                        msg = '%s' %str(e)
                        msg += 'Instead I will use the default boundary: %s\n'\
                            %str(self.default_boundary) 
                        msg += 'Note: Further warnings will be supressed'
                        log.critical(msg)
               
                    # FIXME (Ole): Replace this crude flag with
                    # Python's ability to print warnings only once.
                    # See http://docs.python.org/lib/warning-filter.html
                    self.default_boundary_invoked = True

        return res

class File_boundary(Boundary):
    """The File_boundary reads values for the conserved
    quantities from an sww NetCDF file, and returns interpolated values
    at the midpoints of each associated boundary segment.
    Time dependency is interpolated linearly.

    Assumes that file contains a time series and possibly
    also spatial info. See docstring for File_function in util.py
    for details about admissible file formats

    File boundary must read and interpolate from *smoothed* version
    as stored in sww and cannot work with the discontinuous triangles.

    Example:
    Bf = File_boundary('source_file.sww', domain)


    Note that the resulting solution history is not exactly the same as if
    the models were coupled as there is no feedback into the source model.
    
    Optional keyword argument default_boundary must be either None or 
    an instance of class descending from class Boundary.
    This will be used in case model time exceeds that available in the 
    underlying data.
       
    """

    # FIXME (Ole): I kind of like the name Spatio_Temporal_boundary for this
    # rather than File_boundary

    def __init__(self, filename, domain, 
                 time_thinning=1, 
                 time_limit=None,
                 boundary_polygon=None,    
                 default_boundary=None,
                 use_cache=False, 
                 verbose=False): 

        import time
        from anuga.config import time_format
        from anuga.abstract_2d_finite_volumes.util import file_function

        Boundary.__init__(self)

        # Get x,y vertex coordinates for all triangles
        V = domain.vertex_coordinates

        # Compute midpoint coordinates for all boundary elements
        # Only a subset may be invoked when boundary is evaluated but
        # we don't know which ones at this stage since this object can
        # be attached to
        # any tagged boundary later on.

        if verbose: log.critical('Find midpoint coordinates of entire boundary')
        self.midpoint_coordinates = num.zeros((len(domain.boundary), 2), float)
        boundary_keys = list(domain.boundary.keys())

        xllcorner = domain.geo_reference.get_xllcorner()
        yllcorner = domain.geo_reference.get_yllcorner()        
        

        # Make ordering unique #FIXME: should this happen in domain.py?
        boundary_keys.sort()

        # Record ordering #FIXME: should this also happen in domain.py or general_mesh.py?
        self.boundary_indices = {}
        for i, (vol_id, edge_id) in enumerate(boundary_keys):

            self.midpoint_coordinates[i,:] = domain.get_edge_midpoint_coordinate(vol_id, edge_id, 
                                                                                 absolute=True)

            # Register index of this boundary edge for use with evaluate
            self.boundary_indices[(vol_id, edge_id)] = i

            
            
        if verbose: log.critical('Initialise file_function')
        self.F = file_function(filename,
                               domain,
                               quantities=domain.conserved_quantities,
	                           interpolation_points=self.midpoint_coordinates,
                               time_thinning=time_thinning,
                               time_limit=time_limit,
                               use_cache=use_cache, 
                               verbose=verbose,
                               boundary_polygon=boundary_polygon)
                             
        # Check and store default_boundary
        msg = 'Keyword argument default_boundary must be either None '
        msg += 'or a boundary object.\n I got %s' %(str(default_boundary))
        assert default_boundary is None or\
            isinstance(default_boundary, Boundary), msg
        self.default_boundary = default_boundary
        self.default_boundary_invoked = False    # Flag

        # Store pointer to domain and verbosity
        self.domain = domain
        self.verbose = verbose


        # Here we'll flag indices outside the mesh as a warning
        # as suggested by Joaquim Luis in sourceforge posting
        # November 2007
        # We won't make it an error as it is conceivable that
        # only part of mesh boundary is actually used with a given
        # file boundary sww file. 
        if hasattr(self.F, 'indices_outside_mesh') and\
               len(self.F.indices_outside_mesh) > 0:
            msg = 'WARNING: File_boundary has points outside the mesh '
            msg += 'given in %s. ' %filename
            msg += 'See warning message issued by Interpolation_function '
            msg += 'for details (should appear above somewhere if '
            msg += 'verbose is True).\n'
            msg += 'This is perfectly OK as long as the points that are '
            msg += 'outside aren\'t used on the actual boundary segment.'
            if verbose is True:            
                log.critical(msg)
            #raise Exception(msg)

        # Test that file function can be called
        q = self.F(0, point_id=0)
        d = len(domain.conserved_quantities)
        msg = 'Values specified in file %s must be ' %filename
        msg += ' a list or an array of length %d' %d
        assert len(q) == d, msg


    def __repr__(self):
        return 'File boundary'


    def evaluate(self, vol_id=None, edge_id=None):
        """Return linearly interpolated values based on domain time
        at midpoint of segment defined by vol_id and edge_id.
        """

        t = self.domain.get_time()
        
        if vol_id is not None and edge_id is not None:
            i = self.boundary_indices[ vol_id, edge_id ]
            
            try:
                res = self.F(t, point_id=i)
            except Modeltime_too_early as e:
                raise Modeltime_too_early(e)
            except Modeltime_too_late as e:
                if self.default_boundary is None:
                    raise Exception(e) # Reraise exception
                else:
                    # Pass control to default boundary
                    res = self.default_boundary.evaluate(vol_id, edge_id)
                    
                    # Ensure that result cannot be manipulated
                    # This is a real danger in case the 
                    # default_boundary is a Dirichlet type 
                    # for instance. 
                    res = res.copy() 
                    
                    if self.default_boundary_invoked is False:
                        # Issue warning the first time
                        if self.verbose:
                            msg = '%s' %str(e)
                            msg += 'Instead I will use the default boundary: %s\n'\
                                %str(self.default_boundary) 
                            msg += 'Note: Further warnings will be supressed'
                            log.critical(msg)
                   
                        # FIXME (Ole): Replace this crude flag with
                        # Python's ability to print warnings only once.
                        # See http://docs.python.org/lib/warning-filter.html
                        self.default_boundary_invoked = True
            
            if num.any(res == NAN):
                x,y=self.midpoint_coordinates[i,:]
                msg = 'NAN value found in file_boundary at '
                msg += 'point id #%d: (%.2f, %.2f).\n' %(i, x, y)

                if hasattr(self.F, 'indices_outside_mesh') and\
                       len(self.F.indices_outside_mesh) > 0:
                    # Check if NAN point is due it being outside
                    # boundary defined in sww file.

                    if i in self.F.indices_outside_mesh:
                        msg += 'This point refers to one outside the '
                        msg += 'mesh defined by the file %s.\n'\
                               %self.F.filename
                        msg += 'Make sure that the file covers '
                        msg += 'the boundary segment it is assigned to '
                        msg += 'in set_boundary.'
                    else:
                        msg += 'This point is inside the mesh defined '
                        msg += 'the file %s.\n' %self.F.filename
                        msg += 'Check this file for NANs.'
                raise Exception(msg)
            
            return res 
        else:
            msg = 'Boundary call without point_id not implemented.\n'
            msg += 'vol_id=%s, edge_id=%s' %(str(vol_id), str(edge_id))
            raise Exception(msg)

class AWI_boundary(Boundary):
    """The AWI_boundary reads values for the conserved
    quantities (only STAGE) from an sww NetCDF file, and returns interpolated values
    at the midpoints of each associated boundary segment.
    Time dependency is interpolated linearly.

    Assumes that file contains a time series and possibly
    also spatial info. See docstring for File_function in util.py
    for details about admissible file formats

    AWI boundary must read and interpolate from *smoothed* version
    as stored in sww and cannot work with the discontinuos triangles.

    Example:
    Ba = AWI_boundary('source_file.sww', domain)


    Note that the resulting solution history is not exactly the same as if
    the models were coupled as there is no feedback into the source model.

    This was added by Nils Goseberg et al in April 2009
    """

    def __init__(self, filename, domain, time_thinning=1, 
                 use_cache=False, verbose=False):
        import time
        from anuga.config import time_format
        from anuga.abstract_2d_finite_volumes.util import file_function

        Boundary.__init__(self)

        # Get x,y vertex coordinates for all triangles
        V = domain.vertex_coordinates

        # Compute midpoint coordinates for all boundary elements
        # Only a subset may be invoked when boundary is evaluated but
        # we don't know which ones at this stage since this object can
        # be attached to
        # any tagged boundary later on.

        if verbose: log.critical('Find midpoint coordinates of entire boundary')
        self.midpoint_coordinates = num.zeros((len(domain.boundary), 2), float)
        boundary_keys = list(domain.boundary.keys())

        xllcorner = domain.geo_reference.get_xllcorner()
        yllcorner = domain.geo_reference.get_yllcorner()        

        # Make ordering unique #FIXME: should this happen in domain.py?
        boundary_keys.sort()

        # Record ordering #FIXME: should this also happen in domain.py?
        self.boundary_indices = {}
        for i, (vol_id, edge_id) in enumerate(boundary_keys):
            base_index = 3*vol_id
            x0, y0 = V[base_index, :]
            x1, y1 = V[base_index+1, :]
            x2, y2 = V[base_index+2, :]
            
            # Compute midpoints
            if edge_id == 0: m = num.array([(x1 + x2)/2, (y1 + y2)/2])
            if edge_id == 1: m = num.array([(x0 + x2)/2, (y0 + y2)/2])
            if edge_id == 2: m = num.array([(x1 + x0)/2, (y1 + y0)/2])

            # Convert to absolute UTM coordinates
            m[0] += xllcorner
            m[1] += yllcorner
            
            # Register point and index
            self.midpoint_coordinates[i,:] = m

            # Register index of this boundary edge for use with evaluate
            self.boundary_indices[(vol_id, edge_id)] = i


        if verbose: log.critical('Initialise file_function')
        self.F = file_function(filename, domain,
	                           interpolation_points=self.midpoint_coordinates,
                               time_thinning=time_thinning,
                               use_cache=use_cache, 
                               verbose=verbose)
        self.domain = domain

        # Test

        # Here we'll flag indices outside the mesh as a warning
        # as suggested by Joaquim Luis in sourceforge posting
        # November 2007
        # We won't make it an error as it is conceivable that
        # only part of mesh boundary is actually used with a given
        # file boundary sww file. 
        if (hasattr(self.F, 'indices_outside_mesh') and
               len(self.F.indices_outside_mesh)) > 0:
            msg = 'WARNING: File_boundary has points outside the mesh '
            msg += 'given in %s. ' % filename
            msg += 'See warning message issued by Interpolation_function '
            msg += 'for details (should appear above somewhere if '
            msg += 'verbose is True).\n'
            msg += 'This is perfectly OK as long as the points that are '
            msg += 'outside aren\'t used on the actual boundary segment.'
            if verbose is True:            
                log.critical(msg)
            #raise Exception(msg)

        q = self.F(0, point_id=0)

        d = len(domain.conserved_quantities)
        msg = 'Values specified in file %s must be ' % filename
        msg += 'a list or an array of length %d' % d
        assert len(q) == d, msg


    def __repr__(self):
        return 'File boundary'


    def evaluate(self, vol_id=None, edge_id=None):
        """Return linearly interpolated values based on domain time
	       at midpoint of segment defined by vol_id and edge_id.
        """

        q = self.domain.get_conserved_quantities(vol_id, edge=edge_id)
        t = self.domain.get_time()

        if vol_id is not None and edge_id is not None:
            i = self.boundary_indices[vol_id, edge_id]
            res = self.F(t, point_id=i)

            if res == NAN:
                x,y = self.midpoint_coordinates[i,:]
                msg = 'NAN value found in file_boundary at '
                msg += 'point id #%d: (%.2f, %.2f).\n' % (i, x, y)

                if (hasattr(self.F, 'indices_outside_mesh') and
                       len(self.F.indices_outside_mesh) > 0):
                    # Check if NAN point is due it being outside
                    # boundary defined in sww file.

                    if i in self.F.indices_outside_mesh:
                        msg += 'This point refers to one outside the '
                        msg += 'mesh defined by the file %s.\n' % self.F.filename
                        msg += 'Make sure that the file covers '
                        msg += 'the boundary segment it is assigned to '
                        msg += 'in set_boundary.'
                    else:
                        msg += 'This point is inside the mesh defined '
                        msg += 'the file %s.\n' % self.F.filename
                        msg += 'Check this file for NANs.'
                raise Exception(msg)
            
            q[0] = res[0] # Take stage, leave momentum alone
            return q
            #return res 
        else:
            # raise 'Boundary call without point_id not implemented'
            # FIXME: What should the semantics be?
            return self.F(t)



