"""
Environmental forcing functions, such as wind and rainfall.

Constraints: See GPL license in the user guide
Version: 1.0 ($Revision: 7731 $)
ModifiedBy:
    $Author: hudson $
    $Date: 2010-05-18 14:54:05 +1000 (Tue, 18 May 2010) $
"""


from anuga.abstract_2d_finite_volumes.neighbour_mesh import segment_midpoints
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.fit_interpolate.interpolate import Modeltime_too_early, \
                                              Modeltime_too_late
from anuga.geometry.polygon import is_inside_polygon, inside_polygon, \
                                    polygon_area
from types import IntType, FloatType
from anuga.geospatial_data.geospatial_data import ensure_geospatial

from warnings import warn
import numpy as num



def check_forcefield(f):
    """Check that force object is as expected.
    
    Check that f is either:
    1: a callable object f(t,x,y), where x and y are vectors
       and that it returns an array or a list of same length
       as x and y
    2: a scalar
    """

    if callable(f):
        N = 3
        x = num.ones(3, num.float)
        y = num.ones(3, num.float)
        try:
            q = f(1.0, x=x, y=y)
        except Exception, e:
            msg = 'Function %s could not be executed:\n%s' %(f, e)
            # FIXME: Reconsider this semantics
            raise Exception, msg

        try:
            q = num.array(q, num.float)
        except:
            msg = ('Return value from vector function %s could not '
                   'be converted into a numeric array of floats.\nSpecified '
                   'function should return either list or array.' % f)
            raise Exception, msg

        # Is this really what we want?
        # info is "(func name, filename, defining line)"
        func_info = (f.func_name, f.func_code.co_filename,
                     f.func_code.co_firstlineno)
        func_msg = 'Function %s (defined in %s, line %d)' % func_info
        try:
            result_len = len(q)
        except:
            msg = '%s must return vector' % func_msg
            self.fail(msg)
        msg = '%s must return vector of length %d' % (func_msg, N)
        assert result_len == N, msg
    else:
        try:
            f = float(f)
        except:
            msg = ('Force field %s must be a scalar value coercible to float.'
                   % str(f))
            raise Exception, msg

    return f



class Wind_stress:
    """Apply wind stress to water momentum in terms of
    wind speed [m/s] and wind direction [degrees]
    """
    def __init__(self, *args, **kwargs):
        """Initialise windfield from wind speed s [m/s]
        and wind direction phi [degrees]

        Inputs v and phi can be either scalars or Python functions, e.g.

        W = Wind_stress(10, 178)

        #FIXME - 'normal' degrees are assumed for now, i.e. the
        vector (1,0) has zero degrees.
        We may need to convert from 'compass' degrees later on and also
        map from True north to grid north.

        Arguments can also be Python functions of t,x,y as in

        def speed(t,x,y):
            ...
            return s

        def angle(t,x,y):
            ...
            return phi

        where x and y are vectors.

        and then pass the functions in

        W = Wind_stress(speed, angle)

        The instantiated object W can be appended to the list of
        forcing_terms as in

        Alternatively, one vector valued function for (speed, angle)
        can be applied, providing both quantities simultaneously.
        As in
        W = Wind_stress(F), where returns (speed, angle) for each t.

        domain.forcing_terms.append(W)
        """

        from anuga.config import rho_a, rho_w, eta_w

        if len(args) == 2:
            s = args[0]
            phi = args[1]
        elif len(args) == 1:
            # Assume vector function returning (s, phi)(t,x,y)
            vector_function = args[0]
            s = lambda t,x,y: vector_function(t,x=x,y=y)[0]
            phi = lambda t,x,y: vector_function(t,x=x,y=y)[1]
        else:
           # Assume info is in 2 keyword arguments
           if len(kwargs) == 2:
               s = kwargs['s']
               phi = kwargs['phi']
           else:
               raise Exception, 'Assumes two keyword arguments: s=..., phi=....'

        self.speed = check_forcefield(s)
        self.phi = check_forcefield(phi)

        self.const = eta_w*rho_a/rho_w

    ##
    # @brief 'execute' this class instance.
    # @param domain 
    def __call__(self, domain):
        """Evaluate windfield based on values found in domain"""

        from math import pi, cos, sin, sqrt

        xmom_update = domain.quantities['xmomentum'].explicit_update
        ymom_update = domain.quantities['ymomentum'].explicit_update

        N = len(domain)    # number_of_triangles
        t = domain.time

        if callable(self.speed):
            xc = domain.get_centroid_coordinates()
            s_vec = self.speed(t, xc[:,0], xc[:,1])
        else:
            # Assume s is a scalar
            try:
                s_vec = self.speed * num.ones(N, num.float)
            except:
                msg = 'Speed must be either callable or a scalar: %s' %self.s
                raise msg

        if callable(self.phi):
            xc = domain.get_centroid_coordinates()
            phi_vec = self.phi(t, xc[:,0], xc[:,1])
        else:
            # Assume phi is a scalar

            try:
                phi_vec = self.phi * num.ones(N, num.float)
            except:
                msg = 'Angle must be either callable or a scalar: %s' %self.phi
                raise msg

        assign_windfield_values(xmom_update, ymom_update,
                                s_vec, phi_vec, self.const)


##
# @brief Assign wind field values
# @param xmom_update 
# @param ymom_update
# @param s_vec 
# @param phi_vec 
# @param const 
def assign_windfield_values(xmom_update, ymom_update,
                            s_vec, phi_vec, const):
    """Python version of assigning wind field to update vectors.
    A C version also exists (for speed)
    """

    from math import pi, cos, sin, sqrt

    N = len(s_vec)
    for k in range(N):
        s = s_vec[k]
        phi = phi_vec[k]

        # Convert to radians
        phi = phi*pi/180

        # Compute velocity vector (u, v)
        u = s*cos(phi)
        v = s*sin(phi)

        # Compute wind stress
        S = const * sqrt(u**2 + v**2)
        xmom_update[k] += S*u
        ymom_update[k] += S*v


##
# @brief A class for a general explicit forcing term.
class General_forcing:
    """General explicit forcing term for update of quantity

    This is used by Inflow and Rainfall for instance


    General_forcing(quantity_name, rate, center, radius, polygon)

    domain:     ANUGA computational domain
    quantity_name: Name of quantity to update.
                   It must be a known conserved quantity.

    rate [?/s]: Total rate of change over the specified area.
                This parameter can be either a constant or a
                function of time. Positive values indicate increases,
                negative values indicate decreases.
                Rate can be None at initialisation but must be specified
                before forcing term is applied (i.e. simulation has started).

    center [m]: Coordinates at center of flow point
    radius [m]: Size of circular area
    polygon:    Arbitrary polygon
    default_rate: Rate to be used if rate fails (e.g. if model time exceeds its data)
                  Admissible types: None, constant number or function of t


    Either center, radius or polygon can be specified but not both.
    If neither are specified the entire domain gets updated.
    All coordinates to be specified in absolute UTM coordinates (x, y) assuming the zone of domain.

    Inflow or Rainfall for examples of use
    """


    # FIXME (AnyOne) : Add various methods to allow spatial variations

    ##
    # @brief Create an instance of this forcing term.
    # @param domain 
    # @param quantity_name 
    # @param rate 
    # @param center 
    # @param radius 
    # @param polygon 
    # @param default_rate 
    # @param verbose 
    def __init__(self,
                 domain,
                 quantity_name,
                 rate=0.0,
                 center=None,
                 radius=None,
                 polygon=None,
                 default_rate=None,
                 verbose=False):

        from math import pi, cos, sin

        if center is None:
            msg = 'I got radius but no center.'
            assert radius is None, msg

        if radius is None:
            msg += 'I got center but no radius.'
            assert center is None, msg

        self.domain = domain
        self.quantity_name = quantity_name
        self.rate = rate
        self.center = ensure_numeric(center)
        self.radius = radius
        self.polygon = polygon
        self.verbose = verbose
        self.value = 0.0    # Can be used to remember value at
                            # previous timestep in order to obtain rate

        # Get boundary (in absolute coordinates)
        bounding_polygon = domain.get_boundary_polygon()

        # Update area if applicable
        if center is not None and radius is not None:
            assert len(center) == 2
            msg = 'Polygon cannot be specified when center and radius are'
            assert polygon is None, msg

            # Check that circle center lies within the mesh.
            msg = 'Center %s specified for forcing term did not' % str(center)
            msg += 'fall within the domain boundary.'
            assert is_inside_polygon(center, bounding_polygon), msg

            # Check that circle periphery lies within the mesh.
            N = 100
            periphery_points = []
            for i in range(N):
                theta = 2*pi*i/100

                x = center[0] + radius*cos(theta)
                y = center[1] + radius*sin(theta)

                periphery_points.append([x,y])

            for point in periphery_points:
                msg = 'Point %s on periphery for forcing term' % str(point)
                msg += ' did not fall within the domain boundary.'
                assert is_inside_polygon(point, bounding_polygon), msg

        if polygon is not None:
            # Check that polygon lies within the mesh.
            for point in self.polygon:
                msg = 'Point %s in polygon for forcing term' % str(point)
                msg += ' did not fall within the domain boundary.'
                assert is_inside_polygon(point, bounding_polygon), msg

        # Pointer to update vector
        self.update = domain.quantities[self.quantity_name].explicit_update

        # Determine indices in flow area
        N = len(domain)
        points = domain.get_centroid_coordinates(absolute=True)

        # Calculate indices in exchange area for this forcing term
        self.exchange_indices = None
        if self.center is not None and self.radius is not None:
            # Inlet is circular
            inlet_region = 'center=%s, radius=%s' % (self.center, self.radius)

            self.exchange_indices = []
            for k in range(N):
                x, y = points[k,:]    # Centroid

                c = self.center
                if ((x-c[0])**2+(y-c[1])**2) < self.radius**2:
                    self.exchange_indices.append(k)

        if self.polygon is not None:
            # Inlet is polygon
            inlet_region = 'polygon=%s' % (self.polygon) 
            self.exchange_indices = inside_polygon(points, self.polygon)

        if self.exchange_indices is None:
            self.exchange_area = polygon_area(bounding_polygon)
        else:    
            if len(self.exchange_indices) == 0:
                msg = 'No triangles have been identified in '
                msg += 'specified region: %s' % inlet_region
                raise Exception, msg

            # Compute exchange area as the sum of areas of triangles identified
            # by circle or polygon
            self.exchange_area = 0.0
            for i in self.exchange_indices:
                self.exchange_area += domain.areas[i]
            

        msg = 'Exchange area in forcing term'
        msg += ' has area = %f' %self.exchange_area
        assert self.exchange_area > 0.0            
            
                

            
        # Check and store default_rate
        msg = ('Keyword argument default_rate must be either None '
               'or a function of time.\nI got %s.' % str(default_rate))
        assert (default_rate is None or
                type(default_rate) in [IntType, FloatType] or
                callable(default_rate)), msg

        if default_rate is not None:
            # If it is a constant, make it a function
            if not callable(default_rate):
                tmp = default_rate
                default_rate = lambda t: tmp

            # Check that default_rate is a function of one argument
            try:
                default_rate(0.0)
            except:
                raise Exception, msg

        self.default_rate = default_rate
        self.default_rate_invoked = False    # Flag

    ##
    # @brief Execute this instance.
    # @param domain 
    def __call__(self, domain):
        """Apply inflow function at time specified in domain, update stage"""

        # Call virtual method allowing local modifications
        t = domain.get_time()
        try:
            rate = self.update_rate(t)
        except Modeltime_too_early, e:
            raise Modeltime_too_early, e
        except Modeltime_too_late, e:
            if self.default_rate is None:
                msg = '%s: ANUGA is trying to run longer than specified data.\n' %str(e)
                msg += 'You can specify keyword argument default_rate in the '
                msg += 'forcing function to tell it what to do in the absence of time data.'
                raise Modeltime_too_late, msg    
            else:
                # Pass control to default rate function
                rate = self.default_rate(t)

                if self.default_rate_invoked is False:
                    # Issue warning the first time
                    msg = ('%s\n'
                           'Instead I will use the default rate: %s\n'
                           'Note: Further warnings will be supressed'
                           % (str(e), str(self.default_rate)))
                    warn(msg)

                    # FIXME (Ole): Replace this crude flag with
                    # Python's ability to print warnings only once.
                    # See http://docs.python.org/lib/warning-filter.html
                    self.default_rate_invoked = True

        if rate is None:
            msg = ('Attribute rate must be specified in General_forcing '
                   'or its descendants before attempting to call it')
            raise Exception, msg

        # Now rate is a number
        if self.verbose is True:
            log.critical('Rate of %s at time = %.2f = %f'
                         % (self.quantity_name, domain.get_time(), rate))

        if self.exchange_indices is None:
            self.update[:] += rate
        else:
            # Brute force assignment of restricted rate
            for k in self.exchange_indices:
                self.update[k] += rate

    ##
    # @brief Update the internal rate.
    # @param t A callable or scalar used to set the rate.
    # @return The new rate.
    def update_rate(self, t):
        """Virtual method allowing local modifications by writing an
        overriding version in descendant
        """

        if callable(self.rate):
            rate = self.rate(t)
        else:
            rate = self.rate

        return rate

    ##
    # @brief Get values for the specified quantity.
    # @param quantity_name Name of the quantity of interest.
    # @return The value(s) of the quantity.
    # @note If 'quantity_name' is None, use self.quantity_name.
    def get_quantity_values(self, quantity_name=None):
        """Return values for specified quantity restricted to opening

        Optionally a quantity name can be specified if values from another
        quantity is sought
        """

        if quantity_name is None:
            quantity_name = self.quantity_name

        q = self.domain.quantities[quantity_name]
        return q.get_values(location='centroids',
                            indices=self.exchange_indices)

    ##
    # @brief Set value for the specified quantity.
    # @param val The value object used to set value.
    # @param quantity_name Name of the quantity of interest.
    # @note If 'quantity_name' is None, use self.quantity_name.
    def set_quantity_values(self, val, quantity_name=None):
        """Set values for specified quantity restricted to opening

        Optionally a quantity name can be specified if values from another
        quantity is sought
        """

        if quantity_name is None:
            quantity_name = self.quantity_name

        q = self.domain.quantities[self.quantity_name]
        q.set_values(val,
                     location='centroids',
                     indices=self.exchange_indices)


##
# @brief A class for rainfall forcing function.
# @note Inherits from General_forcing.
class Rainfall(General_forcing):
    """Class Rainfall - general 'rain over entire domain' forcing term.

    Used for implementing Rainfall over the entire domain.

        Current Limited to only One Gauge..

        Need to add Spatial Varying Capability
        (This module came from copying and amending the Inflow Code)

    Rainfall(rain)

    domain
    rain [mm/s]:  Total rain rate over the specified domain.
                  NOTE: Raingauge Data needs to reflect the time step.
                  IE: if Gauge is mm read at a time step, then the input
                  here is as mm/(timeStep) so 10mm in 5minutes becomes
                  10/(5x60) = 0.0333mm/s.

                  This parameter can be either a constant or a
                  function of time. Positive values indicate inflow,
                  negative values indicate outflow.
                  (and be used for Infiltration - Write Seperate Module)
                  The specified flow will be divided by the area of
                  the inflow region and then applied to update the
                  stage quantity.

    polygon: Specifies a polygon to restrict the rainfall.

    Examples
    How to put them in a run File...

    #------------------------------------------------------------------------
    # Setup specialised forcing terms
    #------------------------------------------------------------------------
    # This is the new element implemented by Ole and Rudy to allow direct
    # input of Rainfall in mm/s

    catchmentrainfall = Rainfall(rain=file_function('Q100_2hr_Rain.tms'))
                        # Note need path to File in String.
                        # Else assumed in same directory

    domain.forcing_terms.append(catchmentrainfall)
    """

    ##
    # @brief Create an instance of the class.
    # @param domain Domain of interest.
    # @param rate Total rain rate over the specified domain (mm/s).
    # @param center 
    # @param radius 
    # @param polygon Polygon  to restrict rainfall.
    # @param default_rate 
    # @param verbose True if this instance is to be verbose.
    def __init__(self,
                 domain,
                 rate=0.0,
                 center=None,
                 radius=None,
                 polygon=None,
                 default_rate=None,
                 verbose=False):

        # Converting mm/s to m/s to apply in ANUGA)
        if callable(rate):
            rain = lambda t: rate(t)/1000.0
        else:
            rain = rate/1000.0

        if default_rate is not None:
            if callable(default_rate):
                default_rain = lambda t: default_rate(t)/1000.0
            else:
                default_rain = default_rate/1000.0
        else:
            default_rain = None

            
            
        General_forcing.__init__(self,
                                 domain,
                                 'stage',
                                 rate=rain,
                                 center=center,
                                 radius=radius,
                                 polygon=polygon,
                                 default_rate=default_rain,
                                 verbose=verbose)


##
# @brief A class for inflow (rain and drain) forcing function.
# @note Inherits from General_forcing.
class Inflow(General_forcing):
    """Class Inflow - general 'rain and drain' forcing term.

    Useful for implementing flows in and out of the domain.

    Inflow(flow, center, radius, polygon)

    domain
    rate [m^3/s]: Total flow rate over the specified area.
                  This parameter can be either a constant or a
                  function of time. Positive values indicate inflow,
                  negative values indicate outflow.
                  The specified flow will be divided by the area of
                  the inflow region and then applied to update stage.
    center [m]: Coordinates at center of flow point
    radius [m]: Size of circular area
    polygon:    Arbitrary polygon.

    Either center, radius or polygon must be specified

    Examples

    # Constant drain at 0.003 m^3/s.
    # The outflow area is 0.07**2*pi=0.0154 m^2
    # This corresponds to a rate of change of 0.003/0.0154 = 0.2 m/s
    #
    Inflow((0.7, 0.4), 0.07, -0.003)


    # Tap turning up to a maximum inflow of 0.0142 m^3/s.
    # The inflow area is 0.03**2*pi = 0.00283 m^2
    # This corresponds to a rate of change of 0.0142/0.00283 = 5 m/s
    # over the specified area
    Inflow((0.5, 0.5), 0.03, lambda t: min(0.01*t, 0.0142))


    #------------------------------------------------------------------------
    # Setup specialised forcing terms
    #------------------------------------------------------------------------
    # This is the new element implemented by Ole to allow direct input
    # of Inflow in m^3/s

    hydrograph = Inflow(center=(320, 300), radius=10,
                        rate=file_function('Q/QPMF_Rot_Sub13.tms'))

    domain.forcing_terms.append(hydrograph)
    """

    ##
    # @brief Create an instance of the class.
    # @param domain Domain of interest.
    # @param rate Total rain rate over the specified domain (mm/s).
    # @param center 
    # @param radius 
    # @param polygon Polygon  to restrict rainfall.
    # @param default_rate 
    # @param verbose True if this instance is to be verbose.
    def __init__(self,
                 domain,
                 rate=0.0,
                 center=None,
                 radius=None,
                 polygon=None,
                 default_rate=None,
                 verbose=False):
        # Create object first to make area is available
        General_forcing.__init__(self,
                                 domain,
                                 'stage',
                                 rate=rate,
                                 center=center,
                                 radius=radius,
                                 polygon=polygon,
                                 default_rate=default_rate,
                                 verbose=verbose)

    ##
    # @brief Update the instance rate.
    # @param t New rate object.
    def update_rate(self, t):
        """Virtual method allowing local modifications by writing an
        overriding version in descendant

        This one converts m^3/s to m/s which can be added directly
        to 'stage' in ANUGA
        """

        if callable(self.rate):
            _rate = self.rate(t)/self.exchange_area
        else:
            _rate = self.rate/self.exchange_area

        return _rate


##
# @brief A class for creating cross sections.
# @note Inherits from General_forcing.
class Cross_section:
    """Class Cross_section - a class to setup a cross section from
    which you can then calculate flow and energy through cross section


    Cross_section(domain, polyline)

    domain:
    polyline: Representation of desired cross section - it may contain
              multiple sections allowing for complex shapes. Assume
              absolute UTM coordinates.
              Format [[x0, y0], [x1, y1], ...]
    verbose: 
    """

    ##
    # @brief Create an instance of the class.
    # @param domain Domain of interest.
    # @param polyline Polyline defining cross section
    # @param verbose True if this instance is to be verbose.
    def __init__(self,
                 domain,
                 polyline=None,
                 verbose=False):
        
        self.domain = domain
        self.polyline = polyline
        self.verbose = verbose
        
        # Find all intersections and associated triangles.
        self.segments = self.domain.get_intersecting_segments(self.polyline,
                                                              use_cache=True,
                                                              verbose=self.verbose)
        
        # Get midpoints
        self.midpoints = segment_midpoints(self.segments)

        # Make midpoints Geospatial instances
        self.midpoints = ensure_geospatial(self.midpoints, self.domain.geo_reference)

    ##
    # @brief set verbose mode
    def set_verbose(self,verbose=True):
        """Set verbose mode true or flase
        """

        self.verbose=verbose

    ##
    # @brief calculate current flow through cross section
    def get_flow_through_cross_section(self):
        """ Output: Total flow [m^3/s] across cross section.
        """

        # Get interpolated values
        xmomentum = self.domain.get_quantity('xmomentum')
        ymomentum = self.domain.get_quantity('ymomentum')

        uh = xmomentum.get_values(interpolation_points=self.midpoints,
                                  use_cache=True)
        vh = ymomentum.get_values(interpolation_points=self.midpoints,
                                  use_cache=True)

        # Compute and sum flows across each segment
        total_flow = 0
        for i in range(len(uh)):
            # Inner product of momentum vector with segment normal [m^2/s]
            normal = self.segments[i].normal
            normal_momentum = uh[i]*normal[0] + vh[i]*normal[1]

            # Flow across this segment [m^3/s]
            segment_flow = normal_momentum*self.segments[i].length

            # Accumulate
            total_flow += segment_flow

        return total_flow
 

    ##
    # @brief calculate current energy flow through cross section
    def get_energy_through_cross_section(self, kind='total'):
        """Obtain average energy head [m] across specified cross section.

        Output:
            E: Average energy [m] across given segments for all stored times.

        The average velocity is computed for each triangle intersected by
        the polyline and averaged weighted by segment lengths.

        The typical usage of this function would be to get average energy of
        flow in a channel, and the polyline would then be a cross section
        perpendicular to the flow.

        #FIXME (Ole) - need name for this energy reflecting that its dimension
        is [m].
        """

        from anuga.config import g, epsilon, velocity_protection as h0
        
        # Get interpolated values
        stage = self.domain.get_quantity('stage')
        elevation = self.domain.get_quantity('elevation')
        xmomentum = self.domain.get_quantity('xmomentum')
        ymomentum = self.domain.get_quantity('ymomentum')

        w = stage.get_values(interpolation_points=self.midpoints, use_cache=True)
        z = elevation.get_values(interpolation_points=self.midpoints, use_cache=True)
        uh = xmomentum.get_values(interpolation_points=self.midpoints,
                                  use_cache=True)
        vh = ymomentum.get_values(interpolation_points=self.midpoints,
                                  use_cache=True)
        h = w-z                # Depth

        # Compute total length of polyline for use with weighted averages
        total_line_length = 0.0
        for segment in self.segments:
            total_line_length += segment.length

        # Compute and sum flows across each segment
        average_energy = 0.0
        for i in range(len(w)):
            # Average velocity across this segment
            if h[i] > epsilon:
                # Use protection against degenerate velocities
                u = uh[i]/(h[i] + h0/h[i])
                v = vh[i]/(h[i] + h0/h[i])
            else:
                u = v = 0.0

            speed_squared = u*u + v*v
            kinetic_energy = 0.5*speed_squared/g

            if kind == 'specific':
                segment_energy = h[i] + kinetic_energy
            elif kind == 'total':
                segment_energy = w[i] + kinetic_energy
            else:
                msg = 'Energy kind must be either "specific" or "total".'
                msg += ' I got %s' %kind

            # Add to weighted average
            weigth = self.segments[i].length/total_line_length
            average_energy += segment_energy*weigth

        return average_energy

