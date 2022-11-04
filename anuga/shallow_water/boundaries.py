""" Boundary conditions specific to the shallow water wave equation

Title: ANUGA boundaries with dependencies on shallow_water_domain
Author:
    Ole Nielsen, Ole.Nielsen@ga.gov.au
    Stephen Roberts, Stephen.Roberts@anu.edu.au
    Duncan Gray, Duncan.Gray@ga.gov.au
    Gareth Davies, gareth.davies.ga.code@gmail.com
CreationDate: 2010
Description::
    This module contains boundary functions for ANUGA that are specific to the shallow water Domain class.  
Constraints: See GPL license in the user guide
Version: 1.0 ($Revision: 7731 $)
ModifiedBy::
    Author: hudson
    Date: 2010-05-18 14:54:05 +1000 (Tue, 18 May 2010)

"""

from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Boundary, File_boundary
import numpy as num

import anuga.utilities.log as log
from anuga.fit_interpolate.interpolate import Modeltime_too_late
from anuga.fit_interpolate.interpolate import Modeltime_too_early
from anuga.config import g as gravity
     
from .shallow_water_ext import rotate

#from anuga.utilities import compile
#if compile.can_use_C_extension('shallow_water_ext.c'):
#    # Underlying C implementations can be accessed
#    from shallow_water_ext import rotate
#else:
#    msg = 'C implementations could not be accessed by %s.\n ' % __file__
#    msg += 'Make sure compile_all.py has been run as described in '
#    msg += 'the ANUGA installation guide.'
#    raise Exception, msg


class Reflective_boundary(Boundary):
    """Reflective boundary condition object
    
    Reflective boundary returns same conserved quantities as
    those present in its neighbour volume but with normal momentum reflected.
    """

    def __init__(self, domain=None):
        """Create boundary condition object
        
        :param domain: domain on which to apply BC
        
        Example: 
        
        Set all the tagged boundaries to use the Reflective boundaries
        
        >>> domain = anuga.rectangular_cross_domain(10, 10) 
        >>> BC = anuga.Reflective_boundary(domain)
        >>> domain.set_boundary({'left': BC, 'right': BC, 'top': BC, 'bottom': BC})
        
        """



        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for reflective boundary'
            raise Exception(msg)

        # Handy shorthands
        self.stage = domain.quantities['stage'].edge_values
        self.xmom = domain.quantities['xmomentum'].edge_values
        self.ymom = domain.quantities['ymomentum'].edge_values
        self.normals = domain.normals

        self.conserved_quantities = num.zeros(3, float)

    def __repr__(self):
        return 'Reflective_boundary'

    def evaluate(self, vol_id, edge_id):
        """Calculate BC associated to specified edge

        :param int vol_id: Triangle ID
        :param int edge_id: Edge opposite to Vertex ID  
  
        """

        q = self.conserved_quantities
        q[0] = self.stage[vol_id, edge_id]
        q[1] = self.xmom[vol_id, edge_id]
        q[2] = self.ymom[vol_id, edge_id]

        normal = self.normals[vol_id, 2*edge_id:2*edge_id+2]

        r = rotate(q, normal, direction = 1)
        r[1] = -r[1]
        q = rotate(r, normal, direction = -1)

        return q



    def evaluate_segment(self, domain, segment_edges):
        """Apply BC on the boundary edges defined by segment_edges

        :param domain: Apply BC on this domain
        :param segment_edges: List of boundary cells on which to apply BC

        """

        if segment_edges is None:
            return
        if domain is None:
            return


        ids = segment_edges
        vol_ids  = domain.boundary_cells[ids]
        edge_ids = domain.boundary_edges[ids]

        Stage = domain.quantities['stage']
        Elev  = domain.quantities['elevation']
        Height= domain.quantities['height']
        Xmom  = domain.quantities['xmomentum']
        Ymom  = domain.quantities['ymomentum']
        Xvel  = domain.quantities['xvelocity']
        Yvel  = domain.quantities['yvelocity']

        Normals = domain.normals

        #print vol_ids
        #print edge_ids
        #Normals.reshape((4,3,2))
        #print Normals.shape
        #print Normals[vol_ids, 2*edge_ids]
        #print Normals[vol_ids, 2*edge_ids+1]
        
        n1  = Normals[vol_ids,2*edge_ids]
        n2  = Normals[vol_ids,2*edge_ids+1]

        # Transfer these quantities to the boundary array
        Stage.boundary_values[ids]  = Stage.edge_values[vol_ids,edge_ids]
        Elev.boundary_values[ids]   = Elev.edge_values[vol_ids,edge_ids]
        Height.boundary_values[ids] = Height.edge_values[vol_ids,edge_ids]

        # Rotate and negate Momemtum
        q1 = Xmom.edge_values[vol_ids,edge_ids]
        q2 = Ymom.edge_values[vol_ids,edge_ids]

        r1 = -q1*n1 - q2*n2
        r2 = -q1*n2 + q2*n1

        Xmom.boundary_values[ids] = n1*r1 - n2*r2
        Ymom.boundary_values[ids] = n2*r1 + n1*r2

        # Rotate and negate Velocity
        q1 = Xvel.edge_values[vol_ids,edge_ids]
        q2 = Yvel.edge_values[vol_ids,edge_ids]

        r1 = q1*n1 + q2*n2
        r2 = q1*n2 - q2*n1

        Xvel.boundary_values[ids] = n1*r1 - n2*r2
        Yvel.boundary_values[ids] = n2*r1 + n1*r2



class Transmissive_momentum_set_stage_boundary(Boundary):
    """ Bounday condition object that returns transmissive momentum and sets stage
    
    Returns same momentum conserved quantities as
    those present in its neighbour volume.
    Sets stage by specifying a function f of time which may either be a
    vector function or a scalar function
    """

    def __init__(self, domain=None, function=None):
        """Create boundary condition object.
        
        :param domain: domain on which to apply BC
        :param function: function to set stage
        
        Example: Set all the tagged boundaries to use the 
        
        >>> domain = anuga.rectangular_cross_domain(10, 10) 
        >>> def waveform(t):
        >>>    return sea_level + normalized_amplitude/cosh(t-25)**2
        >>> BC = anuga.Transmissive_momentum_set_stage_boundary(domain, waveform)
        >>> domain.set_boundary({'left': BC, 'right': BC, 'top': BC, 'bottom': BC})
        
        """

        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for this type boundary'
            raise Exception(msg)

        if function is None:
            msg = 'Function must be specified for this type boundary'
            raise Exception(msg)

        self.domain = domain

        if isinstance(function, (int, float)):
            tmp = function
            function = lambda t: tmp
            
        self.function = function


    def __repr__(self):
        """ Return a representation of this object. """

        return 'Transmissive_momentum_set_stage_boundary(%s)' % self.domain

    def evaluate(self, vol_id, edge_id):
        """Transmissive momentum set stage boundaries return the edge momentum
        values of the volume they serve.

        vol_id is volume id
        edge_id is the edge within the volume
        """

        q = self.domain.get_conserved_quantities(vol_id, edge = edge_id)
        t = self.domain.get_time()

        if hasattr(self.function, 'time'):
            # Roll boundary over if time exceeds
            while t > self.function.time[-1]:
                msg = 'WARNING: domain time %.2f has exceeded' % t
                msg += 'time provided in '
                msg += 'transmissive_momentum_set_stage_boundary object.\n'
                msg += 'I will continue, reusing the object from t==0'
                log.critical(msg)
                t -= self.function.time[-1]

        value = self.function(t)
        try:
            x = float(value)
        except:
            x = float(value[0])

        q[0] = x
           
        return q

        # FIXME: Consider this (taken from File_boundary) to allow
        # spatial variation
        # if vol_id is not None and edge_id is not None:
        #     i = self.boundary_indices[ vol_id, edge_id ]
        #     return self.F(t, point_id = i)
        # else:
        #     return self.F(t)


class Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(Boundary):
    """Bounday condition object that returns transmissive normal momentum and sets stage
    
    Returns the same normal momentum as that 
    present in neighbour volume edge. Zero out the tangential momentum. 
    Sets stage by specifying a function f of time which may either be a
    vector function or a scalar function

    """

    def __init__(self, domain=None, function=None, default_boundary=0.0):
        """Create boundary condition object.
        
        :param domain: domain on which to apply BC
        :param function: function to set stage
        :param float default_boundary: 
        
        Example: Set all the tagged boundaries to use the BC
        
        >>> domain = anuga.rectangular_cross_domain(10, 10) 
        >>> def waveform(t):
        >>>    return sea_level + normalized_amplitude/cosh(t-25)**2
        >>> BC = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain, waveform)
        >>> domain.set_boundary({'left': BC, 'right': BC, 'top': BC, 'bottom': BC})
        
        """


        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for this type boundary'
            raise Exception(msg)

        if function is None:
            msg = 'Function must be specified for this type boundary'
            raise Exception(msg)

        self.domain = domain
        self.function = function
        self.default_boundary = default_boundary


    def __repr__(self):
        """ Return a representation of this instance. """
        msg = 'Transmissive_n_momentum_zero_t_momentum_set_stage_boundary'
        msg += '(%s)' % self.domain
        return msg


    def evaluate(self, vol_id, edge_id):
        """Transmissive_n_momentum_zero_t_momentum_set_stage_boundary
        return the edge momentum values of the volume they serve.
        """

        q = self.domain.get_conserved_quantities(vol_id, edge = edge_id)

        normal = self.domain.get_normal(vol_id, edge_id)


        ## t = self.domain.get_time()

        ## if hasattr(self.function, 'time'):
        ##     # Roll boundary over if time exceeds
        ##     while t > self.function.time[-1]:
        ##         msg = 'WARNING: domain time %.2f has exceeded' % t
        ##         msg += 'time provided in '
        ##         msg += 'transmissive_momentum_set_stage_boundary object.\n'
        ##         msg += 'I will continue, reusing the object from t==0'
        ##         log.critical(msg)
        ##         t -= self.function.time[-1]

        ## value = self.function(t)
        ## try:
        ##     x = float(value)
        ## except:
        ##     x = float(value[0])

        value = self.get_boundary_values()
        try:
            x = float(value)
        except:
            x = float(value[0])

        q[0] = x

        ndotq = (normal[0]*q[1] + normal[1]*q[2])
        q[1] = normal[0]*ndotq
        q[2] = normal[1]*ndotq

        return q


    def evaluate_segment(self, domain, segment_edges): 
        """Apply BC on the boundary edges defined by segment_edges

        :param domain: Apply BC on this domain
        :param segment_edges: List of boundary cells on which to apply BC

        Vectorized form for speed. Gareth Davies 14/07/2016
        
        """
        
        Stage = domain.quantities['stage']
        #Elev  = domain.quantities['elevation']
        #Height= domain.quantities['height']
        Xmom  = domain.quantities['xmomentum']
        Ymom  = domain.quantities['ymomentum']

        ids = segment_edges
        vol_ids  = domain.boundary_cells[ids]
        edge_ids = domain.boundary_edges[ids]
        Normals = domain.normals

        n1  = Normals[vol_ids,2*edge_ids]
        n2  = Normals[vol_ids,2*edge_ids+1]
       
        # Call the boundary function which returns stage 
        value = self.get_boundary_values()
        try:
            x = float(value)
        except:
            x = float(value[0])
       
        # Set stage 
        Stage.boundary_values[ids]  = x
       
        # Compute flux normal to edge
        q1 = Xmom.edge_values[vol_ids,edge_ids]
        q2 = Ymom.edge_values[vol_ids,edge_ids]
        ndotq = n1*q1 + n2*q2

        Xmom.boundary_values[ids] = ndotq * n1
        Ymom.boundary_values[ids] = ndotq * n2


class Transmissive_stage_zero_momentum_boundary(Boundary):
    """BC where stage is same as neighbour volume and momentum to zero.

    Underlying domain must be specified when boundary is instantiated
    """

    def __init__(self, domain=None):
        """ Instantiate a Transmissive (zero momentum) boundary. """

        Boundary.__init__(self)

        if domain is None:
            msg = ('Domain must be specified for '
                   'Transmissive_stage_zero_momentum boundary')
            raise Exception(msg)

        self.domain = domain

    def __repr__(self):
        """ Return a representation of this instance. """
        return 'Transmissive_stage_zero_momentum_boundary(%s)' % self.domain


    def evaluate(self, vol_id, edge_id):
        """Calculate transmissive (zero momentum) results. """

        q = self.domain.get_conserved_quantities(vol_id, edge=edge_id)

        q[1] = q[2] = 0.0
        return q



class Time_stage_zero_momentum_boundary(Boundary):
    """Time dependent boundary returns values for stage
    conserved quantities as a function of time.
    Must specify domain to get access to model time and a function of t
    which must return conserved stage quantities as a function time.
    
    Example:
      B = Time_stage_zero_momentum_boundary(domain, 
                        function=lambda t: (60<t<3660)*2)
      
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
            raise Exception('You must specify a domain to Time_stage_zero_momemtum_boundary')
            
        if function is None:
            raise Exception('You must specify a function to Time_stage_zero_momemtum_boundary')
            
        
        try:
            q = function(0.0)
        except Exception as e:
            msg = 'Function for time stage boundary could not be executed:\n%s' %e
            raise Exception(msg)


        try:
            q = float(q)
        except:
            msg = 'Return value from time boundary function could '
            msg += 'not be converted into a float.\n'
            msg += 'I got %s' %str(q)
            raise Exception(msg)


        self.f = function
        self.domain = domain

    def __repr__(self):
        return 'Time_stage_zero_momemtum_boundary'


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

        #-------------------------------------------------
        # Now update boundary values
        #-------------------------------------------------
        domain.quantities['stage'].boundary_values[ids] = q_bdry
        domain.quantities['xmomentum'].boundary_values[ids] = 0.0
        domain.quantities['ymomentum'].boundary_values[ids] = 0.0        



class Characteristic_stage_boundary(Boundary):
    """Sets the stage via a function and the momentum is determined 
    via assumption of simple incoming wave (uses Riemann invariant)


    Example:

    def waveform(t):
        return sea_level + normalized_amplitude/cosh(t-25)**2

    Bcs = Characteristic_stage_boundary(domain, waveform)

    Underlying domain must be specified when boundary is instantiated
    """

    def __init__(self, domain=None, function=None, default_stage = 0.0):
        """ Instantiate a
            Characteristic_stage_boundary.
            domain is the domain containing the boundary
            function is the function to apply the wave
            default_stage is the assumed stage pre the application of wave
        """

        #raise Exception('This boundary type is not implemented yet')
    
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for this type boundary'
            raise Exception(msg)

        if function is None:
            msg = 'Function must be specified for this type boundary'
            raise Exception(msg)

        self.domain = domain
        self.function = function
        self.default_stage = default_stage

        self.elev   = domain.quantities['elevation']
        self.stage  = domain.quantities['stage']
        #self.height = domain.quantities['height']
        self.xmom   = domain.quantities['xmomentum']
        self.ymom   = domain.quantities['ymomentum']

    def __repr__(self):
        """ Return a representation of this instance. """
        msg = 'Characteristic_stage_boundary '
        msg += '(%s) ' % self.domain
        msg += '(%s) ' % self.default_stage
        return msg


    def evaluate(self, vol_id, edge_id):
        """Calculate reflections (reverse outward momentum).

        vol_id   
        edge_id  
        """
        
        t = self.domain.get_time()


        value = self.function(t)
        try:
            w_outside = float(value)
        except:
            w_outside = float(value[0])

        q = num.zeros(len(self.conserved_quantities), float)

        q[0] = self.stage.edge_values[vol_id, edge_id]
        q[1] = self.xmom.edge_values[vol_id, edge_id]
        q[2] = self.ymom.edge_values[vol_id, edge_id]
        elev = self.elev.edge_values[vol_id, edge_id]

        normal = self.normals[vol_id, 2*edge_id:2*edge_id+2]

        uh_inside  = normal[0]*q[1] + normal[1]*q[2]
        vh_inside  = normal[1]*q[1] - normal[0]*q[2]


        # use elev as elev both inside and outside

        h_outside = max(w_outside - elev,0)
        h_inside  = max(q[0] - elev, 0)

        u_inside  = uh_inside/h_inside


        sqrt_g = gravity**0.5
        sqrt_h_inside = h_inside**0.5
        sqrt_h_outside = h_outside**0.5

        h_m = (0.5*(sqrt_h_inside+sqrt_h_outside) + u_inside/4.0/sqrt_g )**2
        u_m = 0.5*u_inside + sqrt_g*(sqrt_h_inside - sqrt_h_outside)

        uh_m = h_m*u_m
        vh_m = vh_inside

        # if uh_inside > 0.0 then outflow
        if uh_inside > 0.0 :
            vh_m = vh_inside
        else:
            vh_m = 0.0

        if h_inside == 0.0:
            q[0] = w_outside
            q[1] = 0.0
            q[2] = 0.0
        else:
            q[0] = h_m + elev       
            q[1] = uh_m*normal[0] + vh_m*normal[1]
            q[2] = uh_m*normal[1] - vh_m*normal[0]

        return q


    def evaluate_segment(self, domain, segment_edges):
        """Apply BC on the boundary edges defined by
        segment_edges
        """

        Stage = domain.quantities['stage']
        Elev  = domain.quantities['elevation']
        Xmom  = domain.quantities['xmomentum']
        Ymom  = domain.quantities['ymomentum']

        ids = segment_edges
        vol_ids  = domain.boundary_cells[ids]
        edge_ids = domain.boundary_edges[ids]
        Normals  = domain.normals

        n1  = Normals[vol_ids,2*edge_ids]
        n2  = Normals[vol_ids,2*edge_ids+1]
   
        # Get stage value 
        t = self.domain.get_time()
        value = self.function(t)
        try:
            w_outside = float(value)
        except:
            w_outside = float(value[0])

        # Transfer these quantities to the boundary array
        Stage.boundary_values[ids] = Stage.edge_values[vol_ids,edge_ids]
        Xmom.boundary_values[ids]  = Xmom.edge_values[vol_ids,edge_ids]
        Ymom.boundary_values[ids]  = Ymom.edge_values[vol_ids,edge_ids]
        Elev.boundary_values[ids]  = Elev.edge_values[vol_ids,edge_ids]



        h_inside = num.maximum(Stage.boundary_values[ids]-Elev.boundary_values[ids], 0.0)
        w_outside = 0.0*Stage.boundary_values[ids] + w_outside 

        # Do vectorized operations here
        #
        # In dry cells, the values will be ....
        q0_dry = num.where(Elev.boundary_values[ids] <= w_outside, w_outside, Elev.boundary_values[ids])
        q1_dry = 0.0 * Xmom.boundary_values[ids]
        q2_dry = 0.0 * Ymom.boundary_values[ids]
        #
        # and in wet cells, the values will be ...
        # (see 'evaluate' method above for more comments on theory,
        # in particular we assume subcritical flow and a zero outside velocity)
        #
        # (note: When cells are dry, this calculation will throw invalid
        # values, but such values will never be selected to be returned)
        sqrt_g = gravity**0.5
        h_inside  = num.maximum(Stage.boundary_values[ids] - Elev.boundary_values[ids], 0)
        uh_inside = n1 * Xmom.boundary_values[ids] + n2 * Ymom.boundary_values[ids]
        vh_inside = n2 * Xmom.boundary_values[ids] - n1 * Ymom.boundary_values[ids]
        u_inside  = num.where(h_inside>0.0, uh_inside/h_inside, 0.0)

        h_outside = num.maximum(w_outside - Elev.boundary_values[ids], 0)

        sqrt_h_inside = h_inside**0.5
        sqrt_h_outside = h_outside**0.5

        h_m = (0.5*(sqrt_h_inside+sqrt_h_outside) + u_inside/4.0/sqrt_g )**2
        u_m = 0.5*u_inside + sqrt_g*(sqrt_h_inside - sqrt_h_outside)

        uh_m = h_m*u_m

        # if uh_inside > 0.0 then outflow
        vh_m = num.where(uh_inside > 0.0, vh_inside, 0.0)

        w_m = h_m + Elev.boundary_values[ids]

        dry_test = num.logical_or(h_inside == 0.0, h_outside == 0.0)

        q1 = uh_m*n1 + vh_m*n2
        q2 = uh_m*n2 - vh_m*n1

        Stage.boundary_values[ids] = num.where(
            dry_test,
            w_outside,
            w_m 
            )

        Xmom.boundary_values[ids] = num.where(
            dry_test, 
            0.0,
            q1 
            )

        Ymom.boundary_values[ids] = num.where(
            dry_test,
            0.0,
            q2)

class Dirichlet_discharge_boundary(Boundary):
    """ Class for a Dirichlet discharge boundary.

    Sets stage (stage0)
    Sets momentum (wh0) in the inward normal direction.

    Underlying domain must be specified when boundary is instantiated
    """

    def __init__(self, domain=None, stage0=None, wh0=None):
        Boundary.__init__(self)
        """ Instantiate a Dirichlet discharge boundary.
            domain underlying domain
            stage0 stag
            wh0 momentum in the inward normal direction.
            """

        if domain is None:
            msg = 'Domain must be specified for this type of boundary'
            raise Exception(msg)

        if stage0 is None:
            raise Exception('Stage must be specified for this type of boundary')

        if wh0 is None:
            wh0 = 0.0

        self.domain = domain
        self.stage0 = stage0
        self.wh0 = wh0


    def __repr__(self):
        """ Return a representation of this instance. """
        return 'Dirichlet_Discharge_boundary(%s)' % self.domain

    def evaluate(self, vol_id, edge_id):
        """Set discharge in the (inward) normal direction"""

        normal = self.domain.get_normal(vol_id,edge_id)
        q = [self.stage0, -self.wh0*normal[0], -self.wh0*normal[1]]
        return q

        # FIXME: Consider this (taken from File_boundary) to allow
        # spatial variation
        # if vol_id is not None and edge_id is not None:
        #     i = self.boundary_indices[ vol_id, edge_id ]
        #     return self.F(t, point_id = i)
        # else:
        #     return self.F(t)


class Inflow_boundary(Boundary):
    """Apply given flow in m^3/s to boundary segment.
    Depth and momentum is derived using Manning's formula.

    Underlying domain must be specified when boundary is instantiated
    """
    
    # FIXME (Ole): This is work in progress and definitely not finished.
    # The associated test has been disabled

    def __init__(self, domain=None, rate=0.0):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for '
            msg += 'Inflow boundary'
            raise Exception(msg)

        self.domain = domain
        
        # FIXME(Ole): Allow rate to be time dependent as well
        self.rate = rate
        self.tag = None # Placeholder for tag associated with this object.

    def __repr__(self):
        return 'Inflow_boundary(%s)' %self.domain

    def evaluate(self, vol_id, edge_id):
        """Apply inflow rate at each edge of this boundary
        """
        
        # First find all segments having the same tag is vol_id, edge_id
        # This will be done the first time evaluate is called.
        if self.tag is None:
            boundary = self.domain.boundary
            self.tag = boundary[(vol_id, edge_id)]        
            
            # Find total length of boundary with this tag
            length = 0.0
            for v_id, e_id in boundary:
                if self.tag == boundary[(v_id, e_id)]:
                    length += self.domain.mesh.get_edgelength(v_id, e_id)            

            self.length = length
            self.average_momentum = self.rate/length
            
            
        # Average momentum has now been established across this boundary
        # Compute momentum in the inward normal direction 
        
        inward_normal = -self.domain.mesh.get_normal(vol_id, edge_id)       
        xmomentum, ymomentum = self.average_momentum * inward_normal
            
        # Compute depth based on Manning's formula v = 1/n h^{2/3} sqrt(S)
        # Where v is velocity, n is manning's coefficient, h is depth
        # and S is the slope into the domain. 
        # Let mu be the momentum (vh), then this equation becomes:
        #            mu = 1/n h^{5/3} sqrt(S) 
        # from which we can isolate depth to get
        #             h = (mu n/sqrt(S) )^{3/5} 
        
        slope = 0 # get gradient for this triangle dot normal
        epsilon = 1.0e-12
        import math
        
        # get manning coef from this triangle
        friction = self.domain.get_quantity('friction').get_values(\
                    location='edges', indices=[vol_id])[0]
        mannings_n = friction[edge_id]

        if slope > epsilon and mannings_n > epsilon:
            depth = pow(self.average_momentum * mannings_n/math.sqrt(slope), \
                        3.0/5) 
        else:
            depth = 1.0
            
        # Elevation on this edge    
        
        z = self.domain.get_quantity('elevation').get_values(\
                    location='edges', indices=[vol_id])[0]
        elevation = z[edge_id]
            
        # Assign conserved quantities and return
        q = num.array([elevation + depth, xmomentum, ymomentum], float)
        return q


        
    
            
        
class Field_boundary(Boundary):
    """Set boundary from given field.
    
    Given field is represented in an sww file containing
    values for stage, xmomentum and ymomentum.

    Optionally, the user can specify mean_stage to offset the stage provided
    in the sww file.

    This function is a thin wrapper around the generic File_boundary. The
    difference between the File_boundary and Field_boundary is only that the
    Field_boundary will allow you to change the level of the stage height when
    you read in the boundary condition. This is very useful when running
    different tide heights in the same area as you need only to convert one
    boundary condition to a SWW file, ideally for tide height of 0 m
    (saving disk space). Then you can use Field_boundary to read this SWW file
    and change the stage height (tide) on the fly depending on the scenario.
    """

    def __init__(self,
                 filename,
                 domain,
                 mean_stage=0.0,
                 time_thinning=1,
                 time_limit=None,
                 boundary_polygon=None,
                 default_boundary=None,
                 use_cache=False,
                 verbose=False):
        """Constructor

        :param filename: Name of sww file containing stage and x/ymomentum

        :param domain: pointer to shallow water domain for which the boundary applies

        :param mean_stage: The mean water level which will be added to stage derived from the boundary condition

        :param time_thinning: Will set how many time steps from the sww file read in will be interpolated to the boundary. 

        :param default_boundary: This will be used in case model time exceeds that available in the underlying data.

        :param time_limit: 

        :param boundary_polygon: 

        :param use_cache:        True if caching is to be used.

        :param verbose:          True if this method is to be verbose.

        For example if
        the sww file has 1 second time steps and is 24 hours
        in length it has 86400 time steps. If you set
        time_thinning to 1 it will read all these steps.
        If you set it to 100 it will read every 100th step eg
        only 864 step. This parameter is very useful to increase
        the speed of a model run that you are setting up
        and testing.

        """

        # Create generic file_boundary object
        self.file_boundary = File_boundary(filename,
                                           domain,
                                           time_thinning=time_thinning,
                                           time_limit=time_limit,
                                           boundary_polygon=boundary_polygon,
                                           default_boundary=default_boundary,
                                           use_cache=use_cache,
                                           verbose=verbose)

        # Record information from File_boundary
        self.F = self.file_boundary.F
        self.domain = self.file_boundary.domain

        # Record mean stage
        self.mean_stage = mean_stage


    def __repr__(self):
        """ Generate a string representation of this instance. """
        return 'Field boundary'


    def evaluate(self, vol_id=None, edge_id=None):
        """ Calculate 'field' boundary results.
            vol_id and edge_id are ignored

            Return linearly interpolated values based on domain.time
        """

        # Evaluate file boundary
        q = self.file_boundary.evaluate(vol_id, edge_id)

        # Adjust stage
        for j, name in enumerate(self.domain.conserved_quantities):
            if name == 'stage':
                q[j] += self.mean_stage
        return q





class Flather_external_stage_zero_velocity_boundary(Boundary):
    """Boundary condition based on a Flather type approach

    Setting the external stage with a function, and a zero external velocity,
    
    The idea is similar (but not identical) to that described on page 239 of
    the following article::

        Article{blayo05,
        Title       = {Revisiting open boundary conditions from the point of view of characteristic variables},
        Author      = {Blayo, E. and Debreu, L.},
        Journal     = {Ocean Modelling},
        Year        = {2005},
        Pages       = {231-252},
        Volume      = {9},
        }

    Approach

     #. The external (outside boundary) stage is set with a function, the
        external velocity is zero, the internal stage and velocity are taken from the
        domain values.
     #. Some 'characteristic like' variables are computed, depending on whether
        the flow is incoming or outgoing. See Blayo and Debreu (2005)
     #. The boundary conserved quantities are computed from these characteristic
        like variables

    This has been useful as a 'weakly reflecting' boundary when the stage should
    be approximately specified but allowed to adapt to outgoing waves.
    
    """

    def __init__(self, domain=None, function=None):
        """Create boundary condition object.

        :param domain: The domain on which to apply boundary condition
        :param function: Function to apply on the boundary

        Example:

        .. code:: python

            def waveform(t):
                return sea_level + normalized_amplitude/cosh(t-25)**2

            Bf = Flather_external_stage_zero_velocity_boundary(domain, waveform)

        """

        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for this type boundary'
            raise Exception(msg)

        if function is None:
            msg = 'Function must be specified for this type boundary'
            raise Exception(msg)

        self.domain = domain
        self.function = function


    def __repr__(self):
        """ Return a representation of this instance. """
        msg = 'Flather_external_stage_zero_velocity_boundary'
        msg += '(%s)' % self.domain
        return msg


    def evaluate(self, vol_id, edge_id):
        """
        """

        q = self.domain.get_conserved_quantities(vol_id, edge = edge_id)
        bed = self.domain.quantities['elevation'].centroid_values[vol_id]
        depth_inside=max(q[0]-bed,0.0)
        dt=self.domain.timestep

        normal = self.domain.get_normal(vol_id, edge_id)


        t = self.domain.get_time()

        value = self.function(t)
        try:
            stage_outside = float(value)
        except:
            stage_outside = float(value[0])

        if(depth_inside==0.):
            q[0] = stage_outside
            q[1] = 0. 
            q[2] = 0. 

        else:

            # Asssume sub-critical flow. Set the values of the characteristics as
            # appropriate, depending on whether we have inflow or outflow

            # These calculations are based on the paper cited above
            sqrt_g_on_depth_inside = (gravity/depth_inside)**0.5
            ndotq_inside = (normal[0]*q[1] + normal[1]*q[2]) # momentum perpendicular to the boundary
            if(ndotq_inside>0.):
                # Outflow (assumed subcritical)
                # Compute characteristics using a particular extrapolation
                #
                # Theory: 2 characteristics coming from inside domain, only
                # need to impose one characteristic from outside
                # 

                # w1 =  u - sqrt(g/depth)*(Stage_outside)  -- uses 'outside' info
                w1 = 0. - sqrt_g_on_depth_inside*stage_outside

                # w2 = v [velocity parallel to boundary] -- uses 'inside' info
                w2 = (+normal[1]*q[1] -normal[0]*q[2])/depth_inside

                # w3 = u + sqrt(g/depth)*(Stage_inside) -- uses 'inside info'
                w3 = ndotq_inside/depth_inside + sqrt_g_on_depth_inside*q[0]
                
            else:
                # Inflow (assumed subcritical)
                # Need to set 2 characteristics from outside information
                
                # w1 =  u - sqrt(g/depth)*(Stage_outside)  -- uses 'outside' info
                w1 = 0. - sqrt_g_on_depth_inside*stage_outside

                # w2 = v [velocity parallel to boundary] -- uses 'outside' info
                w2 = 0.

                # w3 = u + sqrt(g/depth)*(Stage_inside) -- uses 'inside info'
                w3 = ndotq_inside/depth_inside + sqrt_g_on_depth_inside*q[0]


            q[0] = (w3-w1)/(2*sqrt_g_on_depth_inside)
            qperp= (w3+w1)/2.*depth_inside
            qpar=  w2*depth_inside

            # So q[1], q[2] = qperp*(normal[0], normal[1]) + qpar*(-normal[1], normal[0])

            q[1] = qperp*normal[0] + qpar*normal[1]
            q[2] = qperp*normal[1] - qpar*normal[0]

        return q

    
    def evaluate_segment(self, domain, segment_edges): 
        """Applied in vectorized form for speed. Gareth Davies 14/07/2016
        """
 
        Stage = domain.quantities['stage']
        Elev  = domain.quantities['elevation']
        #Height= domain.quantities['height']
        Xmom  = domain.quantities['xmomentum']
        Ymom  = domain.quantities['ymomentum']

        ids = segment_edges
        vol_ids  = domain.boundary_cells[ids]
        edge_ids = domain.boundary_edges[ids]
        Normals = domain.normals

        n1  = Normals[vol_ids,2*edge_ids]
        n2  = Normals[vol_ids,2*edge_ids+1]
   
        # Get stage value 
        t = self.domain.get_time()
        value = self.function(t)
        try:
            stage_outside = float(value)
        except:
            stage_outside = float(value[0])

        # Transfer these quantities to the boundary array
        Stage.boundary_values[ids] = Stage.edge_values[vol_ids,edge_ids]
        Xmom.boundary_values[ids]  = Xmom.edge_values[vol_ids,edge_ids]
        Ymom.boundary_values[ids]  = Ymom.edge_values[vol_ids,edge_ids]
        Elev.boundary_values[ids]  = Elev.edge_values[vol_ids,edge_ids]

        bed = Elev.centroid_values[vol_ids]
        depth_inside = num.maximum(Stage.boundary_values[ids]-bed, 0.0)
        stage_outside = 0.0*Stage.boundary_values[ids] + stage_outside 

        # Do vectorized operations here
        #
        # In dry cells, the values will be ....
        q0_dry = num.where(bed <= stage_outside, stage_outside, Elev.boundary_values[ids])
        q1_dry = 0.0 * Xmom.boundary_values[ids]
        q2_dry = 0.0 * Ymom.boundary_values[ids]
        #
        # and in wet cells, the values will be ...
        # (see 'evaluate' method above for more comments on theory,
        # in particular we assume subcritical flow and a zero outside velocity)
        #
        # (note: When cells are dry, this calculation will throw invalid
        # values, but such values will never be selected to be returned)
        sqrt_g_on_depth_inside = (gravity/depth_inside)**0.5
        ndotq_inside = (n1 * Xmom.boundary_values[ids] + 
            n2 * Ymom.boundary_values[ids])
        # w1 =  u - sqrt(g/depth)*(Stage_outside)  -- uses 'outside' info
        w1 = 0.0 - sqrt_g_on_depth_inside * stage_outside
        # w2 = v [velocity parallel to boundary] -- uses 'inside' or 'outside'
        # info as required
        w2 = num.where(ndotq_inside > 0.0,
            (n2 * Xmom.boundary_values[ids] - n1 * Ymom.boundary_values[ids])/depth_inside, 
            0.0 * ndotq_inside)
        # w3 = u + sqrt(g/depth)*(Stage_inside) -- uses 'inside info'
        w3 = ndotq_inside/depth_inside + sqrt_g_on_depth_inside*Stage.boundary_values[ids]
            
        q0_wet = (w3 - w1)/(2.0 * sqrt_g_on_depth_inside)

        qperp = (w3 + w1)/2.0 * depth_inside
        qpar = w2 * depth_inside

        q1_wet = qperp * n1 + qpar * n2
        q2_wet = qperp * n2 - qpar * n1

        dry_test = num.logical_or(depth_inside == 0.0, stage_outside > bed)

        Stage.boundary_values[ids] = num.where(
            dry_test,
            q0_dry, 
            q0_wet)

        Xmom.boundary_values[ids] = num.where(
            dry_test, 
            q1_dry, 
            q1_wet)

        Ymom.boundary_values[ids] = num.where(
            dry_test,
            q2_dry,
            q2_wet)


        
