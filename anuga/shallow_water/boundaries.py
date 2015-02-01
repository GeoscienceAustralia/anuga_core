"""
Boundary conditions - specific to the shallow water wave equation

Title: ANUGA boundaries with dependencies on shallow_water_domain


Author: Ole Nielsen, Ole.Nielsen@ga.gov.au
        Stephen Roberts, Stephen.Roberts@anu.edu.au
        Duncan Gray, Duncan.Gray@ga.gov.au
        Gareth Davies, gareth.davies.ga.code@gmail.com

CreationDate: 2010

Description:
    This module contains boundary functions for ANUGA that are specific
    to the shallow water Domain class.
    
Constraints: See GPL license in the user guide
Version: 1.0 ($Revision: 7731 $)
ModifiedBy:
    $Author: hudson $
    $Date: 2010-05-18 14:54:05 +1000 (Tue, 18 May 2010) $
"""


from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Boundary, File_boundary
import numpy as num

import anuga.utilities.log as log
from anuga.fit_interpolate.interpolate import Modeltime_too_late
from anuga.fit_interpolate.interpolate import Modeltime_too_early
from anuga.config import g as gravity
     
from shallow_water_ext import rotate

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
    """Reflective boundary returns same conserved quantities as
    those present in its neighbour volume but reflected.

    This class is specific to the shallow water equation as it
    works with the momentum quantities assumed to be the second
    and third conserved quantities.
    """

    def __init__(self, domain=None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for reflective boundary'
            raise Exception, msg

        # Handy shorthands
        self.stage = domain.quantities['stage'].edge_values
        self.xmom = domain.quantities['xmomentum'].edge_values
        self.ymom = domain.quantities['ymomentum'].edge_values
        self.normals = domain.normals

        self.conserved_quantities = num.zeros(3, num.float)

    def __repr__(self):
        return 'Reflective_boundary'

    def evaluate(self, vol_id, edge_id):
        """Calculate reflections (reverse outward momentum).

        vol_id   
        edge_id  
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
        """Apply reflective BC on the boundary edges defined by
        segment_edges
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
    """Returns same momentum conserved quantities as
    those present in its neighbour volume.
    Sets stage by specifying a function f of time which may either be a
    vector function or a scalar function

    Example:

    def waveform(t):
        return sea_level + normalized_amplitude/cosh(t-25)**2

    Bts = Transmissive_momentum_set_stage_boundary(domain, waveform)

    Underlying domain must be specified when boundary is instantiated
    """

    def __init__(self, domain=None, function=None):
        Boundary.__init__(self)
        """ Instantiate a Transmissive_momentum_set_stage_boundary. """


        if domain is None:
            msg = 'Domain must be specified for this type boundary'
            raise Exception, msg

        if function is None:
            msg = 'Function must be specified for this type boundary'
            raise Exception, msg

        self.domain = domain

        if isinstance(function, (int, float)):
            tmp = function
            function = lambda t: tmp
            
        self.function = function


    def __repr__(self):
        """ Return a representation of this instance. """
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
    """Returns the same normal momentum as that 
    present in neighbour volume edge. Zero out the tangential momentum. 
    Sets stage by specifying a function f of time which may either be a
    vector function or a scalar function

    Example:

    def waveform(t):
        return sea_level + normalized_amplitude/cosh(t-25)**2

    Bts = Transmissive_n_momentum_zero_t_momentum_set_stage_boundary\
                            (domain, waveform)

    Underlying domain must be specified when boundary is instantiated
    """

    def __init__(self, domain=None, function=None):
        """ Instantiate a
            Transmissive_n_momentum_zero_t_momentum_set_stage_boundary.
            domain is the domain containing the boundary
            function is the function to apply
        """

        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for this type boundary'
            raise Exception, msg

        if function is None:
            msg = 'Function must be specified for this type boundary'
            raise Exception, msg

        self.domain = domain
        self.function = function


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
        ndotq = (normal[0]*q[1] + normal[1]*q[2])
        q[1] = normal[0]*ndotq
        q[2] = normal[1]*ndotq

        return q

class Transmissive_stage_zero_momentum_boundary(Boundary):
    """Return same stage as those present in its neighbour volume.
    Set momentum to zero.

    Underlying domain must be specified when boundary is instantiated
    """

    def __init__(self, domain=None):
        """ Instantiate a Transmissive (zero momentum) boundary. """

        Boundary.__init__(self)

        if domain is None:
            msg = ('Domain must be specified for '
                   'Transmissive_stage_zero_momentum boundary')
            raise Exception, msg

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
        except Exception, e:
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

        #-------------------------------------------------
        # Now update boundary values
        #-------------------------------------------------
        domain.quantities['stage'].boundary_values[ids] = q_bdry
        domain.quantities['xmomentum'].boundary_values[ids] = 0.0
        domain.quantities['ymomentum'].boundary_values[ids] = 0.0        


    def get_boundary_values(self, t=None):

        if t is None:
            t = self.get_time()
            
        try:
            res = self.f(t)
        except Modeltime_too_early, e:
            raise Modeltime_too_early(e)
        except Modeltime_too_late, e:
            if self.default_boundary is None:
                raise Exception(e) # Reraise exception
            else:
                # Pass control to default boundary
                res = self.default_boundary
                
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
                        msg += 'Note: Further warnings will be suppressed'
                        log.critical(msg)
               
                    # FIXME (Ole): Replace this crude flag with
                    # Python's ability to print warnings only once.
                    # See http://docs.python.org/lib/warning-filter.html
                    self.default_boundary_invoked = True

        return res


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

        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for this type boundary'
            raise Exception, msg

        if function is None:
            msg = 'Function must be specified for this type boundary'
            raise Exception, msg

        self.domain = domain
        self.function = function
        self.default_stage = default_stage

        self.Elev  = domain.quantities['elevation']
        self.Stage = domain.quantities['stage']
        self.Height = domain.quantities['height']

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
            stage = float(value)
        except:
            stage = float(value[0])



        q = self.conserved_quantities
        #q[0] = self.stage[vol_id, edge_id]
        q[0] = stage
        q[1] = self.xmom[vol_id, edge_id]
        q[2] = self.ymom[vol_id, edge_id]

        normal = self.normals[vol_id, 2*edge_id:2*edge_id+2]

        r = rotate(q, normal, direction = 1)
        r[1] = -r[1]
        q = rotate(r, normal, direction = -1)


        return q






    def evaluate_segment(self, domain, segment_edges):
        """Apply reflective BC on the boundary edges defined by
        segment_edges
        """

        if segment_edges is None:
            return
        if domain is None:
            return

        t = self.domain.get_time()


        value = self.function(t)
        try:
            stage = float(value)
        except:
            stage = float(value[0])
                       

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
            raise Exception, msg

        if stage0 is None:
            raise Exception, 'Stage must be specified for this type of boundary'

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
            raise Exception, msg

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
        q = num.array([elevation + depth, xmomentum, ymomentum], num.float)
        return q


        
    
            
        
class Field_boundary(Boundary):
    """Set boundary from given field represented in an sww file containing
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

        filename: Name of sww file containing stage and x/ymomentum
        domain: pointer to shallow water domain for which the boundary applies
        mean_stage: The mean water level which will be added to stage derived
                    from the boundary condition
        time_thinning: Will set how many time steps from the sww file read in
                       will be interpolated to the boundary. For example if
                       the sww file has 1 second time steps and is 24 hours
                       in length it has 86400 time steps. If you set
                       time_thinning to 1 it will read all these steps.
                       If you set it to 100 it will read every 100th step eg
                       only 864 step. This parameter is very useful to increase
                       the speed of a model run that you are setting up
                       and testing.

        default_boundary: Must be either None or an instance of a
                          class descending from class Boundary.
                          This will be used in case model time exceeds
                          that available in the underlying data.

                          Note that mean_stage will also be added to this.

        time_limit: 
        boundary_polygon: 
        use_cache:        True if caching is to be used.
        verbose:          True if this method is to be verbose.

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
    """

    Boundary condition based on a Flather type approach 
    (setting the external stage with a function, and a zero external velocity),
    

    The idea is similar (but not identical) to that described on page 239 of
    the following article:

    Article{blayo05,
      Title                    = {Revisiting open boundary conditions from the point of view of characteristic variables},
      Author                   = {Blayo, E. and Debreu, L.},
      Journal                  = {Ocean Modelling},
      Year                     = {2005},
      Pages                    = {231-252},
      Volume                   = {9},
    }

    Approach
    1) The external (outside boundary) stage is set with a function, the
       external velocity is zero, the internal stage and velocity are taken from the
       domain values.
    2) Some 'characteristic like' variables are computed, depending on whether
       the flow is incoming or outgoing. See Blayo and Debreu (2005)
    3) The boundary conserved quantities are computed from these characteristic
       like variables

    This has been useful as a 'weakly reflecting' boundary when the stage should
    be approximately specified but allowed to adapt to outgoing waves.


    Example:

    def waveform(t):
        return sea_level + normalized_amplitude/cosh(t-25)**2

    Bf = Flather_external_stage_zero_velocity_boundary(domain, waveform)

    Underlying domain must be specified when boundary is instantiated


    
    """

    def __init__(self, domain=None, function=None):
        """ Instantiate a
            Nudge_boundary.
            domain is the domain containing the boundary
            function is the function to apply
        """

        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for this type boundary'
            raise Exception, msg

        if function is None:
            msg = 'Function must be specified for this type boundary'
            raise Exception, msg

        self.domain = domain
        self.function = function


    def __repr__(self):
        """ Return a representation of this instance. """
        msg = 'Nudge_boundary'
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
            q[2] = qperp*normal[1] -qpar*normal[0]

        return q
