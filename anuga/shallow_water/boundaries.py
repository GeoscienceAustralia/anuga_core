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
import numpy as np

import anuga.utilities.log as log
from anuga.fit_interpolate.interpolate import Modeltime_too_late
from anuga.fit_interpolate.interpolate import Modeltime_too_early
from anuga.config import g as gravity

from anuga.shallow_water.sw_domain_openmp_ext import rotate, evaluate_reflective_segment

try:
    from numba import jit
except ImportError:
    def jit(nopython=True):
        """Dummy decorator for numba"""
        def dummy_decorator(func):
            return func
        return dummy_decorator


x = np.arange(100).reshape(10, 10)


class Reflective_boundary(Boundary):
    """Reflective boundary condition.

    Returns the same conserved quantities as the neighbour volume edge but
    with the normal momentum component negated, so the net mass flux through
    the boundary is zero (wall / no-slip analogue for the shallow-water
    equations).

    Parameters
    ----------
    domain : anuga.Domain
        The domain to which this boundary is attached.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> Br = anuga.Reflective_boundary(domain)
    >>> domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    """

    def __init__(self, domain=None):
        """Initialise a reflective boundary.

        Parameters
        ----------
        domain : anuga.Domain
            The domain to which this boundary is attached.

        Raises
        ------
        Exception
            If *domain* is ``None``.
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

        self.conserved_quantities = np.zeros(3, float)


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


    @jit(nopython=True)
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



    # TODO JLGV, reflective boundary condition needs openmp version
    # this one first
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
        ids_array = np.array(ids, dtype=np.int64)

        evaluate_reflective_segment(domain, ids_array, vol_ids, edge_ids)





class Transmissive_momentum_set_stage_boundary(Boundary):
    """Transmissive momentum boundary with prescribed stage.

    Returns the same momentum conserved quantities as the neighbour volume
    edge (transmissive / open condition) while overriding stage with a
    caller-supplied scalar function of time.

    Parameters
    ----------
    domain : anuga.Domain
        The domain to which this boundary is attached.
    function : callable
        A function ``f(t)`` returning the stage value at time *t*.  May
        also be a scalar ``int`` or ``float``, in which case the stage is
        held constant.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> BC = anuga.Transmissive_momentum_set_stage_boundary(domain, lambda t: 0.5)
    >>> domain.set_boundary({'left': BC, 'right': BC, 'top': BC, 'bottom': BC})
    """

    def __init__(self, domain=None, function=None):
        """Initialise a transmissive-momentum / set-stage boundary.

        Parameters
        ----------
        domain : anuga.Domain
            The domain to which this boundary is attached.
        function : callable or float
            Stage function ``f(t)`` or a constant value.

        Raises
        ------
        Exception
            If *domain* or *function* is ``None``.
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
                log.info(msg)
                t -= self.function.time[-1]

        value = self.function(t)
        try:
            x = float(value)
        except (ValueError, TypeError):
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
    """Transmissive normal momentum, zero tangential momentum, prescribed stage.

    Returns the same momentum component normal to the boundary as the
    neighbour volume edge.  The tangential momentum component is zeroed.
    Stage is set by a caller-supplied function of time.

    Parameters
    ----------
    domain : anuga.Domain
        The domain to which this boundary is attached.
    function : callable
        A function ``f(t)`` returning the stage at time *t*.
    default_boundary : float, optional
        Stage value returned when model time exceeds the range of
        *function*.  Default is ``0.0``.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> BC = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
    ...     domain, lambda t: 0.5)
    >>> domain.set_boundary({'left': BC, 'right': BC, 'top': BC, 'bottom': BC})
    """

    def __init__(self, domain=None, function=None, default_boundary=0.0):
        """Initialise the boundary.

        Parameters
        ----------
        domain : anuga.Domain
            The domain to which this boundary is attached.
        function : callable
            Stage function ``f(t)``.
        default_boundary : float, optional
            Fallback stage when model time is out of range.  Default ``0.0``.

        Raises
        ------
        Exception
            If *domain* or *function* is ``None``.
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


        value = self.get_boundary_values()
        try:
            x = float(value)
        except (ValueError, TypeError):
            x = float(value[0])

        q[0] = x

        ndotq = (normal[0]*q[1] + normal[1]*q[2])
        q[1] = normal[0]*ndotq
        q[2] = normal[1]*ndotq

        return q

    # TODO JLGV, needs openmp version
    def evaluate_segment(self, domain, segment_edges):
        """Apply BC on the boundary edges defined by segment_edges

        :param domain: Apply BC on this domain
        :param segment_edges: List of boundary cells on which to apply BC

        Vectorized form for speed. Gareth Davies 14/07/2016

        """

        Stage = domain.quantities['stage']
        Elev  = domain.quantities['elevation']
        Height= domain.quantities['height']
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
        except (ValueError, TypeError):
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
    """Transmissive stage with zero momentum boundary.

    Copies the stage from the neighbour volume edge while forcing both
    momentum components to zero.  Useful as a simple open boundary where
    water level is inherited from the interior but no momentum is imposed.

    Parameters
    ----------
    domain : anuga.Domain
        The domain to which this boundary is attached.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> Bt = anuga.Transmissive_stage_zero_momentum_boundary(domain)
    >>> domain.set_boundary({'left': Bt, 'right': Bt, 'top': Bt, 'bottom': Bt})
    """

    def __init__(self, domain=None):
        """Initialise a transmissive-stage / zero-momentum boundary.

        Parameters
        ----------
        domain : anuga.Domain
            The domain to which this boundary is attached.

        Raises
        ------
        Exception
            If *domain* is ``None``.
        """

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
    """Time-varying stage boundary with zero momentum.

    Sets stage as a scalar function of time while holding both momentum
    components at zero.

    Parameters
    ----------
    domain : anuga.Domain
        The domain to which this boundary is attached.
    function : callable
        A function ``f(t)`` returning the stage value at model time *t*.
        Must be convertible to a scalar float.
    default_boundary : float or None, optional
        Stage value returned when model time is out of the function's
        range.  ``None`` means raise an exception on out-of-range.
    verbose : bool, optional
        If ``True``, emit informational messages.  Default ``False``.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> # 2 m square wave from t=60 s to t=3660 s, zero elsewhere
    >>> Bt = anuga.Time_stage_zero_momentum_boundary(
    ...     domain, function=lambda t: (60 < t < 3660) * 2.0)
    >>> domain.set_boundary({'left': Bt, 'right': Bt, 'top': Bt, 'bottom': Bt})
    """

    def __init__(self, domain=None,
                 #f=None, # Should be removed and replaced by function below
                 function=None,
                 default_boundary=None,
                 verbose=False):
        """Initialise a time-stage / zero-momentum boundary.

        Parameters
        ----------
        domain : anuga.Domain
            The domain to which this boundary is attached.
        function : callable
            Stage function ``f(t)``.
        default_boundary : float or None, optional
            Fallback stage when model time is out of range.
        verbose : bool, optional
            Verbosity flag.  Default ``False``.

        Raises
        ------
        Exception
            If *domain* or *function* is ``None``, or if *function* cannot
            be evaluated at ``t=0`` or does not return a scalar float.
        """
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
        except (ValueError, TypeError):
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
    """Stage boundary using characteristic (Riemann-invariant) extrapolation.

    Sets the exterior stage via a function of time.  Momentum at the boundary
    is determined from a characteristic decomposition, giving a weakly
    reflecting open boundary: outgoing waves leave cleanly while incoming
    waves are set by the prescribed stage.

    Parameters
    ----------
    domain : anuga.Domain
        The domain to which this boundary is attached.
    function : callable
        A function ``f(t)`` returning the exterior stage at model time *t*.
    default_stage : float, optional
        Ambient stage assumed before the wave arrives.  Default ``0.0``.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> Bcs = anuga.Characteristic_stage_boundary(domain, lambda t: 0.1)
    >>> domain.set_boundary({'left': Bcs, 'right': Bcs, 'top': Bcs, 'bottom': Bcs})
    """

    def __init__(self, domain=None, function=None, default_stage=0.0):
        """Initialise a characteristic-stage boundary.

        Parameters
        ----------
        domain : anuga.Domain
            The domain to which this boundary is attached.
        function : callable
            Exterior stage function ``f(t)``.
        default_stage : float, optional
            Ambient stage assumed before the wave.  Default ``0.0``.

        Raises
        ------
        Exception
            If *domain* or *function* is ``None``.
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
        except (ValueError, TypeError):
            w_outside = float(value[0])

        q = np.zeros(len(self.conserved_quantities), float)

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
        except (ValueError, TypeError):
            w_outside = float(value[0])

        # Transfer these quantities to the boundary array
        Stage.boundary_values[ids] = Stage.edge_values[vol_ids,edge_ids]
        Xmom.boundary_values[ids]  = Xmom.edge_values[vol_ids,edge_ids]
        Ymom.boundary_values[ids]  = Ymom.edge_values[vol_ids,edge_ids]
        Elev.boundary_values[ids]  = Elev.edge_values[vol_ids,edge_ids]



        h_inside = np.maximum(Stage.boundary_values[ids]-Elev.boundary_values[ids], 0.0)
        w_outside = 0.0*Stage.boundary_values[ids] + w_outside

        # Do vectorized operations here
        #
        # In dry cells, the values will be ....
        q0_dry = np.where(Elev.boundary_values[ids] <= w_outside, w_outside, Elev.boundary_values[ids])
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
        h_inside  = np.maximum(Stage.boundary_values[ids] - Elev.boundary_values[ids], 0)
        uh_inside = n1 * Xmom.boundary_values[ids] + n2 * Ymom.boundary_values[ids]
        vh_inside = n2 * Xmom.boundary_values[ids] - n1 * Ymom.boundary_values[ids]
        u_inside  = np.divide(uh_inside, h_inside,
                              out=np.zeros_like(uh_inside),
                              where=h_inside > 0.0)

        h_outside = np.maximum(w_outside - Elev.boundary_values[ids], 0)

        sqrt_h_inside = h_inside**0.5
        sqrt_h_outside = h_outside**0.5

        h_m = (0.5*(sqrt_h_inside+sqrt_h_outside) + u_inside/4.0/sqrt_g )**2
        u_m = 0.5*u_inside + sqrt_g*(sqrt_h_inside - sqrt_h_outside)

        uh_m = h_m*u_m

        # if uh_inside > 0.0 then outflow
        vh_m = np.where(uh_inside > 0.0, vh_inside, 0.0)

        w_m = h_m + Elev.boundary_values[ids]

        dry_test = np.logical_or(h_inside == 0.0, h_outside == 0.0)

        q1 = uh_m*n1 + vh_m*n2
        q2 = uh_m*n2 - vh_m*n1

        Stage.boundary_values[ids] = np.where(
            dry_test,
            w_outside,
            w_m
            )

        Xmom.boundary_values[ids] = np.where(
            dry_test,
            0.0,
            q1
            )

        Ymom.boundary_values[ids] = np.where(
            dry_test,
            0.0,
            q2)

class Dirichlet_discharge_boundary(Boundary):
    """Dirichlet boundary with prescribed stage and inward-normal discharge.

    Sets stage to a constant *stage0* and momentum in the inward-normal
    direction to *wh0*.  The tangential momentum is always zero.

    Parameters
    ----------
    domain : anuga.Domain
        The domain to which this boundary is attached.
    stage0 : float
        Prescribed water stage (m).
    wh0 : float, optional
        Momentum magnitude in the inward-normal direction (m^2/s).
        Default is ``0.0`` (stage only, no imposed flow).

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> Bd = anuga.Dirichlet_discharge_boundary(domain, stage0=1.0, wh0=0.5)
    >>> domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})
    """

    def __init__(self, domain=None, stage0=None, wh0=None):
        """Initialise a Dirichlet discharge boundary.

        Parameters
        ----------
        domain : anuga.Domain
            The domain to which this boundary is attached.
        stage0 : float
            Prescribed stage.
        wh0 : float, optional
            Inward-normal momentum.  Default ``0.0``.

        Raises
        ------
        Exception
            If *domain* or *stage0* is ``None``.
        """
        Boundary.__init__(self)

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
    """Inflow boundary that imposes a volumetric flow rate.

    Distributes the prescribed flow (m^3/s) uniformly along the
    boundary segment.  Depth and momentum are derived from Manning's
    formula using the local bed gradient and friction coefficient.

    .. note::
        This class is work in progress and the associated test is disabled.

    Parameters
    ----------
    domain : anuga.Domain
        The domain to which this boundary is attached.
    rate : float, optional
        Total volumetric inflow rate in m^3/s.  Default ``0.0``.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> Bi = anuga.Inflow_boundary(domain, rate=1.0)
    >>> domain.set_boundary({'left': Bi, 'right': Bi, 'top': Bi, 'bottom': Bi})
    """

    # FIXME (Ole): This is work in progress and definitely not finished.
    # The associated test has been disabled

    def __init__(self, domain=None, rate=0.0):
        """Initialise an inflow boundary.

        Parameters
        ----------
        domain : anuga.Domain
            The domain to which this boundary is attached.
        rate : float, optional
            Volumetric inflow rate (m^3/s).  Default ``0.0``.

        Raises
        ------
        Exception
            If *domain* is ``None``.
        """
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
        q = np.array([elevation + depth, xmomentum, ymomentum], float)
        return q






class Field_boundary(Boundary):
    """Boundary condition driven by an SWW field file with optional stage offset.

    Reads stage, x-momentum and y-momentum time series from an SWW file and
    applies them as a boundary condition, linearly interpolating in time.
    An optional *mean_stage* offset can be added to the stage values, which
    avoids regenerating the SWW file when running at different tide levels.

    This is a thin wrapper around :class:`File_boundary`; the only difference
    is the *mean_stage* offset capability.

    Parameters
    ----------
    filename : str
        Path to the SWW file containing stage and momentum time series.
    domain : anuga.Domain
        The domain to which this boundary is attached.
    mean_stage : float, optional
        Constant offset added to the stage read from the file.  Useful for
        running at different tidal datums without recreating the SWW file.
        Default ``0.0``.
    time_thinning : int, optional
        Read every *time_thinning*-th time step from the file.  Larger values
        speed up model setup at the cost of temporal resolution.  Default
        ``1`` (all steps).
    time_limit : float or None, optional
        Stop reading the file after this time (seconds).  ``None`` means
        read to the end.
    boundary_polygon : list or None, optional
        Clip the SWW points to this polygon.  ``None`` means use all points.
    default_boundary : float or None, optional
        Stage returned when model time exceeds the file's time range.
        ``None`` raises an exception on out-of-range.
    use_cache : bool, optional
        Cache the interpolated field function.  Default ``False``.
    verbose : bool, optional
        Emit progress messages.  Default ``False``.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> Bf = anuga.Field_boundary('boundary.sww', domain, mean_stage=0.5)
    >>> domain.set_boundary({'left': Bf, 'right': Bf, 'top': Bf, 'bottom': Bf})
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
        """Initialise a field boundary.

        Parameters
        ----------
        filename : str
            Path to the SWW file.
        domain : anuga.Domain
            The domain to which this boundary is attached.
        mean_stage : float, optional
            Stage offset (m).  Default ``0.0``.
        time_thinning : int, optional
            Step stride when reading time steps.  Default ``1``.
        time_limit : float or None, optional
            Maximum time to read from file.
        boundary_polygon : list or None, optional
            Polygon used to clip SWW points.
        default_boundary : float or None, optional
            Fallback stage when model time is out of range.
        use_cache : bool, optional
            Enable caching.  Default ``False``.
        verbose : bool, optional
            Verbosity flag.  Default ``False``.
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
    """Weakly-reflecting open boundary using a Flather-type characteristic approach.

    Sets the exterior stage via a function of time and assumes zero exterior
    velocity.  Interior values are taken from the domain.  The boundary
    conserved quantities are then computed from characteristic-like variables,
    making this boundary weakly reflecting — outgoing waves leave with minimal
    spurious reflection while incoming wave forcing is prescribed.

    The approach is similar (but not identical) to that described on page 239
    of:

    .. code-block:: bibtex

        @Article{blayo05,
          title   = {Revisiting open boundary conditions from the point of
                     view of characteristic variables},
          author  = {Blayo, E. and Debreu, L.},
          journal = {Ocean Modelling},
          year    = {2005},
          volume  = {9},
          pages   = {231--252},
        }

    Algorithm
    ---------
    1. The exterior stage is set from *function(t)*; exterior velocity is zero;
       interior stage and velocity are taken from the domain edge values.
    2. Characteristic-like variables are computed depending on whether flow is
       incoming or outgoing (see Blayo & Debreu 2005).
    3. The boundary conserved quantities (stage, x-momentum, y-momentum) are
       recovered from these characteristic variables.

    Parameters
    ----------
    domain : anuga.Domain
        The domain to which this boundary is attached.
    function : callable
        A function ``f(t)`` returning the exterior stage at model time *t*.
        Typically a :func:`~anuga.file_function` time series.
    default_boundary : float, optional
        Stage value returned when model time exceeds the range of *function*
        (e.g. when a file-function time series ends).  ``0.0`` corresponds to
        ambient sea level / no wave forcing.  Default ``0.0``.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> Bf = anuga.Flather_external_stage_zero_velocity_boundary(
    ...     domain, lambda t: 0.1, default_boundary=0.0)
    >>> domain.set_boundary({'left': Bf, 'right': Bf, 'top': Bf, 'bottom': Bf})
    """

    def __init__(self, domain=None, function=None, default_boundary=0.0):
        """Initialise a Flather-type open boundary.

        Parameters
        ----------
        domain : anuga.Domain
            The domain to which this boundary is attached.
        function : callable
            Exterior stage function ``f(t)``.
        default_boundary : float, optional
            Fallback stage when model time is out of range.  Default ``0.0``.

        Raises
        ------
        Exception
            If *domain* or *function* is ``None``.
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

        value = self.get_boundary_values(t)
        try:
            stage_outside = float(value)
        except (ValueError, TypeError):
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
        value = self.get_boundary_values(t)
        try:
            stage_outside = float(value)
        except (ValueError, TypeError):
            stage_outside = float(value[0])

        # Transfer these quantities to the boundary array
        Stage.boundary_values[ids] = Stage.edge_values[vol_ids,edge_ids]
        Xmom.boundary_values[ids]  = Xmom.edge_values[vol_ids,edge_ids]
        Ymom.boundary_values[ids]  = Ymom.edge_values[vol_ids,edge_ids]
        Elev.boundary_values[ids]  = Elev.edge_values[vol_ids,edge_ids]

        bed = Elev.centroid_values[vol_ids]
        depth_inside = np.maximum(Stage.boundary_values[ids]-bed, 0.0)
        stage_outside = 0.0*Stage.boundary_values[ids] + stage_outside

        # Do vectorized operations here
        #
        # In dry cells, the values will be ....
        q0_dry = np.where(bed <= stage_outside, stage_outside, Elev.boundary_values[ids])
        q1_dry = 0.0 * Xmom.boundary_values[ids]
        q2_dry = 0.0 * Ymom.boundary_values[ids]
        #
        # and in wet cells, the values will be ...
        # (see 'evaluate' method above for more comments on theory,
        # in particular we assume subcritical flow and a zero outside velocity)
        #
        # (note: When cells are dry, this calculation will throw invalid
        # values, but such values will never be selected to be returned)
        with np.errstate(invalid='ignore', divide='ignore'):
            sqrt_g_on_depth_inside = (gravity/depth_inside)**0.5
            ndotq_inside = (n1 * Xmom.boundary_values[ids] +
                n2 * Ymom.boundary_values[ids])
            # w1 =  u - sqrt(g/depth)*(Stage_outside)  -- uses 'outside' info
            w1 = 0.0 - sqrt_g_on_depth_inside * stage_outside
            # w2 = v [velocity parallel to boundary] -- uses 'inside' or 'outside'
            # info as required
            w2 = np.where(ndotq_inside > 0.0,
                (n2 * Xmom.boundary_values[ids] - n1 * Ymom.boundary_values[ids])/depth_inside,
                0.0 * ndotq_inside)
            # w3 = u + sqrt(g/depth)*(Stage_inside) -- uses 'inside info'
            w3 = ndotq_inside/depth_inside + sqrt_g_on_depth_inside*Stage.boundary_values[ids]

        q0_wet = (w3 - w1)/(2.0 * sqrt_g_on_depth_inside)

        qperp = (w3 + w1)/2.0 * depth_inside
        qpar = w2 * depth_inside

        q1_wet = qperp * n1 + qpar * n2
        q2_wet = qperp * n2 - qpar * n1

        dry_test = np.logical_or(depth_inside == 0.0, stage_outside > bed)

        Stage.boundary_values[ids] = np.where(
            dry_test,
            q0_dry,
            q0_wet)

        Xmom.boundary_values[ids] = np.where(
            dry_test,
            q1_dry,
            q1_wet)

        Ymom.boundary_values[ids] = np.where(
            dry_test,
            q2_dry,
            q2_wet)


class Absorbing_wave_boundary(Boundary):
    """Active-absorption open boundary with prescribed incoming wave.

    Simultaneously prescribes an incoming wave and absorbs outgoing
    (reflected) waves by setting the ghost-cell stage to::

        stage_ghost = 2 * wave(t) - stage_interior

    so that the boundary face always sees exactly ``wave(t)`` regardless
    of what is propagating back from inside the domain.  The ghost normal
    momentum is set to ``ghost_depth × interior_velocity`` (velocity-
    preserving transmissive): this keeps the ghost velocity bounded when
    ghost depth is small (e.g. after bed-clamping) and gives the correct
    incoming-wave energy flux.  Zero ghost momentum would halve the
    incoming wave amplitude.  Tangential momentum is zeroed.  When the
    ghost cell is dry (bed-clamped to zero depth), both momentum components
    are set to zero to avoid a pathological dry-ghost state in the Riemann
    solver.

    This is the numerical equivalent of an *active-absorption wavemaker* used
    in physical wave-tank experiments: the paddle continuously adjusts to
    cancel any reflected wave arriving at the boundary while still generating
    the desired incoming signal.

    .. note::
        Absorption efficiency depends on the phase of the standing-wave cycle.
        At a standing-wave antinode (``u ≈ 0`` at the boundary) the net
        momentum flux is zero and energy drains slowly over multiple wave
        periods.  For highly reflective short-domain problems (e.g. the
        Okushiri benchmark) this is still preferable to a pure Flather
        condition, which can exhibit 2× stage amplification at the boundary
        gauge.

    Parameters
    ----------
    domain : anuga.Domain
        The domain to which this boundary is attached.
    function : callable
        A function ``f(t)`` returning the prescribed stage at model time *t*.
        Typically a :func:`~anuga.file_function` time series.
    default_boundary : float, optional
        Stage returned when model time exceeds the range of *function*.
        Defaults to ``0.0``.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> Ba = anuga.Absorbing_wave_boundary(domain, lambda t: 0.1 * (t < 5))
    >>> domain.set_boundary({'left': Ba, 'right': Ba, 'top': Ba, 'bottom': Ba})
    """

    def __init__(self, domain=None, function=None, default_boundary=0.0):
        """Initialise an active-absorption open boundary.

        Parameters
        ----------
        domain : anuga.Domain
            The domain to which this boundary is attached.
        function : callable
            Prescribed stage function ``f(t)``.
        default_boundary : float, optional
            Fallback stage when model time is out of range.  Default ``0.0``.

        Raises
        ------
        Exception
            If *domain* or *function* is ``None``.
        """

        Boundary.__init__(self)

        if domain is None:
            raise Exception('Domain must be specified for this type boundary')

        if function is None:
            raise Exception('Function must be specified for this type boundary')

        self.domain = domain
        self.function = function
        self.default_boundary = default_boundary

    def __repr__(self):
        return 'Absorbing_wave_boundary(%s)' % self.domain

    def evaluate(self, vol_id, edge_id):
        """Return ghost-cell state for active-absorption boundary."""

        q = self.domain.get_conserved_quantities(vol_id, edge=edge_id)
        bed = self.domain.quantities['elevation'].centroid_values[vol_id]
        depth_inside = max(q[0] - bed, 0.0)

        t = self.domain.get_time()
        value = self.get_boundary_values(t)
        try:
            wave = float(value)
        except (ValueError, TypeError):
            wave = float(value[0])

        normal = self.domain.get_normal(vol_id, edge_id)

        if depth_inside == 0.0:
            # Dry interior: prescribe wave stage directly
            q[0] = wave
            q[1] = 0.0
            q[2] = 0.0
        else:
            # Active absorption: ghost stage ensures face = wave(t).
            # Clamp to bed so the ghost depth is never negative — this
            # can occur when a large reflected wave makes stage_interior
            # exceed 2*wave(t).
            ghost_stage = max(2.0 * wave - q[0], bed)
            ghost_depth = ghost_stage - bed

            # Interior velocity (bounded by depth_inside > 0 guard above)
            u_int = q[1] / depth_inside
            v_int = q[2] / depth_inside

            q[0] = ghost_stage

            if ghost_depth > 0.0:
                # Ghost momentum = ghost_depth × interior_velocity (normal
                # component only, tangential zeroed).  Preserving velocity
                # (not raw momentum) keeps the ghost velocity bounded when
                # ghost_depth << interior_depth; raw transmissive momentum
                # (ghost_hu = interior_hu) causes ghost velocity → ∞ as
                # ghost_depth → 0.  For deep boundaries (e.g. Okushiri)
                # ghost_depth ≈ interior_depth and the result is identical.
                ndotu = normal[0] * u_int + normal[1] * v_int
                q[1] = ghost_depth * ndotu * normal[0]
                q[2] = ghost_depth * ndotu * normal[1]
            else:
                q[1] = 0.0
                q[2] = 0.0

        return q

    def evaluate_segment(self, domain, segment_edges):
        """Vectorised form for speed."""

        Stage = domain.quantities['stage']
        Elev  = domain.quantities['elevation']
        Xmom  = domain.quantities['xmomentum']
        Ymom  = domain.quantities['ymomentum']

        ids      = segment_edges
        vol_ids  = domain.boundary_cells[ids]
        edge_ids = domain.boundary_edges[ids]
        Normals  = domain.normals

        n1 = Normals[vol_ids, 2 * edge_ids]
        n2 = Normals[vol_ids, 2 * edge_ids + 1]

        t = self.domain.get_time()
        value = self.get_boundary_values(t)
        try:
            wave_val = float(value)
        except (ValueError, TypeError):
            wave_val = float(value[0])

        stage_interior = Stage.edge_values[vol_ids, edge_ids]
        bed            = Elev.centroid_values[vol_ids]
        depth_inside   = np.maximum(stage_interior - bed, 0.0)

        # Active absorption ghost stage; clamp to bed to prevent negative
        # ghost depth when a large reflected wave makes stage_interior > 2*wave
        raw_ghost = np.where(depth_inside > 0.0, 2.0 * wave_val - stage_interior, wave_val)
        stage_ghost = np.maximum(raw_ghost, bed)
        ghost_depth = stage_ghost - bed
        Stage.boundary_values[ids] = stage_ghost

        # Ghost momentum = ghost_depth × interior_velocity (normal component
        # only, tangential zeroed).  Preserving velocity (not raw momentum)
        # keeps ghost velocity bounded when ghost_depth → 0 (bed-clamped);
        # raw transmissive momentum (ghost_hu = interior_hu) would cause
        # ghost velocity → ∞.  Zero ghost velocity when either side is dry.
        both_wet   = (depth_inside > 0.0) & (ghost_depth > 0.0)
        safe_depth = np.where(depth_inside > 0.0, depth_inside, 1.0)
        q1    = Xmom.edge_values[vol_ids, edge_ids]
        q2    = Ymom.edge_values[vol_ids, edge_ids]
        u_int = q1 / safe_depth
        v_int = q2 / safe_depth
        ndotu = n1 * u_int + n2 * v_int
        Xmom.boundary_values[ids] = np.where(both_wet, ghost_depth * ndotu * n1, 0.0)
        Ymom.boundary_values[ids] = np.where(both_wet, ghost_depth * ndotu * n2, 0.0)


class Characteristic_wave_boundary(Boundary):
    """Nonlinear characteristic open boundary with prescribed incoming wave.

    Prescribes the incoming Riemann invariant from a wave perturbation and
    extrapolates the outgoing Riemann invariant from the interior, solving
    the characteristic equations exactly (no linearisation).

    The ghost state is derived from::

        c_ghost     = (v_n_int + 2*c_int - 2*c0 + 4*c_wave) / 4
        v_n_ghost   = c0 - 2*c_wave + v_n_int/2 + c_int
        h_ghost     = c_ghost**2 / g
        stage_ghost = h_ghost + bed

    where *c0* = sqrt(g*h0) is the background wave speed (computed from
    ``background_stage`` at each cell), *c_wave* = sqrt(g*h_wave) is the
    prescribed wave speed, and *v_n* is the outward-normal velocity.

    The wave function returns a **perturbation** from ``background_stage``
    (not an absolute stage).  This is important: the background depth h0
    = ``background_stage`` − bed is used to set the incoming Riemann
    invariant assuming no outgoing wave at the exterior.

    Compared to :class:`Absorbing_wave_boundary`:

    * Better for large-amplitude waves (η ~ h) where linearisation error
      in the Flather form is significant.
    * The face stage is not guaranteed to equal wave(t) exactly; the
      characteristic condition controls energy transport rather than
      stage directly.
    * For small-amplitude waves (η ≪ h) both classes give similar results.

    Parameters
    ----------
    domain : anuga.Domain
        The domain to which this boundary is attached.
    function : callable
        ``f(t)`` returning the stage *perturbation* at time *t*.
    background_stage : float, optional
        Still-water stage around which the perturbation is measured.
        Defaults to ``0.0`` (sea level).
    default_boundary : float, optional
        Perturbation returned when model time is out of range.  Default ``0.0``.
    """

    def __init__(self, domain=None, function=None,
                 background_stage=0.0, default_boundary=0.0):
        Boundary.__init__(self)
        if domain is None:
            raise Exception('Domain must be specified for this type boundary')
        if function is None:
            raise Exception('Function must be specified for this type boundary')
        self.domain = domain
        self.function = function
        self.background_stage = float(background_stage)
        self.default_boundary = default_boundary

    def __repr__(self):
        return 'Characteristic_wave_boundary(%s)' % self.domain

    def evaluate(self, vol_id, edge_id):
        """Return ghost-cell state from nonlinear characteristic decomposition."""

        q = self.domain.get_conserved_quantities(vol_id, edge=edge_id)
        bed         = self.domain.quantities['elevation'].centroid_values[vol_id]
        depth_inside = max(q[0] - bed, 0.0)
        normal      = self.domain.get_normal(vol_id, edge_id)

        t     = self.domain.get_time()
        value = self.get_boundary_values(t)
        try:
            perturb = float(value)
        except (ValueError, TypeError):
            perturb = float(value[0])

        stage_wave = self.background_stage + perturb
        h0     = max(self.background_stage - bed, 0.0)
        c0     = (gravity * max(h0, 1e-10)) ** 0.5
        h_wave = max(stage_wave - bed, 0.0)
        c_wave = (gravity * max(h_wave, 1e-10)) ** 0.5

        if depth_inside == 0.0:
            q[0] = stage_wave
            q[1] = 0.0
            q[2] = 0.0
        else:
            u_int   = q[1] / depth_inside
            v_int   = q[2] / depth_inside
            c_int   = (gravity * depth_inside) ** 0.5
            v_n_int = normal[0] * u_int + normal[1] * v_int

            c_ghost   = max((v_n_int + 2.0*c_int - 2.0*c0 + 4.0*c_wave) / 4.0, 0.0)
            v_n_ghost = c0 - 2.0*c_wave + 0.5*v_n_int + c_int
            h_ghost   = c_ghost**2 / gravity

            q[0] = h_ghost + bed
            if c_ghost > 0.0:
                q[1] = h_ghost * v_n_ghost * normal[0]
                q[2] = h_ghost * v_n_ghost * normal[1]
            else:
                q[1] = 0.0
                q[2] = 0.0

        return q

    def evaluate_segment(self, domain, segment_edges):
        """Vectorised form for speed."""

        Stage = domain.quantities['stage']
        Elev  = domain.quantities['elevation']
        Xmom  = domain.quantities['xmomentum']
        Ymom  = domain.quantities['ymomentum']

        ids      = segment_edges
        vol_ids  = domain.boundary_cells[ids]
        edge_ids = domain.boundary_edges[ids]
        Normals  = domain.normals

        n1 = Normals[vol_ids, 2 * edge_ids]
        n2 = Normals[vol_ids, 2 * edge_ids + 1]

        t     = self.domain.get_time()
        value = self.get_boundary_values(t)
        try:
            perturb = float(value)
        except (ValueError, TypeError):
            perturb = float(value[0])

        bed        = Elev.centroid_values[vol_ids]
        stage_wave = self.background_stage + perturb

        # Background wave speed from still-water depth at each cell
        h0 = np.maximum(self.background_stage - bed, 0.0)
        c0 = np.sqrt(gravity * np.maximum(h0, 1e-10))

        # Prescribed wave
        h_wave = np.maximum(stage_wave - bed, 0.0)
        c_wave = np.sqrt(gravity * np.maximum(h_wave, 1e-10))

        # Interior state
        stage_interior = Stage.edge_values[vol_ids, edge_ids]
        depth_inside   = np.maximum(stage_interior - bed, 0.0)
        safe_depth     = np.where(depth_inside > 0.0, depth_inside, 1.0)
        q1     = Xmom.edge_values[vol_ids, edge_ids]
        q2     = Ymom.edge_values[vol_ids, edge_ids]
        c_int   = np.sqrt(gravity * safe_depth)
        v_n_int = (n1 * q1 + n2 * q2) / safe_depth

        # Nonlinear characteristic solution
        c_ghost   = np.maximum((v_n_int + 2.0*c_int - 2.0*c0 + 4.0*c_wave) / 4.0, 0.0)
        v_n_ghost = c0 - 2.0*c_wave + 0.5*v_n_int + c_int
        h_ghost   = c_ghost**2 / gravity
        stage_ghost = h_ghost + bed

        wet = depth_inside > 0.0
        Stage.boundary_values[ids] = np.where(wet, stage_ghost, stage_wave)
        Xmom.boundary_values[ids]  = np.where(wet & (c_ghost > 0.0),
                                               h_ghost * v_n_ghost * n1, 0.0)
        Ymom.boundary_values[ids]  = np.where(wet & (c_ghost > 0.0),
                                               h_ghost * v_n_ghost * n2, 0.0)

