from __future__ import annotations

import anuga
import numpy as num
import math
from . import inlet_enquiry

from anuga.utilities.system_tools import log_to_file
from anuga.utilities.numerical_tools import ensure_numeric


_c_culvert_funcs = None


def _get_c_culvert_funcs():
    """Lazily import the shared C culvert kernels (apply + inlet gather).

    Returns (apply_fn, gather_fn), or (None, None) if the extension is
    unavailable. These are the single source of truth for the Boyd/weir
    discharge + water-transfer physics: mode-1 (here) and mode-2 (the GPU
    culvert batch) call the same C code, so the two compute modes agree.
    """
    global _c_culvert_funcs
    if _c_culvert_funcs is None:
        try:
            from anuga.shallow_water.sw_domain_gpu_ext import (
                culvert_apply_one_host_py, culvert_gather_inlet_host_py)
            _c_culvert_funcs = (culvert_apply_one_host_py,
                                culvert_gather_inlet_host_py)
        except Exception:
            _c_culvert_funcs = (None, None)
    return _c_culvert_funcs


def _c_culvert_type_of(op):
    """CULVERT_TYPE_* code (box=0, pipe=1, weir_trapezoid=2) for the shared C
    kernel, or None if the operator has no C implementation. Handles both the
    serial and parallel variants (same predicates the GPU culvert manager uses)."""
    try:
        from anuga.structures.gpu_culvert_manager import GPUCulvertManager
    except Exception:
        return None
    if GPUCulvertManager._is_boyd_box(op):
        return 0
    if GPUCulvertManager._is_boyd_pipe(op):
        return 1
    if GPUCulvertManager._is_weir_trapezoid(op):
        return 2
    return None


def _can_use_c_culvert(op):
    """True if this operator's per-step update can run through the shared C
    kernel: a supported Boyd/weir culvert that is fully local to this rank (so no
    MPI is needed) with a C-covered outflow model. Cross-boundary parallel
    culverts keep the Python + MPI path."""
    if _c_culvert_type_of(op) is None:
        return False
    if op.inlets[0] is None or op.inlets[1] is None:
        return False
    # The C kernel models jet OR zero outflow momentum; Boyd/weir always satisfy
    # this (zero_outflow_momentum == not use_momentum_jet).
    if not (op.use_momentum_jet or op.zero_outflow_momentum):
        return False
    # Parallel operator: only take the C path when this rank fully owns the
    # culvert (master, both enquiry points, both inlets). Otherwise MPI is needed.
    if hasattr(op, 'myid'):
        myid = op.myid
        if myid != op.master_proc:
            return False
        if op.enquiry_proc[0] != myid or op.enquiry_proc[1] != myid:
            return False
        if op.inlet_master_proc[0] != myid or op.inlet_master_proc[1] != myid:
            return False
    return _get_c_culvert_funcs()[0] is not None


def _call_c_culvert(op):
    """Per-step culvert update via the shared C kernel — bit-identical to the
    Python discharge_routine + transfer, and to the mode-2 GPU culvert batch.
    This is the single source of truth for the Boyd/weir physics on the mode-1
    (host) path; only reached for fully-local culverts (see _can_use_c_culvert)."""
    apply_fn, gather_fn = _get_c_culvert_funcs()
    ctype = _c_culvert_type_of(op)
    d = op.domain
    stage_c = d.quantities['stage'].centroid_values
    xmom_c = d.quantities['xmomentum'].centroid_values
    ymom_c = d.quantities['ymomentum'].centroid_values
    bed_c = d.quantities['elevation'].centroid_values
    i0, i1 = op.inlets[0], op.inlets[1]
    ov0 = i0.outward_culvert_vector
    ov1 = i1.outward_culvert_vector
    inv0, inv1 = i0.invert_elevation, i1.invert_elevation
    h0 = 1 if inv0 is not None else 0
    h1 = 1 if inv1 is not None else 0

    def gather(inl):
        idx = num.ascontiguousarray(inl.triangle_indices, dtype=num.intc)
        areas = num.ascontiguousarray(d.areas[inl.triangle_indices], dtype=float)
        a_stage, a_depth, a_xmom, a_ymom = gather_fn(
            idx, areas, stage_c, xmom_c, ymom_c, bed_c, inl.get_area())
        return (inl.get_enquiry_stage(), inl.get_enquiry_xmom(),
                inl.get_enquiry_ymom(), inl.get_enquiry_elevation(),
                a_stage, a_depth, a_xmom, a_ymom, inl.get_area())

    g0 = gather(i0)
    g1 = gather(i1)
    # Index 5 of a gather tuple is that inlet's average depth — the value the C
    # kernel derived its new depths from, and so the one the write-back levels
    # the surface from (see Inlet.set_average_depth).
    average_depths = (g0[5], g1[5])
    res = apply_fn(
        ctype, d.g,
        float(getattr(op, 'culvert_width', 0.0)),
        float(getattr(op, 'culvert_height', 0.0)),
        float(getattr(op, 'culvert_diameter', 0.0)),
        float(getattr(op, 'culvert_z1', 0.0)),
        float(getattr(op, 'culvert_z2', 0.0)),
        op.culvert_length, op.manning, op.sum_loss,
        float(getattr(op, 'culvert_blockage', getattr(op, 'blockage', 0.0))),
        float(getattr(op, 'culvert_barrels', getattr(op, 'barrels', 1.0))),
        1 if op.use_velocity_head else 0,
        1 if op.use_momentum_jet else 0,
        1 if getattr(op, 'use_old_momentum_method', True) else 0,
        1 if getattr(op, 'always_use_Q_wetdry_adjustment', True) else 0,
        op.max_velocity, op.smoothing_timescale,
        ov0[0], ov0[1], ov1[0], ov1[1],
        float(inv0) if h0 else 0.0, float(inv1) if h1 else 0.0, h0, h1,
        op.smooth_delta_total_energy, op.smooth_Q, d.get_timestep(),
        *g0, *g1)
    (sdte, sQ, inflow_idx, nid, nix, niy, nod, nox, noy,
     rg, rd, rv, rde, rdte, ocd) = res

    op.smooth_delta_total_energy = sdte
    op.smooth_Q = sQ
    inflow = op.inlets[inflow_idx]
    outflow = op.inlets[1 - inflow_idx]
    op.inflow = inflow
    op.outflow = outflow

    inflow.set_average_depth(nid, average_depths[inflow_idx])
    inflow.set_average_momenta(nix, niy)
    outflow.set_average_depth(nod, average_depths[1 - inflow_idx])
    outflow.set_average_momenta(nox, noy)

    # Reporting (matches the Python discharge routine and the mode-2 readback)
    op.accumulated_flow += rg
    if d.yieldstep:
        op.discharge_abs_timemean += rg / d.yieldstep
    op.discharge = rd
    op.velocity = rv
    op.driving_energy = rde
    op.delta_total_energy = rdte
    op.outlet_depth = ocd


class Structure_operator(anuga.Operator):
    """Structure Operator - transfer water from one rectangular box to another.
    Sets up the geometry of problem

    This is the base class for structures (culverts, pipes, bridges etc). Inherit from this class (and overwrite
    discharge_routine method for specific subclasses)

    Input: Two points, pipe_size (either diameter or width, depth),
    mannings_rougness,
    """

    counter = 0

    def __init__(self,
                 domain,
                 end_points: list | num.ndarray | None = None,
                 exchange_lines: list | num.ndarray | None = None,
                 enquiry_points: list | num.ndarray | None = None,
                 invert_elevations: list | num.ndarray | None = None,
                 width: float | None = None,
                 height: float | None = None,
                 diameter: float | None = None,
                 z1: float | None = None,
                 z2: float | None = None,
                 blockage: float | None = None,
                 barrels: float | None = None,
                 #culvert_slope=None,
                 apron: float | None = None,
                 manning: float | None = None,
                 enquiry_gap: float | None = None,
                 use_momentum_jet: bool = False,
                 zero_outflow_momentum: bool = True,
                 use_old_momentum_method: bool = True,
                 always_use_Q_wetdry_adjustment: bool = True,
                 force_constant_inlet_elevations: bool = False,
                 description: str | None = None,
                 label: str | None = None,
                 structure_type: str | None = None,
                 logging: bool | None = None,
                 verbose: bool | None = None) -> None:

        """
        exchange_lines define the input lines for each inlet.

        If end_points = None, then the culvert_vector is calculated in the
        directions from the centre of echange_line[0] to centre of exchange_line[1}

        If end_points != None, then culvert_vector is unit vector in direction
        end_point[1] - end_point[0]
        """

        anuga.Operator.__init__(self,domain)

        self.master_proc = 0
        self.end_points = ensure_numeric(end_points, float)
        self.exchange_lines = ensure_numeric(exchange_lines, float)
        self.enquiry_points = ensure_numeric(enquiry_points, float)
        self.invert_elevations = ensure_numeric(invert_elevations, float)

        assert self.end_points is None or self.exchange_lines is None


        if height is None:
            height = width

        if width is None:
            width = diameter

        if apron is None:
            apron = width


        assert width is not None


        self.width  = width
        self.height = height
        self.diameter = diameter
        self.z1 = z1
        self.z2 = z2
        self.blockage = blockage
        self.barrels = barrels
        self.apron  = apron
        self.manning = manning
        self.enquiry_gap = enquiry_gap
        self.use_momentum_jet = use_momentum_jet
        self.zero_outflow_momentum = zero_outflow_momentum
        if use_momentum_jet and zero_outflow_momentum:
            msg = "Can't have use_momentum_jet and zero_outflow_momentum both True"
            raise Exception(msg)
        self.use_old_momentum_method = use_old_momentum_method
        self.always_use_Q_wetdry_adjustment = always_use_Q_wetdry_adjustment


        if description is None:
            self.description = ' '
        else:
            self.description = description

        if label is None:
            self.label = "structure_%g" % Structure_operator.counter
        else:
            self.label = label + '_%g' % Structure_operator.counter

        if structure_type is None:
            self.structure_type = 'generic structure'
        else:
            self.structure_type = structure_type

        self.verbose = verbose

        # Keep count of structures
        Structure_operator.counter += 1

        # Slots for recording current statistics
        self.accumulated_flow = 0.0
        self.discharge = 0.0
        self.discharge_abs_timemean = 0.0
        self.velocity = 0.0
        self.outlet_depth = 0.0
        self.delta_total_energy = 0.0
        self.driving_energy = 0.0

        if exchange_lines is not None:
            self.__process_skew_culvert()
        elif end_points is not None:
            self.__process_non_skew_culvert()
        else:
            raise Exception('Define either exchange_lines or end_points')


        self.inlets = []
        line0 = self.exchange_lines[0] #self.inlet_lines[0]
        if self.apron is None:
            poly0 = line0
        else:
            offset = -self.apron*self.outward_vector_0
            #print line0
            #print offset
            poly0 = num.array([ line0[0], line0[1], line0[1]+offset, line0[0]+offset])
            #print poly0
        if self.invert_elevations is None:
            invert_elevation0 = None
        else:
            invert_elevation0 = self.invert_elevations[0]

        enquiry_point0 = self.enquiry_points[0]

        #outward_vector0 = - self.culvert_vector
        self.inlets.append(inlet_enquiry.Inlet_enquiry(
                           self.domain,
                           poly0,
                           enquiry_point0,
                           invert_elevation = invert_elevation0,
                           outward_culvert_vector = self.outward_vector_0,
                           verbose = self.verbose))

        if force_constant_inlet_elevations:
            # Try to enforce a constant inlet elevation
            inlet_global_elevation = self.inlets[-1].get_average_elevation()
            self.inlets[-1].set_elevations(inlet_global_elevation)

        tris_0 = self.inlets[0].triangle_indices
        #print tris_0
        #print self.domain.centroid_coordinates[tris_0]

        line1 = self.exchange_lines[1]
        if self.apron is None:
            poly1 = line1
        else:
            offset = -self.apron*self.outward_vector_1
            #print line1
            #print offset
            poly1 = num.array([ line1[0], line1[1], line1[1]+offset, line1[0]+offset])
            #print poly1

        if self.invert_elevations is None:
            invert_elevation1 = None
        else:
            invert_elevation1 = self.invert_elevations[1]
        enquiry_point1 = self.enquiry_points[1]

        self.inlets.append(inlet_enquiry.Inlet_enquiry(
                           self.domain,
                           poly1,
                           enquiry_point1,
                           invert_elevation = invert_elevation1,
                           outward_culvert_vector = self.outward_vector_1,
                           verbose = self.verbose))

        if force_constant_inlet_elevations:
            # Try to enforce a constant inlet elevation
            inlet_global_elevation = self.inlets[-1].get_average_elevation()
            self.inlets[-1].set_elevations(inlet_global_elevation)

        tris_1 = self.inlets[1].triangle_indices

        self.set_logging(logging)





    def __call__(self):
        # Tracers ride with the water (anuga/structures/inlet_tracers.py):
        # capture both regions, move the water, then move the tracer.
        tracers = getattr(self, '_structure_tracers', None)
        if tracers is None:
            from anuga.structures.inlet_tracers import StructureTracers
            tracers = self._structure_tracers = StructureTracers(self.domain)
        state = tracers.capture(self.inlets) if tracers.active() else None
        self._transfer_water()
        if state is not None:
            tracers.apply(state)

    def _transfer_water(self):
        """One step of the structure's water transfer (no tracers)."""

        if _can_use_c_culvert(self):
            _call_c_culvert(self)
            return

        timestep = self.domain.get_timestep()

        # implement different types of structures by redefining
        # the discharge_routine
        Q, barrel_speed, outlet_depth = self.discharge_routine()

        old_inflow_depth = self.inflow.get_average_depth()
        old_inflow_stage = self.inflow.get_average_stage()
        old_inflow_xmom = self.inflow.get_average_xmom()
        old_inflow_ymom = self.inflow.get_average_ymom()

        # Implement the update of flow over a timestep by
        # using a semi-implict update. This ensures that
        # the update does not create a negative depth
        if old_inflow_depth > 0.0 :
            dt_Q_on_d = timestep*Q/old_inflow_depth
        else:
            dt_Q_on_d = 0.0

        always_use_Q_wetdry_adjustment = self.always_use_Q_wetdry_adjustment
        use_Q_wetdry_adjustment = ((always_use_Q_wetdry_adjustment) |\
            (old_inflow_depth*self.inflow.get_area() <= Q*timestep))
        # If we use the 'Q_adjustment', then the discharge is rescaled so that
        # the depth update is:
        #    new_inflow_depth*inflow_area =
        #    old_inflow_depth*inflow_area -
        #    timestep*Q*(new_inflow_depth/old_inflow_depth)
        # The last term in () is a wet-dry improvement trick (which rescales Q)
        #
        # Before Feb 2015 this rescaling was always done (even if the flow was
        # not near wet-dry). Now if use_old_Q_adjustment=False, the rescaling
        # is only done if required to avoid drying
        #
        #
        factor = 1.0/(1.0 + dt_Q_on_d/self.inflow.get_area())


        if use_Q_wetdry_adjustment:
            # FIXME:
            new_inflow_depth = old_inflow_depth*factor

            if(old_inflow_depth>0.):
                timestep_star = timestep*new_inflow_depth/old_inflow_depth
            else:
                timestep_star = 0.

        else:
            new_inflow_depth = old_inflow_depth - timestep*Q/self.inflow.get_area()
            timestep_star = timestep

        if(self.use_old_momentum_method):
            # This method is here for consistency with the old version of the
            # routine
            new_inflow_xmom = old_inflow_xmom*factor
            new_inflow_ymom = old_inflow_ymom*factor

        else:
            # For the momentum balance, note that Q also advects the momentum,
            # The volumetric momentum flux should be Q*momentum/depth, where
            # momentum has an average value of new_inflow_mom (or
            # old_inflow_mom).  We use old_inflow_depth for depth
            #
            #     new_inflow_xmom*inflow_area =
            #     old_inflow_xmom*inflow_area -
            #     [timestep*Q]*new_inflow_xmom/old_inflow_depth
            # and:
            #     new_inflow_ymom*inflow_area =
            #     old_inflow_ymom*inflow_area -
            #     [timestep*Q]*new_inflow_ymom/old_inflow_depth
            #
            # The units balance: m^2/s*m^2 = m^2/s*m^2 - s*m^3/s*m^2/s *m^(-1)
            #
            if old_inflow_depth > 0.:
                if use_Q_wetdry_adjustment:
                    # Replace dt*Q with dt*Q*new_inflow_depth/old_inflow_depth = dt_Q_on_d*new_inflow_depth
                    factor2 = 1.0/(1.0 + dt_Q_on_d*new_inflow_depth/(old_inflow_depth*self.inflow.get_area()))
                else:
                    factor2 = 1.0/(1.0 + timestep*Q/(old_inflow_depth*self.inflow.get_area()))
            else:
                factor2 = 0.

            new_inflow_xmom = old_inflow_xmom*factor2
            new_inflow_ymom = old_inflow_ymom*factor2

        self.inflow.set_average_depth(new_inflow_depth, old_inflow_depth)

        #inflow.set_xmoms(Q/inflow.get_area())
        #inflow.set_ymoms(0.0)

        self.inflow.set_average_momenta(new_inflow_xmom, new_inflow_ymom)

        loss = (old_inflow_depth - new_inflow_depth)*self.inflow.get_area()
        xmom_loss = (old_inflow_xmom - new_inflow_xmom)*self.inflow.get_area()
        ymom_loss = (old_inflow_ymom - new_inflow_ymom)*self.inflow.get_area()

        # set outflow
        outflow_extra_depth = Q*timestep_star/self.outflow.get_area()
        outflow_direction = - self.outflow.outward_culvert_vector
        #outflow_extra_momentum = outflow_extra_depth*barrel_speed*outflow_direction

        gain = outflow_extra_depth*self.outflow.get_area()

        #print gain, loss
        assert num.allclose(gain-loss, 0.0)

        # Stats

        #print('gain', gain)

        self.accumulated_flow += gain
        self.discharge  = Q*timestep_star/timestep
        self.discharge_abs_timemean += gain/self.domain.yieldstep
        self.velocity =   barrel_speed
        self.outlet_depth = outlet_depth


        old_outflow_depth = self.outflow.get_average_depth()
        new_outflow_depth = old_outflow_depth + outflow_extra_depth

        self.outflow.set_average_depth(new_outflow_depth, old_outflow_depth)

        if self.use_momentum_jet:
            # FIXME (SR) Review momentum to account for possible hydraulic jumps at outlet
            # FIXME (GD) Depending on barrel speed I think this will be either
            # a source or sink of momentum (considering the momentum losses
            # above). Might not always be reasonable.
            #new_outflow_xmom = self.outflow.get_average_xmom() + outflow_extra_momentum[0]
            #new_outflow_ymom = self.outflow.get_average_ymom() + outflow_extra_momentum[1]
            new_outflow_xmom = barrel_speed*new_outflow_depth*outflow_direction[0]
            new_outflow_ymom = barrel_speed*new_outflow_depth*outflow_direction[1]

        elif self.zero_outflow_momentum:
            new_outflow_xmom = 0.0
            new_outflow_ymom = 0.0
            #new_outflow_xmom = outflow.get_average_xmom()
            #new_outflow_ymom = outflow.get_average_ymom()

        else:
            # Add the momentum lost from the inflow to the outflow. For
            # structures where barrel_speed is unknown + direction doesn't
            # change from inflow to outflow
            new_outflow_xmom = self.outflow.get_average_xmom() + xmom_loss/self.outflow.get_area()
            new_outflow_ymom = self.outflow.get_average_ymom() + ymom_loss/self.outflow.get_area()

        self.outflow.set_average_momenta(new_outflow_xmom, new_outflow_ymom)



    def set_culvert_height(self, height: float) -> None:

        self.culvert_height = height

    def set_culvert_width(self, width: float) -> None:

        self.culvert_width = width

    def set_culvert_z1(self, z1: float) -> None:

        self.culvert_z1 = z1

    def set_culvert_z2(self, z2: float) -> None:

        self.culvert_z2 = z2

    def set_culvert_blockage(self, blockage: float) -> None:

        self.culvert_blockage = blockage

    def set_culvert_barrels(self, barrels: float) -> None:

        self.culvert_barrels = barrels


    def __process_non_skew_culvert(self):

        """Create lines at the end of a culvert inlet and outlet.
        At either end two lines will be created; one for the actual flow to pass through and one a little further away
        for enquiring the total energy at both ends of the culvert and transferring flow.
        """

        self.culvert_vector = self.end_points[1] - self.end_points[0]
        self.culvert_length = math.sqrt(num.sum(self.culvert_vector**2))
        assert self.culvert_length > 0.0, 'The length of culvert is less than 0'

        self.culvert_vector /= self.culvert_length
        self.outward_vector_0 =   self.culvert_vector
        self.outward_vector_1 = - self.culvert_vector


        culvert_normal = num.array([-self.culvert_vector[1], self.culvert_vector[0]])  # Normal vector
        w = 0.5*self.width*culvert_normal # Perpendicular vector of 1/2 width

        self.exchange_lines = []

        # Build exchange polyline and enquiry point
        if self.enquiry_points is None:

            gap = (self.apron + self.enquiry_gap)*self.culvert_vector
            self.enquiry_points = []

            for i in [0, 1]:
                p0 = self.end_points[i] + w
                p1 = self.end_points[i] - w
                self.exchange_lines.append(num.array([p0, p1]))
                ep = self.end_points[i] + (2*i - 1)*gap #(2*i - 1) determines the sign of the points
                self.enquiry_points.append(ep)

        else:
            for i in [0, 1]:
                p0 = self.end_points[i] + w
                p1 = self.end_points[i] - w
                self.exchange_lines.append(num.array([p0, p1]))


    def __process_skew_culvert(self):

        """Compute skew culvert.
        If exchange lines are given, the enquiry points are determined. This is for enquiring
        the total energy at both ends of the culvert and transferring flow.
        """

        centre_point0 = 0.5*(self.exchange_lines[0][0] + self.exchange_lines[0][1])
        centre_point1 = 0.5*(self.exchange_lines[1][0] + self.exchange_lines[1][1])

        n_exchange_0 = len(self.exchange_lines[0])
        n_exchange_1 = len(self.exchange_lines[1])

        assert n_exchange_0 == n_exchange_1, 'There should be the same number of points in both exchange_lines'

        if n_exchange_0 == 2:

            if self.end_points is None:
                self.culvert_vector = centre_point1 - centre_point0
            else:
                self.culvert_vector = self.end_points[1] - self.end_points[0]

            self.outward_vector_0 =   self.culvert_vector
            self.outward_vector_1 = - self.culvert_vector


        elif n_exchange_0 == 4:

            self.outward_vector_0 = self.exchange_lines[0][3] - self.exchange_lines[0][2]
            self.outward_vector_1 = self.exchange_lines[1][3] - self.exchange_lines[1][2]

            self.culvert_vector = centre_point1 - centre_point0

        else:
            raise Exception('n_exchange_0 != 2 or 4')


        self.culvert_length = math.sqrt(num.sum(self.culvert_vector**2))
        assert self.culvert_length > 0.0, 'The length of culvert is less than 0'
        self.culvert_vector /= self.culvert_length

        outward_vector_0_length = math.sqrt(num.sum(self.outward_vector_0**2))
        assert outward_vector_0_length > 0.0, 'The length of outlet_vector_0 is less than 0'
        self.outward_vector_0 /= outward_vector_0_length

        outward_vector_1_length = math.sqrt(num.sum(self.outward_vector_1**2))
        assert outward_vector_1_length > 0.0, 'The length of outlet_vector_1 is less than 0'
        self.outward_vector_1 /= outward_vector_1_length


        if self.enquiry_points is None:

            gap = (self.apron + self.enquiry_gap)*self.culvert_vector

            self.enquiry_points = []

            self.enquiry_points.append(centre_point0 - gap)
            self.enquiry_points.append(centre_point1 + gap)


    def discharge_routine(self):

        msg = 'Need to implement '
        raise


    def statistics(self):


        message  = '=====================================\n'
        message += 'Structure Operator: %s\n' % self.label
        message += '=====================================\n'

        message += 'Structure Type: %s\n' % self.structure_type

        message += 'Description\n'
        message += '%s' % self.description

        #add the culvert dimensions, blockage factor here
        if self.structure_type == 'boyd_pipe':
            message += 'Culvert Diameter: %s\n'% self.diameter
            message += 'Culvert Blockage: %s\n'% self.blockage
            message += 'No.  of  barrels: %s\n'% self.barrels
        elif self.structure_type == 'boyd_box':
            message += 'Culvert  Height: %s\n'% self.height
            message += 'Culvert    Width: %s\n'% self.width
            message += 'Culvert Blockage: %s\n'% self.blockage
            message += 'No.  of  barrels: %s\n'% self.barrels
        else:
            message += 'Culvert Height: %s\n'% self.height
            message += 'Culvert  Width: %s\n'% self.width
            message += 'Batter Slope 1: %s\n'% self.z1
            message += 'Batter Slope 2: %s\n'% self.z2
            message += 'Culvert Blockage: %s\n'% self.blockage
            message += 'No.  of  barrels: %s\n'% self.barrels

        message += '\n'

        for i, inlet in enumerate(self.inlets):
            message += '-------------------------------------\n'
            message +=  'Inlet %i\n' % i
            message += '-------------------------------------\n'

            message += 'inlet triangle indices and centres and elevations\n'
            message += '%s' % inlet.triangle_indices
            message += '\n'

            message += '%s' % self.domain.get_centroid_coordinates()[inlet.triangle_indices]
            message += '\n'

            elev = self.domain.quantities['elevation'].centroid_values[inlet.triangle_indices]
            message += '%s' % elev
            message += '\n'

            elevation_range = elev.max() - elev.min()
            if not num.allclose(elevation_range, 0.):
                message += 'Warning: non-constant inlet elevation can cause well-balancing problems'

            message += 'region\n'
            message += '%s' % inlet.region
            message += '\n'

        message += '=====================================\n'

        return message


    def print_statistics(self):

        print(self.statistics())


    def print_timestepping_statistics(self):

        message = '---------------------------\n'
        message += 'Structure report for %s:\n' % self.label
        message += '--------------------------\n'
        message += 'Type: %s\n' % self.structure_type

        message += 'inlets[0]_enquiry_depth [m]:  %.2f\n' %self.inlets[0].get_enquiry_depth()
        message += 'inlets[0]_enquiry_speed [m/s]:  %.2f\n' %self.inlets[0].get_enquiry_speed()
        message += 'inlets[0]_enquiry_stage [m]:  %.2f\n' %self.inlets[0].get_enquiry_stage()
        message += 'inlets[0]_enquiry_elevation [m]:  %.2f\n' %self.inlets[0].get_enquiry_elevation()
        message += 'inlets[0]_average_depth [m]:  %.2f\n' %self.inlets[0].get_average_depth()
        message += 'inlets[0]_average_speed [m/s]:  %.2f\n' %self.inlets[0].get_average_speed()
        message += 'inlets[0]_average_stage [m]:  %.2f\n' %self.inlets[0].get_average_stage()
        message += 'inlets[0]_average_elevation [m]:  %.2f\n' %self.inlets[0].get_average_elevation()

        message += '\n'

        message += 'inlets[1]_enquiry_depth [m]:  %.2f\n' %self.inlets[1].get_enquiry_depth()
        message += 'inlets[1]_enquiry_speed [m/s]:  %.2f\n' %self.inlets[1].get_enquiry_speed()
        message += 'inlets[1]_enquiry_stage [m]:  %.2f\n' %self.inlets[1].get_enquiry_stage()
        message += 'inlets[1]_enquiry_elevation [m]:  %.2f\n' %self.inlets[1].get_enquiry_elevation()

        message += 'inlets[1]_average_depth [m]:  %.2f\n' %self.inlets[1].get_average_depth()
        message += 'inlets[1]_average_speed [m/s]:  %.2f\n' %self.inlets[1].get_average_speed()
        message += 'inlets[1]_average_stage [m]:  %.2f\n' %self.inlets[1].get_average_stage()
        message += 'inlets[1]_average_elevation [m]:  %.2f\n' %self.inlets[1].get_average_elevation()


        message += 'Discharge [m^3/s]: %.2f\n' % self.discharge
        message += 'Discharge_function_value [m^3/s]: %.2f\n' % self.discharge_abs_timemean
        message += 'Velocity  [m/s]: %.2f\n' % self.velocity
        message += 'Outlet Depth  [m]: %.2f\n' % self.outlet_depth
        message += 'Accumulated Flow [m^3]: %.2f\n' % self.accumulated_flow
        message += 'Inlet Driving Energy %.2f\n' % self.driving_energy
        message += 'Delta Total Energy %.2f\n' % self.delta_total_energy
        message += 'Control at this instant: %s\n' % self.case




        print(message)


    def set_logging(self, flag=True):

        self.logging = flag

        # If flag is true open file with mode = "w" to form a clean file for logging
        if self.logging:
            self.log_filename = self.domain.get_datadir() + '/' + self.label + '.log'
            log_to_file(self.log_filename, self.statistics(), mode='w')
            log_to_file(self.log_filename, 'time, discharge_instantaneous, discharge_abs_timemean, velocity, accumulated_flow, driving_energy_instantaneous, delta_total_energy_instantaneous')

            #log_to_file(self.log_filename, self.culvert_type)


    def timestepping_statistics(self):

        message  = '%.5f, ' % self.domain.get_time()
        message += '%.5f, ' % self.discharge
        message += '%.5f, ' % self.discharge_abs_timemean
        message += '%.5f, ' % self.velocity
        message += '%.5f, ' % self.accumulated_flow
        message += '%.5f, ' % self.driving_energy
        message += '%.5f' % self.delta_total_energy

        # Reset discharge_abs_timemean since last time this function was called
        # (FIXME: This assumes that the function is called only just after a
        # yield step)
        self.discharge_abs_timemean = 0.

        return message


    def get_inlets(self) -> list:

        return self.inlets


    def get_culvert_length(self) -> float:

        return self.culvert_length



    def get_culvert_slope(self) -> float:

        inlet0 = self.inlets[0]
        inlet1 = self.inlets[1]

        elev0 = inlet0.get_enquiry_invert_elevation()
        elev1 = inlet1.get_enquiry_invert_elevation()

        return (elev1-elev0)/self.get_culvert_length()



    def get_culvert_width(self) -> float | None:

        return self.width


    def get_culvert_diameter(self) -> float | None:

            return self.diameter


    def get_culvert_height(self) -> float | None:

        return self.height

    def get_culvert_z1(self) -> float | None:

        return self.z1

    def get_culvert_z2(self) -> float | None:

        return self.z2

    def get_culvert_blockage(self) -> float | None:

        return self.blockage

    def get_culvert_barrels(self) -> float | None:

        return self.barrels

    def get_culvert_apron(self) -> float | None:

        return self.apron


    def get_master_proc(self) -> int:

        return 0



    #--------------------------------------------------------
    # Set of enquiry functions so that in the sequential and paralle case
    # we can get equiry info fron the master Proc
    #---------------------------------------------------------

    def get_enquiry_stages(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_stage()
        enq1 = self.inlets[1].get_enquiry_stage()

        return [enq0, enq1]


    def get_enquiry_depths(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_depth()
        enq1 = self.inlets[1].get_enquiry_depth()

        return [enq0, enq1]


    def get_enquiry_positions(self) -> list[num.ndarray]:

        enq0 = self.inlets[0].get_enquiry_position()
        enq1 = self.inlets[1].get_enquiry_position()

        return [enq0, enq1]


    def get_enquiry_xmoms(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_xmom()
        enq1 = self.inlets[1].get_enquiry_xmom()

        return [enq0, enq1]

    def get_enquiry_ymoms(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_ymom()
        enq1 = self.inlets[1].get_enquiry_ymom()

        return [enq0, enq1]


    def get_enquiry_elevations(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_elevation()
        enq1 = self.inlets[1].get_enquiry_elevation()

        return [enq0, enq1]



    def get_enquiry_water_depths(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_water_depth()
        enq1 = self.inlets[1].get_enquiry_water_depth()

        return [enq0, enq1]


    def get_enquiry_invert_elevations(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_invert_elevation()
        enq1 = self.inlets[1].get_enquiry_invert_elevation()

        return [enq0, enq1]


    def get_enquiry_velocitys(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_velocity()
        enq1 = self.inlets[1].get_enquiry_velocity()

        return [enq0, enq1]


    def get_enquiry_xvelocitys(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_xvelocity()
        enq1 = self.inlets[1].get_enquiry_xvelocity()

        return [enq0, enq1]

    def get_enquiry_yvelocitys(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_yvelocity()
        enq1 = self.inlets[1].get_enquiry_yvelocity()

        return [enq0, enq1]


    def get_enquiry_speeds(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_speed()
        enq1 = self.inlets[1].get_enquiry_speed()

        return [enq0, enq1]


    def get_enquiry_velocity_heads(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_velocity_head()
        enq1 = self.inlets[1].get_enquiry_velocity_head()

        return [enq0, enq1]


    def get_enquiry_total_energys(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_total_energy()
        enq1 = self.inlets[1].get_enquiry_total_energy()

        return [enq0, enq1]


    def get_enquiry_specific_energys(self) -> list[float]:

        enq0 = self.inlets[0].get_enquiry_specific_energy()
        enq1 = self.inlets[1].get_enquiry_specific_energy()

        return [enq0, enq1]

