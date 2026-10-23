
import anuga
import numpy
from . import inlet

import warnings

class Inlet_operator(anuga.Operator):
    """Inlet Operator - add water to an inlet.
    Sets up the geometry of problem

    Inherit from this class (and overwrite
    discharge_routine method for specific subclasses)

    """


    def __init__(self,
                 domain,
                 region,
                 Q = 0.0,
                 velocity = None,
                 zero_velocity = False,
                 default = 0.0,
                 tracer_concentrations = None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
        """Inlet Operator - add water to a domain via an inlet.

        :param domain: Specify domain
        :param region: Apply Inlet flow over a region (which can be a Region, Polygon or line)
        :param Q: function(t) or scalar discharge (m^3/s)
        :param velocity: Optional [u,v] to set velocity of applied discharge
        :param zero_velocity: If set to True, velocity of inlet region set to 0
        :param default: If outside time domain of the Q function, use this default discharge
        :param tracer_concentrations: dict of tracer name -> concentration (float
            or function of t) of the water the inlet adds. Tracers not named come
            in at zero concentration. Water the inlet removes always takes its
            tracer with it, at the inlet's mean concentration.
        :param description: Describe the Inlet_operator
        :param label: Give Inlet_operator a label (name)
        :param verbose: Provide verbose output



        Example:

        >>> inflow_region  = anuga.Region(domain, center=[0.0,0.0], radius=1.0)
        >>> inflow = anuga.Inlet_operator(domain, inflow_region, Q = lambda t : 1 + 0.5*math.sin(t/60))

        """

        anuga.Operator.__init__(self, domain, description, label, logging, verbose)

        self.inlet = inlet.Inlet(self.domain, region, verbose= verbose)

        # constant or function of time, m^3/s
        self.Q = Q

        if velocity is not None:
            assert len(velocity)==2

        self.velocity = velocity

        self.zero_velocity = zero_velocity

        self.applied_Q = 0.0

        self.total_applied_volume = 0.0
        self.total_requested_volume = 0.0

        self.set_default(default)

        from anuga.structures.inlet_tracers import InletTracers
        self.tracers = InletTracers(self.domain, tracer_concentrations,
                                    label=getattr(self, 'label', ''))

        self.activate_logging()

        # GPU state (lazy init on first __call__)
        self._gpu_op_id = None
        self._gpu_initialized = False

    def _init_gpu(self):
        """Initialize GPU inlet operator (lazy, called on first __call__ in GPU mode).

        On any failure — a missing/partial GPU interface, or the
        ``MAX_INLET_OPERATORS`` slot limit being exceeded — leave
        ``self._gpu_initialized`` False and return so ``__call__`` falls through
        to the Python path rather than crashing. Mirrors the graceful-fallback
        contract of ``Rate_operator._init_gpu``.
        """
        if self._gpu_initialized:
            return

        # If the GPU interface isn't ready, fall back to the Python path silently.
        gpu_interface = getattr(self.domain, 'gpu_interface', None)
        if gpu_interface is None or getattr(gpu_interface, 'gpu_dom', None) is None:
            return

        try:
            from anuga.shallow_water import sw_domain_gpu_ext as gpu_ext

            tri_indices = numpy.ascontiguousarray(
                self.inlet.triangle_indices, dtype=numpy.intc)
            areas = numpy.ascontiguousarray(
                self.inlet.get_areas(), dtype=numpy.float64)

            op_id = gpu_ext.init_inlet_operator(gpu_interface.gpu_dom, tri_indices, areas)
        except Exception as e:
            warnings.warn(
                f"GPU inlet operator init failed for "
                f"'{getattr(self, 'label', repr(self))}' ({e}); "
                f"falling back to the Python path.",
                stacklevel=2,
            )
            return

        if op_id < 0:
            warnings.warn(
                f"Failed to register GPU inlet operator "
                f"'{getattr(self, 'label', repr(self))}' (the C init "
                f"returned {op_id}: slot limit exceeded, out-of-range "
                f"triangle indices, or allocation failure; see stderr). "
                f"Falling back to the Python path. If it is the slot limit, "
                f"reduce the number of Inlet_operator instances or increase "
                f"MAX_INLET_OPERATORS in gpu_domain.h.",
                stacklevel=2,
            )
            return

        self._gpu_op_id = op_id
        self._gpu_initialized = True

    def _call_gpu(self):
        """GPU path for __call__ - transfers only inlet data (~6KB)."""
        from anuga.shallow_water import sw_domain_gpu_ext as gpu_ext

        gpu_dom = self.domain.gpu_interface.gpu_dom
        op_id = self._gpu_op_id
        timestep = self.domain.get_timestep()
        t = self.domain.get_time()

        # Get current volume from GPU (small reduction, returns scalar)
        current_volume = gpu_ext.inlet_get_volume_gpu(gpu_dom, op_id)
        total_area = self.inlet.area

        assert current_volume >= 0.0

        # Compute Q on CPU (scalar, cheap)
        Q1 = self.update_Q(t)
        Q2 = self.update_Q(t + timestep)
        Q = 0.5 * (Q1 + Q2)
        # update_Q may return a length-1 array (time-varying / file Q); the GPU
        # kernel (and the scalar volume tracking below) need a Python float —
        # numpy 2.x rejects float() on a >0-d array.
        try:
            Q = float(Q)
        except (TypeError, ValueError):
            Q = float(numpy.asarray(Q).ravel()[0])
        volume = Q * timestep

        self.applied_Q = Q

        # Get velocities from GPU (small D2H)
        vel_u, vel_v = gpu_ext.inlet_get_velocities_gpu(gpu_dom, op_id)

        has_velocity = 1 if self.velocity is not None else 0
        ext_vel_u = self.velocity[0] if has_velocity else 0.0
        ext_vel_v = self.velocity[1] if has_velocity else 0.0
        zero_vel = 1 if self.zero_velocity else 0

        # Tracers the water carries (see inlet_tracers.py): the kernel moves
        # them with the water; here only the inflow concentrations and the
        # bookkeeping.
        c_in = dmass = None
        if self.tracers.active():
            c_in = self.tracers.inflow_concentrations(t, timestep)
            dmass = numpy.zeros(self.domain.number_of_tracers)

        # Apply on GPU (handles all 3 cases)
        actual_volume = gpu_ext.inlet_apply_gpu(
            gpu_dom, op_id, volume, current_volume, total_area,
            vel_u, vel_v, has_velocity, ext_vel_u, ext_vel_v, zero_vel,
            c_in, dmass)
        if dmass is not None:
            self.tracers.record(dmass)

        # Update tracking variables
        self.total_requested_volume += volume
        if volume >= 0.0:
            self.domain.fractional_step_volume_integral += volume
        elif current_volume + volume >= 0.0:
            self.domain.fractional_step_volume_integral += volume
        else:
            self.applied_Q = -current_volume / timestep
            self.domain.fractional_step_volume_integral -= current_volume

        self.total_applied_volume += actual_volume

    def __call__(self):

        # GPU path: skip full domain sync, transfer only inlet data.
        # NOT taken when _gpu_host_writes_suppressed is set: that flag means we
        # are inside apply_fractional_steps' sync_from_device()/sync_to_device()
        # bracket (a CPU-only fractional operator is present). The GPU path
        # writes the *device* inlet cells, but the trailing sync_to_device()
        # (host->device) would then overwrite them with host data that never
        # received the inflow — silently dropping it. Fall through to the host
        # path so the batch sync_to_device() carries the inflow to the device.
        if (getattr(self.domain, 'multiprocessor_mode', 0) == 2
                and not getattr(self.domain, '_gpu_host_writes_suppressed', False)):
            if not self._gpu_initialized:
                self._init_gpu()
            if self._gpu_initialized:
                return self._call_gpu()

        timestep = self.domain.get_timestep()

        #print('Timestep', timestep)

        t = self.domain.get_time()


        # Need to run global command on all processors
        current_volume = self.inlet.get_total_water_volume()
        total_area = self.inlet.get_area()


        #print(current_volume)
        #print(total_area)

        assert current_volume >= 0.0

        Q1 = self.update_Q(t)
        Q2 = self.update_Q(t + timestep)


        #print Q1,Q2
        Q = 0.5*(Q1+Q2)
        volume = Q*timestep

        #print(volume)
        #print volume

        #print Q, volume

        # store last discharge
        self.applied_Q = Q

        #print(self.domain.fractional_step_volume_integral)

        u,v = self.inlet.get_velocities()

        # Tracer state before the water moves (see inlet_tracers.py).
        carry_tracers = self.tracers.active()
        if carry_tracers:
            tr_idx = self.inlet.triangle_indices
            tr_areas = self.inlet.get_areas()
            tr_h_old, tr_m_old = self.tracers.capture(tr_idx)
        requested_volume = volume

        # Distribute positive volume so as to obtain flat surface otherwise
        # just pull water off to have a uniform depth.
        if volume >= 0.0 :
            #print('volume>=0.0')
            self.inlet.set_stages_evenly(volume)
            self.domain.fractional_step_volume_integral+=volume
            self.total_requested_volume += volume

            if self.velocity is not None:
                depths = self.inlet.get_depths()
                self.inlet.set_xmoms(depths*self.velocity[0])
                self.inlet.set_ymoms(depths*self.velocity[1])
            else:
                depths = self.inlet.get_depths()

                self.inlet.set_xmoms(depths*u)
                self.inlet.set_ymoms(depths*v)

            if self.zero_velocity:
                self.inlet.set_xmoms(0.0)
                self.inlet.set_ymoms(0.0)

        elif current_volume + volume >= 0.0 :
            depth = (current_volume + volume)/total_area
            self.inlet.set_depths(depth)
            self.total_requested_volume += volume
            self.domain.fractional_step_volume_integral+=volume

            if self.velocity is not None:
                depths = self.inlet.get_depths()
                self.inlet.set_xmoms(depths*self.velocity[0])
                self.inlet.set_ymoms(depths*self.velocity[1])
            else:
                depths = self.inlet.get_depths()
                self.inlet.set_xmoms(depths*u)
                self.inlet.set_ymoms(depths*v)

            if self.zero_velocity:
                self.inlet.set_xmoms(0.0)
                self.inlet.set_ymoms(0.0)

        else: #extracting too much water!
            self.inlet.set_depths(0.0)
            self.total_requested_volume += volume
            volume = -current_volume
            self.applied_Q = -current_volume/timestep
            self.domain.fractional_step_volume_integral-=current_volume
            self.applied_Q = - current_volume/timestep
            self.inlet.set_xmoms(0.0)
            self.inlet.set_ymoms(0.0)


        self.total_applied_volume += volume

        if carry_tracers:
            self.tracers.apply(tr_idx, tr_areas, tr_h_old, tr_m_old,
                               requested_volume, t, timestep)


    def update_Q(self, t):
        """Allowing local modifications of Q
        """
        from anuga.fit_interpolate.interpolate import Modeltime_too_early, Modeltime_too_late

        if callable(self.Q):
            try:
                Q = self.Q(t)
            except Modeltime_too_early as e:
                Q = self.get_default(t,err_msg=e)
            except Modeltime_too_late as e:
                Q = self.get_default(t,err_msg=e)
        else:
            Q = self.Q

        return Q

    def statistics(self):


        message  = '=====================================\n'
        message += 'Inlet Operator: %s\n' % self.label
        message += '=====================================\n'

        message += 'Description\n'
        message += '%s' % self.description
        message += '\n'

        inlet = self.inlet

        message += '-------------------------------------\n'
        message +=  'Inlet\n'
        message += '-------------------------------------\n'

        message += 'inlet triangle indices and centres\n'
        message += '%s' % inlet.triangle_indices
        message += '\n'

        message += '%s' % self.domain.get_centroid_coordinates()[inlet.triangle_indices]
        message += '\n'

        message += 'region\n'
        message += '%s' % inlet
        message += '\n'

        message += '=====================================\n'

        return message


    def timestepping_statistics(self):

        message = '---------------------------\n'
        message += 'Inlet report for %s:\n' % self.label
        message += '--------------------------\n'
        message += 'Q [m^3/s]: %.2f\n' % self.applied_Q
        message += 'Total volume [m^3]: %.2f\n' % self.total_applied_volume

        return message

    def print_timestepping_statisitics(self):

        message = self.timestepping_statistics()
        print(message)


    def set_default(self, default=None):

        """ Either leave default as None or change it into a function"""

        if default is not None:
            # If it is a constant, make it a function
            if not callable(default):
                tmp = default
                default = lambda t: tmp

            # Check that default_rate is a function of one argument
            try:
                default(0.0)
            except TypeError:
                msg = "could not call default"
                raise Exception(msg)

        self.default = default
        self.default_invoked = False



    def get_default(self,t, err_msg=' '):
        """ Call get_default only if exception
        Modeltime_too_late(msg) has been raised
        """


        # Pass control to default rate function
        value = self.default(t)

        if self.default_invoked is False:
            # Issue warning the first time
            msg = ('%s\n'
                   'Instead I will use the default rate: %s\n'
                   'Note: Further warnings will be suppressed'
                   % (str(err_msg), str(self.default(t))))

            warnings.warn(msg)

            # FIXME (Ole): Replace this crude flag with
            # Python's ability to print warnings only once.
            # See http://docs.python.org/lib/warning-filter.html
            self.default_invoked = True

        return value

    def set_Q(self, Q):

        self.Q = Q

    def get_Q(self):

        return self.Q

    def get_inlet(self):

        return self.inlet

    def get_applied_Q(self):

        return self.applied_Q

    def get_total_applied_volume(self):

        return self.total_applied_volume
