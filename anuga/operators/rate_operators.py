"""
Rate operators (such as rain)

Constraints: See GPL license in the user guide
Version: 1.0 ($Revision: 7731 $)
"""
__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"



from anuga.config import indent
from anuga.config import MULTIPROCESSOR_OPENMP, MULTIPROCESSOR_GPU
import warnings
import numpy as num
import anuga.utilities.log as log
from anuga.utilities.function_utils import evaluate_temporal_function


from anuga import Quantity
from anuga.operators.base_operator import Operator
from anuga import Region

class Rate_operator(Operator):


    def __init__(self,
                 domain,
                 rate=0.0,
                 factor=1.0,
                 region=None,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 default_rate=0.0,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False,
                 monitor = False,
                 tracer_concentrations = None,
                 tracer_extraction = 'retain'):
        """Create a Rate_operator that adds water over a region at a specified rate.

        The applied rate is ``rate * factor`` in m/s (depth per second).
        For common use-cases prefer the factory constructors:

        * ``Rate_operator.rainfall(domain, rate_mm_hr)`` — rainfall in mm/hr
        * ``Rate_operator.inflow(domain, rate_m3_s)``    — inflow in m³/s

        Parameters
        ----------
        domain : anuga.Domain
            The simulation domain.
        rate : scalar, callable, Quantity, or ndarray
            Rate in m/s (after multiplication by *factor*).  May be a scalar,
            a function of ``t``, ``(x, y)``, or ``(x, y, t)``, a Quantity,
            a numpy array of shape ``(n_triangles,)``, or an xarray DataArray.
        factor : scalar or callable(t)
            Multiplier applied to *rate* before adding to stage.  Use to
            convert units (e.g. ``1/(1000*3600)`` for mm/hr → m/s).
        region : Region, optional
            Pre-built Region.  Cannot be combined with *polygon*, *center*,
            *radius*, or *indices*.
        indices : list, optional
            Triangle indices where the rate is applied.
        polygon : list of [x, y], optional
            Polygon bounding the application area.
        center : [x, y], optional
            Centre of a circular application area.
        radius : float, optional
            Radius of the circular application area.
        default_rate : scalar or callable(t), optional
            Rate to use outside the time interval of the rate function/xarray.
        description, label, logging, verbose, monitor :
            Passed to the base ``Operator``.
        tracer_concentrations : dict, optional
            Tracer name -> concentration (float or function of t) of the water
            a positive rate adds. Tracers not named come in at zero: clean rain.
        tracer_extraction : {'retain', 'carry'} or dict, optional
            What a negative rate does with each tracer: ``'retain'`` (default,
            the old behaviour) leaves it behind, as salt stays under
            evaporation; ``'carry'`` removes it with the water at the cell's
            concentration, as for pumping or drains. A dict gives one choice
            per tracer name, the rest retaining.
        """



        # --------------------------------------------------
        # Input validation
        # --------------------------------------------------
        if isinstance(region, Region) and any(
                arg is not None for arg in (polygon, center, radius, indices)):
            raise ValueError(
                'Rate_operator: cannot specify both a Region object and '
                'polygon/center/radius/indices — the Region already encodes location')

        if not (rate is None
                or isinstance(rate, (int, float, num.ndarray))
                or callable(rate)
                or hasattr(rate, 'centroid_values')):  # Quantity duck-type
            raise TypeError(
                'Rate_operator: rate must be a number, callable, Quantity, or '
                'numpy array; got %s' % type(rate).__name__)

        # Validate the tracer settings BEFORE Operator.__init__ registers the
        # operator with the domain: a constructor that raised after that would
        # leave a half-built operator in the fractional-step list.
        from anuga.structures.inlet_tracers import RateTracers
        tracers = RateTracers(domain, tracer_concentrations, tracer_extraction,
                              label=label if label else '')

        Operator.__init__(self, domain, description, label, logging, verbose)
        self.tracers = tracers
        if getattr(domain, 'use_active_set', False):
            log.warning('Rate operator attached while active-set stepping is '
                        'enabled: results stay correct (wetted cells activate '
                        'on the next step), but widespread rain activates the '
                        'whole mesh and the active-set speedup degrades.')

        #-----------------------------------------------------
        # Make sure region is actually an instance of a region
        # Otherwise create a new region based on the other
        # input arguments
        #-----------------------------------------------------
        if isinstance(region,Region):
            region.set_verbose(verbose)
            self.region = region

        else:
            self.region = Region(domain,
                        indices=indices,
                        polygon=polygon,
                        center=center,
                        radius=radius,
                        verbose=verbose)


        #------------------------------------------
        # Local variables
        #------------------------------------------
        self.indices = self.region.indices
        self.set_areas()
        self.set_full_indices()

        #--------------------------------
        # Setting up rate
        #--------------------------------
        self.rate_callable = False
        self.rate_spatial = False
        self.rate_xarray = False


        #-------------------------------
        # Check if rate is actually an xarray.
        # Need xarray package installed
        #-------------------------------
        try:
            import xarray
        except ImportError as e:
            log.debug('xarray not available, xarray rate inputs disabled: %s' % e)
            xarray = None
        if xarray is not None:
            if type(rate) is xarray.core.dataarray.DataArray:
                self.rate_xarray = True
                xa = rate
                rate = 0.0
                self._prepare_xarray_rate(xa)


        self.set_rate(rate)
        self.set_default_rate(default_rate)
        self.default_rate_invoked = False    # Flag

        #------------------------------
        # Setting up factor, can be a scalar
        # or a function of time.
        # Limitation: factor does not currently support time-series files
        # or arrays.  Only scalars and callables of the form f(t) are
        # accepted.  Extend set_factor() if file/array support is needed.
        #------------------------------

        self.factor_callable = False
        self.set_factor(factor)



        # ----------------
        # Mass tracking
        #-----------------
        self.monitor = monitor
        self.cumulative_influx = 0.0
        self.local_influx=0.0
        self.local_max = 0.0
        self.local_min = 0.0

        # ----------------
        # GPU support
        #-----------------
        self._gpu_op_id = None  # GPU operator ID (set on first GPU call)
        self._gpu_initialized = False
        self._gpu_rate_array_cache = None  # Cached rate array for GPU (avoids recreating every call)
        self._gpu_rate_min_cache = None    # Cached min value for statistics
        self._gpu_rate_max_cache = None    # Cached max value for statistics
        self._gpu_rate_changed = True      # Flag to indicate rate data needs to be transferred to GPU

    def _init_gpu(self):
        """Initialize GPU operator for this rate operator."""
        if self._gpu_initialized:
            return

        # Check if domain is in GPU mode
        if not hasattr(self.domain, 'multiprocessor_mode') or self.domain.multiprocessor_mode != MULTIPROCESSOR_GPU:
            return

        # Check if we have a GPU interface
        if not hasattr(self.domain, 'gpu_interface') or self.domain.gpu_interface is None:
            return

        gpu_interface = self.domain.gpu_interface
        if not hasattr(gpu_interface, 'gpu_dom') or gpu_interface.gpu_dom is None:
            return

        # Only support non-spatial, non-xarray rates for GPU
        # Spatial rates (x,y or x,y,t functions) and xarray rates need more complex handling
        if self.rate_spatial or self.rate_xarray:
            return

        # Supported rate types: scalar, t, quantity, centroid_array
        if self.rate_type not in ('scalar', 't', 'quantity', 'centroid_array'):
            return

        # Get indices - if None, apply to all elements
        if self.indices is None:
            indices = num.arange(self.domain.number_of_elements, dtype=num.intc)
            areas = self.domain.areas.copy()
        elif hasattr(self.indices, '__len__') and len(self.indices) == 0:
            # No local elements on this rank (e.g. rainfall polygon doesn't
            # overlap this rank's partition).  This operator is a no-op here.
            # Mark as GPU-initialized so _has_cpu_only_fractional_operators
            # doesn't trigger an unnecessary GPU<->CPU sync every timestep.
            # __call__ still returns immediately at the empty-indices guard.
            self._gpu_initialized = True
            return
        else:
            indices = num.asarray(self.indices, dtype=num.intc)
            areas = self.domain.areas[indices].copy()

        # Get full indices for mass tracking
        if self.full_indices is not None and len(self.full_indices) > 0:
            full_indices = num.asarray(self.full_indices, dtype=num.intc)
        else:
            full_indices = None

        # On any failure — kernel import/init error, or the MAX_RATE_OPERATORS
        # slot limit being exceeded — leave _gpu_initialized False and return so
        # __call__ falls through to the Python path rather than crashing.
        try:
            from anuga.shallow_water.sw_domain_gpu_ext import init_rate_operator
            self._gpu_op_id = init_rate_operator(
                gpu_interface.gpu_dom,
                indices,
                areas.astype(num.float64),
                full_indices
            )
        except Exception as e:
            warnings.warn(
                f"GPU rate operator init failed for "
                f"'{getattr(self, 'label', repr(self))}' ({e}); "
                f"falling back to the Python path.",
                stacklevel=2,
            )
            self._gpu_op_id = None
            return

        if self._gpu_op_id < 0:
            warnings.warn(
                f"Failed to register GPU rate operator "
                f"'{getattr(self, 'label', repr(self))}': slot limit exceeded "
                f"(MAX_RATE_OPERATORS=64). Falling back to the Python path; "
                f"reduce the number of Rate_operator instances or increase "
                f"MAX_RATE_OPERATORS in gpu_domain.h to use the GPU kernel.",
                stacklevel=2,
            )
            self._gpu_op_id = None
            return

        self._gpu_initialized = True

    def __call__(self):
        """Apply rate operator to the domain for one timestep.

        Adds water (or removes it when the rate is negative) to the stage
        quantity of each triangle selected by ``indices``, scaled by the
        current timestep and factor.

        - If ``indices`` is an empty list, no triangles are modified.
        - If ``indices`` is None, all triangles are modified.
        - Otherwise only the triangles at the given indices are modified.

        Returns
        -------
        None
            Modifies ``domain.quantities['stage']`` (and momentum quantities
            when the rate is negative) in place.
        """

        if self.indices is not None and len(self.indices) == 0:
            return

        # Check for GPU execution path. Skip it when _gpu_host_writes_suppressed
        # is set: that flag means we are inside apply_fractional_steps'
        # sync_from_device()/sync_to_device() bracket (a CPU-only fractional
        # operator is present). The GPU path writes the *device* stage/momentum,
        # but the trailing sync_to_device() (host->device) would overwrite them
        # with host data that never received the rate — silently dropping it.
        # Fall through to the host path so the batch sync carries the change.
        if (hasattr(self.domain, 'multiprocessor_mode') and
            self.domain.multiprocessor_mode == MULTIPROCESSOR_GPU and
            not getattr(self.domain, '_gpu_host_writes_suppressed', False) and
            not self.rate_spatial and
            not self.rate_xarray):

            # Lazy initialization of GPU operator
            if not self._gpu_initialized:
                self._init_gpu()

            if self._gpu_initialized and self._gpu_op_id is not None and self._gpu_op_id >= 0:
                t = self.domain.get_time()
                timestep = self.domain.get_timestep()
                factor = self.get_factor(t)
                # Tracers ride with the water (inlet_tracers.RateTracers).
                gpu_tracers = self.tracers.needed()
                if gpu_tracers:
                    from anuga.shallow_water.sw_domain_gpu_ext import rate_operator_tracers_capture_gpu
                    rate_operator_tracers_capture_gpu(self.domain.gpu_interface.gpu_dom, self._gpu_op_id)

                # DEBUG: Confirm GPU path taken
                #print(f"DEBUG Rate_operator: GPU path, op_id={self._gpu_op_id}, rate_type={self.rate_type}, t={t}, timestep={timestep}, factor={factor}")
                # DEBUG: Check rate quantity info
                #if self.rate_type == 'quantity':
                #    print(f"DEBUG Rate_operator: rate quantity name={self.rate.name}, rate object id={id(self.rate)}")

                if self.rate_type == 'quantity':
                    # Quantity type - use array-based GPU kernel
                    from anuga.shallow_water.sw_domain_gpu_ext import apply_rate_operator_array_gpu
                    # Use cached rate array if available (avoids expensive array copy every RK2 step)
                    if self._gpu_rate_array_cache is None:
                        self._gpu_rate_array_cache = num.ascontiguousarray(self.rate.centroid_values, dtype=num.float64)
                        # Cache min/max too (avoids iterating 5M elements every call)
                        if self.indices is None:
                            self._gpu_rate_min_cache = self._gpu_rate_array_cache.min()
                            self._gpu_rate_max_cache = self._gpu_rate_array_cache.max()
                        else:
                            self._gpu_rate_min_cache = self._gpu_rate_array_cache[self.indices].min()
                            self._gpu_rate_max_cache = self._gpu_rate_array_cache[self.indices].max()
                        self._gpu_rate_changed = True  # New cache, needs transfer
                    rate_array = self._gpu_rate_array_cache
                    # rate_array is full domain size, use_indices_into_rate=1
                    # Pass rate_changed flag to avoid unnecessary H2D transfer
                    rate_changed = 1 if self._gpu_rate_changed else 0
                    self.local_influx = apply_rate_operator_array_gpu(
                        self.domain.gpu_interface.gpu_dom,
                        self._gpu_op_id,
                        rate_array,
                        1,  # use_indices_into_rate
                        rate_changed,
                        float(factor),
                        float(timestep)
                    )
                    self._gpu_rate_changed = False  # Data transferred, don't repeat
                    # Use cached min/max for statistics
                    self.local_max = self._gpu_rate_max_cache * factor
                    self.local_min = self._gpu_rate_min_cache * factor

                elif self.rate_type == 'centroid_array':
                    # Centroid array type - use array-based GPU kernel
                    from anuga.shallow_water.sw_domain_gpu_ext import apply_rate_operator_array_gpu
                    # Use cached rate array if available
                    if self._gpu_rate_array_cache is None:
                        self._gpu_rate_array_cache = num.ascontiguousarray(self.rate, dtype=num.float64)
                        # Cache min/max too
                        if self.indices is None:
                            self._gpu_rate_min_cache = self._gpu_rate_array_cache.min()
                            self._gpu_rate_max_cache = self._gpu_rate_array_cache.max()
                        else:
                            self._gpu_rate_min_cache = self._gpu_rate_array_cache[self.indices].min()
                            self._gpu_rate_max_cache = self._gpu_rate_array_cache[self.indices].max()
                        self._gpu_rate_changed = True  # New cache, needs transfer
                    rate_array = self._gpu_rate_array_cache
                    # rate_array is full domain size, use_indices_into_rate=1
                    # Pass rate_changed flag to avoid unnecessary H2D transfer
                    rate_changed = 1 if self._gpu_rate_changed else 0
                    self.local_influx = apply_rate_operator_array_gpu(
                        self.domain.gpu_interface.gpu_dom,
                        self._gpu_op_id,
                        rate_array,
                        1,  # use_indices_into_rate
                        rate_changed,
                        float(factor),
                        float(timestep)
                    )
                    self._gpu_rate_changed = False  # Data transferred, don't repeat
                    # Use cached min/max for statistics
                    self.local_max = self._gpu_rate_max_cache * factor
                    self.local_min = self._gpu_rate_min_cache * factor

                else:
                    # Scalar or time-dependent rate - use scalar GPU kernel
                    from anuga.shallow_water.sw_domain_gpu_ext import apply_rate_operator_gpu
                    rate = self.get_non_spatial_rate(t)
                    # Handle case where rate is an array (from file_function)
                    try:
                        rate_scalar = float(rate)
                    except (TypeError, ValueError):
                        rate_scalar = float(rate[0])
                    self.local_influx = apply_rate_operator_gpu(
                        self.domain.gpu_interface.gpu_dom,
                        self._gpu_op_id,
                        rate_scalar,
                        float(factor),
                        float(timestep)
                    )
                    # Estimate min/max rate for statistics
                    self.local_max = rate_scalar * factor if rate_scalar >= 0 else 0.0
                    self.local_min = rate_scalar * factor if rate_scalar < 0 else 0.0

                # Update tracking
                if gpu_tracers:
                    from anuga.shallow_water.sw_domain_gpu_ext import rate_operator_tracers_apply_gpu
                    c_in, carry = self.tracers.settings(t, timestep)
                    rate_operator_tracers_apply_gpu(self.domain.gpu_interface.gpu_dom,
                                                    self._gpu_op_id, c_in, carry)
                self.cumulative_influx += self.local_influx
                self.domain.fractional_step_volume_integral += self.local_influx

                return

        # Fall back to CPU path
        if self.rate_xarray:
            # setup centroid_array from xarray corresponding to current time
            self._update_Q_xarray()

        t = self.domain.get_time()
        timestep = self.domain.get_timestep()
        factor = self.get_factor()
        indices = self.indices


        if self.rate_spatial:
            if indices is None:
                x = self.coord_c[:,0]
                y = self.coord_c[:,1]
            else:
                x = self.coord_c[indices,0]
                y = self.coord_c[indices,1]

            rate = self.get_spatial_rate(x,y,t)
        elif self.rate_type == 'quantity':
            if indices is None:
                rate  = self.rate.centroid_values
            else:
                rate = self.rate.centroid_values[indices]
        elif self.rate_type == 'centroid_array':
            if indices is None:
                rate  = self.rate
            else:
                rate = self.rate[indices]
        else:
            rate = self.get_non_spatial_rate(t)


        factor = self.get_factor(t)

        # We need to adjust the momentums if rate < 0 since otherwise
        # the xmom and ymom stay the same but height -> 0 which leads to xvel, yvel -> infty


        fid = self.full_indices
        host_tracers = self.tracers.needed()
        if host_tracers:
            tr_h_old = self.tracers.capture(indices)
        if num.all(rate >= 0.0):
            # Record the local flux for mass conservation tracking
            if indices is None:
                local_rates = factor*timestep*rate
                self.local_influx = (local_rates*self.areas)[fid].sum()

                self.stage_c[:] = self.stage_c[:] + local_rates
            else:
                local_rates = factor*timestep*rate
                self.local_influx=(local_rates*self.areas)[fid].sum()

                self.stage_c[indices] = self.stage_c[indices] + local_rates
        else: # Be more careful if rate < 0
            if indices is None:
                #self.local_influx=(num.minimum(factor*timestep*rate, self.stage_c[:]-self.elev_c[:])*self.areas)[fid].sum()
                #self.stage_c[:] = num.maximum(self.stage_c  \
                #       + factor*rate*timestep, self.elev_c )
                self.height_c[:] = self.stage_c[:] - self.elev_c[:]
                local_rates = num.maximum(factor*timestep*rate, -self.height_c)
                local_factors = num.where(local_rates < 0.0, (local_rates+self.height_c)/(self.height_c+1.0e-10), 1.0)

                #print(local_factors, local_rates)
                self.local_influx = (local_rates*self.areas)[fid].sum()

                self.stage_c[:] = self.stage_c + local_rates
                self.xmom_c[:] = self.xmom_c[:]*local_factors
                self.ymom_c[:] = self.ymom_c[:]*local_factors
            else:
                #self.local_influx=(num.minimum(factor*timestep*rate, self.stage_c[indices]-self.elev_c[indices])*self.areas)[fid].sum()
                #self.stage_c[indices] = num.maximum(self.stage_c[indices] \
                #       + factor*rate*timestep, self.elev_c[indices])

                #local_rates = num.maximum(factor*timestep*rate, self.elev_c[indices]-self.stage_c[indices])

                heights = self.stage_c[indices] - self.elev_c[indices]
                local_rates = num.maximum(factor*timestep*rate, -heights)
                local_factors = num.where(local_rates < 0.0, (local_rates+heights)/(heights+1.0e-10), 1.0)

                #print(local_factors, local_rates, fid)

                self.local_influx = (local_rates*self.areas)[fid].sum()

                self.stage_c[indices] = self.stage_c[indices] + local_rates
                self.xmom_c[indices] = self.xmom_c[indices]*local_factors
                self.ymom_c[indices] = self.ymom_c[indices]*local_factors


        if host_tracers:
            self.tracers.apply(indices, tr_h_old, t, timestep)

        try:
            arr = local_rates[fid]
            if isinstance(arr, num.ndarray) and arr.size == 0:
                self.local_max = 0.0
                self.local_min = 0.0
            else:
                self.local_max = arr.max()/timestep if timestep != 0.0 else 0.0
                self.local_min = arr.min()/timestep if timestep != 0.0 else 0.0
        except (TypeError, IndexError):
            self.local_max = local_rates/timestep if timestep != 0.0 else 0.0
            self.local_min = local_rates/timestep if timestep != 0.0 else 0.0

        # print(self.local_min, self.local_max)

        self.cumulative_influx += self.local_influx

        # Update mass inflows from fractional steps
        self.domain.fractional_step_volume_integral+=self.local_influx

        if self.monitor:
            log.info('Local Flux at time %.2f = %f'
                         % (self.domain.get_time(), self.local_influx))



        return

    def get_non_spatial_rate(self, t=None):
        """Provide a rate to calculate added volume
        """

        if t is None:
            t = self.get_time()

        assert not self.rate_spatial

        rate = evaluate_temporal_function(self.rate, t,
                                          default_right_value=self.default_rate,
                                          default_left_value=self.default_rate)

        if rate is None:
            msg = ('Attribute rate must be specified in '+self.__name__+
                   ' before attempting to call it')
            raise Exception(msg)

        return rate

    def get_spatial_rate(self, x=None, y=None, t=None):
        """Provide a rate to calculate added volume
        only call if self.rate_spatial = True
        """

        assert self.rate_spatial

        if t is None:
            t = self.get_time()

        if x is None:
            assert y is None
            if self.indices is None:
                x = self.coord_c[:,0]
                y = self.coord_c[:,1]
            else:
                x = self.coord_c[self.indices,0]
                y = self.coord_c[self.indices,1]

        assert x is not None
        assert y is not None

        assert isinstance(t, (int, float))
        assert len(x) == len(y)

        #print xy
        #print t

        #print self.rate_type, self.rate_type == 'x,y,t'
        if self.rate_type == 'x,y,t':
            rate = self.rate(x,y,t)
        else:
            rate = self.rate(x,y)

        return rate


    def set_rate(self, rate):
        """Set rate. Can change rate while running


        Can be a scalar, numpy array, or a function of t or x,y or x,y,t or a quantity
        """

        # Test if rate is a quantity
        if isinstance(rate, Quantity):
            self.rate_type = 'quantity'
        elif isinstance(rate, num.ndarray):
            rate_shape = rate.shape
            msg =  f"The shape {rate_shape} of the input rate "
            msg += f"should match (number of triangles,) i.e. ({self.domain.number_of_triangles},)"
            assert rate_shape == (self.domain.number_of_triangles,) \
                or rate_shape == (self.domain.number_of_triangles, 1), msg
            self.rate_type = 'centroid_array'
            rate = rate.reshape((-1,))
        else:
            # Possible types are 'scalar', 't', 'x,y' and 'x,y,t'
            from anuga.utilities.function_utils import determine_function_type
            self.rate_type = determine_function_type(rate)


        self.rate = rate

        # Invalidate GPU caches (will be recreated on next GPU call)
        # Use hasattr since set_rate() can be called from __init__ before cache is initialized
        if hasattr(self, '_gpu_rate_array_cache'):
            self._gpu_rate_array_cache = None
            self._gpu_rate_min_cache = None
            self._gpu_rate_max_cache = None
            self._gpu_rate_changed = True  # Signal that rate data needs to be transferred to GPU

        if self.rate_type == 'scalar':
            self.rate_callable = False
            self.rate_spatial = False
        elif self.rate_type == 'quantity':
            self.rate_callable = False
            self.rate_spatial = False
        elif self.rate_type == 'centroid_array':
            self.rate_callable = False
            self.rate_spatial = False
        elif self.rate_type == 't':
            self.rate_callable = True
            self.rate_spatial = False
        else:
            self.rate_callable = True
            self.rate_spatial = True


    def set_factor(self, factor):
        """Set factor. Can change factor while running


        Can be a scalar, a function of t, or an n by 2 numpy array defining a time sequence
        """

        if isinstance(factor, num.ndarray):
            factor_shape = factor.shape
            msg =  f"The shape {factor_shape} of the input factor "
            msg += "should be (2,n) so that a time function can be constructed"
            assert factor_shape[0] == 2, msg
            self.factor_type = 'time_sequence'
            from scipy.interpolate import interp1d
            factor = interp1d(factor[0,:], factor[1,:], kind='zero', bounds_error=False,  fill_value = (0.0, 0.0))


        from anuga.utilities.function_utils import determine_function_type
        self.factor_type = determine_function_type(factor)


        self.factor = factor

        if self.factor_type == 'scalar':
            self.factor_callable = False
        elif self.factor_type == 't':
            self.factor_callable = True
        else:
            msg = f'factor must be a scalar or a function of t. It was determined to be a function of {self.factor_type}'
            raise Exception(msg)

    def get_factor(self, t=None):
        """Provide a factor to calculate added volume
        """

        if t is None:
            t = self.get_time()

        assert isinstance(t, (int, float))


        if self.factor_type == 'scalar':
            factor = self.factor
        else:
            factor = self.factor(t)

        return factor

    def set_areas(self):

        if self.indices is None:
            self.areas = self.domain.areas
            return

        if self.indices is not None and len(self.indices) == 0:
            self.areas = num.array([])
            return

        self.areas = self.domain.areas[self.indices]

    def set_full_indices(self):

        if self.indices is None:
            self.full_indices = num.where(self.domain.tri_full_flag ==1)[0]
            return

        if self.indices is not None and len(self.indices) == 0:
            self.full_indices = num.array([], dtype=int)
            return

        self.full_indices = num.where(self.domain.tri_full_flag[self.indices] == 1)[0]

    def get_Q(self, full_only=True):
        """ Calculate current overall discharge
        """

        # FIXME SR: this does not take into account the zeroing of large negative rates

        factor = self.get_factor()

        if full_only:
            if self.rate_spatial:
                rate = self.get_spatial_rate() # rate is an array
                fid = self.full_indices
                return num.sum(self.areas[fid]*rate[fid])*factor
            elif self.rate_type == 'quantity':
                rate = self.rate.centroid_values # rate is a quantity
                fid = self.full_indices
                return num.sum(self.areas[fid]*rate[fid])*factor
            elif self.rate_type == 'centroid_array':
                rate = self.rate # rate is already a centroid sized array
                fid = self.full_indices
                return num.sum(self.areas[fid]*rate[fid])*factor
            else:
                rate = self.get_non_spatial_rate() # rate is a scalar
                fid = self.full_indices
                return num.sum(self.areas[fid]*rate)*factor
        else:
            if self.rate_spatial:
                rate = self.get_spatial_rate() # rate is an array
                return num.sum(self.areas*rate)*factor
            elif self.rate_type == 'quantity':
                rate = self.rate.centroid_values # rate is a quantity
                return num.sum(self.areas*rate)*factor
            elif self.rate_type == 'centroid_array':
                rate = self.rate # rate is already a centroid sized array
                return num.sum(self.areas*rate)*factor
            else:
                rate = self.get_non_spatial_rate() # rate is a scalar
                return num.sum(self.areas*rate)*factor

    def set_default_rate(self, default_rate):
        """
        Check and store default_rate
        """
        msg = ('Default_rate must be either None '
               'a scalar, or a function of time.\nI got %s.' % str(default_rate))
        assert (default_rate is None or
                isinstance(default_rate, (int, float)) or
                callable(default_rate)), msg


        #------------------------------------------
        # Allow longer than data
        #------------------------------------------
        if default_rate is not None:
            # If it is a constant, make it a function
            if not callable(default_rate):
                tmp = default_rate
                default_rate = lambda t: tmp

            # Check that default_rate is a function of one argument
            try:
                default_rate(0.0)
            except TypeError:
                raise Exception(msg)

        self.default_rate = default_rate

    # ------------------------------------------------------------------
    # Factory constructors
    # ------------------------------------------------------------------

    @classmethod
    def rainfall(cls, domain, rate, polygon=None, region=None,
                 center=None, radius=None, indices=None,
                 default_rate=0.0, label=None, description=None,
                 logging=False, verbose=False, monitor=False,
                 tracer_concentrations=None, tracer_extraction='retain'):
        """Create a Rate_operator for rainfall.

        Parameters
        ----------
        domain : anuga.Domain
        rate : scalar, callable(t), or array
            Rainfall intensity in **mm/hr**.  All other rate forms accepted
            by ``Rate_operator`` (callables, arrays) are also supported and
            are interpreted as mm/hr.
        polygon, region, center, radius, indices :
            Location arguments — same as ``Rate_operator.__init__``.

        Returns
        -------
        Rate_operator
            Operator with ``factor = 1 / (1000 * 3600)`` so that mm/hr is
            converted to m/s automatically.

        Examples
        --------
        >>> rain = Rate_operator.rainfall(domain, rate=10.0)          # 10 mm/hr
        >>> rain = Rate_operator.rainfall(domain, rate=lambda t: 5.0) # time-varying
        """
        MM_HR_TO_M_S = 1.0 / (1000.0 * 3600.0)
        return cls(domain, rate=rate, factor=MM_HR_TO_M_S,
                   polygon=polygon, region=region,
                   center=center, radius=radius, indices=indices,
                   default_rate=default_rate, label=label,
                   description=description, logging=logging,
                   verbose=verbose, monitor=monitor,
                   tracer_concentrations=tracer_concentrations,
                   tracer_extraction=tracer_extraction)

    @classmethod
    def inflow(cls, domain, rate, polygon=None, region=None,
               center=None, radius=None, indices=None,
               default_rate=0.0, label=None, description=None,
               logging=False, verbose=False, monitor=False,
               tracer_concentrations=None, tracer_extraction='retain'):
        """Create a Rate_operator for a volumetric inflow.

        Parameters
        ----------
        domain : anuga.Domain
        rate : scalar or callable(t)
            Volumetric flow rate in **m³/s**.  The operator divides by the
            total region area so that the net inflow equals *rate* m³/s.

        Returns
        -------
        Rate_operator

        Raises
        ------
        ValueError
            If the specified region has zero area.

        Examples
        --------
        >>> op = Rate_operator.inflow(domain, rate=0.5, polygon=poly)
        >>> op = Rate_operator.inflow(domain, rate=lambda t: 0.1*t)
        """
        # Build the operator with rate=0 first to resolve the region/area.
        op = cls(domain, rate=0.0, factor=1.0,
                 polygon=polygon, region=region,
                 center=center, radius=radius, indices=indices,
                 default_rate=default_rate, label=label,
                 description=description, logging=logging,
                 verbose=verbose, monitor=monitor,
                 tracer_concentrations=tracer_concentrations,
                 tracer_extraction=tracer_extraction)

        total_area = float(op.areas.sum()) if op.areas is not None and len(op.areas) > 0 else 0.0
        if total_area <= 0.0:
            raise ValueError(
                'Rate_operator.inflow: region has zero area — '
                'check that the polygon/center/radius overlaps the domain')

        # Convert m³/s → m/s by dividing by region area.
        # Use a closure factory so the wrapper has exactly one argument (t),
        # which determine_function_type correctly classifies as type 't'.
        if callable(rate):
            def _make_scaled(fn, area):
                def scaled(t):
                    return fn(t) / area
                return scaled
            op.set_rate(_make_scaled(rate, total_area))
        else:
            op.set_rate(rate / total_area)

        return op

    def _prepare_xarray_rate(self, xa):

        import numpy as np
        import pandas

        # to speed up parallel code it helps to load the xarray
        self.xa = xa.load()

        self.xa['time'] = pandas.to_datetime(xa['time'], utc=True)

        # these are absolute coord (since we haven't implemented offsets)
        # Convert to relative coords to domain xllcorner and yllcorner

        xllcorner = self.domain.geo_reference.xllcorner
        yllcorner = self.domain.geo_reference.yllcorner
        self.xy = np.array([self.xa['eastings']-xllcorner, self.xa['northings']-yllcorner]).T


        # Determine data timestep from xarray. We assume the timestep is constant, so just test first 2 timeslices.
        try:
            data_dt = (self.xa['time'][1].values.astype('int64')-self.xa['time'][0].values.astype('int64'))/1.0e9
            self.domain.set_evolve_max_timestep(min(data_dt, self.domain.get_evolve_max_timestep()))
        except Exception:  # if we can't determine the timestep probably means there is just one timeslice so just
            pass

        from scipy.spatial import KDTree
        tree = KDTree(self.xy)
        if self.verbose:
            print(tree.size, self.xy.shape)

        #dd, ii = tree.query(self.domain.centroid_coordinates)
        dd, ii = tree.query(self.domain.get_centroid_coordinates(absolute=True))

        self.ii = ii

        self.previous_Q_ref_time = None
        self.previous_Q_numpy = None


    def _update_Q_xarray(self):

        import pandas
        current_utc_datetime64 = pandas.to_datetime(self.domain.get_datetime()).tz_convert('UTC')#.replace(tzinfo=None)

        if self.verbose:
            print(f"{self.domain.get_time()} {self.domain.get_datetime()}")
            print(f"UTC time {current_utc_datetime64} type {type(current_utc_datetime64)} ")
            print(self.xa.sel(time=current_utc_datetime64, method='nearest'))

        try:
            Q_ref = self.xa.sel(time=current_utc_datetime64, method="ffill", tolerance='5m')

            Q_ref_time = Q_ref['time'].values

            if self.verbose:
                print(f"UTC time {current_utc_datetime64} Q_ref time {Q_ref_time} {Q_ref_time == self.previous_Q_ref_time} ")

            optimize = True
            if optimize:
                if Q_ref_time == self.previous_Q_ref_time :
                    Q_numpy = self.previous_Q_numpy
                else:
                    Q_numpy = Q_ref[self.ii].to_numpy()
                    self.previous_Q_numpy = Q_numpy
                    self.previous_Q_ref_time = Q_ref_time
            else:
                Q_numpy = Q_ref[self.ii].to_numpy()

        except Exception:
            Q_numpy = self.default_rate
            if self.verbose:
                print(f"UTC time {current_utc_datetime64} Using default rate Q = {Q_numpy(self.get_time())}")

        self.set_rate(rate=Q_numpy)

    def parallel_safe(self):
        """Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True

    def statistics(self):

        message = 'You need to implement operator statistics for your operator'
        return message


    def timestepping_statistics(self):

        # retrieve data from last __call__ call
        message  = indent + self.label + f': At time {self.domain.get_time()} Min rate = {self.local_min:.2e} m/s, Max rate = {self.local_max:.2e} m/s, Total Q = {self.cumulative_influx:.2e} m^3'


        # if self.rate_spatial:
        #     rate = self.get_spatial_rate()
        #     try:
        #         min_rate = num.min(rate)
        #     except ValueError:
        #         min_rate = 0.0
        #     try:
        #         max_rate = num.max(rate)
        #     except ValueError:
        #         max_rate = 0.0

        #     Q = self.get_Q()
        #     message  = indent + self.label + ': Min rate = %g m/s, Max rate = %g m/s, Total Q = %g m^3/s'% (min_rate,max_rate, Q)

        # elif self.rate_type == 'quantity':
        #     rate = self.get_non_spatial_rate() # return quantity
        #     min_rate = rate.get_minimum_value()
        #     max_rate = rate.get_maximum_value()
        #     Q = self.get_Q()
        #     message  = indent + self.label + ': Min rate = %g m/s, Max rate = %g m/s, Total Q = %g m^3/s'% (min_rate,max_rate, Q)

        # elif self.rate_type == 'centroid_array':
        #     rate = self.get_non_spatial_rate() # return centroid_array
        #     min_rate = rate.min()
        #     max_rate = rate.max()
        #     Q = self.get_Q()
        #     message  = indent + self.label + ': Min rate = %g m/s, Max rate = %g m/s, Total Q = %g m^3/s'% (min_rate,max_rate, Q)

        # else:
        #     rate = self.get_non_spatial_rate()
        #     Q = self.get_Q()
        #     message  = indent + self.label + ': Rate = %g m/s, Total Q = %g m^3/s' % (rate, Q)


        #print(message)

        return message

# ===============================================================================
# Specific Rate Operators for circular region.
# ===============================================================================
class Circular_rate_operator(Rate_operator):
    """
    Add water at certain rate (ms^{-1} = vol/Area/sec) over a
    circular region

    rate can be a function of time.

    Other units can be used by using the factor argument.

    """

    def __init__(self, domain,
                 rate=0.0,
                 factor=1.0,
                 center=None,
                 radius=None,
                 default_rate=None,
                 verbose=False):


        Rate_operator.__init__(self,
                               domain,
                               rate=rate,
                               factor=factor,
                               center=center,
                               radius=radius,
                               default_rate=default_rate,
                               verbose=verbose)



#===============================================================================
# Specific Rate Operators for polygonal region.
#===============================================================================
class Polygonal_rate_operator(Rate_operator):
    """
    Add water at certain rate (ms^{-1} = vol/Area/sec) over a
    polygonal region

    rate can be a function of time.

    Other units can be used by using the factor argument.

    """

    def __init__(self, domain,
                 rate=0.0,
                 factor=1.0,
                 polygon=None,
                 default_rate=None,
                 verbose=False):


        Rate_operator.__init__(self,
                               domain,
                               rate=rate,
                               factor=factor,
                               polygon=polygon,
                               default_rate=default_rate,
                               verbose=verbose)
