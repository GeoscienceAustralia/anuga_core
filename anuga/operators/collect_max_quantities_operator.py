"""
Collect max stage quantity info
"""

__author__ = "gareth"
__date__ = "$01/08/2014 4:46:39 PM$"


import numpy as num
from anuga.operators.base_operator import Operator
from anuga import Quantity
from anuga.parallel import myid
from anuga.config import max_float


class Collect_max_quantities_operator(Operator):
    """
    Operator to collect the max stage, depth, speed, and momentum magnitude
    (speed * depth, i.e. ||(uh, vh)||) during a run.

    Maxima are updated every update_frequency timesteps [any integer >= 1],
    after t exceeds collection_start_time.

    NOTE: the INITIAL condition is not sampled. Collection happens on the
    fractional-step call, which first runs after the first update, so a cell
    whose maximum is its starting value -- the reservoir behind a dam break,
    the source cells of a plume that then washes out -- reports the first
    post-step value instead. Long-standing behaviour for stage/depth/speed/uh;
    the tracer maxima follow it so that every max_* variable in an sww means
    the same thing.

    If store_to_sww=True the running maxima are written to the SWW file every
    yield step as centroid quantities max_stage_c, max_depth_c, max_speed_c,
    max_uh_c.  This avoids the need to compute maxima by scanning all timeslices
    of the SWW stage output after the run.

    Every registered tracer -- and so every suspended sediment class, a
    sediment class being a tracer -- is tracked the same way and stored as
    max_<name>_c.  Register the tracers before creating this operator: the
    names are fixed here, and an sww variable cannot appear after the file is
    created.

    Velocity is zeroed below velocity_zero_height (defaults to
    minimum_allowed_height) to suppress spurious spikes in nearly-dry cells.

    When the domain is in GPU mode (MULTIPROCESSOR_GPU) the kernel runs
    entirely on the device; a host<->device transfer for the four max arrays
    only occurs at yield steps when store_to_sww=True, or on an explicit call
    to export_max_quantities_to_csv().
    """

    def __init__(self,
                 domain,
                 update_frequency=1,
                 collection_start_time=0.,
                 velocity_zero_height=None,
                 store_to_sww=False,
                 description=None,
                 label=None,
                 logging=False,
                 verbose=False):

        Operator.__init__(self, domain, description, label, logging, verbose)

        self.domain = domain

        n = len(domain.centroid_coordinates[:, 0])
        self.max_stage = num.full(n, -max_float)
        self.max_depth = num.zeros(n)
        self.max_speed = num.zeros(n)
        self.max_uh    = num.zeros(n)

        # Running max of each tracer's concentration -- and therefore of each
        # suspended sediment class, which is a tracer. Zero rather than
        # -max_float: a concentration is non-negative, so zero is the floor
        # and not a sentinel.
        #
        # The names are fixed HERE. add_tracer reallocates the tracer block at
        # a new width, and an sww variable cannot appear after the file is
        # created, so a tracer registered later cannot be collected; __call__
        # says so rather than quietly dropping it.
        self._tracer_names = list(getattr(domain, 'get_tracer_names', list)())
        self.max_tracer = num.zeros((len(self._tracer_names), n))
        self._minimum_allowed_height = domain.minimum_allowed_height

        self.xy    = domain.centroid_coordinates
        self.stage = domain.quantities['stage']
        self.elev  = domain.quantities['elevation']
        self.xmom  = domain.quantities['xmomentum']
        self.ymom  = domain.quantities['ymomentum']

        self.counter = 0
        assert update_frequency > 0, 'Update frequency must be >= 1'
        self.update_frequency = update_frequency
        self.collection_start_time = collection_start_time

        if velocity_zero_height is not None:
            self.velocity_zero_height = velocity_zero_height
        else:
            self.velocity_zero_height = domain.minimum_allowed_height

        self.store_to_sww = store_to_sww
        if store_to_sww:
            names = ['max_stage', 'max_depth', 'max_speed', 'max_uh']
            names += ['max_' + t for t in self._tracer_names]
            for qname in names:
                Quantity(domain, name=qname, register=True, qty_type='centroid_only')
            # flag=4: centroid-only static variable overwritten each yield step
            domain.quantities_to_be_stored.update({qname: 4 for qname in names})

        self._gpu_initialized = False

    # ------------------------------------------------------------------
    # GPU initialisation (lazy, called on first __call__ in GPU mode)
    # ------------------------------------------------------------------

    def _init_gpu(self):
        if self._gpu_initialized:
            return

        domain = self.domain
        try:
            from anuga.config import MULTIPROCESSOR_GPU
        except ImportError:
            return

        if not hasattr(domain, 'multiprocessor_mode') or domain.multiprocessor_mode != MULTIPROCESSOR_GPU:
            return
        if not hasattr(domain, 'gpu_interface') or domain.gpu_interface is None:
            return
        gpu_interface = domain.gpu_interface
        if not hasattr(gpu_interface, 'gpu_dom') or gpu_interface.gpu_dom is None:
            return

        from anuga.shallow_water.sw_domain_gpu_ext import init_max_quantities_gpu
        n   = len(domain.centroid_coordinates[:, 0])
        ret = init_max_quantities_gpu(gpu_interface.gpu_dom, n, self.velocity_zero_height)
        if ret != 0:
            raise RuntimeError(
                "gpu_max_quantities_init failed — check device memory availability"
            )
        self._gpu_initialized = True

    # ------------------------------------------------------------------

    def __call__(self):
        """
        Update max quantities every update_frequency timesteps once
        t > collection_start_time.
        """

        t = self.domain.get_time()

        if t >= self.collection_start_time:
            self.counter += 1

            if self.counter == self.update_frequency:
                self.counter = 0

                # ---- GPU path ----------------------------------------
                domain = self.domain
                try:
                    from anuga.config import MULTIPROCESSOR_GPU
                    _gpu_mode = (
                        hasattr(domain, 'multiprocessor_mode') and
                        domain.multiprocessor_mode == MULTIPROCESSOR_GPU and
                        hasattr(domain, 'gpu_interface') and
                        domain.gpu_interface is not None
                    )
                except ImportError:
                    _gpu_mode = False

                if _gpu_mode:
                    self._init_gpu()
                    if self._gpu_initialized:
                        from anuga.shallow_water.sw_domain_gpu_ext import (
                            update_max_quantities_gpu, get_max_quantities_gpu)

                        update_max_quantities_gpu(domain.gpu_interface.gpu_dom)

                        if self.store_to_sww:
                            # D2H sync: pull device maxima to host numpy arrays,
                            # then push to domain quantities for the SWW writer.
                            get_max_quantities_gpu(
                                domain.gpu_interface.gpu_dom,
                                self.max_stage, self.max_depth,
                                self.max_speed, self.max_uh)
                            d = domain
                            d.quantities['max_stage'].centroid_values[:] = self.max_stage
                            d.quantities['max_depth'].centroid_values[:] = self.max_depth
                            d.quantities['max_speed'].centroid_values[:] = self.max_speed
                            d.quantities['max_uh'].centroid_values[:]    = self.max_uh
                            self._store_tracer_maxima_from_device(domain)
                        return

                # ---- CPU / NumPy path --------------------------------
                stage_c = self.stage.centroid_values
                elev_c  = self.elev.centroid_values
                xmom_c  = self.xmom.centroid_values
                ymom_c  = self.ymom.centroid_values

                self.max_stage = num.maximum(self.max_stage, stage_c)

                mom_norm = num.hypot(xmom_c, ymom_c)
                self.max_uh = num.maximum(self.max_uh, mom_norm)

                depth = num.maximum(stage_c - elev_c, 0.0)
                self.max_depth = num.maximum(self.max_depth, depth)

                vel = (mom_norm / num.maximum(depth, self.velocity_zero_height)) \
                      * (depth > self.velocity_zero_height)
                self.max_speed = num.maximum(self.max_speed, vel)

                if self._tracer_names:
                    self._check_tracers_unchanged(domain)
                    # DERIVE c from the conserved m rather than reading
                    # domain.tracer_centroid_values. c is only recomputed
                    # during extrapolation, which runs at the START of a step,
                    # so by the time a fractional-step operator is called it is
                    # one step behind the m that update_conserved_quantities
                    # has just written. Reading it would quietly miss every
                    # new maximum reached on the final step -- max_stage does
                    # not suffer this because stage IS a conserved variable.
                    #
                    # Same rule the kernel uses: dry cells carry no
                    # concentration.
                    inv_h = num.where(depth > self._minimum_allowed_height,
                                      1.0 / num.maximum(depth, 1e-300), 0.0)
                    num.maximum(self.max_tracer,
                                domain.tracer_conserved_values * inv_h,
                                out=self.max_tracer)

                if self.store_to_sww:
                    d = self.domain
                    d.quantities['max_stage'].centroid_values[:] = self.max_stage
                    d.quantities['max_depth'].centroid_values[:] = self.max_depth
                    d.quantities['max_speed'].centroid_values[:] = self.max_speed
                    d.quantities['max_uh'].centroid_values[:]    = self.max_uh
                    for s_idx, tname in enumerate(self._tracer_names):
                        d.quantities['max_' + tname].centroid_values[:] = \
                            self.max_tracer[s_idx]

    def _check_tracers_unchanged(self, domain):
        """A tracer added after this operator was built cannot be collected.

        add_tracer reallocates the tracer block at a new width and the sww
        file's variables are fixed when it is created, so there is nowhere for
        a late tracer to go. Fail rather than silently record maxima for the
        wrong tracers -- the arrays would still line up by index, so the
        numbers would look entirely plausible.
        """
        current = list(domain.get_tracer_names())
        if current != self._tracer_names:
            raise RuntimeError(
                'tracers changed after Collect_max_quantities_operator was '
                'created (%r -> %r); register every tracer before collecting '
                'maxima' % (self._tracer_names, current))

    def _store_tracer_maxima_from_device(self, domain):
        """Pull the device-side tracer maxima and hand them to the sww writer.

        Only at yield steps, and only when storing: the maxima themselves are
        updated on the device by gpu_max_quantities_update, so the running
        collection stays sync-free.
        """
        if not (self._tracer_names and self.store_to_sww):
            return
        self._check_tracers_unchanged(domain)
        from anuga.shallow_water.sw_domain_gpu_ext import get_max_tracers_gpu
        get_max_tracers_gpu(domain.gpu_interface.gpu_dom, self.max_tracer)
        for s_idx, tname in enumerate(self._tracer_names):
            domain.quantities['max_' + tname].centroid_values[:] = \
                self.max_tracer[s_idx]

    def parallel_safe(self):
        """Operator is applied independently on each cell and so is parallel safe."""
        return True

    def statistics(self):
        return self.label + ': Collect_max_quantity operator'

    def timestepping_statistics(self):
        from anuga import indent
        return indent + self.label + ': Collecting_max_quantity'

    def export_max_quantities_to_csv(self, filename_start='Max_Quantities_'):
        """Export max quantities to a CSV file."""

        # For GPU mode, ensure host arrays are up to date before export.
        domain = self.domain
        if self._gpu_initialized:
            try:
                from anuga.config import MULTIPROCESSOR_GPU
                from anuga.shallow_water.sw_domain_gpu_ext import get_max_quantities_gpu
                if (hasattr(domain, 'multiprocessor_mode') and
                        domain.multiprocessor_mode == MULTIPROCESSOR_GPU and
                        hasattr(domain, 'gpu_interface') and
                        domain.gpu_interface is not None):
                    get_max_quantities_gpu(
                        domain.gpu_interface.gpu_dom,
                        self.max_stage, self.max_depth,
                        self.max_speed, self.max_uh)
            except (ImportError, Exception):
                pass

        full_ids = domain.tri_full_flag.nonzero()[0]
        geo = domain.geo_reference
        out = num.vstack([
            self.xy[full_ids, 0] + geo.xllcorner,
            self.xy[full_ids, 1] + geo.yllcorner,
            self.max_stage[full_ids],
            self.max_depth[full_ids],
            self.max_speed[full_ids],
            self.max_uh[full_ids],
        ]).T
        outname = filename_start + 'P' + str(myid) + '_X_Y_Stage_Depth_Speed_UH_MAX.csv'
        num.savetxt(outname, out, delimiter=',')
