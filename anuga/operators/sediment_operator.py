"""Fractional-step operator for suspended sediment transport and bed evolution.

Applies, once per timestep and after the hydrodynamic step:

    m_s  <-  m_s + dt (E_s - D_s)          [G-3] source term
    z    <-  z   + dt (D - E)/(1 - lambda) [G-4] Exner, suspended contribution

Both use the same limited source, so the sediment volume leaving suspension is
exactly `(1 - lambda) dz` and the mass budget closes by construction, whatever
timestepping method the flow algorithm uses.

The advective transport of `m` is NOT here: that lives in the flux kernel and
runs every substep with the rest of the shallow-water system. This operator
carries only the bed exchange.

GPU safety
----------
The work is done by `core_apply_sediment_source`, the same kernel in both
compute modes. In mode 2 it is called through the device path, so the tracer and
bed arrays are updated in place on the GPU with no host round trip. That matters
more than it looks: a fractional-step operator that writes on the host forces
`_gpu_host_writes_suppressed` on, which reactivates the sync bracket and drops
the whole run onto the host path -- a large, silent performance loss. This
operator must therefore stay on the GPU-safe list in
`Domain._has_cpu_only_fractional_operators`.
"""

import math
import warnings

import anuga.utilities.log as log
from anuga.operators.base_operator import Operator


class Sediment_transport_operator(Operator):
    """Suspended sediment transport and bed evolution, as a fractional step.

    Not normally constructed directly. The entry point is
    :meth:`~anuga.shallow_water.shallow_water_domain.Domain.initialize_sediment_operator`,
    which creates this operator and returns it, with
    :meth:`~anuga.shallow_water.shallow_water_domain.Domain.add_sediment_fraction`
    for each grain size:

    .. code-block:: python

        domain.initialize_sediment_operator(porosity=0.28)
        domain.add_sediment_fraction('sand', diameter=2.0e-4)
        domain.add_sediment_fraction('silt', diameter=2.0e-5)

    Constructing it directly is equivalent to
    `domain.initialize_sediment_operator()` with no domain-wide parameters, and
    is there for code that wants the handle without setting anything.

    ONE OPERATOR PER DOMAIN. Constructing a second returns the first: the
    kernel makes a single pass over every registered grain size, so two
    operators in the fractional-step list would apply the bed exchange twice
    per timestep.

    Parameters
    ----------
    domain : Domain
        The domain to transport sediment on.
    description, label, logging, verbose
        See :class:`~anuga.operators.base_operator.Operator`. Ignored when this
        domain already has a sediment operator.
    """

    def __new__(cls, domain, *args, **kwargs):
        # One per domain: return the existing operator rather than registering
        # a duplicate fractional step.
        for op in getattr(domain, 'fractional_step_operators', ()):
            if isinstance(op, cls):
                return op
        return super().__new__(cls)

    def __init__(self, domain, description=None, label=None, logging=False,
                 verbose=False):
        # __new__ may have handed back the operator this domain already has, in
        # which case Python still calls __init__ on it. Re-running
        # Operator.__init__ would register a second fractional step.
        if getattr(self, '_sediment_initialised', False):
            return

        Operator.__init__(self, domain, description, label, logging, verbose)
        self._sediment_initialised = True

        self._warned_no_grain_size = False

        # Reporting for the repose relaxation (spec 7).
        self.repose_sweeps = 0
        self.repose_sweeps_total = 0
        self.repose_cap_hits = 0

    def add_sediment_fraction(self, name, diameter, **kwargs):
        """Register a grain size. Delegates to `domain.add_sediment_fraction`."""
        return self.domain.add_sediment_fraction(name, diameter, **kwargs)

    def __call__(self):
        # Creating the operator without a grain size is legal -- it is how
        # operator ORDER is controlled (see initialize_sediment_operator) --
        # but evolving that way transports nothing, and the run would otherwise
        # complete silently with the bed untouched.
        if self.domain.n_sediment_classes == 0:
            if not self._warned_no_grain_size:
                self._warned_no_grain_size = True
                msg = ('sediment transport is on but no grain size is '
                       'registered, so nothing will be transported and the '
                       'bed will not change. This includes bedload, which '
                       'takes the grain diameter and R from a registered '
                       'grain size: for bed evolution with no suspended '
                       'sediment, register the grain size and select '
                       "set_bedload('engelund_hansen'), which turns the "
                       'suspended exchange off. Call '
                       'domain.add_sediment_fraction(name, diameter)')
                log.critical(msg)
                warnings.warn(msg, UserWarning, stacklevel=2)
            return

        timestep = self.domain.get_timestep()
        if timestep <= 0.0:
            return

        domain = self.domain
        suspended = getattr(domain, '_sediment_suspended_enabled', True)
        # _gpu_host_writes_suppressed marks a HOST-COHERENT window: because a
        # CPU-only operator is registered, apply_fractional_steps has already
        # pulled the device state to the host and will push the host state back
        # when the loop finishes. Running on the device inside that window is
        # not merely redundant -- the trailing push overwrites whatever the
        # device computed with the host copy, so the work is silently lost.
        # That is what made an external source contribute nothing in mode 2
        # whenever any CPU-only operator was present (#288).
        on_gpu = (domain.multiprocessor_mode == 2
                  and domain.gpu_interface is not None
                  and not getattr(domain, '_gpu_host_writes_suppressed', False))

        # Mode 1 in parallel: the fractional steps run before the step's ghost
        # exchange, so ghost cells still hold their own, unexchanged state.
        # The bedload flux on an edge shared with another rank is formed from
        # the ghost's transport vector, so refresh the ghosts first (#424).
        if (domain.parallel and not on_gpu and domain.sediment_bedload_mode
                and getattr(domain, 'sediment_bed_evolution', False)):
            domain.update_ghosts()

        if on_gpu:
            from anuga.shallow_water.sw_domain_gpu_ext import (
                apply_sediment_source_gpu, apply_bedload_gpu)
            if suspended:
                apply_sediment_source_gpu(domain.gpu_interface.gpu_dom, timestep)
            if domain.sediment_bedload_mode:
                apply_bedload_gpu(domain.gpu_interface.gpu_dom, timestep)
        else:
            from anuga.shallow_water.sw_domain_openmp_ext import (
                apply_sediment_source, apply_bedload)
            if suspended:
                apply_sediment_source(domain, timestep)
            if domain.sediment_bedload_mode:
                apply_bedload(domain, timestep)

        # Angle-of-repose relaxation LAST, so it relaxes the bed this step
        # actually produced rather than the one it started from. It is the only
        # non-cell-local sediment kernel, and the only one that can iterate.
        if domain.sediment_repose_tan > 0.0:
            if on_gpu:
                from anuga.shallow_water.sw_domain_gpu_ext import (
                    apply_repose_gpu)
                sweeps = apply_repose_gpu(domain.gpu_interface.gpu_dom)
            else:
                from anuga.shallow_water.sw_domain_openmp_ext import (
                    apply_repose)
                sweeps = apply_repose(domain)

            self.repose_sweeps = sweeps
            self.repose_sweeps_total += sweeps
            if sweeps >= domain.sediment_repose_max_sweeps:
                # Spec 7 requires this be reported, not swallowed: hitting the
                # cap means the bed may still be over-steep where the whole
                # point of the kernel was that it should not be.
                self.repose_cap_hits += 1
                if self.repose_cap_hits == 1 or self.verbose:
                    log.critical(
                        '%s: angle-of-repose relaxation hit its %d-sweep cap '
                        'at t = %g s; the bed may still exceed %.1f degrees. '
                        'Raise max_sweeps, or relax less aggressively.'
                        % (self.label, domain.sediment_repose_max_sweeps,
                           domain.get_time(),
                           math.degrees(math.atan(domain.sediment_repose_tan))))

            # Spec 7 requires a halo exchange per sweep. This exchanges once
            # per timestep instead: the sweep loop lives inside the kernel so
            # that it stays on the device, and pulling it into Python to
            # exchange between sweeps would put a host round trip in the middle
            # of every sweep and drop the run off the GPU path. In parallel the
            # consequence is that relaxation propagates across a subdomain
            # boundary one sweep per TIMESTEP rather than one per sweep, so a
            # slump spanning a boundary relaxes more slowly there. Recorded in
            # PHYSICS_SPEC 7.1; serial results are unaffected.
            # (Mode 1: the exchange below.)

        # Mode 1 in parallel: give every ghost cell its owner's bed (#423). The
        # sediment kernels also update ghost cells, but from their state BEFORE
        # this step's ghost exchange (and bedload from an incomplete
        # neighbourhood), so a ghost's bed drifted from its owner's. Elevation
        # is not in the per-step exchange, so the drift accumulated and the
        # fluxes across the partition boundary saw the wrong depth: sediment
        # mass was lost and the run differed from serial by ~1%. Mode 2
        # exchanges before its fractional steps, so its ghosts compute their
        # owners' suspended bed change exactly.
        if (domain.parallel and not on_gpu
                and (getattr(domain, 'sediment_bed_evolution', False)
                     or domain.sediment_repose_tan > 0.0)):
            _exchange_bed_ghosts(domain)
        # Mode 2 exchanges before its fractional steps, so its ghosts already
        # compute their owners' suspended bed change exactly. Bedload and
        # repose read neighbours, and an outer ghost lacks some, so with either
        # the device halo (which does not carry elevation) needs the bed too
        # (#424).
        elif (domain.parallel and on_gpu
              and getattr(domain, 'sediment_bed_evolution', False)
              and (domain.sediment_bedload_mode or domain.sediment_repose_tan > 0.0)):
            from anuga.shallow_water.sw_domain_gpu_ext import exchange_bed_ghosts_gpu
            exchange_bed_ghosts_gpu(domain.gpu_interface.gpu_dom)

    def parallel_safe(self):
        """Safe in parallel.

        The suspended exchange is cell-local. The bedload divergence DOES read
        neighbours, but only their transport vector, which is computed from
        stage, momentum and friction -- all quantities the halo exchange already
        keeps current -- so ghost cells carry a valid q_b and the divergence is
        correct on owned cells.

        Angle-of-repose relaxation also reads neighbours, and exchanges
        elevation after its sweeps. See the note at the call site for what it
        does NOT do, which is exchange between them.
        """
        return True

    def statistics(self):
        return '%s: %d sediment class(es)' % (self.label,
                                              self.domain.n_sediment_classes)

    def timestepping_statistics(self):
        if self.domain.sediment_repose_tan > 0.0:
            return ', repose sweeps %d' % self.repose_sweeps
        return ''


def _exchange_bed_ghosts(domain):
    """Copy the owner's bed into every ghost cell (mode 1, parallel).

    The bed is discontinuous in the DE algorithms: the sediment kernels shift a
    cell's three edge values together with its centroid. Copying only the
    centroid would leave a ghost's edges at its own, stale bed, so each ghost's
    edges are shifted by the same amount its centroid moves."""
    import numpy as num
    ghosts = [domain.ghost_recv_dict[p][0] for p in domain.ghost_recv_dict]
    if not ghosts:
        return
    ghosts = num.concatenate(ghosts)
    bed = domain.quantities['elevation']
    before = bed.centroid_values[ghosts].copy()
    domain.update_ghosts(['elevation'])
    shift = bed.centroid_values[ghosts] - before
    bed.edge_values[ghosts, :] += shift[:, None]
