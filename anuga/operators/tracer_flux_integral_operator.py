"""Integrate the tracer fluxes through the domain boundaries.

The tracer counterpart of `boundary_flux_integral_operator`, and it applies the
same timestepping weights for the same reason: the flux kernel records one total
per substep, and each substep enters the time integral with the weight its
method gives it.

Relies on `tracer_boundary_flux_sum` being filled in compute_fluxes.
"""

import numpy as num

from anuga.operators.base_operator import Operator


class tracer_flux_integral_operator(Operator):
    """Collect the time integral of the tracer boundary fluxes during a run.

    Registered automatically by `Domain.add_tracer`; there is nothing to set
    up. A domain with no tracers never creates one.
    """

    def __init__(self, domain, description=None, label=None,
                 logging=False, verbose=False):
        Operator.__init__(self, domain, description, label, logging, verbose)
        self.domain = domain

    def __call__(self):
        """Accumulate the tracer boundary flux for this timestep."""
        domain = self.domain
        ns = domain.number_of_tracers
        if ns == 0 or domain.tracer_boundary_flux_sum is None:
            return

        dt = domain.timestep
        ts_method = domain.timestepping_method
        # (substep, tracer), matching the kernel's [substep*ns + s] indexing.
        fs = domain.tracer_boundary_flux_sum.reshape(-1, ns)

        if ts_method == 'euler':
            contribution = dt * fs[0]
        elif ts_method == 'rk2':
            contribution = 0.5 * dt * fs[0:2].sum(axis=0)
        elif ts_method == 'rk3':
            contribution = (1.0 / 6.0) * dt * (fs[0] + fs[1] + 4.0 * fs[2])
        elif ts_method == 'ader2':
            # fs[0] holds the flux from the Q^{n+1/2} midpoint state.
            contribution = dt * fs[0]
        else:
            raise Exception(
                'Cannot compute tracer boundary flux integral with '
                'timestepping method %r' % ts_method)

        domain._tracer_flux_integral[:ns] += contribution
        # Zero for the next step, exactly as the water operator does.
        domain.tracer_boundary_flux_sum[:] = 0.0

    def parallel_safe(self):
        """Applied independently on each parallel sub-domain.

        Each rank counts only the boundary edges of cells it owns -- the kernel
        skips ghosts -- so the per-rank integrals sum to the whole-domain one,
        and `get_tracer_mass` reduces the same way.
        """
        return True

    def statistics(self):
        return self.label + ': tracer_flux_integral_operator'

    def timestepping_statistics(self):
        return ''
