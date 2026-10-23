"""A time-varying tracer boundary is re-evaluated EVERY STEP on the GPU path.

``set_tracer_boundary`` accepts a callable, re-evaluated each timestep by
``update_tracer_boundary_values`` -- which is hooked into ``update_boundary``.
The C-driven mode-2 evolve loops (``_evolve_one_*_step_c``) evaluate the
hydrodynamic boundaries on the device and never call ``update_boundary``
themselves; the only calls left are the ones the generic evolve loop makes
at yield points.  So the callable was evaluated once per YIELDSTEP, not once
per timestep, and every inflow edge injected a concentration up to a whole
yieldstep stale -- silently, and by an amount that depends on the caller's
yieldstep.

A linear ramp c(t) = t makes that plain: held at its yield-point value the
inflow carries about half the tracer it should over the first yield
interval.  Mode 1 (one evaluation per step, at its start, under DE0) is the
reference and the two must agree to roundoff.  Same per-process isolation
rule as test_tracers_gpu.py.
"""

import os
import warnings

import numpy as np
import pytest

import anuga
from anuga import Dirichlet_boundary, Reflective_boundary, rectangular_cross_domain

pytest.importorskip("anuga.shallow_water.sw_domain_gpu_ext")

if anuga.gpu_offload_supported() and not os.environ.get("ANUGA_GPU_TESTS_ISOLATED"):
    _skip_reason = (
        "GPU-offload build: run this file via anuga_run_isolated_tests (one fresh "
        "process per test); running it in one process aborts the NVHPC runtime."
    )
    warnings.warn(_skip_reason, stacklevel=1)
    pytest.skip(_skip_reason, allow_module_level=True)


def _ramp(t):
    return float(t)


def _run(mode):
    d = rectangular_cross_domain(10, 10, len1=100.0, len2=100.0)
    d.set_flow_algorithm("DE0")
    d.store = False
    d.set_quantity("elevation", 0.0)
    d.set_quantity("stage", 1.0)
    Br = Reflective_boundary(d)
    d.set_boundary({"left": Dirichlet_boundary([1.5, 0.0, 0.0]), "right": Br, "top": Br, "bottom": Br})
    d.add_tracer("c", beta=1.0)
    d.set_tracer("c", 0.0)
    d.set_tracer_boundary("c", "left", _ramp)
    d.set_multiprocessor_mode(mode)
    assert d.multiprocessor_mode == mode, "mode %d did not engage" % mode
    m0 = d.get_tracer_mass("c")
    # One yield interval that is many timesteps long: a per-yield refresh
    # holds c = 0 for the whole of it.
    for _ in d.evolve(yieldstep=1.0, finaltime=1.0):
        pass
    return d, d.get_tracer_mass("c") - m0


def test_a_callable_tracer_boundary_is_refreshed_every_step_on_the_device():
    cpu, gained_cpu = _run(1)
    gpu, gained_gpu = _run(2)
    assert gained_cpu > 0.0, "reference: no tracer entered"
    # Per-yield refresh: c stays at its t = 0 value, 0.0, and nothing enters.
    assert gained_gpu > 0.5 * gained_cpu, "stale inflow concentration on the device: %r vs %r" % (
        gained_gpu,
        gained_cpu,
    )
    ma = cpu.tracer_conserved_values[0]
    mb = gpu.tracer_conserved_values[0]
    assert np.abs(ma - mb).max() < 1e-8
