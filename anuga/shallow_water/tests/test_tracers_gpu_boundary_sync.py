"""A tracer boundary value set AFTER the GPU interface exists reaches the device.

The device copy of ``tracer_boundary_values`` is made once, when the arrays
are mapped, and the per-step boundary push carries only the hydrodynamic
values.  ``set_tracer_boundary`` wrote the host array and nothing else, so a
concentration prescribed after ``set_multiprocessor_mode(2)`` was never seen
by the flux kernel: every inflow edge kept injecting the zero-filled mapped
copy, and mode 2 silently disagreed with mode 1 depending purely on call
order.

Only an offload build can show it -- on a CPU build the "device" array is the
host array -- so this is a mode-1 vs mode-2 oracle in the style of
test_tracers_gpu.py, with the same per-process isolation rule.
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
    # The interface is built here, with the boundary array still zero-filled.
    d.set_multiprocessor_mode(mode)
    assert d.multiprocessor_mode == mode, "mode %d did not engage" % mode
    # ...and the concentration is prescribed only afterwards.
    d.set_tracer_boundary("c", "left", 1.0)
    m0 = d.get_tracer_mass("c")
    for _ in d.evolve(yieldstep=0.5, finaltime=2.0):
        pass
    return d, d.get_tracer_mass("c") - m0


def test_a_boundary_value_set_after_the_interface_exists_reaches_the_device():
    cpu, gained_cpu = _run(1)
    gpu, gained_gpu = _run(2)
    assert gained_cpu > 0.0, "no tracer entered in mode 1"
    # With the bug the device keeps the zero-filled copy and nothing enters.
    assert gained_gpu > 0.5 * gained_cpu, "tracer inflow lost on the device: %r vs %r" % (gained_gpu, gained_cpu)
    ma = cpu.tracer_conserved_values[0]
    mb = gpu.tracer_conserved_values[0]
    assert np.abs(ma - mb).max() < 1e-10
