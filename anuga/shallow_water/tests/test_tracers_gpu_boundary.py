"""A tracer carried in across a boundary must not fault the device in mode 2.

Regression test for a host-struct dereference inside the offloaded flux loop:
the tracer advection block read ``D->boundary_length`` inside the
``omp target`` element loop, and ``D`` is a host pointer that is never mapped,
so on a real GPU-offload build the read faulted with
``CUDA_ERROR_ILLEGAL_ADDRESS`` as soon as a tracer was registered, and the
NVHPC runtime aborted the whole process.

That abort cannot be caught in-process, which is also why the fault was never
seen: the mode-1 vs mode-2 oracles were only ever run on the CPU build, where
the target region executes on the host and a ``D->`` load is harmless.  So
the evolve runs in a CHILD process and the assertion is on its exit status.

The child drives water (and tracer) IN through a Dirichlet boundary with a
prescribed concentration, so the inflow branch that indexes the boundary
array with ``boundary_length`` is exercised on every step, and it checks the
mode-2 tracer field against a mode-1 run of the same problem: a wrong hoist
(the right type, the wrong member) would index the wrong slot and disagree,
where a mere "did not crash" check would pass.

Only meaningful on an offload build with a device in use, so it skips
elsewhere rather than passing vacuously.
"""

import os
import subprocess
import sys
import textwrap

import pytest

import anuga

pytest.importorskip("anuga.shallow_water.sw_domain_gpu_ext")

CHILD = textwrap.dedent("""
    import numpy as np
    import anuga
    from anuga import Dirichlet_boundary, Reflective_boundary
    from anuga import rectangular_cross_domain

    def run(mode):
        d = rectangular_cross_domain(10, 10, len1=100.0, len2=100.0)
        d.set_flow_algorithm('DE0')
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 1.0)
        # Higher stage on the left: water flows IN there for the whole run,
        # carrying the prescribed concentration; walls elsewhere.
        Br = Reflective_boundary(d)
        d.set_boundary({'left': Dirichlet_boundary([1.5, 0.0, 0.0]),
                        'right': Br, 'top': Br, 'bottom': Br})
        d.add_tracer('c', beta=1.0)
        d.set_tracer('c', 0.0)
        d.set_tracer_boundary('c', 'left', 1.0)
        d.set_multiprocessor_mode(mode)
        if d.multiprocessor_mode != mode:
            raise SystemExit('mode %d did not engage' % mode)
        m0 = d.get_tracer_mass('c')
        for _ in d.evolve(yieldstep=0.5, finaltime=2.0):
            pass
        return d, d.get_tracer_mass('c') - m0

    cpu, gained_cpu = run(1)
    gpu, gained_gpu = run(2)
    if not gained_cpu > 0.0:
        raise SystemExit('no tracer entered in mode 1: %r' % gained_cpu)
    ma = cpu.tracer_conserved_values[0]
    mb = gpu.tracer_conserved_values[0]
    err = float(np.abs(ma - mb).max())
    if not err < 1e-10:
        raise SystemExit('mode 2 tracer disagrees with mode 1: max |dm| = %g'
                         % err)
    print('OK')
""")


@pytest.mark.skipif(not anuga.gpu_offload_supported(), reason="needs a GPU-offload build with a device")
def test_a_tracer_carried_in_across_a_boundary_does_not_fault_the_device():
    env = dict(os.environ, OMP_NUM_THREADS="1")
    # An inherited 'disabled' would run the target regions on the host and
    # make the whole check vacuous.
    env.pop("OMP_TARGET_OFFLOAD", None)
    proc = subprocess.run([sys.executable, "-c", CHILD], capture_output=True, text=True, timeout=180, env=env)
    assert proc.returncode == 0, "mode-2 evolve with a tracer inflow failed (rc=%d):\n%s" % (
        proc.returncode,
        (proc.stdout + proc.stderr)[-2000:],
    )
