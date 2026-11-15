"""The implicit viscosity solve must reach the device in mode 2.

Regression test for two sync faults around the host-side solve. They cannot show
on a CPU build, where host and device are the same memory:

* The operator's last call, ``distribute_to_vertices_and_edges()``, pulled the
  centroids from the device and discarded the momentum just computed, so viscosity
  had exactly no effect.
* The solve read host boundary values that were never pulled from the device, so
  every open edge looked dry and still. Reflective walls hid this: they carry no
  stress, so the problems below need a Dirichlet edge.

The oracle is the same build on the host, with ``OMP_TARGET_OFFLOAD=disabled`` set
at launch. It has to be the same build: comparing with another install would also
compare any other difference between the two sources. Runs are in CHILD processes
because that variable is only read at launch.

Only meaningful on an offload build with a device in use, so it skips elsewhere
rather than passing vacuously.
"""

import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

import anuga

CHILD = textwrap.dedent("""
    import sys

    import numpy as np

    from anuga import Dirichlet_boundary, Reflective_boundary, rectangular_cross_domain
    from anuga.operators.kinematic_viscosity_operator import Kinematic_viscosity_operator

    out, form = sys.argv[1], sys.argv[2]
    domain = rectangular_cross_domain(10, 10)
    domain.set_name('kv_gpu_host_solve')
    domain.set_store(False)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', 1.0)
    wall = Reflective_boundary(domain)
    domain.set_boundary({'left': Dirichlet_boundary([1.1, 0.5, 0.0]),
                         'right': wall, 'top': wall, 'bottom': wall})
    domain.set_compute_mode('unified')
    if form != 'none':
        Kinematic_viscosity_operator(domain, diffusivity=0.5, form=form)
    for _ in domain.evolve(yieldstep=0.1, finaltime=0.5):
        pass
    np.save(out, np.array([domain.quantities[q].centroid_values
                           for q in ('stage', 'xmomentum', 'ymomentum')]))
""")


def run_child(tmp_path, name, form, on_device):
    out = tmp_path / (name + '.npy')
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("OMP_TARGET_OFFLOAD", None)
    if not on_device:
        env["OMP_TARGET_OFFLOAD"] = "disabled"
    proc = subprocess.run([sys.executable, "-c", CHILD, str(out), form],
                          capture_output=True, text=True, timeout=180, env=env, cwd=tmp_path)
    assert proc.returncode == 0, "mode-2 evolve failed (rc=%d):\n%s" % (
        proc.returncode,
        (proc.stdout + proc.stderr)[-2000:],
    )
    return np.load(out)


pytestmark = pytest.mark.skipif(not anuga.gpu_offload_supported(),
                                reason="needs a GPU-offload build with a device")


@pytest.mark.parametrize('form', ['velocity', 'momentum'])
def test_viscosity_changes_the_flow_on_the_device(tmp_path, form):
    viscous = run_child(tmp_path, 'viscous', form, on_device=True)
    inviscid = run_child(tmp_path, 'inviscid', 'none', on_device=True)
    assert np.abs(viscous[1:] - inviscid[1:]).max() > 1e-3


@pytest.mark.parametrize('form', ['velocity', 'momentum'])
def test_device_matches_host_with_an_open_boundary(tmp_path, form):
    device = run_child(tmp_path, 'device', form, on_device=True)
    host = run_child(tmp_path, 'host', form, on_device=False)
    np.testing.assert_allclose(device, host, rtol=0.0, atol=1e-11)
