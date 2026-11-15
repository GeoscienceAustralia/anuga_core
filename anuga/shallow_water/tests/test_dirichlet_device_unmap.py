"""A domain must not pick up the Dirichlet values of a domain freed before it.

Regression test for an incomplete unmap in mode 2. ``gpu_domain_map_arrays``
maps six Dirichlet arrays to the device, but ``gpu_domain_unmap_arrays`` deleted
only the three index arrays and then cleared ``mapped``, so
``gpu_dirichlet_finalize`` skipped its own delete as well. The stage and momentum
arrays were freed with their device copies still registered under the old host
addresses.

When ``malloc`` handed those addresses to the next domain, the OpenMP runtime
found them already present and skipped the copy, so that domain ran with device
values left over from the earlier one. The result was silently wrong boundary
values, or, when the sizes differed, an abort of the whole process
("partially present on the device").

The domains run in CHILD processes because that abort cannot be caught in
process, and because the process-wide offload setting left behind by another
test must not decide what this one exercises. Each result is compared with the
same domain built first in a fresh process, where nothing can be stale.

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
    import gc
    import sys

    import numpy as np

    from anuga import Dirichlet_boundary, Domain, Reflective_boundary, rectangular_cross


    def run(value):
        points, vertices, boundary = rectangular_cross(8, 8, len1=8.0, len2=8.0)
        domain = Domain(points, vertices, boundary, verbose=False)
        domain.set_name('dirichlet_unmap')
        domain.set_store(False)
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        wall = Reflective_boundary(domain)
        domain.set_boundary({'left': Dirichlet_boundary([value, 0.0, 0.0]),
                             'right': wall, 'top': wall, 'bottom': wall})
        domain.set_compute_mode('unified')
        for _ in domain.evolve(yieldstep=0.2, finaltime=0.2):
            pass
        return domain.get_quantity('stage').centroid_values.copy()


    out = sys.argv[1]
    results = []
    for value in sys.argv[2:]:
        results.append(run(float(value)))
        gc.collect()
    np.save(out, np.array(results))
""")


def run_child(tmp_path, name, values):
    out = tmp_path / (name + '.npy')
    env = dict(os.environ, OMP_NUM_THREADS="1")
    # An inherited 'disabled' would run the target regions on the host and make
    # the whole check vacuous.
    env.pop("OMP_TARGET_OFFLOAD", None)
    proc = subprocess.run([sys.executable, "-c", CHILD, str(out)] + [str(v) for v in values],
                          capture_output=True, text=True, timeout=180, env=env, cwd=tmp_path)
    assert proc.returncode == 0, "mode-2 evolve failed (rc=%d):\n%s" % (
        proc.returncode,
        (proc.stdout + proc.stderr)[-2000:],
    )
    return np.load(out)


@pytest.mark.skipif(not anuga.gpu_offload_supported(), reason="needs a GPU-offload build with a device")
def test_a_domain_does_not_inherit_the_dirichlet_values_of_a_freed_domain(tmp_path):
    high, low = 1.5, 1.1
    fresh = {high: run_child(tmp_path, 'fresh_high', [high])[0],
             low: run_child(tmp_path, 'fresh_low', [low])[0]}
    # The two values must give visibly different flows, or the check is vacuous.
    assert np.abs(fresh[high] - fresh[low]).max() > 0.01

    sequence = [high, low, high, low]
    results = run_child(tmp_path, 'sequence', sequence)
    for i, (value, stage) in enumerate(zip(sequence, results)):
        np.testing.assert_allclose(stage, fresh[value], rtol=0.0, atol=1e-12,
                                   err_msg=f'domain {i} in the sequence, Dirichlet stage {value}')
