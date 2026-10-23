"""The other two distribution paths carry tracers as well (#278).

distribute() is covered end-to-end by run_parallel_tracers.py. These two share
partition_mesh with it, so the tracer values ride along automatically -- but
each has its own place where the reserved entries have to be taken back out
before they are mistaken for quantities, and that is what this exercises:

  distribute_collaborative()      pops them after its own Scatterv
  sequential_distribute_dump/load pops them when the partition file is read

Asserts internally and exits non-zero on failure; the wrapper just runs it.
"""

import os
import shutil
import tempfile

import numpy as num

from anuga import rectangular_cross_domain, Reflective_boundary
from anuga.parallel import myid, numprocs, finalize, barrier
from anuga.parallel.parallel_api import distribute_collaborative
from anuga.parallel.sequential_distribute import (
    sequential_distribute_dump, sequential_distribute_load)

from anuga.utilities.parallel_abstraction import global_except_hook
import sys
sys.excepthook = global_except_hook

BETA = 1.5


def field(x, y):
    return 0.1 + 0.6 * x + 0.3 * y


def build():
    d = rectangular_cross_domain(12, 12)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    b = Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    d.add_tracer('salinity', beta=BETA)
    d.add_tracer('dye', beta=BETA)
    xy = d.centroid_coordinates
    d.set_tracer('salinity', field(xy[:, 0], xy[:, 1]))
    d.set_tracer('dye', 0.25)
    return d


def check(domain, what):
    assert domain.get_tracer_names() == ['salinity', 'dye'], \
        '%s: rank %d has tracers %r' % (what, myid, domain.get_tracer_names())
    assert domain.beta_tracer == BETA, \
        '%s: rank %d has beta %r' % (what, myid, domain.beta_tracer)
    xy = domain.centroid_coordinates
    err = num.abs(domain.get_tracer('salinity') - field(xy[:, 0], xy[:, 1])).max()
    assert err < 1e-12, \
        '%s: rank %d salinity misplaced by %g' % (what, myid, err)
    assert num.allclose(domain.get_tracer('dye'), 0.25), \
        '%s: rank %d lost the uniform tracer' % (what, myid)


# --- collaborative ---------------------------------------------------------
domain = build() if myid == 0 else None
check(distribute_collaborative(domain), 'distribute_collaborative')

# --- dump / load -----------------------------------------------------------
# A partition written to disk and read back on each rank.
workdir = tempfile.mkdtemp(prefix='tracer_partition_') if myid == 0 else None
if numprocs > 1:
    from mpi4py import MPI
    workdir = MPI.COMM_WORLD.bcast(workdir, root=0)

if myid == 0:
    sequential_distribute_dump(build(), numprocs, verbose=False,
                               partition_dir=workdir)
barrier()

check(sequential_distribute_load(filename='domain', partition_dir=workdir,
                                 verbose=False),
      'sequential_distribute_load')

barrier()
if myid == 0:
    shutil.rmtree(workdir, ignore_errors=True)
    print('OK: tracers carried by distribute_collaborative and by dump/load')

finalize()
