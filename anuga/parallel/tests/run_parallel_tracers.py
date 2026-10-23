"""Tracers survive distribute() and give the same answer in parallel (#278).

Run sequentially and under mpiexec; each writes the final tracer field, sorted
by centroid so the two orderings can be compared cell for cell, plus the mass
budget. test_parallel_tracers.py runs both and diffs them.

The tracers are added BEFORE distribute() -- that is the whole point of #278.
Boundary concentrations are set after, since tracer_boundary_values is sized by
each sub-domain's own boundary.
"""

import numpy as num

import anuga
from anuga import rectangular_cross_domain
from anuga import Reflective_boundary, Transmissive_boundary
from anuga.parallel import distribute, myid, numprocs, finalize

from anuga.utilities.parallel_abstraction import global_except_hook
import sys
sys.excepthook = global_except_hook

verbose = False

NXY = 20
LEN = 1.0
FINALTIME = 2.0
BETA = 1.0


def salinity_field(x, y):
    """Varies in both directions, so a mis-partition cannot hide in it."""
    return 0.1 + 0.6 * x / LEN + 0.3 * y / LEN


def build_domain():
    d = rectangular_cross_domain(NXY, NXY, len1=LEN, len2=LEN)
    d.set_flow_algorithm('DE0')
    d.set_quantity('elevation', lambda x, y: -0.2 * x)
    d.set_quantity('stage', lambda x, y: num.maximum(0.4 - 0.2 * x, 0.05))
    d.set_quantity('friction', 0.0)
    d.store = False

    # Two tracers: one structured, one uniform. Two of them, so a partition
    # that collapsed the slots would show up.
    d.add_tracer('salinity', beta=BETA)
    d.add_tracer('dye', beta=BETA)
    x = d.centroid_coordinates[:, 0]
    y = d.centroid_coordinates[:, 1]
    d.set_tracer('salinity', salinity_field(x, y))
    d.set_tracer('dye', 0.5)
    return d


# --------------------------------------------------------------------------
# Build on rank 0 only, exactly as a real parallel script does
# --------------------------------------------------------------------------
if numprocs > 1:
    domain = build_domain() if myid == 0 else None
    domain = distribute(domain, verbose=verbose)
else:
    domain = build_domain()

# Every rank must end up with the same tracers, in the same slots.
assert domain.get_tracer_names() == ['salinity', 'dye'], \
    'rank %d has tracers %r' % (myid, domain.get_tracer_names())
assert domain.beta_tracer == BETA, \
    'rank %d has beta_tracer %r' % (myid, domain.beta_tracer)

# The initial condition must have survived the partition: c is a known
# function of position, so check it against this rank's own centroids --
# including ghost cells, which the partition carries too.
xy = domain.centroid_coordinates
expected = salinity_field(xy[:, 0], xy[:, 1])
err = num.abs(domain.get_tracer('salinity') - expected).max()
assert err < 1e-12, \
    'rank %d: initial salinity misplaced by %g after distribute' % (myid, err)
assert num.allclose(domain.get_tracer('dye'), 0.5), \
    'rank %d: uniform tracer did not survive distribute' % myid

# --------------------------------------------------------------------------
# Boundaries, then evolve
# --------------------------------------------------------------------------
Br = Reflective_boundary(domain)
Bt = Transmissive_boundary(domain)
domain.set_boundary({'left': Br, 'top': Br, 'bottom': Br, 'right': Bt,
                     'ghost': None})
domain.set_tracer_boundary('salinity', 'left', 0.9)

for t in domain.evolve(yieldstep=0.5, finaltime=FINALTIME):
    pass

# --------------------------------------------------------------------------
# Collect the owned cells and their tracer values on rank 0
# --------------------------------------------------------------------------
owned = num.flatnonzero(domain.tri_full_flag == 1)
local = num.column_stack((domain.centroid_coordinates[owned, 0],
                          domain.centroid_coordinates[owned, 1],
                          domain.get_tracer('salinity')[owned],
                          domain.get_tracer('dye')[owned]))

# Collective, so every rank calls them.
mass_salinity = domain.get_tracer_mass('salinity')
change, flux, disc = domain.check_tracer_conservation('salinity')

if numprocs > 1:
    from mpi4py import MPI
    gathered = MPI.COMM_WORLD.gather(local, root=0)
    rows = num.vstack(gathered) if myid == 0 else None
else:
    rows = local

if myid == 0:
    order = num.lexsort((rows[:, 1], rows[:, 0]))
    rows = rows[order]
    name = ('tracers_parallel.txt' if numprocs > 1
            else 'tracers_sequential.txt')
    with open(name, 'w') as fid:
        fid.write('%d\n' % rows.shape[0])
        fid.write('%.16e %.16e %.16e\n' % (mass_salinity, change, flux))
        for r in rows:
            fid.write('%.16e %.16e %.16e %.16e\n' % (r[0], r[1], r[2], r[3]))

finalize()
