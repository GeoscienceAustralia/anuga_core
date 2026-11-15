"""A bed changed mid-run, on owned cells only, across ranks.

Run sequentially and under mpiexec by test_parallel_bed_change.py. A sloping
channel fed by an inlet; at t = 1 s a strip of the bed is raised, in one of
three ways, on the cells each rank owns (not its ghosts), so the ghosts are
right only if the domain copies the owners' new bed into them. Writes the
number of ranks, the water volume and the final stage over the owned cells,
ordered by centroid so a parallel run lines up with the sequential one.

    python run_parallel_bed_change.py OUTFILE {none,set_quantity,direct,operator}
                                      [--unified]
"""
import sys
import numpy as num

import anuga
from anuga.parallel import distribute, myid, numprocs, finalize
from anuga.utilities.parallel_abstraction import global_except_hook
sys.excepthook = global_except_hook

outfile, change = sys.argv[1], sys.argv[2]
mode = 'unified' if '--unified' in sys.argv else 'legacy'

if myid == 0:
    d = anuga.rectangular_cross_domain(20, 10, len1=200.0, len2=100.0)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.set_quantity('elevation', lambda x, y: -x / 100.0)
    d.set_quantity('friction', 0.01)
    d.set_quantity('stage', expression='elevation + 0.2')
else:
    d = None
d = distribute(d)
d.set_store(False)
d.set_boundary({'left': anuga.Reflective_boundary(d),
                'right': anuga.Dirichlet_boundary([-3.0, 0, 0]),
                'top': anuga.Reflective_boundary(d),
                'bottom': anuga.Reflective_boundary(d)})
anuga.Inlet_operator(d, [[10.0, 20.0], [10.0, 80.0]], 5.0, verbose=False)

own = d.tri_full_flag == 1 if hasattr(d, 'tri_full_flag') else num.ones(d.number_of_elements, bool)
x = d.centroid_coordinates[:, 0]
strip = num.flatnonzero((x > 80.0) & (x < 120.0) & own)

if change == 'operator':
    from anuga.operators.set_elevation_operator import Set_elevation_operator
    # Keeps the depth: stage moves with the bed.
    Set_elevation_operator(d, elevation=lambda t: -1.0 if t < 1.0 else -0.9,
                           indices=strip)
d.set_compute_mode(mode)

for t in d.evolve(yieldstep=0.5, finaltime=3.0):
    if abs(t - 1.0) < 1e-9:
        if change == 'set_quantity':
            d.set_quantity('elevation', -0.9, indices=strip, location='centroids')
        elif change == 'direct':
            bed = d.quantities['elevation']
            bed.centroid_values[strip] += 0.15
            bed.edge_values[strip] += 0.15
            d.sync_elevation_to_device()

volume = d.get_water_volume()
if d.get_compute_mode() == 'unified' and d.gpu_interface is not None:
    d.gpu_interface.sync_from_device()
mine = num.vstack([d.centroid_coordinates[own].T,
                   d.quantities['stage'].centroid_values[own]])
if numprocs > 1:
    from mpi4py import MPI
    mine = num.hstack(MPI.COMM_WORLD.gather(mine, root=0) or [mine])
if myid == 0:
    order = num.lexsort((mine[1], mine[0]))
    with open(outfile, 'w') as f:
        f.write('%d\n' % numprocs)
        f.write('%.17g\n' % volume)
        f.write(' '.join('%.17g' % v for v in mine[2, order]) + '\n')
finalize()
