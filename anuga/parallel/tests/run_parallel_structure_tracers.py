"""Structures carry tracers in parallel (inlet_tracers.StructureTracers).

Run sequentially and under mpiexec by test_parallel_structure_tracers.py. Two
basins separated by a ridge are joined by a Boyd box culvert; the upstream
basin is dyed. Distributed over three ranks the culvert's two ends land on
different sub-domains, so the tracer given up upstream has to be summed over
the structure's ranks and delivered to the rank holding the outflow. Writes the
global tracer mass before and after, the upstream and downstream tracer, and
the water the culvert moved.
"""
import sys
import numpy as num

import anuga
from anuga import rectangular_cross_domain, Reflective_boundary
from anuga.parallel import distribute, myid, numprocs, finalize
from anuga.utilities.parallel_abstraction import global_except_hook
sys.excepthook = global_except_hook

L, W = 60.0, 6.0      # a strip: the partitioner cuts it across its length


def ridge(x, y):
    return num.where(num.abs(x - L / 2) < 6.0, 5.0, 0.0)


if myid == 0:
    d = rectangular_cross_domain(60, 6, len1=L, len2=W)
    d.set_flow_algorithm('DE1')
    d.set_quantity('elevation', ridge)
    d.set_quantity('stage', lambda x, y: num.maximum(num.where(x < L / 2, 2.0, 0.5), ridge(x, y)))
    d.set_quantity('friction', 0.01)
    d.add_tracer('dye')
    x = d.centroid_coordinates[:, 0]
    d.set_tracer('dye', num.where(x < L / 2, 0.02, 0.0))
else:
    d = None
d = distribute(d)
d.store = False
d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})

op = anuga.Boyd_box_operator(d, losses=1.5, width=1.0, height=1.0, manning=0.013,
                             end_points=[[12.0, W / 2], [L - 12.0, W / 2]])   # ends on different ranks


def side_mass(side):
    x = d.centroid_coordinates[:, 0] + d.geo_reference.xllcorner
    full = d.tri_full_flag == 1
    sel = full & ((x < L / 2) if side == 'up' else (x > L / 2))
    local = float((d.tracer_conserved_values[0] * d.areas)[sel].sum())
    if numprocs == 1:
        return local
    from mpi4py import MPI
    return MPI.COMM_WORLD.allreduce(local, op=MPI.SUM)


m0, up0, dn0 = d.get_tracer_mass('dye'), side_mass('up'), side_mass('down')
d.evolve_to_end(finaltime=20.0)
m1, up1, dn1 = d.get_tracer_mass('dye'), side_mass('up'), side_mass('down')

flow = op.accumulated_flow if (op is not None and getattr(op, 'myid', 0) == getattr(op, 'master_proc', 0)) else 0.0
if numprocs > 1:
    from mpi4py import MPI
    flow = MPI.COMM_WORLD.allreduce(flow, op=MPI.SUM)

if myid == 0:
    out = sys.argv[1] if len(sys.argv) > 1 else 'structure_tracers.txt'
    with open(out, 'w') as f:
        f.write('%d\n' % numprocs)
        f.write('%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n' % (m0, m1, up0, up1, dn0, dn1, flow))
finalize()
