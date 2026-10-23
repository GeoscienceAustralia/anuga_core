"""Inlet operators carry tracers in parallel (inlet_tracers.py).

Run sequentially and under mpiexec by test_parallel_inlet_tracers.py. A closed
basin at a uniform concentration c0 has one inlet adding water at c_in and one
pumping water out; the extraction region straddles sub-domains, so its pool
spans ranks and water moves between them inside it. The dye the inflow brings
cannot reach the outlet in the run, so the outlet removes c0 per unit of water
and the budget is exact; a second, mixed case (--mixed) puts the outlet where
the concentration varies, which only a global pool concentration conserves.
Writes the global water volume and tracer mass, the tracer each operator moved
(summed over ranks) and the water it moved.
"""
import sys
import numpy as num

import anuga
from anuga import rectangular_cross_domain, Reflective_boundary
from anuga.parallel import distribute, myid, numprocs, finalize
from anuga.utilities.parallel_abstraction import global_except_hook
sys.excepthook = global_except_hook

C0, C_IN = 0.01, 0.05
Q_IN, Q_OUT = 20.0, -30.0
FINALTIME = 60.0

if myid == 0:
    d = rectangular_cross_domain(24, 24, len1=100.0, len2=100.0)
    d.set_flow_algorithm('DE1')
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_quantity('friction', 0.0)
    d.add_tracer('dye', initial_value=C0)
else:
    d = None
d = distribute(d)
d.store = False
d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})

# 50 m edge to edge from the outlet: the inflow's dye cannot reach it in the
# run, so the outlet sees c0 throughout, which is what makes the budget exact.
inflow = anuga.Inlet_operator(d, anuga.Region(d, center=(10.0, 10.0), radius=6.0),
                              Q=Q_IN, tracer_concentrations={'dye': C_IN})
# centred on the middle so a 3-way partition cuts it
MIXED = '--mixed' in sys.argv
if MIXED:
    # an outlet straddling a sharp concentration front
    xc = d.centroid_coordinates[:, 0] + d.geo_reference.xllcorner
    d.set_tracer('dye', num.where(xc < 60.0, C0, 5 * C0))
outflow = anuga.Inlet_operator(d, anuga.Region(d, center=(60.0, 60.0), radius=15.0),
                               Q=Q_OUT)

m0 = d.get_tracer_mass('dye')
w0 = d.get_water_volume()
d.evolve_to_end(finaltime=FINALTIME)
m1 = d.get_tracer_mass('dye')
w1 = d.get_water_volume()


def total(op, attr):
    local = 0.0 if op is None else getattr(op.tracers, attr).get('dye', 0.0)
    if numprocs == 1:
        return local
    from mpi4py import MPI
    return MPI.COMM_WORLD.allreduce(local, op=MPI.SUM)


# net transfers: across ranks one rank's cells can gain tracer from an extraction
t_in = total(inflow, 'total_in') - total(inflow, 'total_out')
t_out = total(outflow, 'total_out') - total(outflow, 'total_in')

if myid == 0:
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    out = args[0] if args else 'inlet_tracers.txt'
    with open(out, 'w') as f:
        f.write('%d\n' % numprocs)
        f.write('%.17g %.17g %.17g %.17g %.17g %.17g\n' % (w0, w1, m0, m1, t_in, t_out))

finalize()
