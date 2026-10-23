"""Phase 2: tracer halo exchange across MPI ranks.

Run the same dam-break with a tracer on N ranks and compare global integrals
against the 1-rank answer. If the tracer halo is not exchanged, ghost cells
carry stale m, the fluxes at partition boundaries are wrong, and the first
moment (centre of mass) diverges even when total mass happens to survive.

    OMP_NUM_THREADS=1 mpiexec -n 1 python run_parallel_tracer.py --mode 1
    OMP_NUM_THREADS=1 mpiexec -n 2 python run_parallel_tracer.py --mode 1

NOTE mode 2 needs ONE RANK PER GPU. With more ranks than GPUs the unified path
oversubscribes the device and deadlocks (HANDOVER.md), so -n 2 --mode 2 is only
meaningful on a multi-GPU machine.
"""
import argparse
import numpy as np
import anuga
from anuga import rectangular_cross_domain, Reflective_boundary
from anuga import myid, distribute, numprocs, finalize

ap = argparse.ArgumentParser()
ap.add_argument('--mode', type=int, default=1)
ap.add_argument('--nxy', type=int, default=40)
ap.add_argument('--finaltime', type=float, default=6.0)
args = ap.parse_args()

LEN = 500.0

domain = rectangular_cross_domain(args.nxy, args.nxy, len1=LEN, len2=LEN)
domain.set_flow_algorithm('DE0')
domain.set_low_froude(0)
domain.set_datadir('.')
domain.store = False
domain.set_quantity('elevation', 0.0)
domain.set_quantity('stage', lambda x, y: np.where(x < LEN / 2, 2.0, 0.5))
domain.set_quantity('friction', 0.0)

if numprocs > 1:
    domain = distribute(domain, verbose=False)

domain.set_boundary({t: Reflective_boundary(domain)
                     for t in domain.get_boundary_tags()})

# Tracers are registered AFTER distribute: the distributed domain is a
# different object with its own element count, and add_tracer sizes the arrays
# from it.
domain.add_tracer('wedge', beta=1.0)
# A uniform tracer is a self-check that needs no reference run: with c == 1
# everywhere, m == h, so the tracer's integral must equal the water volume
# EXACTLY, on any number of ranks. If the tracer halo is not exchanged while
# the hydro halo is, the two drift apart and the gap is the halo error.
domain.add_tracer('uniform', beta=1.0)
x = domain.centroid_coordinates[:, 0]
domain.set_tracer('wedge', np.where(x < LEN / 2, 1.0, 0.0))
domain.set_tracer('uniform', 1.0)

domain.set_multiprocessor_mode(args.mode)

domain.evolve_to_end(finaltime=args.finaltime)

# Integrate over OWNED cells only -- ghosts are duplicates of another rank's
# owned cells and would be double counted.
full = domain.tri_full_flag == 1
m = domain.tracer_conserved_values[0]
areas = domain.areas
xc = domain.centroid_coordinates[:, 0]

mu = domain.tracer_conserved_values[1]
h = np.maximum(domain.quantities['stage'].centroid_values
               - domain.quantities['elevation'].centroid_values, 0.0)

local_mass = float((m[full] * areas[full]).sum())
local_momx = float((m[full] * areas[full] * xc[full]).sum())
local_umass = float((mu[full] * areas[full]).sum())
local_vol = float((h[full] * areas[full]).sum())

if numprocs > 1:
    from anuga import pypar_available
    import mpi4py.MPI as MPI
    comm = MPI.COMM_WORLD
    mass = comm.allreduce(local_mass, op=MPI.SUM)
    momx = comm.allreduce(local_momx, op=MPI.SUM)
    umass = comm.allreduce(local_umass, op=MPI.SUM)
    vol = comm.allreduce(local_vol, op=MPI.SUM)
else:
    mass, momx = local_mass, local_momx
    umass, vol = local_umass, local_vol

if myid == 0:
    print('RESULT nprocs=%d mode=%d mass=%.12e momx=%.12e com=%.12e'
          % (numprocs, args.mode, mass, momx, momx / mass))
    print('SELFCHECK uniform_tracer_mass=%.12e water_volume=%.12e rel_gap=%.3e'
          % (umass, vol, abs(umass - vol) / max(abs(vol), 1.0)))

finalize()
