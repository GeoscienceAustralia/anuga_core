"""
Run a kinematic viscosity simulation serially and in parallel, saving the
xvelocity centroid values (in global triangle order) for comparison.

Usage:
  python run_parallel_kv_operator.py           # serial  → kv_serial_xvel.npy
  mpiexec -np 3 python run_parallel_kv_operator.py   # parallel → kv_parallel_xvel.npy
"""

import os
import numpy as num
import anuga
from anuga import rectangular_cross_domain, Reflective_boundary
from anuga import distribute, myid, numprocs, finalize
from anuga.operators.kinematic_viscosity_operator import Kinematic_viscosity_operator
from anuga.utilities.parallel_abstraction import global_except_hook
import sys
sys.excepthook = global_except_hook

verbose = False

# -----------------------------------------------------------------------
# Build domain (identical mesh for serial and parallel)
# -----------------------------------------------------------------------
m, n = 8, 8
domain = rectangular_cross_domain(m, n)
domain.set_quantity('elevation', 0.0)
domain.set_quantity('stage',     expression='1.0 + 0.5*x')
domain.set_quantity('xmomentum', expression='x - x*x')
domain.set_quantity('ymomentum', expression='y - y*y')

out_dir = os.path.dirname(os.path.abspath(__file__))
domain.set_datadir(out_dir)
domain.store = False  # no SWW output needed — we only compare xvelocity arrays

if numprocs > 1:
    domain = distribute(domain)

Br = Reflective_boundary(domain)
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

kv = Kinematic_viscosity_operator(domain, diffusivity='height',
                                   use_triangle_areas=True, verbose=verbose)

for t in domain.evolve(yieldstep=0.05, finaltime=0.15):
    pass

# -----------------------------------------------------------------------
# Gather xvelocity centroid values in global triangle order
# -----------------------------------------------------------------------
xvel   = domain.quantities['xvelocity'].centroid_values
n_full = domain.number_of_full_triangles

if numprocs == 1:
    # Serial: global index == local index, all triangles are full
    num.save(os.path.join(out_dir, 'kv_serial_xvel.npy'), xvel.copy())
else:
    import anuga.utilities.parallel_abstraction as pypar
    from mpi4py.MPI import SUM, MAX

    # Number of global triangles = sum of each rank's full triangle count
    n_full_arr = num.array([n_full], dtype=num.int64)
    n_global   = num.zeros(1, dtype=num.int64)
    pypar.comm.Allreduce(n_full_arr, n_global, op=SUM)
    n_g = int(n_global[0])

    # Each rank writes its full-triangle values into a global-sized array
    full_global = domain.tri_l2g[:n_full]   # global indices of full tris
    local_buf   = num.zeros(n_g, dtype=float)
    local_buf[full_global] = xvel[:n_full]

    # Reduce (sum) to rank 0 — no overlap between full sets, so SUM = correct value
    global_buf = num.zeros(n_g, dtype=float)
    pypar.comm.Reduce(local_buf, global_buf, op=SUM, root=0)

    if myid == 0:
        num.save(os.path.join(out_dir, 'kv_parallel_xvel.npy'), global_buf)

finalize()
