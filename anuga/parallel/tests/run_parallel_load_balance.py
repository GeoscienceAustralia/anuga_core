"""
Demonstrate load-balance statistics on a runup-style test case where one side
of the domain starts dry and gets inundated during the simulation.

Run:
  python run_parallel_load_balance.py              # serial (baseline)
  mpiexec -np 4 python run_parallel_load_balance.py
"""

import sys
import numpy as num
import anuga
from anuga import rectangular_cross_domain, Reflective_boundary, Dirichlet_boundary
from anuga import distribute, myid, numprocs, finalize
from anuga.utilities.parallel_abstraction import global_except_hook

sys.excepthook = global_except_hook

verbose = (myid == 0)

# ---------------------------------------------------------------------------
# Build domain: 200x50 mesh, left half initially dry (inundation scenario)
# ---------------------------------------------------------------------------
m, n = 200, 50
Lx, Ly = 1000.0, 250.0

domain = rectangular_cross_domain(m, n, len1=Lx, len2=Ly)
domain.set_flow_algorithm('DE0')
domain.store = False

def elevation(x, y):
    # Sloped bed: rises from -2 m at x=0 to +2 m at x=Lx
    return -2.0 + 4.0 * x / Lx

domain.set_quantity('elevation', elevation)
domain.set_quantity('stage', lambda x, y: num.maximum(elevation(x, y), -0.5))
domain.set_quantity('friction', 0.02)

if numprocs > 1:
    domain = distribute(domain)

Br = Reflective_boundary(domain)
Bd = Dirichlet_boundary([-0.5, 0.0, 0.0])
domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

# ---------------------------------------------------------------------------
# Evolve and print load balance at each yieldstep
# ---------------------------------------------------------------------------
yieldstep  = 10.0
finaltime  = 60.0

for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
    if myid == 0:
        domain.print_timestepping_statistics()
    domain.print_load_balance_statistics()

finalize()
