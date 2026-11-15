"""Run a storm with rainfall losses on a tilted plane (sequentially or in
parallel) to support test_parallel_rainfall_losses.py
"""

import sys

import numpy as np

import anuga
from anuga import distribute, finalize, myid, numprocs

# Set up MPI to abort on error
from anuga.utilities.parallel_abstraction import global_except_hook
sys.excepthook = global_except_hook

if myid == 0:
    base_domain = anuga.rectangular_cross_domain(20, 10, len1=200., len2=100.)
    base_domain.set_flow_algorithm('DE0')
    base_domain.set_datadir('.')
    base_domain.set_quantity('elevation', lambda x, y: 0.01*x)
    base_domain.set_quantity('friction', 0.03)
    base_domain.set_quantity('stage', expression='elevation')
else:
    base_domain = None

if numprocs == 1:
    domain = base_domain
    domain.set_name('s_rainfall_losses')
else:
    domain = distribute(base_domain)
    domain.set_name('p_rainfall_losses')

Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

# 60 mm in three 5-minute bursts, on curve numbers that vary across the plane
losses = anuga.SCS_curve_number(
    curve_number=lambda x, y: np.where(y < 50.0, 95.0, 70.0))
rain = anuga.Hyetograph_rate_operator(
    domain, times=[0.0, 300.0, 600.0, 900.0], depths=[10.0, 30.0, 20.0],
    losses=losses)

for t in domain.evolve(yieldstep=300.0, finaltime=1200.0):
    pass

# Volumes over the whole domain: summed over the processes
from anuga.utilities import parallel_abstraction as pypar   # noqa: E402
totals = np.array([rain.cumulative_influx, rain.cumulative_loss])
if numprocs > 1:
    for p in range(1, numprocs):
        if myid == 0:
            totals = totals + pypar.receive(p)
        elif myid == p:
            pypar.send(totals, 0)
if myid == 0:
    print('TOTALS %.12e %.12e' % tuple(totals))

domain.sww_merge(delete_old=True)
finalize()
