"""Sediment mass across ranks (#423, #424).

Run sequentially and under mpiexec by test_parallel_sediment_mass.py. A closed
box with a current over a sloping, erodible bed and two suspended fractions; in
a closed box the solid balance -- suspended plus (1 - porosity) times the bed
change, over owned cells -- must stay at zero. Writes that balance and the
fraction totals, summed over ranks.

    python run_parallel_sediment_mass.py OUTFILE [--composition] [--bedload]
                                         [--thin] [--unified]
"""
import sys
import numpy as num

import anuga
from anuga.parallel import distribute, myid, numprocs, finalize
from anuga.utilities.parallel_abstraction import global_except_hook
sys.excepthook = global_except_hook

outfile = sys.argv[1]
composition = '--composition' in sys.argv
bedload = '--bedload' in sys.argv
thin = '--thin' in sys.argv          # the bed exhausts: [L-5] binds
mode = 'unified' if '--unified' in sys.argv else 'legacy'
POROSITY = 0.3

if myid == 0:
    d = anuga.rectangular_cross_domain(40, 12, len1=40.0, len2=12.0)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.set_quantity('elevation', lambda x, y: -0.01 * x)
    d.set_quantity('stage', lambda x, y: -0.01 * x + 0.6)
    d.set_quantity('xmomentum', 0.7)
    d.set_quantity('ymomentum', 0.0)
    d.set_quantity('friction', 0.03)
else:
    d = None
d = distribute(d)
d.set_store(False)
d.set_compute_mode(mode)
d.set_boundary({t: anuga.Reflective_boundary(d) for t in d.get_boundary_tags()})
d.set_sediment_parameters(porosity=POROSITY)
if bedload:
    d.set_bedload('wong_parker_eq24')
d.add_sediment_fraction(name='fine', diameter=1.0e-4, tau_c_star=0.04,
                        initial_concentration=0.0)
d.add_sediment_fraction(name='coarse', diameter=2.0e-3, tau_c_star=0.04,
                        initial_concentration=0.0)
d.set_erodible_base(depth=0.01 if thin else 0.5)
if composition:
    d.set_bed_composition({'fine': lambda x, y: 0.2 + 0.6 * (x / 40.0),
                           'coarse': lambda x, y: 0.8 - 0.6 * (x / 40.0)},
                          active_layer=0.005)
z0 = d.quantities['elevation'].centroid_values.copy()

for t in d.evolve(yieldstep=10.0, finaltime=30.0):
    pass

own = d.tri_full_flag == 1 if hasattr(d, 'tri_full_flag') else num.ones(d.number_of_elements, bool)


def total(v):
    s = float((v * d.areas)[own].sum())
    if numprocs > 1:
        from mpi4py import MPI
        s = MPI.COMM_WORLD.allreduce(s)
    return s


dz = d.quantities['elevation'].centroid_values - z0
m = d.tracer_conserved_values
values = [total(m[0] + m[1] + (1.0 - POROSITY) * dz),     # solid balance
          total(m[0]), total(m[1]), total(num.minimum(dz, 0.0))]
if myid == 0:
    with open(outfile, 'w') as f:
        f.write('%d\n' % numprocs)
        f.write(' '.join('%.17g' % v for v in values) + '\n')
finalize()
