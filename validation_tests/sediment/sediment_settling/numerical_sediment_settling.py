"""Settling of suspended sediment in a tank of still water.

Uniform initial concentration, no flow, reflective walls. The only process
is deposition onto the bed, so the run checks the deposition term [D-1], its
coupling into the conserved sediment mass [G-3], and the bed update [G-4],
against the reference solution in analytical_sediment_settling.py.
"""
import json
import numpy as np
import anuga
from anuga import Domain, myid, finalize, distribute

output_dir = '.'
output_file = 'sediment_settling'

L, W = 20.0, 4.0          # tank size, m
dx = 1.0
h0 = 1.0                  # water depth, m
c0 = 0.01                 # initial volumetric concentration (about 26 g/L)
diameter = 1.0e-4         # 100 um: v_s of a few mm/s, so the e-folding time
                          # h0/v_s is a few minutes and a short run resolves it
d_star = 1.0              # constant near-bed profile factor: D = v_s c
porosity = 0.30
finaltime = 300.0
yieldstep = 5.0

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

if myid == 0:
    points, vertices, boundary = anuga.rectangular_cross(int(L / dx), int(W / dx), L, W)
    domain = Domain(points, vertices, boundary)
    domain.set_name(output_file)
    domain.set_datadir(output_dir)
    domain.set_flow_algorithm(alg)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', h0)
else:
    domain = None

domain = distribute(domain)

Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

# The sediment: one fraction, a constant profile factor, bed evolution on.
# Still water has zero shear, so the erosion law never fires whatever its
# threshold; deposition is the only exchange.
domain.initialize_sediment_operator(porosity=porosity)
domain.set_deposition(law='d_star', near_bed='constant')
domain.add_sediment_fraction('silt', diameter=diameter, d_star=d_star,
                             initial_concentration=c0)

if myid == 0:
    v_s = float(domain.sediment_settling_velocity[0])
    with open('settling_parameters.json', 'w') as f:
        json.dump({'h0': h0, 'c0': c0, 'v_s': v_s, 'd_star': d_star,
                   'porosity': porosity, 'finaltime': finaltime,
                   'yieldstep': yieldstep}, f, indent=1)
    if verbose:
        print(domain.sediment_summary())
        print('settling velocity %.4g m/s, e-folding time %.1f s'
              % (v_s, h0 / (d_star * v_s)))

from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)

for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
    if myid == 0 and verbose:
        print(domain.timestepping_statistics())

domain.sww_merge(delete_old=True)
finalize()
