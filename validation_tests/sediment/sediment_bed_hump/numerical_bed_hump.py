"""A migrating bed hump under Grass bedload.

Uniform subcritical flow along a frictionless channel, 10 m deep at 1 m/s,
held by Dirichlet boundaries at both ends, over a flat bed with a low
sin^2 hump. Bedload follows Grass's law [K-6] and moves the bed by its
divergence [K-3]; the inflow and outflow are open to bedload so the flat
reaches at either end are steady. The hump migrates downstream and its
front steepens; before the front breaks the bed follows the characteristic
solution in analytical_bed_hump.py. The run checks the Grass law, the
bedload divergence and its bed update, the open bedload boundaries, and
that the bed volume is conserved.
"""
import json
import numpy as np
import anuga
from anuga import Domain, myid, finalize, distribute
import analytical_bed_hump as analytic

output_dir = '.'
output_file = 'bed_hump'

L, W = 1000.0, 10.0       # channel length and width, m
dx = 5.0
w0 = 10.0                 # stage, m: depth 9.9 m over the flat bed
q = 10.0                  # discharge per unit width, m^2/s: u = 1.01 m/s
A_g = 0.01                # Grass coefficient; Hudson & Sweby's larger value
m = 3.0
porosity = 0.4
hump = dict(z_flat=0.1, x0=300.0, width=200.0, height=1.0)
t_shock = analytic.shock_time(w0, q, A_g, m, porosity, **hump)
finaltime = round(0.5 * t_shock, -2)   # well before the front breaks
yieldstep = finaltime / 20.0

args = anuga.get_args()
alg = args.alg
verbose = args.verbose


def bed(x, y):
    return analytic.initial_bed(x, **hump)


if myid == 0:
    points, vertices, boundary = anuga.rectangular_cross(int(L / dx), int(W / dx), L, W)
    domain = Domain(points, vertices, boundary)
    domain.set_name(output_file)
    domain.set_datadir(output_dir)
    domain.set_flow_algorithm(alg)
    domain.set_quantity('elevation', bed)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', w0)
    domain.set_quantity('xmomentum', q)
else:
    domain = None

domain = distribute(domain)

Bd = anuga.Dirichlet_boundary([w0, q, 0.0])
Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

# The sediment: Grass bedload (total load, so nothing is suspended) on an
# evolving bed, with the inflow and outflow open to it. A fraction must be
# registered for the bedload kernel to run; its grain size does not enter.
domain.initialize_sediment_operator(porosity=porosity, bed_evolution=True)
domain.add_sediment_fraction('sand', diameter=1.0e-3, initial_concentration=0.0)
domain.set_bedload('grass', K=A_g, m=m, open_boundaries=('left', 'right'))

if myid == 0:
    with open('bed_hump_parameters.json', 'w') as f:
        json.dump({'L': L, 'W': W, 'dx': dx, 'w0': w0, 'q': q, 'A_g': A_g, 'm': m,
                   'porosity': porosity, 'hump': hump, 't_shock': t_shock,
                   'finaltime': finaltime, 'yieldstep': yieldstep}, f, indent=1)
    if verbose:
        print(domain.sediment_summary())
        print('shock at %.0f s; running to %.0f s; crest speed %.2e m/s'
              % (t_shock, finaltime,
                 analytic.wave_speed(hump['z_flat'] + hump['height'], w0, q, A_g, m, porosity)))

from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)

for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
    if myid == 0 and verbose:
        print(domain.timestepping_statistics())

domain.sww_merge(delete_old=True)
finalize()
