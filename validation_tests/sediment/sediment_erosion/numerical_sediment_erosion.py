"""Entrainment of bed sediment into a tank of still water.

The bed is a plane of gentle slope. With the depth-slope shear closure
[T-7], its slope FROZEN at the setup value, the bed shear stress of a cell
is a function of its depth alone, so the sand fraction is above its Shields
threshold and is entrained until deposition balances it, while a second
fraction with a critical stress far above the imposed one must stay at
zero: that is the threshold control. The bed evolves: it lowers by the net
entrained volume over the packing fraction, the depth grows, and the stress
with it. There is no flow, so every cell is a closed system and follows the
per-cell reference in analytical_sediment_erosion.py. The run checks the
entrainment law [E-1]/[E-2] and its threshold, the depth-slope closure in
every cell (walls and corners included), the balance with deposition [D-1]
in the conserved sediment mass [G-3], and the bed update [G-4] under
erosion. See the analytical module for why the slope is frozen.
"""
import json
import numpy as np
import anuga
from anuga import Domain, myid, finalize, distribute

output_dir = '.'
output_file = 'sediment_erosion'

L, W = 20.0, 4.0          # tank size, m
dx = 1.0
S_bed = 1.0e-3            # bed slope
h0 = 1.0                  # water depth at x = 0, m
diameter = 5.0e-4         # 0.5 mm sand: v_s ~ 0.1 m/s, so the relaxation
                          # time h0/v_s is ~10 s and a short run resolves it;
                          # tau* = h S/(R d) ~ 1.2 on the slope, ~30 tau_c*
d_star = 1.0              # well-mixed near-bed profile factor: D = v_s c
tau_c_control = 10.0      # critical Shields stress of the control fraction,
                          # ~8x the tau* the slope gives: must not entrain
porosity = 0.30
finaltime = 120.0
yieldstep = 2.0


def bed(x, y):
    return -S_bed * x


args = anuga.get_args()
alg = args.alg
verbose = args.verbose

if myid == 0:
    points, vertices, boundary = anuga.rectangular_cross(int(L / dx), int(W / dx), L, W)
    domain = Domain(points, vertices, boundary)
    domain.set_name(output_file)
    domain.set_datadir(output_dir)
    domain.set_flow_algorithm(alg)
    domain.set_quantity('elevation', bed)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', h0)
else:
    domain = None

domain = distribute(domain)

Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

# The sediment: one fraction, no initial load, the Shields (non-cohesive)
# entrainment law with its default critical stress, the depth-slope shear
# closure so that still water carries a bed stress -- with its slope frozen
# at the setup bed, so the closure does not feed on the roughness that
# uneven erosion would otherwise create -- and an evolving bed.
domain.initialize_sediment_operator(porosity=porosity, bed_evolution=True)
domain.set_shear_closure('depth_slope', freeze_slope=True)
domain.set_bed_material('noncohesive')
domain.set_deposition(law='d_star', near_bed='constant')
domain.add_sediment_fraction('sand', diameter=diameter, d_star=d_star,
                             initial_concentration=0.0)
# The control: same grain, but a critical stress the imposed slope cannot
# reach. Nothing may be entrained into it.
domain.add_sediment_fraction('control', diameter=diameter, d_star=d_star,
                             tau_c_star=tau_c_control,
                             initial_concentration=0.0)

if myid == 0:
    v_s = float(domain.sediment_settling_velocity[0])
    h_cell = (domain.quantities['stage'].centroid_values
              - domain.quantities['elevation'].centroid_values)
    prm = {'h0': h0, 'S_bed': S_bed, 'diameter': diameter, 'v_s': v_s,
           'R': float(domain.sediment_R[0]),
           'tau_c_star': float(domain.sediment_tau_c_star[0]),
           'tau_c_control': tau_c_control,
           'gamma0': float(domain.sediment_gamma0),
           'd_star': d_star, 'porosity': porosity,
           'finaltime': finaltime, 'yieldstep': yieldstep,
           'h_cell': h_cell.tolist(),
           'x_cell': domain.centroid_coordinates[:, 0].tolist(),
           'y_cell': domain.centroid_coordinates[:, 1].tolist()}
    with open('erosion_parameters.json', 'w') as f:
        json.dump(prm, f, indent=1)
    if verbose:
        print(domain.sediment_summary())
        print('settling velocity %.4g m/s, relaxation time %.1f s'
              % (v_s, h0 / (d_star * v_s)))

from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)

for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
    if myid == 0 and verbose:
        print(domain.timestepping_statistics())

domain.sww_merge(delete_old=True)
finalize()
