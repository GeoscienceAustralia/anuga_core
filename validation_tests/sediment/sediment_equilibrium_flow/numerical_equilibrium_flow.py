"""Equilibrium suspended load in normal flow.

Uniform flow at normal depth down a plane channel with Manning roughness,
held at both ends by Dirichlet boundaries carrying the normal state. Under
the quadratic-drag shear closure the bed stress in normal flow is exactly
rho g h S, so the Shields stress, the entrainment rate and the equilibrium
concentration are closed-form. Clear water enters at the inflow and the
suspended load relaxes to the equilibrium along the channel over the
settling length: see analytical_equilibrium_flow.py. The bed is held fixed
so that the flow, and with it the reference, stay exact. The run checks
the quadratic-drag closure [T-1] with the Manning friction factor [T-6],
the entrainment law [E-1]/[E-2] in flowing water, its balance with
deposition [D-1] in the sediment mass balance [G-3] under advection, and
that a fixed bed stays fixed.
"""
import json
import numpy as np
import anuga
from anuga import Domain, myid, finalize, distribute
import analytical_equilibrium_flow as analytic

output_dir = '.'
output_file = 'equilibrium_flow'

L, W = 300.0, 10.0        # channel length and width, m
dx = 2.5
S = 2.0e-4                # bed slope
n = 0.03                  # Manning roughness
h0 = 1.0                  # normal depth, m: u = 0.47 m/s, Fr = 0.15
diameter = 1.0e-4         # 100 um: v_s ~ 8 mm/s, settling length q/v_s ~ 59 m;
                          # tau* = h S/(R d) = 1.2, 30 times tau_c*, c_eq ~ 4%
d_star = 1.0
porosity = 0.30
finaltime = 1500.0        # steady after L/u = 640 s
yieldstep = 50.0

u0 = analytic.normal_velocity(h0, S, n)
q = u0 * h0

args = anuga.get_args()
alg = args.alg
verbose = args.verbose


def bed(x, y):
    return S * (L - x)


if myid == 0:
    points, vertices, boundary = anuga.rectangular_cross(int(L / dx), int(W / dx), L, W)
    domain = Domain(points, vertices, boundary)
    domain.set_name(output_file)
    domain.set_datadir(output_dir)
    domain.set_flow_algorithm(alg)
    domain.set_quantity('elevation', bed)
    domain.set_quantity('friction', n)
    domain.set_quantity('stage', lambda x, y: bed(x, y) + h0)
    domain.set_quantity('xmomentum', q)
else:
    domain = None

domain = distribute(domain)

# Normal flow is the steady state of the sloping rough channel, so the
# normal state at both ends holds it.
Bin = anuga.Dirichlet_boundary([S * L + h0, q, 0.0])
Bout = anuga.Dirichlet_boundary([h0, q, 0.0])
Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'left': Bin, 'right': Bout, 'top': Br, 'bottom': Br})

# The sediment: one fraction, no initial load, the Shields entrainment law
# with its default critical stress under the quadratic-drag closure (the
# default), well-mixed deposition, a fixed bed, and clear water at the inflow.
domain.initialize_sediment_operator(porosity=porosity, bed_evolution=False)
domain.set_shear_closure('quadratic_drag')
domain.set_bed_material('noncohesive')
domain.set_deposition(law='d_star', near_bed='constant')
domain.add_sediment_fraction('sand', diameter=diameter, d_star=d_star,
                             initial_concentration=0.0)
domain.set_tracer_boundary('sand', 'left', 0.0)

if myid == 0:
    v_s = float(domain.sediment_settling_velocity[0])
    R = float(domain.sediment_R[0])
    tau_c_star = float(domain.sediment_tau_c_star[0])
    gamma0 = float(domain.sediment_gamma0)
    c_eq = analytic.equilibrium_concentration(h0, S, v_s, R, diameter,
                                              tau_c_star, gamma0, d_star)
    with open('equilibrium_flow_parameters.json', 'w') as f:
        json.dump({'L': L, 'W': W, 'S': S, 'n': n, 'h0': h0, 'u0': u0, 'q': q,
                   'diameter': diameter, 'v_s': v_s, 'R': R,
                   'tau_c_star': tau_c_star, 'gamma0': gamma0,
                   'd_star': d_star, 'porosity': porosity, 'c_eq': c_eq,
                   'finaltime': finaltime, 'yieldstep': yieldstep}, f, indent=1)
    if verbose:
        print(domain.sediment_summary())
        print('normal velocity %.4g m/s, settling velocity %.4g m/s, '
              'settling length %.1f m, tau* %.3g, c_eq %.4g'
              % (u0, v_s, analytic.settling_length(q, v_s, d_star),
                 analytic.shields_stress(h0, S, R, diameter), c_eq))

from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)

for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
    if myid == 0 and verbose:
        print(domain.timestepping_statistics())

domain.sww_merge(delete_old=True)
finalize()
