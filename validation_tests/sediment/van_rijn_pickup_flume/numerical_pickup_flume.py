"""Van Rijn's pick-up flume: clear water entraining sand into suspension.

Delft Hydraulics flume experiment reported by van Rijn (1986b, Fig. 15;
1986c, Table 2 and Figs. 2-3): a uniform flow 0.25 m deep at 0.67 m/s,
initially free of sediment, runs onto a bed of 230 um sand and picks up
sediment until the suspended load reaches equilibrium. The measured
depth-integrated suspended load rises from half its equilibrium value at
x = 4 d to 1.3 times the computed equilibrium (0.036 kg/s/m) at x = 40 d.

Here the normal flow is held by Dirichlet boundaries carrying the normal
state, as in the equilibrium-flow analytical case, with the bed slope and
Manning roughness chosen to reproduce the measured bed-shear velocity
u* = 0.0477 m/s (k_s = 0.01 m). The single sand fraction has the settling
velocity of the suspended sediment (0.022 m/s, equivalent Ferguson-Church
diameter 193 um with natural-grain constants C1 = 18, C2 = 1). Clear
water enters at the inflow and the bed is fixed (the experiment kept the
scour to a minimum and van Rijn's own simulation ignored bed change).
Entrainment is the de Leeuw et al. (2020) law, whose near-bed
concentration is defined at 0.1 h, so the Rouse near-bed ratio is
referenced there; the Rouse number w_s/(kappa u*) is 1.15, the suspension
is strongly stratified, and the well-mixed d* = 1 would be wrong by an
order of magnitude.

The comparison is the depth-integrated suspended load rho_s c q along the
flume against the four measured stations.
"""
import json
import numpy as np
import anuga
from anuga import Domain, myid, finalize, distribute

output_file = 'pickup_flume'

L, W = 12.0, 0.5          # modelled reach (m); stations to x = 40 d = 10 m
dx = 0.05
h0 = 0.25                 # depth (m)
u0 = 0.67                 # velocity (m/s)
u_star = 0.0477           # measured bed-shear velocity (m/s)
w_s = 0.022               # settling velocity of the suspended sand (m/s)
diameter = 1.928e-4       # Ferguson-Church (C1 = 18, C2 = 1) inverse of w_s
d50 = 2.3e-4              # bed material
rho_s = 2650.0
qs_inf = 0.036            # kg/s/m, van Rijn's computed equilibrium load
porosity = 0.4
finaltime = 400.0         # L/u = 18 s; settling time h/w_s = 11 s
yieldstep = 20.0
g = anuga.g

# Normal flow reproducing the measured u*: tau_b = rho u*^2 = rho g h S
S = u_star ** 2 / (g * h0)
n = h0 ** (2.0 / 3.0) * np.sqrt(S) / u0
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
    domain.set_datadir('.')
    domain.set_flow_algorithm(alg)
    domain.set_quantity('elevation', bed)
    domain.set_quantity('friction', n)
    domain.set_quantity('stage', lambda x, y: bed(x, y) + h0)
    domain.set_quantity('xmomentum', q)
else:
    domain = None

domain = distribute(domain)

Bin = anuga.Dirichlet_boundary([S * L + h0, q, 0.0])
Bout = anuga.Dirichlet_boundary([h0, q, 0.0])
Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'left': Bin, 'right': Bout, 'top': Br, 'bottom': Br})

domain.initialize_sediment_operator(porosity=porosity, bed_evolution=False)
domain.set_shear_closure('quadratic_drag')
# de Leeuw et al. (2020) entrainment [E-6], with the near-bed concentration
# it is defined at, 0.1 h, as the Rouse reference level. The Smith-McLean
# law [E-1] with the default 0.01 h reference over-predicts this flume's
# equilibrium load seven-fold (see the report).
domain.set_bed_material('noncohesive', entrainment='de_leeuw',
                        de_leeuw_fit='de_leeuw_2020')
domain.set_deposition(law='d_star', near_bed='rouse', reference_height_floor=0.1)
domain.add_sediment_fraction('sand', diameter=diameter, rho_s=rho_s,
                             initial_concentration=0.0, C1=18.0, C2=1.0)
domain.set_tracer_boundary('sand', 'left', 0.0)

if myid == 0:
    v_s = float(domain.sediment_settling_velocity[0])
    with open('pickup_flume_parameters.json', 'w') as f:
        json.dump({'L': L, 'W': W, 'dx': dx, 'S': S, 'n': n, 'h0': h0, 'u0': u0,
                   'q': q, 'u_star': u_star, 'w_s': w_s, 'v_s': v_s,
                   'diameter': diameter, 'd50': d50, 'rho_s': rho_s,
                   'qs_inf': qs_inf, 'porosity': porosity,
                   'rouse': w_s / (0.4 * u_star),
                   'finaltime': finaltime, 'yieldstep': yieldstep}, f, indent=1)
    if verbose:
        print(domain.sediment_summary())
        print('S = %.3e, n = %.4f, u* = %.4f m/s, Rouse number %.2f, v_s = %.4f m/s'
              % (S, n, u_star, w_s / (0.4 * u_star), v_s))

from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)

for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
    if myid == 0 and verbose:
        print(domain.timestepping_statistics())

domain.sww_merge(delete_old=True)
finalize()
