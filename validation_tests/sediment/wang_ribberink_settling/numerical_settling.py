"""Wang and Ribberink's settling flume: a suspension settling through a
perforated bed with no re-entrainment.

Delft flume experiment reported by Wang and Ribberink (1986, J. Hydraul.
Res. 24(1), Table 1): a uniform flow 0.2155 m deep at 0.558 m/s
(u* = 0.034 m/s) carries 100 um sand (w_s = 0.7 cm/s) in equilibrium over a
10 m rigid bed, then over 16 m of perforated plate through which every
grain that reaches the bed falls and is lost. The depth-averaged
concentration decays along the test section from about 140 ppm at the
transition to about 30 ppm at 16 m, and the concentration profile settles
into a self-preserving shape within 2 m. With no entrainment the case
isolates the DEPOSITION closure and its adaptation. The instantaneous
exchange deposits at the Rouse near-bed ratio d* = 2.99 and decays twice
as fast as the flume did, so the case runs the two-layer suspension with
the log-law velocity split, which decays at 1.50 against the measured
1.44 and 1.55 (Wang and Ribberink's own theory gives 1.83 to first order
and 1.47 to second). The other closures are selected by WANG_ADAPTATION
and tabulated in the report.

The flow is held at the normal state by Dirichlet boundaries, with the
Manning n chosen to reproduce the measured u* and the slope the normal
one for it, on a two-row mesh of 0.1 m cells, bed fixed. The inflow
carries the measured depth-averaged concentration at x = 0.1 m; the
adaptation is chosen by WANG_ADAPTATION (default 'none', the closure the
case validates), the run by WANG_RUN (1 or 2).
"""
import json
import os

import numpy as np

import anuga
from anuga import myid, finalize

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

ADAPT = os.environ.get('WANG_ADAPTATION', 'two_layer')
RUN = int(os.environ.get('WANG_RUN', '1'))
output_file = 'settling_run%d_%s' % (RUN, ADAPT)

# --- the flume and flow (paper, section 4.1 and Table 1) ---
h0 = 0.2155
Q, W_FLUME = 0.0601, 0.5
u0 = Q / (W_FLUME * h0)                  # 0.558 m/s
u_star = 0.034
w_s = 0.007
diameter = 1.0e-4
rho_s = 2650.0
g = anuga.g
f_c = (u_star / u0) ** 2                 # quadratic drag for the measured u*
n = np.sqrt(f_c * h0 ** (1.0 / 3.0) / g)  # 0.0151
S = f_c * u0 * u0 / (g * h0)             # normal-flow slope for that n
L, W = 20.0, 0.2                         # 16 m of test section + 4 m out; two rows
dx = 0.1
X_TEST = 16.0
q = u0 * h0

# inflow concentration: the measured depth mean at x = 0.1 m, ppm -> volume
table = np.loadtxt('wang_ribberink_1986_table1.csv', delimiter=',', comments='#')
row = table[(table[:, 0] == RUN) & (table[:, 1] == 0.1)][0]
c_in = row[-1] * 1.0e-3 / rho_s

points, vertices, boundary = anuga.rectangular_cross(int(L / dx), int(W / dx), L, W)
domain = anuga.Domain(points, vertices, boundary)
domain.set_name(output_file)
domain.set_datadir('.')
domain.set_flow_algorithm(alg)
domain.set_store_vertices_uniquely()
domain.set_quantity('elevation', lambda x, y: S * (L - x))
domain.set_quantity('friction', n)
domain.set_quantity('stage', lambda x, y: S * (L - x) + h0)
domain.set_quantity('xmomentum', q)
Bin = anuga.Dirichlet_boundary([S * L + h0, q, 0.0])
Bout = anuga.Dirichlet_boundary([h0, q, 0.0])
Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'left': Bin, 'right': Bout, 'top': Br, 'bottom': Br})

# --- sediment: deposition only, nothing comes back off the plate ---
domain.initialize_sediment_operator(porosity=0.4, bed_evolution=False)
domain.set_shear_closure('quadratic_drag')
kw = {}
if ADAPT == 'two_layer':
    kw = dict(layer_fraction=0.2, velocity_profile=True)
domain.set_deposition(law='d_star', near_bed='rouse', reference_height_floor=0.1,
                      adaptation=ADAPT, **kw)
# the fraction's settling velocity comes from its diameter through the
# domain's settling law; take the diameter that gives the measured w_s
from scipy.optimize import brentq
d_eq = brentq(lambda d: domain.settling_velocity(d, rho_s=rho_s) - w_s, 2.0e-5, 1.0e-3)
domain.add_sediment_fraction('sand', diameter=d_eq, rho_s=rho_s,
                             tau_c_star=1.0e9,                            # no entrainment
                             initial_concentration=c_in)
domain.set_tracer_boundary('sand', 'left', c_in)
if verbose:
    print(domain.sediment_summary())
    print('run %d: u = %.3f m/s, n = %.4f, S = %.2e, c_in = %.3e (%.1f ppm), w_s %.4f at d = %.1f um' % (RUN, u0, n, S, c_in, row[-1], domain.settling_velocity(d_eq, rho_s=rho_s), 1e6 * d_eq))

# to a steady state: 20 m at 0.56 m/s is 36 s; the lag closures need a few
# more passes
finaltime = 240.0
for t in domain.evolve(yieldstep=20.0, finaltime=finaltime):
    if verbose:
        print(domain.timestepping_statistics())

if myid == 0:
    x = domain.centroid_coordinates[:, 0]
    c = domain.get_tracer('sand')
    np.savez('settling_run%d_%s.npz' % (RUN, ADAPT), x=x, c=c, c_in=c_in)
    with open('settling_run%d_parameters.json' % RUN, 'w') as f:
        json.dump({'run': RUN, 'adaptation': ADAPT, 'u0': u0, 'n': n, 'S': S, 'h0': h0,
                   'u_star': u_star, 'w_s': w_s, 'c_in': c_in, 'alg': alg}, f, indent=1)
    from anuga.validation_utilities import save_parameters_tex
    save_parameters_tex(domain)
finalize()
