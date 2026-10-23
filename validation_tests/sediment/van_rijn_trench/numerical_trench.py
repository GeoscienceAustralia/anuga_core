"""Van Rijn's migrating trench: suspended sand settling into a dredged trench.

Delft Hydraulics flume experiments reported by van Rijn (1986b, Figs. 16-17):
a 0.39 m deep flow at 0.51 m/s carrying 160 um sand in equilibrium
(0.04 kg/s/m in total, 0.03 suspended and 0.01 as bedload) crosses a
trench cut across the flume. The flow slows over the trench, the suspended
sand settles onto the downstream side slope and the floor, and the trench
fills and migrates downstream. The bed profiles measured after 15 hours
for three side slopes (1:3, 1:7, 1:10) are the data.

ANUGA runs the approach flow as normal flow on a plane (Manning n from
k_s = 0.025 m), held by Dirichlet boundaries, with the trench cut into the
plane 5 m from the inflow. One sand fraction has the settling velocity of
the suspended sand (0.013 m/s; Ferguson-Church equivalent 138 um), the
Rouse near-bed ratio, an evolving bed, and the inflow carries the measured
suspended supply. As van Rijn did for his own model ("the constant of Eq.
26 was adjusted somewhat to give a suspended load of 0.03 kg/sm at the
inlet"), the entrainment constant is calibrated so that the APPROACH FLOW
is in equilibrium at the measured supply: two short pre-runs on the plane
without the trench fix it. What the trench then tests is the response of
the transport to the non-uniform flow -- the settling lag, the deposition
on the downstream slope and floor, and the Exner coupling -- against
measured bed levels.

    python numerical_trench.py [-alg DE1] [-v]

The test number (1, 2, 3) is TEST below; validate_trench.py runs all three.
The bed evolves under a morphological acceleration factor (MORFAC, default
10), so the 15 h of bed evolution take 1.5 h of simulated flow; the sww time
axis is simulated time, and the morphological time is MORFAC times it.
"""
import json
import os
import sys
import numpy as np
import anuga
from anuga import Domain, myid, finalize, distribute
import flume_tools as F

TEST = int(os.environ.get('VAN_RIJN_TRENCH_TEST', '1'))
ENTRAIN = os.environ.get('VAN_RIJN_TRENCH_ENTRAIN', 'de_leeuw')   # or 'smith_mclean'
# Morphological acceleration: the bed change per step is multiplied by MORFAC
# and the run covers 15 h of bed evolution in 15/MORFAC h of flow. The
# suspension adapts over q/(alpha v_s) ~ 20 m ~ 40 s, the bed over hours,
# so MORFAC = 10 leaves the two well separated; MORFAC = 1 is the reference.
MORFAC = float(os.environ.get('VAN_RIJN_TRENCH_MORFAC', '10'))

output_file = 'trench_test%d' % TEST

X_OFF = 5.0               # station 1 of the figure sits 5 m from the inflow
L, W = 21.0, 0.5          # x = -5 .. 16 m in the figure's coordinate
dx = 0.1
h0, u0 = 0.39, 0.51       # approach flow
k_s = 0.025               # bed roughness (m), from the velocity profiles
w_s = 0.013               # settling velocity of the suspended sand (m/s)
diameter = 1.378e-4       # Ferguson-Church (C1 = 18, C2 = 1) inverse of w_s
rho_s = 2650.0
qs_in = 0.03              # measured suspended supply, kg/s/m
qb_in = 0.01              # measured bedload supply, kg/s/m (reported only)
porosity = 0.4
hours = float(os.environ.get('VAN_RIJN_TRENCH_HOURS', '15'))   # morphological hours
finaltime = hours * 3600.0 / MORFAC                                # simulated seconds
yieldstep = 1800.0 / MORFAC
g = anuga.g

n = k_s ** (1.0 / 6.0) / (8.1 * np.sqrt(g))      # Manning-Strickler
S = n * n * u0 * u0 / h0 ** (4.0 / 3.0)          # normal-flow slope
q = u0 * h0
c_in = qs_in / (q * rho_s)                        # volumetric

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

init = np.loadtxt('van_rijn_1986_fig16_initial.csv', delimiter=',', comments='#')
init = init[init[:, 0] == TEST][:, 1:]


def plane(x):
    return S * (L - x)


def trench_depth(x):
    """depth of the trench below the plane, at figure coordinate x - X_OFF"""
    return np.interp(x - X_OFF, init[:, 0], init[:, 1], left=0.0, right=0.0)


def build(with_trench, name):
    if myid == 0:
        # the flow is one-dimensional: two rows of cells across the flume
        points, vertices, boundary = anuga.rectangular_cross(int(L / dx), 2, L, W)
        domain = Domain(points, vertices, boundary)
        domain.set_name(name)
        domain.set_datadir('.')
        domain.set_flow_algorithm(alg)
        if with_trench:
            domain.set_quantity('elevation', lambda x, y: plane(x) - trench_depth(x))
        else:
            domain.set_quantity('elevation', lambda x, y: plane(x))
        domain.set_quantity('friction', n)
        domain.set_quantity('stage', lambda x, y: plane(x) + h0)
        domain.set_quantity('xmomentum', q)
    else:
        domain = None
    domain = distribute(domain)
    Bin = anuga.Dirichlet_boundary([S * L + h0, q, 0.0])
    Bout = anuga.Dirichlet_boundary([h0, q, 0.0])
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({'left': Bin, 'right': Bout, 'top': Br, 'bottom': Br})
    return domain


def configure(domain, scale, bed_evolution):
    domain.initialize_sediment_operator(porosity=porosity, bed_evolution=bed_evolution,
                                        morphological_factor=MORFAC if bed_evolution else None)
    domain.set_shear_closure('quadratic_drag')
    if ENTRAIN == 'de_leeuw':
        A0 = domain.DE_LEEUW_FITS['de_leeuw_2020'][0]
        domain.set_bed_material('noncohesive', entrainment='de_leeuw',
                                de_leeuw_fit='de_leeuw_2020', A=A0 * scale)
        domain.set_deposition(law='d_star', near_bed='rouse', reference_height_floor=0.1)
    else:
        domain.set_bed_material('noncohesive')
        domain.sediment_gamma0 = 0.0024 * scale
        domain.set_deposition(law='d_star', near_bed='rouse')
    domain.add_sediment_fraction('sand', diameter=diameter, rho_s=rho_s,
                                 initial_concentration=c_in, C1=18.0, C2=1.0)
    domain.set_tracer_boundary('sand', 'left', c_in)


def equilibrium_load(scale):
    """Suspended load the approach flow settles to on the plane, kg/s/m."""
    domain = build(False, 'trench_prerun')
    domain.set_store(False)
    configure(domain, scale, bed_evolution=False)
    for t in domain.evolve(yieldstep=100.0, finaltime=300.0):
        pass
    c = domain.get_tracer('sand')
    uh = domain.quantities['xmomentum'].centroid_values
    x = domain.centroid_coordinates[:, 0]
    sel = (x > 0.5 * L) & (x < 0.9 * L)
    # [D-5v] under the default closure the layers move at their own log-law
    # speeds, so the transport is the profile integral of u c, which is
    # beta u h c with beta the fraction's speed factor, not u h c.
    sf = getattr(domain, 'tracer_speed_factor', None)
    global beta_in
    beta_in = float(np.mean(sf[0][sel])) if sf is not None else 1.0
    return rho_s * float(np.mean(c[sel] * uh[sel])) * beta_in


# Calibration: scale the entrainment constant so the approach flow carries
# the measured supply. E* is close to linear in the constant at these
# transport stages, so two secant steps land within a percent.
scale = 1.0
beta_in = 1.0
history = []
for it in range(3):
    qs = equilibrium_load(scale)
    history.append((scale, qs))
    if myid == 0 and verbose:
        print('calibration %d: scale %.4g -> q_s %.4g kg/s/m (target %.3g)' % (it, scale, qs, qs_in))
    if abs(qs / qs_in - 1.0) < 0.01:
        break
    scale *= qs_in / qs

# [D-5v] the inflow supply crosses the boundary at the profile-integrated
# speed too, so the concentration that carries qs_in is c_in / beta.
c_in = c_in / beta_in

# Bedload: Wong & Parker (2006) Eq 24, q* = K (tau* - tau_c*)^1.5, gives
# 0.037 kg/s/m for this fine sand in the approach flow, against the
# measured 0.01 (the relation was fitted to gravel-bed data). As van Rijn
# did ("the constant of the bed load formula was adjusted to give
# s_b = 0.01 kg/sm at the inlet"), K is scaled so that the approach flow
# carries the measured bedload; the flow is uniform there, so the closure
# is closed-form under the quadratic-drag Manning shear [T-1], [T-6].
R = rho_s / 1000.0 - 1.0
f_c = g * n * n * h0 ** (-1.0 / 3.0)
tau_star = f_c * u0 * u0 / (R * g * diameter)
K0, m0, tau_c0 = 3.97, 1.5, 0.0495
qb_wp = rho_s * K0 * max(tau_star - tau_c0, 0.0) ** m0 * np.sqrt(R * g) * diameter ** 1.5
K_cal = K0 * qb_in / qb_wp

domain = build(True, output_file)
configure(domain, scale, bed_evolution=True)
# The inflow carries the measured bedload supply as a PRESCRIBED flux, as
# van Rijn's flume was fed at a set rate. With the zero-gradient import
# (the open boundary's default) the inflow cell's import equals its own
# export, so a cell that aggrades under the fixed inflow stage sees its
# transport and its import rise, and after a few hours the boundary is an
# unbounded sediment source (a checkerboard in the inflow cells, the whole
# reach aggrading, the time step collapsing). A rigid approach strip was
# tried instead and is worse: the first erodible cell behind it receives
# no bedload, scours a 25-30 cm hole in 15 h, and feeds the trench. The
# prescribed supply removes the feedback without either.
domain.set_bedload('wong_parker_eq24', K=K_cal, open_boundaries=['right'],
                   supply={'left': qb_in / rho_s})
LOCK = 0.0

if myid == 0:
    v_s = float(domain.sediment_settling_velocity[0])
    with open('%s_parameters.json' % output_file, 'w') as f:
        json.dump({'test': TEST, 'entrain': ENTRAIN, 'X_OFF': X_OFF, 'L': L, 'W': W,
                   'dx': dx, 'S': S, 'n': n, 'h0': h0, 'u0': u0, 'q': q,
                   'k_s': k_s, 'w_s': w_s, 'v_s': v_s, 'diameter': diameter,
                   'rho_s': rho_s, 'qs_in': qs_in, 'qb_in': qb_in, 'c_in': c_in,
                   'beta_in': beta_in,
                   'porosity': porosity, 'scale': scale, 'calibration': history,
                   'morfac': MORFAC,
                   'tau_star': tau_star, 'qb_wong_parker': qb_wp, 'K_bedload': K_cal,
                   'lock': LOCK,
                   'finaltime': finaltime, 'yieldstep': yieldstep}, f, indent=1)
    if verbose:
        print(domain.sediment_summary())
        print('test %d: S = %.3e, n = %.4f, c_in = %.3e, entrainment scale %.3f, '
              'tau* %.3f, Wong-Parker bedload %.4f kg/s/m -> K = %.3f'
              % (TEST, S, n, c_in, scale, tau_star, qb_wp, K_cal))

from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)

for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
    if myid == 0 and verbose:
        print(domain.timestepping_statistics())

domain.sww_merge(delete_old=True)
finalize()
