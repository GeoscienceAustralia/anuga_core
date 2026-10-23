"""Dam-break over an erodible sand bed around a pier: the Zaragoza Kinect
experiments of Segovia-Burillo et al. (2026), J. Hydraulic Research 64(2),
data DOI 10.5281/zenodo.17777387.

A 6 m x 0.24 m PVC flume (Manning n = 0.010) fed by a pneumatic gate from an
81 cm x 157 cm reservoir holding 8 cm of water; the flume bed is dry. From
1.0 to 2.5 m past the gate the bed is a 5 cm layer of sand (0.7 and 1.3 mm,
50/50, rho_s 2650, porosity 0.34) levelled with the PVC; a pier stands on
the centreline 1.6 m from the gate, a 3 cm cylinder (P1) or a 3 x 6.5 cm
rounded rectangle (P2). The flat reach ends at 3.0 m and the flume then
drops at 4 % to the sediment traps. The 8 cm dam-break is released three
times over the evolving bed; a Kinect sensor records the bed after each and
the water surface during each.

Selected by environment: ZARAGOZA_CASE (P1 default, or P2) and
ZARAGOZA_EVENTS (1 default, up to 3). Each event runs T_EVENT seconds from
the gate opening; between events the reservoir is refilled to 8 cm, the
flume dried, and the bed kept.
"""
import json
import os

import numpy as np

import anuga
from anuga import Reflective_boundary, Transmissive_boundary, myid, finalize

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

CASE = os.environ.get('ZARAGOZA_CASE', 'P1')
EVENTS = int(os.environ.get('ZARAGOZA_EVENTS', '1'))
T_EVENT = float(os.environ.get('ZARAGOZA_T_EVENT', '5.0'))     # s per dam-break
# Gate: 0 = instantaneous removal; otherwise the time (s) at which the gate
# is fully open. The paper's gate lifts from the bottom, under 1 cm in the
# first 50 ms and the rest in under 100 ms; here the opening a(t) is 1 cm
# at GATE_TIME/3 and the full depth at GATE_TIME, and the gate is a sill
# across the flume mouth whose crest is set each step so that the weir
# flow over it equals the orifice flow under the real gate.
GATE_TIME = float(os.environ.get('ZARAGOZA_GATE_TIME', '0.0'))
GATE_DT = 0.005
# Extra refinement of the contraction at the gate (m^2 per triangle), 0 = none
A_GATE = float(os.environ.get('ZARAGOZA_A_GATE', '0.0'))
# [T-16] pier-scour correction: the bed shear the sediment sees is amplified
# by 1 + PIER_AMP exp(-(r - R)/(PIER_LEN R)) around the pier, r the distance
# from its centreline (the P2 pier's, from its axis), R its half-width and
# PIER_LEN the decay length in half-widths, standing in for the horseshoe
# vortex. PIER_AMP = 0 is off.
PIER_AMP = float(os.environ.get('ZARAGOZA_PIER_AMP', '4.0'))
PIER_LEN = float(os.environ.get('ZARAGOZA_PIER_LEN', '2.0'))
# 1: the amplification acts on the front and flanks only (the horseshoe
# vortex), fading to none in the wake; 0: all round.
PIER_FRONT = int(os.environ.get('ZARAGOZA_PIER_FRONT', '0'))
# Suspended load as well as bedload (off by default: the case is bedload
# only). With it, the sand can be entrained (Smith-McLean, critical Shields
# number 0.047) and settle with the Rouse near-bed profile, for testing whether
# sand lifted at the pier is carried past the near wake.
SUSPENSION = int(os.environ.get('ZARAGOZA_SUSPENSION', '0'))
# A downstream tail on the shear amplification (on by default, A = 2 over
# 5 pier diameters; ZARAGOZA_WAKE_AMP=0 turns it off). Behind the
# pier the depth-averaged flow is sheltered, so its shear drops and the bedload
# stops within a diameter; in the flume the horseshoe-vortex legs and the shed
# wake vortices keep the bed shear raised for several diameters. The tail adds
# WAKE_AMP * exp(-s / (WAKE_LEN * D)) * exp(-y'^2 / (2 D^2)) to the factor, with
# s the distance downstream of the pier's rear face and D the pier diameter.
WAKE_AMP = float(os.environ.get('ZARAGOZA_WAKE_AMP', '2.0'))
WAKE_LEN = float(os.environ.get('ZARAGOZA_WAKE_LEN', '5.0'))
# ZARAGOZA_TAG suffixes the output names, so experiments can run side by side
output_file = 'pier_%s%s' % (CASE, os.environ.get('ZARAGOZA_TAG', ''))

# --- geometry (m); x from the gate along the flume, y from the right wall ---
W = 0.24
X_RES, W_RES = -1.57, 0.81            # reservoir, centred on the flume
Y_RES0 = 0.5 * (W - W_RES)
X_SAND0, X_SAND1 = 1.0, 2.5           # erodible reach
SAND_DEPTH = 0.05
X_FLAT, X_END, S_DOWN = 3.0, 3.5, 0.04
X_PIER, Y_PIER = 1.6, 0.5 * W
H0 = 0.08                             # reservoir depth
N_PVC, N_SAND = 0.010, 0.015          # Manning; sand from the paper's Table 1
POROSITY, D_SAND = 0.34, 1.0e-3       # one fraction at the mix's geometric mean
REPOSE = 32.0                         # degrees

# Mesh resolution, m^2 per triangle. A_COARSE is the overall cap and applies
# to the reservoir, where the flow is a slow drawdown; the whole 24 cm flume
# is held at A_FLUME so the channel is resolved across its width from the
# gate to the exit and not only over the sand, which is finer again, and the
# pier finer still.
A_COARSE, A_FLUME, A_SAND, A_PIER = 4.0e-4, 1.5e-4, 1.0e-4, 2.0e-5
# Wake-resolution experiments (defaults leave the case unchanged): the pier
# box's triangle area, how far downstream the box reaches, and the output
# interval, fine enough to see vortex shedding (period about 0.15 s).
A_PIER = float(os.environ.get('ZARAGOZA_A_PIER', str(A_PIER)))
BOX_X1 = float(os.environ.get('ZARAGOZA_BOX_X1', '1.75'))
YIELD = float(os.environ.get('ZARAGOZA_YIELD', '0.1'))


def pier_polygon(case, n=16):
    r = 0.015
    if case == 'P1':
        th = np.linspace(0.0, 2 * np.pi, 2 * n, endpoint=False)
        return [(X_PIER + r * np.cos(a), Y_PIER + r * np.sin(a)) for a in th]
    # P2: a stadium 6.5 cm long along the flow, 3 cm wide
    half = 0.5 * (0.065 - 0.03)
    pts = []
    for a in np.linspace(-0.5 * np.pi, 0.5 * np.pi, n):          # downstream cap
        pts.append((X_PIER + half + r * np.cos(a), Y_PIER + r * np.sin(a)))
    for a in np.linspace(0.5 * np.pi, 1.5 * np.pi, n):           # upstream cap
        pts.append((X_PIER - half + r * np.cos(a), Y_PIER + r * np.sin(a)))
    return pts


bounding = [(X_RES, Y_RES0), (0.0, Y_RES0), (0.0, 0.0), (X_END, 0.0),
            (X_END, W), (0.0, W), (0.0, Y_RES0 + W_RES), (X_RES, Y_RES0 + W_RES)]
boundary_tags = {'wall': [0, 1, 2, 4, 5, 6, 7], 'outflow': [3]}
hole = pier_polygon(CASE)
# The reaches that span the full width of the flume are delimited by
# breaklines across it and their resolution set by a point inside each,
# NOT by an interior-region polygon: a polygon inset from the walls leaves
# a sliver between itself and the wall that the mesher fills with a band of
# tiny triangles. The pier box is a genuine interior polygon, well clear of
# both walls.
X_IN, X_OUT = X_SAND0 - 0.05, X_SAND1 + 0.05      # ends of the sand reach
pier_box = [(X_PIER - 0.12, 0.03), (BOX_X1, 0.03),
            (BOX_X1, W - 0.03), (X_PIER - 0.12, W - 0.03)]
regions = [(pier_box, A_PIER)]
region_points = [[0.5 * X_IN, 0.5 * W, A_FLUME],                 # gate to sand
                 [0.5 * (X_IN + X_PIER) - 0.1, 0.06, A_SAND],   # the sand reach
                 [0.5 * (X_OUT + X_END), 0.5 * W, A_FLUME]]      # sand to exit
# The gate line across the mouth of the flume. It is interior to the mesh --
# the reservoir and the flume are one domain -- so without a breakline the
# triangles straddle it, the contraction is resolved raggedly and the sill
# that ZARAGOZA_GATE_TIME raises has a ragged crest. As a breakline the mesh
# generator puts cell edges along it.
gate_line = [(0.0, 0.0), (0.0, W)]
# and one across each end of the sand reach, so the three flume compartments
# are closed and each takes the resolution of its own marker point
sand_in_line = [(X_IN, 0.0), (X_IN, W)]
sand_out_line = [(X_OUT, 0.0), (X_OUT, W)]
if A_GATE > 0.0:
    regions += [([(-0.2, -0.06), (-0.005, -0.06), (-0.005, W + 0.06), (-0.2, W + 0.06)], A_GATE),
                ([(0.005, 0.005), (0.4, 0.005), (0.4, W - 0.005), (0.005, W - 0.005)], A_GATE)]

if myid == 0:
    domain = anuga.create_domain_from_regions(
        bounding, boundary_tags=boundary_tags, maximum_triangle_area=A_COARSE,
        interior_regions=regions,
        interior_holes=[hole], hole_tags=[{'wall': list(range(len(hole)))}],
        breaklines=[gate_line, sand_in_line, sand_out_line],
        regionPtArea=region_points,
        mesh_geo_reference=anuga.Geo_reference(xllcorner=0.0, yllcorner=0.0),
        use_cache=False, verbose=verbose)
    domain.set_name(output_file)
    domain.set_datadir('.')
    domain.set_flow_algorithm(alg)
    domain.set_store_vertices_uniquely()

    def elevation(x, y):
        return np.where(x > X_FLAT, -S_DOWN * (x - X_FLAT), 0.0)

    def friction(x, y):
        return np.where((x >= X_SAND0) & (x <= X_SAND1), N_SAND, N_PVC)

    domain.set_quantity('elevation', elevation)
    domain.set_quantity('friction', friction)
    domain.set_quantity('stage', expression='elevation')

    # --- sediment: bedload only, over a 5 cm erodible layer in the sand reach ---
    domain.initialize_sediment_operator(porosity=POROSITY, bed_evolution=True)
    if SUSPENSION:
        domain.set_deposition(near_bed='rouse')
        domain.add_sediment_fraction('sand', diameter=D_SAND, tau_c_star=0.047,
                                     initial_concentration=0.0)
    else:
        domain.set_deposition(law='threshold', tau_d=0.0)          # no settling exchange
        domain.add_sediment_fraction('sand', diameter=D_SAND, tau_c_star=1.0e9,   # no entrainment
                                     initial_concentration=0.0)
    domain.set_bedload('wong_parker_eq24', open_boundaries=['outflow'], min_depth=2.0 * D_SAND)
    xc = domain.centroid_coordinates[:, 0]
    z = domain.quantities['elevation'].centroid_values
    base = np.where((xc >= X_SAND0) & (xc <= X_SAND1), z - SAND_DEPTH, z)
    domain.set_erodible_base(elevation=base)
    domain.set_angle_of_repose(REPOSE)
    if PIER_AMP > 0.0:
        R = 0.015
        half = 0.0 if CASE == 'P1' else 0.5 * (0.065 - 0.03)

        def amplification(x, y):
            dx = np.maximum(np.abs(x - X_PIER) - half, 0.0)
            r = np.hypot(dx, y - Y_PIER)
            amp = PIER_AMP * np.exp(-np.maximum(r - R, 0.0) / (PIER_LEN * R))
            if PIER_FRONT:
                # weight 1 upstream of the pier's axis, falling to 0 over the
                # downstream quadrant: cos of the angle from upstream, clipped
                ang = np.arctan2(np.abs(y - Y_PIER), -(x - X_PIER))
                amp = amp * np.clip(1.5 - ang / (0.5 * np.pi), 0.0, 1.0)
            if WAKE_AMP > 0.0:
                D = 2.0 * R
                s_down = x - (X_PIER + half + R)
                tail = (WAKE_AMP * np.exp(-np.maximum(s_down, 0.0) / (WAKE_LEN * D))
                        * np.exp(-0.5 * ((y - Y_PIER) / D) ** 2))
                amp = amp + np.where(s_down > 0.0, tail, 0.0)
            return 1.0 + amp
        domain.set_shear_amplification(amplification)
    if verbose:
        print(domain.sediment_summary())
        print('mesh: %d triangles' % domain.number_of_elements)
else:
    domain = None

domain = anuga.distribute(domain)

Br = Reflective_boundary(domain)
Bt = Transmissive_boundary(domain)
domain.set_boundary({'wall': Br, 'outflow': Bt})


def release_reservoir():
    x = domain.centroid_coordinates[:, 0]
    z = domain.quantities['elevation'].centroid_values
    stage = np.where(x < 0.0, H0, z)
    domain.set_quantity('stage', stage, location='centroids')
    domain.set_quantity('xmomentum', 0.0)
    domain.set_quantity('ymomentum', 0.0)


if myid == 0:
    from anuga.validation_utilities import save_parameters_tex
    save_parameters_tex(domain)

z0 = domain.quantities['elevation'].centroid_values.copy()
xc_all = domain.centroid_coordinates[:, 0]
gate_cells = (xc_all > -0.03) & (xc_all <= 0.0) & (domain.centroid_coordinates[:, 1] > 0.0) \
             & (domain.centroid_coordinates[:, 1] < W)


def gate_crest(tau):
    """Sill crest (m) at time tau after the release: the weir head that
    passes the orifice flow under a gate open by a(tau)."""
    if tau >= GATE_TIME:
        return 0.0
    t1 = GATE_TIME / 3.0
    a = 0.01 * tau / t1 if tau < t1 else 0.01 + (H0 - 0.01) * (tau - t1) / (GATE_TIME - t1)
    q = 0.6 * a * np.sqrt(2.0 * 9.8 * H0)                    # orifice
    head = (q / (0.385 * np.sqrt(2.0 * 9.8))) ** (2.0 / 3.0)  # broad-crested weir
    return max(H0 - head, 0.0)


def set_gate(crest):
    z = domain.quantities['elevation'].centroid_values.copy()
    z[gate_cells] = crest
    domain.set_quantity('elevation', z, location='centroids')
    base = domain.sediment_z_base.copy()
    base[gate_cells] = crest
    domain.set_erodible_base(elevation=base)


beds = []
for event in range(EVENTS):
    t_start = event * T_EVENT
    release_reservoir()
    if GATE_TIME > 0.0:
        set_gate(gate_crest(0.0))
        n_gate = int(round(GATE_TIME / GATE_DT))
        for k in range(1, n_gate + 1):
            for t in domain.evolve(yieldstep=GATE_DT, finaltime=t_start + k * GATE_DT):
                pass
            set_gate(gate_crest(k * GATE_DT))
    for t in domain.evolve(yieldstep=YIELD, finaltime=(event + 1) * T_EVENT):
        if myid == 0 and verbose:
            print(domain.timestepping_statistics())
    beds.append(domain.quantities['elevation'].centroid_values.copy())
    if myid == 0 and verbose:
        dz = beds[-1] - z0
        print('event %d: bed change min %.1f mm max %.1f mm' % (event + 1, 1000 * dz.min(), 1000 * dz.max()))

domain.sww_merge(delete_old=True)
if myid == 0:
    np.savez('%s_bed.npz' % output_file, x=domain.centroid_coordinates[:, 0],
             y=domain.centroid_coordinates[:, 1], z0=z0, beds=np.array(beds),
             areas=domain.areas)
    with open('%s_parameters.json' % output_file, 'w') as f:
        json.dump({'case': CASE, 'events': EVENTS, 't_event': T_EVENT, 'alg': alg,
                   'porosity': POROSITY, 'd_sand': D_SAND, 'n_pvc': N_PVC, 'n_sand': N_SAND,
                   'repose': REPOSE, 'triangles': int(domain.number_of_elements),
                   'gate_time': GATE_TIME, 'a_gate': A_GATE, 'pier_amp': PIER_AMP, 'pier_len': PIER_LEN, 'pier_front': PIER_FRONT}, f, indent=1)
finalize()
