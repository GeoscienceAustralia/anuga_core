"""A settling basin with flow.

Uniform steady flow along a flat frictionless channel, fixed at both ends
by Dirichlet boundaries carrying the same state, so the flow is exact.
Sediment-laden water enters at the inflow (a tracer boundary concentration)
with entrainment switched off, and the suspended load settles along the
channel at the well-mixed deposition rate; the bed rises where it lands.
The steady concentration profile is a closed-form exponential in the
distance from the inflow, and the bed-rise rate follows from it: see
analytical_settling_basin.py. The run checks the tracer advection scheme
with the inflow boundary concentration, the deposition term [D-1] in the
sediment mass balance [G-3], the bed update [G-4] under deposition, and
conservation of sediment between the boundaries, the water column and the
bed.
"""
import json
import numpy as np
import anuga
from anuga import Domain, myid, finalize, distribute

output_dir = '.'
output_file = 'settling_basin'

L, W = 300.0, 10.0        # channel length and width, m
dx = 2.5
h0 = 1.0                  # depth, m
u0 = 0.5                  # velocity, m/s; q = u0 h0 = 0.5 m^2/s
c0 = 2.0e-4               # inflow concentration (volumetric), small enough
                          # that the bed rise stays a fraction of a percent
                          # of the depth over the run
diameter = 1.0e-4         # 100 um: v_s ~ 8 mm/s, settling length q/v_s ~ 62 m,
                          # so the load falls by e^-4.8 along the channel
d_star = 1.0
porosity = 0.30
finaltime = 1500.0        # steady after L/u0 = 600 s; the last 100 s give
yieldstep = 50.0          # the bed-rise rate

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
    domain.set_quantity('xmomentum', u0 * h0)
else:
    domain = None

domain = distribute(domain)

# The uniform flow is an exact steady state of the flat frictionless
# channel, so the same state at both ends holds it.
Bd = anuga.Dirichlet_boundary([h0, u0 * h0, 0.0])
Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

# The sediment: one fraction, no initial load, no entrainment (critical
# Shields stress zero), well-mixed deposition, and the inflow carries c0.
domain.initialize_sediment_operator(porosity=porosity, bed_evolution=True)
domain.set_deposition(law='d_star', near_bed='constant')
domain.add_sediment_fraction('silt', diameter=diameter, d_star=d_star,
                             tau_c_star=0.0, initial_concentration=0.0)
domain.set_tracer_boundary('silt', 'left', c0)

if myid == 0:
    v_s = float(domain.sediment_settling_velocity[0])
    if verbose:
        print(domain.sediment_summary())
        print('settling velocity %.4g m/s, settling length %.1f m'
              % (v_s, u0 * h0 / (d_star * v_s)))

from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)

z0 = domain.quantities['elevation'].centroid_values.copy()
for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
    if myid == 0 and verbose:
        print(domain.timestepping_statistics())

# Sediment budget at the end: boundary flux in, water column, bed. The
# conservation check is collective, so every rank calls it.
change, flux, discrepancy = domain.check_tracer_conservation('silt')
bed_volume = (1.0 - porosity) * float(
    ((domain.quantities['elevation'].centroid_values - z0)
     * domain.areas)[domain.tri_full_flag == 1].sum())
if domain.numproc > 1:
    from mpi4py import MPI
    bed_volume = MPI.COMM_WORLD.allreduce(bed_volume, op=MPI.SUM)

if myid == 0:
    with open('settling_basin_parameters.json', 'w') as f:
        json.dump({'L': L, 'W': W, 'h0': h0, 'u0': u0, 'q': u0 * h0, 'c0': c0,
                   'v_s': v_s, 'd_star': d_star, 'porosity': porosity,
                   'finaltime': finaltime, 'yieldstep': yieldstep,
                   'water_column_change': float(change),
                   'boundary_flux_in': float(flux),
                   'bed_sediment_volume': float(bed_volume)}, f, indent=1)

domain.sww_merge(delete_old=True)
finalize()
