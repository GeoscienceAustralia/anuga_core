"""run_smpl_rectangular_load_evolve  (smpl = Sequential Mesh, Parallel Load)

Phase 2 of the two-phase parallel workflow.  Run after the partition files
have been created by run_smpl_rectangular_create_dump.py:

    mpiexec -np N python -u run_smpl_rectangular_load_evolve.py [options]

Each rank loads its own partition file from disk, then sets quantities
independently (no data is transferred between ranks for initial conditions).
Boundary conditions and operators are attached in the usual way.

See run_smpl_rectangular_create_dump.py for the full workflow description.

Authors: Steve Roberts, 2026
"""

import argparse
import sys
import time

import anuga
from anuga import (
    Transmissive_boundary, Reflective_boundary,
    myid, numprocs, finalize, barrier,
    Set_stage,
)


# ---------------------------------------------------------------------------
# Command-line arguments  (must match create-dump script for domain_name)
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Load mesh partitions and evolve in parallel')
parser.add_argument('-sn', '--sqrtN',     type=int,   default=100,
                    help='Must match value used in create_dump script')
parser.add_argument('-gl', '--ghost_layer', type=int, default=2,
                    help='Must match value used in create_dump script')
parser.add_argument('-ys', '--yieldstep', type=float, default=0.005,
                    help='Output timestep (s)')
parser.add_argument('-ft', '--finaltime', type=float, default=0.015,
                    help='Final simulation time (s)')
parser.add_argument('-pd', '--partition_dir', type=str, default='Partitions',
                    help='Directory containing partition files')
parser.add_argument('-n',  '--name',       type=str, default=None,
                    help='Partition file base name (default: derived from -sn/-gl/-np)')
parser.add_argument('-fdt', '--fixed_dt',  type=float, default=0.0,
                    help='Fixed flux timestep (0 = adaptive CFL)')
parser.add_argument('-v',  '--verbose',    action='store_true',
                    help='Verbose output')
parser.add_argument('-ve', '--evolve_verbose', action='store_true',
                    help='Per-rank evolve timing')
parser.add_argument('-sww', '--store_sww', action='store_true',
                    help='Write SWW output file')
args = parser.parse_args()

sqrtN         = args.sqrtN
yieldstep     = args.yieldstep
finaltime     = args.finaltime
partition_dir = args.partition_dir
verbose       = args.verbose
fixed_flux_timestep = args.fixed_dt if args.fixed_dt != 0.0 else None

if args.name is not None:
    domain_name = args.name
else:
    domain_name = f'rect_smpl_gl{args.ghost_layer}_sn{sqrtN}_np{numprocs}'

if myid == 0:
    print(args)
    sys.stdout.flush()

# ---------------------------------------------------------------------------
# Load this rank's mesh partition
# ---------------------------------------------------------------------------

barrier()
t0 = time.time()

domain = anuga.sequential_mesh_load(name=domain_name,
                                    partition_dir=partition_dir,
                                    verbose=verbose)

barrier()
if myid == 0:
    print(f'Mesh loaded in {time.time() - t0:.2f} s  '
          f'({domain.number_of_global_triangles} global triangles, '
          f'{numprocs} ranks)')

# ---------------------------------------------------------------------------
# Set quantities  — each rank evaluates functions on its own coordinates;
# no quantity data is read from the partition files or sent between ranks.
# ---------------------------------------------------------------------------

domain.set_name(domain_name)
domain.set_store(args.store_sww)
domain.set_flow_algorithm('DE0')
domain.set_cfl(1.0)
domain.set_fixed_flux_timestep(fixed_flux_timestep)

domain.set_quantity('elevation', lambda x, y: -1.0 - x)
domain.set_quantity('stage',     1.0)
domain.set_quantity('friction',  0.01)

# ---------------------------------------------------------------------------
# Boundary conditions
# ---------------------------------------------------------------------------

R = Reflective_boundary(domain)
T = Transmissive_boundary(domain)
domain.set_boundary({'left': R, 'right': R, 'bottom': R, 'top': R})

# ---------------------------------------------------------------------------
# Operators  (optional)
# ---------------------------------------------------------------------------

setter = Set_stage(domain, center=(0.0, 0.0), radius=0.5, stage=2.0)
setter()

# ---------------------------------------------------------------------------
# Evolve
# ---------------------------------------------------------------------------

barrier()
t_evolve = time.time()

for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
    if myid == 0:
        domain.write_time()
        sys.stdout.flush()

evolve_time = time.time() - t_evolve

if args.evolve_verbose:
    for p in range(numprocs):
        barrier()
        if myid == p:
            print(f'Rank {p}: evolve {evolve_time:.2f} s  '
                  f'comm {domain.communication_time:.2f} s  '
                  f'Allreduce {domain.communication_reduce_time:.2f} s')
            sys.stdout.flush()

# ---------------------------------------------------------------------------
# Merge per-rank SWW files  (rank 0 only)
# ---------------------------------------------------------------------------

if args.store_sww:
    domain.sww_merge(delete_old=True)

total_time = time.time() - domain.initial_walltime

if myid == 0:
    print(80 * '=')
    print('np, ntri, total_time, evolve_time')
    print(f'{numprocs}, {domain.number_of_global_triangles}, '
          f'{total_time:.3f}, {evolve_time:.3f}')

finalize()
