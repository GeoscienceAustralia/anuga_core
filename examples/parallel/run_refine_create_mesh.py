"""run_refine_create_mesh  —  coarse-partition + offline refinement

Preprocessing step for the refined-mesh parallel workflow:

  Phase 1 (this script) — run once on a workstation or login node:

      python run_refine_create_mesh.py -np N -rl 2 [options]

  Phase 2 — run on the HPC cluster with the same -np and -sn arguments:

      mpiexec -np N python -u run_smpl_rectangular_load_evolve.py -sn <sqrtN>

How it works
------------
:func:`anuga.create_parallel_mesh` partitions the *coarse* mesh (sqrtN × sqrtN
cells) into *numprocs* pieces and then refines each piece *refinement_levels*
times.  Because the refinement works partition-by-partition, the peak RAM on
this preprocessing machine is determined by the *coarse* mesh size, not the
final refined size.

Triangle counts
---------------
  coarse triangles  = 4 × sqrtN²   (rectangular_cross gives 4 per cell)
  refined triangles = 4^levels × coarse_triangles

  Examples (sqrtN=100, numprocs=8):
    levels=0 →   40 000 tri  (pure partitioning, no refinement)
    levels=1 →  160 000 tri
    levels=2 →  640 000 tri
    levels=3 → 2 560 000 tri

Authors: Steve Roberts, 2026
"""

import argparse
import sys
import time

import anuga
from anuga import rectangular_cross_domain


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description='Partition a coarse mesh and refine offline')
parser.add_argument('-np',  '--numprocs',    type=int,   default=4,
                    help='Number of MPI ranks / partition files to create')
parser.add_argument('-sn',  '--sqrtN',       type=int,   default=100,
                    help='Coarse grid: sqrtN cells per side  '
                         '(creates 2×sqrtN² coarse triangles)')
parser.add_argument('-rl',  '--refine_levels', type=int, default=1,
                    help='Refinement levels (each level 4× the triangles)')
parser.add_argument('-gl',  '--ghost_layer', type=int,   default=2,
                    help='Ghost layer width (default 2)')
parser.add_argument('-ps',  '--partition_scheme', type=str, default='metis',
                    choices=['metis', 'morton', 'hilbert'],
                    help='Partitioning algorithm')
parser.add_argument('-pd',  '--partition_dir', type=str, default='Partitions',
                    help='Output directory for partition files')
parser.add_argument('-nw',  '--num_workers',  type=int, default=1,
                    help='Worker processes for parallel partition refinement '
                         '(Pass 2 only; default 1 = single-process)')
parser.add_argument('-v',   '--verbose',     action='store_true',
                    help='Verbose output')
args = parser.parse_args()

sqrtN            = args.sqrtN
numprocs         = args.numprocs
refinement_levels = args.refine_levels
partition_dir    = args.partition_dir
num_workers      = args.num_workers
verbose          = args.verbose

coarse_tri = 4 * sqrtN ** 2
fine_tri   = (4 ** refinement_levels) * coarse_tri

print(f'Coarse mesh       : {sqrtN}×{sqrtN} cells → {coarse_tri:,} triangles')
print(f'Refinement levels : {refinement_levels}  (4^{refinement_levels} = {4**refinement_levels}×)')
print(f'Fine mesh         : {fine_tri:,} triangles')
print(f'Partitions        : {numprocs} ranks')
print(f'Partition scheme  : {args.partition_scheme}')
print(f'Output directory  : {partition_dir}')
print(f'Refine workers    : {num_workers}')
sys.stdout.flush()

# ---------------------------------------------------------------------------
# Build the coarse domain (topology only — quantities are not needed)
# ---------------------------------------------------------------------------

t0 = time.time()

coarse = rectangular_cross_domain(
    sqrtN, sqrtN,
    len1=2.0, len2=2.0,
    origin=(-1.0, -1.0),
    verbose=verbose,
)
coarse.set_name(f'rect_refined_sn{sqrtN}_rl{refinement_levels}_np{numprocs}')

print(f'\nCoarse mesh built in {time.time() - t0:.2f} s  '
      f'({coarse.number_of_triangles} triangles, '
      f'{coarse.number_of_nodes} nodes)')
sys.stdout.flush()

# ---------------------------------------------------------------------------
# Partition the coarse mesh and refine each partition offline
# ---------------------------------------------------------------------------

t1 = time.time()

dist_params = {
    'ghost_layer_width': args.ghost_layer,
    'partition_scheme':  args.partition_scheme,
}

anuga.create_parallel_mesh(
    coarse,
    numprocs=numprocs,
    refinement_levels=refinement_levels,
    partition_dir=partition_dir,
    verbose=verbose,
    parameters=dist_params,
    num_workers=num_workers,
)

elapsed = time.time() - t1
name = coarse.get_name()
print(f'\nDone in {elapsed:.2f} s')
print(f'Files  : {partition_dir}/{name}_mesh_P{numprocs}_<rank>.nc')
print(f'\nTo evolve in parallel:')
print(f'  mpiexec -np {numprocs} python -u run_smpl_rectangular_load_evolve.py '
      f'-n {name} -gl {args.ghost_layer} -pd {partition_dir}')
print(f'\nTip: use -nw N to refine {numprocs} partitions across N worker processes.')
