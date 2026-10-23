"""run_smpl_rectangular_create_dump  (smpl = Sequential Mesh, Parallel Load)

Offline preprocessing step for the two-phase parallel workflow:

  Phase 1 (this script) — run once on any machine with enough RAM:

      python run_smpl_rectangular_create_dump.py -np N [options]

  Phase 2 — run on the HPC cluster as many times as needed:

      mpiexec -np N python -u run_smpl_rectangular_load_evolve.py [options]

Difference from the ``sdpl`` (sequential_distribute) workflow
-------------------------------------------------------------
``sequential_mesh_dump`` saves only mesh topology and halo structure —
no quantity data.  The partition files are therefore much smaller and do
not need to be regenerated when initial conditions change.  Each MPI rank
sets its own quantities independently after loading the mesh.

This is the recommended approach when:
  * the full mesh + quantities would exceed available RAM on a single node, OR
  * you need to run many ensemble/scenario variants on the same mesh.

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

parser = argparse.ArgumentParser(description='Create and dump mesh partitions (no quantities)')
parser.add_argument('-np', '--numprocs',  type=int,   default=4,
                    help='Number of partitions to create')
parser.add_argument('-sn', '--sqrtN',     type=int,   default=100,
                    help='sqrt(N) grid cells per side  (100 → 80 000 triangles)')
parser.add_argument('-gl', '--ghost_layer', type=int, default=2,
                    help='Ghost layer width (default 2)')
parser.add_argument('-ps', '--partition_scheme', type=str, default='metis',
                    choices=['metis', 'morton', 'hilbert'],
                    help='METIS / Morton / Hilbert partitioning')
parser.add_argument('-pd', '--partition_dir', type=str, default='Partitions',
                    help='Directory for partition files')
parser.add_argument('-v',  '--verbose',   action='store_true',
                    help='Verbose output')
args = parser.parse_args()

sqrtN         = args.sqrtN
ncpus         = args.numprocs
partition_dir = args.partition_dir
verbose       = args.verbose
length = width = 2.0

domain_name = f'rect_smpl_gl{args.ghost_layer}_sn{sqrtN}_np{ncpus}'

dist_params = {
    'ghost_layer_width': args.ghost_layer,
    'partition_scheme':  args.partition_scheme,
}

print(f'Creating mesh: {2 * sqrtN**2} triangles, partitioning into {ncpus} ranks')
print(f'Partition scheme : {args.partition_scheme}')
print(f'Output directory : {partition_dir}')
sys.stdout.flush()

# ---------------------------------------------------------------------------
# Build mesh (topology only — quantities are NOT needed for mesh dump)
# ---------------------------------------------------------------------------

t0 = time.time()

domain = rectangular_cross_domain(sqrtN, sqrtN,
                                  len1=length, len2=width,
                                  origin=(-length / 2, -width / 2),
                                  verbose=verbose)
domain.set_name(domain_name)

print(f'Mesh created in {time.time() - t0:.2f} s  '
      f'({domain.number_of_triangles} triangles, '
      f'{domain.number_of_nodes} nodes)')

# ---------------------------------------------------------------------------
# Partition and write one .nc file per rank  (mesh + halos, no quantities)
# ---------------------------------------------------------------------------

t1 = time.time()

anuga.sequential_mesh_dump(domain,
                           numprocs=ncpus,
                           partition_dir=partition_dir,
                           verbose=verbose,
                           parameters=dist_params)

print(f'Mesh partitioned and written in {time.time() - t1:.2f} s')
print(f'Files: {partition_dir}/{domain_name}_mesh_P{ncpus}_<rank>.nc')
