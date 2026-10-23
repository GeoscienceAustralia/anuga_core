#!/usr/bin/env python3
"""
ANUGA GPU Weak Scaling Benchmark
----------------------------------
Measures wall time for MPI+GPU runs where the mesh grows proportionally with
the number of ranks, keeping elements-per-rank constant (weak scaling).

Must be launched via mpiexec/srun with MPI.  All ranks execute this script;
rank 0 writes the JSON result file.

Usage
-----
Single-node sweep (CPU_ONLY_MODE, testing structure):

    for N in 1 2 4; do
        mpiexec -np $N python benchmarks/run_weak_scaling.py --base-nx 230
    done

HPC GPU run (see scripts/hpc/submit_weak_scaling.sh):

    OMP_NUM_THREADS=1 OMP_TARGET_OFFLOAD=mandatory \\
        srun --ntasks=$SLURM_NTASKS python benchmarks/run_weak_scaling.py \\
        --base-nx 500 --mode 2

Analyse collected results:

    python benchmarks/run_weak_scaling.py --analyse benchmarks/results/weak_scaling/

Mesh scaling (constant cell size — CFL-invariant)
-------------------------------------------------
Both the mesh resolution AND the physical domain extent are scaled by
sqrt(N) so the triangle edge length (base_len / base_nx) is invariant:

  nx_total  = base_nx  × sqrt(N)     ny_total  = base_ny  × sqrt(N)
  len1      = base_len × sqrt(N)     len2      = base_len × sqrt(N)
  cell_size = base_len / base_nx     (constant — CFL dt is invariant)
  tris/rank ≈ 2 × base_nx × base_ny (constant)
  n_steps   = finaltime / CFL_dt    (constant — equal work per rank)

Default base_nx = base_ny = 500, base_len = 1000 m  →  ~500 K tris/rank,
cell size ≈ 2 m.
"""

import argparse
import glob
import json
import math
import os
import platform
import subprocess
import sys
import tempfile
import time
from datetime import datetime


def _git_info(repo_root):
    def _run(args):
        try:
            return subprocess.check_output(
                args, cwd=repo_root,
                stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return 'unknown'
    return {
        'commit': _run(['git', 'rev-parse', '--short', 'HEAD']),
        'branch': _run(['git', 'rev-parse', '--abbrev-ref', 'HEAD']),
    }


def _env_info():
    try:
        import anuga
        anuga_version = anuga.__version__
    except Exception:
        anuga_version = 'unknown'
    return {
        'python_version': sys.version.split()[0],
        'platform': platform.system(),
        'hostname': platform.node().split('.')[0],
        'omp_num_threads': os.environ.get('OMP_NUM_THREADS', 'unset'),
        'omp_target_offload': os.environ.get('OMP_TARGET_OFFLOAD', 'unset'),
        'anuga_version': anuga_version,
    }


# ---------------------------------------------------------------------------
# Single weak-scaling run (all ranks call this)
# ---------------------------------------------------------------------------

def run_weak_scaling(base_nx, base_ny, base_len, mode, finaltime, yieldstep, verbose):
    """
    Run one weak-scaling experiment.

    Both the mesh resolution (nx, ny) AND the physical domain extent (len1,
    len2) are scaled by sqrt(nranks) so that the triangle edge length —
    and therefore the CFL-limited inner timestep — remains constant across
    all rank counts.  Only the number of triangles per rank stays fixed.

    Scaling invariants (all independent of nranks):
      cell size   = base_len / base_nx  (constant)
      CFL dt      ∝ cell size           (constant)
      n_steps     = finaltime / CFL_dt  (constant)
      tris/rank   ≈ 2 * base_nx * base_ny (constant)

    Returns a result dict on rank 0; None on other ranks.
    """
    from mpi4py import MPI
    import numpy as np
    import anuga
    from anuga.parallel.parallel_api import distribute

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.Get_size()

    # Scale mesh AND physical extent together so cell size is constant.
    # cell_size = len1_total / nx_total = base_len / base_nx  (invariant)
    scale = math.sqrt(nranks)
    nx_total = max(1, round(base_nx * scale))
    ny_total = max(1, round(base_ny * scale))
    len1_total = base_len * scale
    len2_total = base_len * scale
    cell_size = base_len / base_nx   # invariant across N
    dam_x = len1_total / 2.0         # dam at physical centre

    if rank == 0 and verbose:
        print(f'[weak_scaling] nranks={nranks}  nx={nx_total}  ny={ny_total}  '
              f'len={len1_total:.0f}m  cell_size={cell_size:.2f}m  '
              f'total_tris≈{2*nx_total*ny_total:,}  '
              f'tris/rank≈{2*nx_total*ny_total//nranks:,}')
        sys.stdout.flush()

    # --- Create and distribute domain ---
    tmpdir = tempfile.mkdtemp()
    try:
        comm.Barrier()
        if rank == 0:
            domain = anuga.rectangular_cross_domain(
                nx_total, ny_total, len1=len1_total, len2=len2_total)
            domain.set_flow_algorithm('DE0')
            domain.set_low_froude(0)
            domain.set_name('weak_scaling')
            domain.set_datadir(tmpdir)
            domain.store = False
            domain.set_quantity('elevation', 0.0)
            # Dam at physical centre; capture dam_x in default arg to avoid
            # closure over a variable that doesn't change, but be explicit.
            domain.set_quantity('stage',
                                lambda x, y, d=dam_x: np.where(x < d, 2.0, 0.5))
            domain.set_quantity('xmomentum', 0.0)
            domain.set_quantity('ymomentum', 0.0)
            domain.set_boundary(
                {t: anuga.Reflective_boundary(domain)
                 for t in domain.get_boundary_tags()})
        else:
            domain = None

        comm.Barrier()
        t_dist0 = MPI.Wtime()
        domain = distribute(domain, verbose=False)
        t_dist1 = MPI.Wtime()
        dist_time_s = t_dist1 - t_dist0

        domain.set_flow_algorithm('DE0')
        domain.set_low_froude(0)
        domain.store = False

        if mode >= 1:
            domain.set_multiprocessor_mode(mode)

        local_n_tris = domain.number_of_triangles

        # --- Timed evolve loop ---
        comm.Barrier()
        t0 = MPI.Wtime()
        for _ in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
            pass
        comm.Barrier()
        wall_time_s = MPI.Wtime() - t0

        n_steps = domain.number_of_steps

        # Collect timing stats across ranks
        all_wall_times = comm.gather(wall_time_s, root=0)
        all_local_tris = comm.gather(local_n_tris, root=0)

        if rank == 0:
            wall_max = max(all_wall_times)
            wall_min = min(all_wall_times)
            wall_mean = sum(all_wall_times) / len(all_wall_times)
            total_tris = sum(all_local_tris)
            tris_per_rank_mean = total_tris / nranks
            # cells/s = total_tris * n_steps / wall_time (max across ranks,
            # since the barrier makes all ranks wait for the slowest)
            cells_per_s = total_tris * n_steps / wall_max if wall_max > 0 else 0.0

            return {
                'nranks': nranks,
                'base_nx': base_nx,
                'base_ny': base_ny,
                'base_len': base_len,
                'nx_total': nx_total,
                'ny_total': ny_total,
                'len1_total': round(len1_total, 2),
                'cell_size': round(cell_size, 4),
                'total_triangles': total_tris,
                'tris_per_rank_mean': round(tris_per_rank_mean),
                'mode': mode,
                'finaltime': finaltime,
                'n_steps': n_steps,
                'wall_time_max_s': round(wall_max, 4),
                'wall_time_min_s': round(wall_min, 4),
                'wall_time_mean_s': round(wall_mean, 4),
                'cells_per_s': round(cells_per_s),
                'dist_time_s': round(dist_time_s, 3),
            }
        return None

    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Analysis mode: load multiple JSON files and print efficiency table
# ---------------------------------------------------------------------------

def analyse(results_dir):
    """
    Read all weak_scaling_*.json files from results_dir and print a
    parallel efficiency table.

    Parallel efficiency (weak scaling):
        eta(N) = T(1) / T(N)  × 100 %

    Perfect weak scaling: eta = 100% for all N.
    Target for SC26: eta >= 80%.
    """
    pattern = os.path.join(results_dir, 'weak_scaling_*.json')
    paths = sorted(glob.glob(pattern))
    if not paths:
        print(f'No weak_scaling_*.json files found in {results_dir}')
        return

    records = []
    for p in paths:
        with open(p) as fh:
            data = json.load(fh)
        for r in data.get('results', []):
            records.append(r)

    if not records:
        print('No results found in JSON files.')
        return

    # Sort by nranks
    records.sort(key=lambda r: r['nranks'])

    # Find T(1) reference
    ref1 = next((r for r in records if r['nranks'] == 1), None)
    t1 = ref1['wall_time_max_s'] if ref1 else None

    header = (f"{'Ranks':>6}  {'Tris/rank':>10}  {'Total tris':>12}  "
              f"{'cell_m':>7}  {'steps':>6}  "
              f"{'wall(s)':>8}  {'min(s)':>7}  {'cells/s':>12}  "
              f"{'efficiency':>11}  {'dist(s)':>7}")
    print()
    print('Weak scaling results')
    if t1 is not None:
        print(f'  Reference T(1) = {t1:.4f} s')
    print('  cell_m = triangle edge length (should be constant across ranks)')
    print('  steps  = inner CFL timesteps (should be constant across ranks)')
    print()
    print(header)
    print('-' * len(header))

    for r in records:
        t_wall = r['wall_time_max_s']
        eta = (t1 / t_wall * 100.0) if (t1 and t_wall > 0) else float('nan')
        cell_size = r.get('cell_size', float('nan'))
        n_steps = r.get('n_steps', 0)
        print(f"  {r['nranks']:>4}  {r['tris_per_rank_mean']:>10,}  "
              f"{r['total_triangles']:>12,}  "
              f"{cell_size:>7.2f}  {n_steps:>6}  "
              f"{t_wall:>8.3f}  {r['wall_time_min_s']:>7.3f}  "
              f"{r['cells_per_s']:>12,.0f}  "
              f"{'--' if math.isnan(eta) else f'{eta:.1f}%':>11}  "
              f"{r['dist_time_s']:>7.2f}")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='ANUGA GPU weak scaling benchmark (run via mpiexec).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--base-nx', type=int, default=500,
        help='Base nx per rank (default 500, giving ~500K tris/rank)',
    )
    parser.add_argument(
        '--base-ny', type=int, default=None,
        help='Base ny per rank (default: same as --base-nx)',
    )
    parser.add_argument(
        '--base-len', type=float, default=1000.0,
        help='Physical domain length per rank in metres (default 1000.0).  '
             'Scaled by sqrt(nranks) so cell size = base_len/base_nx stays '
             'constant and the CFL timestep is invariant with rank count.',
    )
    parser.add_argument(
        '--mode', type=int, default=2,
        help='multiprocessor_mode: 1=CPU OpenMP, 2=GPU/C (default 2)',
    )
    parser.add_argument(
        '--finaltime', type=float, default=10.0,
        help='Simulation end time in seconds (default 10.0)',
    )
    parser.add_argument(
        '--yieldstep', type=float, default=5.0,
        help='Yieldstep in seconds (default 5.0)',
    )
    parser.add_argument(
        '--output', default=None,
        help='JSON output path (default: benchmarks/results/weak_scaling/weak_scaling_NP<N>_*.json)',
    )
    parser.add_argument(
        '--analyse', metavar='DIR', default=None,
        help='Analysis mode: read all weak_scaling_*.json from DIR and print efficiency table',
    )
    parser.add_argument(
        '--verbose', action='store_true',
        help='Print progress messages from rank 0',
    )
    args = parser.parse_args()

    # --- Analysis mode (single-process, no MPI needed) ---
    if args.analyse:
        analyse(args.analyse)
        return 0

    # --- Benchmark mode (requires MPI) ---
    try:
        from mpi4py import MPI
    except ImportError:
        print('ERROR: mpi4py is required for benchmark mode.  '
              'Install with: conda install mpi4py', file=sys.stderr)
        return 1

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.Get_size()

    base_ny = args.base_ny if args.base_ny is not None else args.base_nx
    base_len = args.base_len

    if rank == 0:
        cell_size = base_len / args.base_nx
        print(f'ANUGA weak scaling benchmark  nranks={nranks}  mode={args.mode}')
        print(f'  base_nx={args.base_nx}  base_ny={base_ny}  '
              f'base_len={base_len:.0f}m  cell_size={cell_size:.2f}m  '
              f'finaltime={args.finaltime}s  yieldstep={args.yieldstep}s')
        env = _env_info()
        print(f'  OMP_NUM_THREADS={env["omp_num_threads"]}  '
              f'OMP_TARGET_OFFLOAD={env["omp_target_offload"]}')
        sys.stdout.flush()

    result = run_weak_scaling(
        base_nx=args.base_nx,
        base_ny=base_ny,
        base_len=base_len,
        mode=args.mode,
        finaltime=args.finaltime,
        yieldstep=args.yieldstep,
        verbose=args.verbose or (rank == 0),
    )

    # Rank 0 writes output
    if rank == 0 and result is not None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        git = _git_info(repo_root)
        env = _env_info()

        if args.output:
            outpath = args.output
        else:
            outdir = os.path.join(repo_root, 'benchmarks', 'results',
                                  'weak_scaling')
            os.makedirs(outdir, exist_ok=True)
            branch_slug = git['branch'].replace('/', '_')
            outpath = os.path.join(
                outdir,
                f"weak_scaling_NP{nranks:04d}_{branch_slug}_{git['commit']}_{timestamp}.json")

        payload = {
            'timestamp': timestamp,
            'git': git,
            'env': env,
            'results': [result],
        }
        with open(outpath, 'w') as fh:
            json.dump(payload, fh, indent=2)

        # Print summary line
        t_max = result['wall_time_max_s']
        t_min = result['wall_time_min_s']
        c_s = result['cells_per_s']
        print(f'\n  NP={nranks:3d}  tris/rank={result["tris_per_rank_mean"]:>10,}  '
              f'total={result["total_triangles"]:>12,}  '
              f'wall_max={t_max:.3f}s  wall_min={t_min:.3f}s  '
              f'cells/s={c_s:,.0f}')
        print(f'  Saved → {outpath}')
        sys.stdout.flush()

    return 0


if __name__ == '__main__':
    sys.exit(main())
