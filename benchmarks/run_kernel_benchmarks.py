#!/usr/bin/env python
"""
Per-kernel microbenchmarks for the shallow-water step (audit issue #337).

run_benchmarks.py measures whole-evolve throughput, which is the right headline
number but cannot say *which* kernel regressed. This script times the kernels
of one timestep individually, on a fixed dam-break mesh, in either compute
mode, and reports microseconds per call and cells per second per kernel.

    python benchmarks/run_kernel_benchmarks.py                    # legacy + unified, 40k tris
    python benchmarks/run_kernel_benchmarks.py --modes legacy     # one mode
    python benchmarks/run_kernel_benchmarks.py --nx 300 --repeats 50
    python benchmarks/run_kernel_benchmarks.py --threads 8     # OpenMP threads
    python benchmarks/run_kernel_benchmarks.py --output /tmp/k.json
    python benchmarks/run_kernel_benchmarks.py --compare a.json b.json

Each kernel runs on the state left by the previous one, in the order a real
step applies them, and the whole sequence is repeated `--repeats` times after
a warm-up pass, so timings include realistic cache state but not JIT/first-
touch effects. The mesh does not evolve between repeats (no timestep is
applied), so every repeat sees the same wet/dry pattern.

On a GPU-offload build the 'unified' mode runs the kernels on the device;
each timing then includes the kernel launch and a device synchronisation, as
it does in a real evolve.

Threads: --threads N sets the process-wide OpenMP thread count before the
domain is built (anuga.set_omp_num_threads). Without it the count comes from
OMP_NUM_THREADS, and from ANUGA's default of 1 when that is unset too, so an
unadorned run is single-threaded rather than "all cores". The count actually
applied is printed in the table header and recorded in the JSON as
omp_num_threads, so two result files are comparable after the fact.

Noise: the cheap kernels (protect, backup, saxpy) take a few microseconds at
the default size and vary by 5-10% between runs. Compare with a larger --nx
and more --repeats before reading anything into a small change; the
--compare mode flags 10% by default.
"""

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone


# ---------------------------------------------------------------------------
# Domain
# ---------------------------------------------------------------------------

def _create_domain(nx, ny, mode, tmpdir):
    """Dam-break on a rectangular cross mesh; same setup as run_benchmarks.py."""
    import numpy as np
    from anuga import rectangular_cross_domain, Reflective_boundary

    domain = rectangular_cross_domain(nx, ny, len1=1000.0, len2=1000.0)
    domain.set_flow_algorithm('DE0')
    domain.set_low_froude(0)
    domain.set_name('kernel_bench')
    domain.set_datadir(tmpdir)
    domain.store = False
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('friction', 0.03)
    domain.set_quantity('stage', lambda x, y: np.where(x < 500.0, 2.0, 0.5))
    domain.set_boundary({t: Reflective_boundary(domain)
                         for t in domain.get_boundary_tags()})
    domain.set_compute_mode(mode)
    # Make sure the C struct / device arrays exist before timing anything.
    domain.update_boundary()
    if mode == 'unified':
        domain._ensure_gpu_interface()
    return domain


# ---------------------------------------------------------------------------
# Kernel sequences
# ---------------------------------------------------------------------------

def _legacy_kernels(domain):
    """(name, callable) pairs for mode 1, in the order of one Euler step."""
    from anuga.shallow_water import sw_domain_openmp_ext as ext
    dt = 0.1
    return [
        ('extrapolate_second_order_edge_sw',
         lambda: ext.extrapolate_second_order_edge_sw(domain)),
        ('distribute_edges_to_vertices',
         lambda: ext.distribute_edges_to_vertices(domain)),
        ('compute_fluxes_ext_central',
         lambda: ext.compute_fluxes_ext_central(domain, domain.evolve_max_timestep)),
        ('manning_friction_flat_semi_implicit',
         lambda: ext.manning_friction_flat_semi_implicit(domain)),
        ('protect_new',
         lambda: ext.protect_new(domain)),
        ('update_conserved_quantities',
         lambda: ext.update_conserved_quantities(domain, dt)),
        ('backup_conserved_quantities',
         lambda: ext.backup_conserved_quantities(domain)),
        ('saxpy_conserved_quantities',
         lambda: ext.saxpy_conserved_quantities(domain, 0.5, 0.5, 1.0)),
    ]


def _unified_kernels(domain):
    """(name, callable) pairs for mode 2 through the GPU/CPU-unified interface."""
    gi = domain.gpu_interface
    dt = 0.1
    return [
        ('extrapolate_second_order_edge_sw',
         lambda: gi.extrapolate_second_order_edge_sw_kernel(domain)),
        ('compute_fluxes_ext_central',
         lambda: gi.compute_fluxes_ext_central_kernel(domain, domain.evolve_max_timestep)),
        ('manning_friction',
         lambda: gi.manning_friction_kernel(domain)),
        ('protect',
         lambda: gi.protect_against_infinitesimal_and_negative_heights_kernel(domain)),
        ('update_conserved_quantities',
         lambda: gi.update_conserved_quantities_kernel(domain, dt)),
        ('backup_conserved_quantities',
         lambda: gi.backup_conserved_quantities_kernel(domain)),
        ('saxpy_conserved_quantities',
         lambda: gi.saxpy_conserved_quantities_kernel(domain, 0.5, 0.5)),
    ]


# ---------------------------------------------------------------------------
# Timing
# ---------------------------------------------------------------------------

def _time_kernels(kernels, repeats, warmup=2):
    """Run the sequence `repeats` times; return {name: [seconds per call]}."""
    samples = {name: [] for name, _ in kernels}
    for r in range(warmup + repeats):
        for name, fn in kernels:
            t0 = time.perf_counter()
            fn()
            dt = time.perf_counter() - t0
            if r >= warmup:
                samples[name].append(dt)
    return samples


def run_mode(mode, nx, ny, repeats):
    tmpdir = tempfile.mkdtemp()
    try:
        domain = _create_domain(nx, ny, mode, tmpdir)
        n_tris = domain.number_of_triangles
        kernels = _legacy_kernels(domain) if mode == 'legacy' else _unified_kernels(domain)
        samples = _time_kernels(kernels, repeats)
    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)

    rows = []
    for name, _ in kernels:
        s = samples[name]
        med = statistics.median(s)
        rows.append({
            'mode': mode,
            'kernel': name,
            'n_triangles': n_tris,
            'repeats': repeats,
            'median_us': round(med * 1e6, 1),
            'min_us': round(min(s) * 1e6, 1),
            'cells_per_s': round(n_tris / med, 0) if med > 0 else 0.0,
        })
    return rows


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _print_table(rows, threads):
    print(f"OpenMP threads: {threads}")
    header = f"{'mode':<8} {'kernel':<38} {'tris':>8} {'median us':>11} {'min us':>10} {'cells/s':>13}"
    print(header)
    print('-' * len(header))
    for r in rows:
        print(f"{r['mode']:<8} {r['kernel']:<38} {r['n_triangles']:>8} "
              f"{r['median_us']:>11.1f} {r['min_us']:>10.1f} {r['cells_per_s']:>13,.0f}")


def _compare(path_a, path_b, threshold=0.05):
    a = json.load(open(path_a))
    b = json.load(open(path_b))
    key = lambda r: (r['mode'], r['kernel'])
    ra = {key(r): r for r in a['results']}
    rb = {key(r): r for r in b['results']}
    ta, tb = a.get('omp_num_threads'), b.get('omp_num_threads')
    if ta != tb:
        print(f"NOTE: thread counts differ (before {ta}, after {tb}); the comparison is not like for like")
    header = f"{'mode':<8} {'kernel':<38} {'before us':>10} {'after us':>10} {'change':>8}"
    print(header)
    print('-' * len(header))
    worst = 0.0
    for k in sorted(set(ra) & set(rb)):
        x, y = ra[k]['median_us'], rb[k]['median_us']
        change = (y - x) / x if x else 0.0
        worst = max(worst, change)
        flag = '  <-- slower' if change > threshold else ('  faster' if change < -threshold else '')
        print(f"{k[0]:<8} {k[1]:<38} {x:>10.1f} {y:>10.1f} {change:>+7.1%}{flag}")
    return worst


def _git_info():
    try:
        sha = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], text=True).strip()
        branch = subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'], text=True).strip()
        return {'sha': sha, 'branch': branch}
    except Exception:
        return {}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--modes', default='legacy,unified',
                        help="comma-separated: legacy, unified (default both)")
    parser.add_argument('--nx', type=int, default=100, help='cells per side (nx*ny*4 triangles)')
    parser.add_argument('--ny', type=int, default=None)
    parser.add_argument('--repeats', type=int, default=20)
    parser.add_argument('--threads', type=int, default=None,
                        help='OpenMP thread count (default: OMP_NUM_THREADS, else 1)')
    parser.add_argument('--output', help='write results as JSON')
    parser.add_argument('--compare', nargs=2, metavar=('BEFORE', 'AFTER'),
                        help='compare two JSON result files instead of running')
    parser.add_argument('--threshold', type=float, default=0.10,
                        help='fractional slowdown flagged by --compare (default 0.10; the '
                             'few-microsecond kernels are noisy at the 5%% level, so use '
                             '--repeats 50 and a larger --nx before trusting a small change)')
    args = parser.parse_args()

    if args.compare:
        worst = _compare(*args.compare, threshold=args.threshold)
        sys.exit(1 if worst > args.threshold else 0)

    import anuga
    # Apply the thread count to the OpenMP runtime before any domain exists.
    # None means "OMP_NUM_THREADS, else 1", which is also what a domain would
    # get; doing it here makes the value we record the value that was used.
    threads = anuga.set_omp_num_threads(args.threads, verbose=False)

    ny = args.ny or args.nx
    rows = []
    for mode in [m.strip() for m in args.modes.split(',') if m.strip()]:
        rows.extend(run_mode(mode, args.nx, ny, args.repeats))

    _print_table(rows, threads)

    if args.output:
        record = {
            'timestamp': datetime.now(timezone.utc).isoformat(),
            'anuga_version': anuga.__version__,
            'git': _git_info(),
            'platform': platform.platform(),
            'omp_num_threads': threads,
            'gpu_offload': anuga.gpu_offload_enabled(),
            'nx': args.nx, 'ny': ny, 'repeats': args.repeats,
            'results': rows,
        }
        with open(args.output, 'w') as f:
            json.dump(record, f, indent=2)
        print(f"\nWrote {args.output}")


if __name__ == '__main__':
    main()
