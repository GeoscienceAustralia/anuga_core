#!/usr/bin/env python3
"""
ANUGA GPU Performance Benchmark Suite
--------------------------------------
Measures wall time, cells/s, and GFLOP/s for GPU-scale meshes using the
mode=2 C/GPU RK loop compared with mode=1 (CPU OpenMP).

Uses the FLOP counter infrastructure built into sw_domain_gpu_ext
(gpu_flop.c / flop_counters_* Cython bindings).  When sw_domain_gpu_ext
is not available or has no real GPU target, the script still runs in
CPU_ONLY_MODE and reports GFLOP/s as measured by the C-side counters.

Results are written to a JSON file in benchmarks/results/ and are
compatible with compare_benchmarks.py (the 'name' key is the join key).

Usage
-----
    # Default: small and medium GPU sizes, modes 1 and 2
    python benchmarks/run_gpu_benchmarks.py

    # All sizes (large ~20 M tris — needs ~80 GB RAM or real GPU VRAM)
    python benchmarks/run_gpu_benchmarks.py --sizes small,medium,large

    # Only GPU mode
    python benchmarks/run_gpu_benchmarks.py --modes 2

    # Custom output file
    python benchmarks/run_gpu_benchmarks.py --output my_gpu_results.json

    # Print FLOP counter per-kernel breakdown
    python benchmarks/run_gpu_benchmarks.py --flop-detail

Mesh sizes
----------
    small   ~100 000 triangles (nx=ny=230)
    medium  ~ 2 000 000 triangles (nx=ny=1000)
    large   ~20 000 000 triangles (nx=ny=3162)

Metrics
-------
- wall_time_s      : wall-clock seconds for the evolve loop only
- n_steps          : number of internal (CFL) timesteps taken
- cells_per_s      : n_triangles * n_steps / wall_time_s
- gflops_per_second: from FLOP counters (0 if counters not available)
- total_gflops     : total floating-point operations (billions)
- setup_rss_mb     : RSS after domain creation (before any evolve)
- peak_rss_mb      : peak RSS sampled during evolve (100 ms polling)
- mb_per_ktri      : peak_rss_mb / (n_triangles / 1000)
- flop_detail      : per-kernel FLOP breakdown dict (if --flop-detail)
- gpu_available    : True if a real GPU offload target was detected
- cpu_only_mode    : True if running in CPU_ONLY_MODE (no real GPU)
"""

import argparse
import json
import os
import platform
import subprocess
import sys
import tempfile
import threading
import time
from datetime import datetime


# ---------------------------------------------------------------------------
# FLOP counter helpers (graceful fallback when extension unavailable)
# ---------------------------------------------------------------------------

def _import_flop_api():
    """
    Try to import FLOP counter functions from sw_domain_gpu_ext.
    Returns a dict of callables or None if unavailable.
    """
    try:
        from anuga.shallow_water.sw_domain_gpu_ext import (
            flop_counters_reset,
            flop_counters_enable,
            flop_counters_start_timer,
            flop_counters_stop_timer,
            flop_counters_get_stats,
        )
        return dict(
            reset=flop_counters_reset,
            enable=flop_counters_enable,
            start=flop_counters_start_timer,
            stop=flop_counters_stop_timer,
            stats=flop_counters_get_stats,
        )
    except (ImportError, AttributeError):
        return None


def _check_gpu_available():
    """Return (gpu_ext_present, gpu_device_available)."""
    try:
        from anuga.shallow_water.sw_domain_gpu_ext import gpu_available
        return True, bool(gpu_available())
    except (ImportError, AttributeError):
        return False, False


# ---------------------------------------------------------------------------
# Memory sampler (background thread) — identical to run_benchmarks.py
# ---------------------------------------------------------------------------

class _MemSampler:
    """Poll process RSS every `interval_s` seconds in a daemon thread."""

    def __init__(self, interval_s=0.1):
        self._interval = interval_s
        self.peak_mb = 0.0
        self._running = False
        self._thread = None

    def start(self):
        self._running = True
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self):
        self._running = False
        if self._thread:
            self._thread.join(timeout=2.0)

    def current_mb(self):
        try:
            import psutil
            return psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2
        except ImportError:
            pass
        try:
            with open('/proc/self/status') as fh:
                for line in fh:
                    if line.startswith('VmRSS:'):
                        return int(line.split()[1]) / 1024
        except OSError:
            pass
        return 0.0

    def _run(self):
        while self._running:
            rss = self.current_mb()
            if rss > self.peak_mb:
                self.peak_mb = rss
            time.sleep(self._interval)

    def reset_peak(self):
        self.peak_mb = self.current_mb()


# ---------------------------------------------------------------------------
# Domain factory
# ---------------------------------------------------------------------------

def _create_domain(nx, ny, mode, tmpdir):
    """
    Create a rectangular dam-break domain with ~2*nx*ny triangles.

    Uses DE0 flow algorithm (required for mode=2 GPU path), reflective
    boundaries, and a simple left/right stage initial condition.
    store=False to avoid I/O overhead.
    """
    import numpy as np
    import anuga
    from anuga import rectangular_cross_domain, Reflective_boundary

    domain = rectangular_cross_domain(nx, ny, len1=1000.0, len2=1000.0)
    domain.set_flow_algorithm('DE0')
    domain.set_low_froude(0)
    domain.set_name('gpu_bench')
    domain.set_datadir(tmpdir)
    domain.store = False

    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', lambda x, y: np.where(x < 500.0, 2.0, 0.5))
    domain.set_quantity('xmomentum', 0.0)
    domain.set_quantity('ymomentum', 0.0)
    domain.set_boundary({t: Reflective_boundary(domain)
                         for t in domain.get_boundary_tags()})

    if mode >= 1:
        domain.set_multiprocessor_mode(mode)

    return domain


# ---------------------------------------------------------------------------
# Scenario definitions
# GPU benchmarks use much larger meshes than the CPU benchmark suite.
# rectangular_cross_domain(nx, ny) creates 2*nx*ny triangles.
# ---------------------------------------------------------------------------

SCENARIOS = {
    #                    ~triangles       run time (short — just enough for timing)
    'small':  dict(nx=230,  ny=230,  finaltime=20.0,  yieldstep=5.0),   # ~106 K tris
    'medium': dict(nx=1000, ny=1000, finaltime=10.0,  yieldstep=5.0),   # ~2 M tris
    'large':  dict(nx=3162, ny=3162, finaltime=5.0,   yieldstep=5.0),   # ~20 M tris
}


# ---------------------------------------------------------------------------
# Run one scenario
# ---------------------------------------------------------------------------

def run_one(size, mode, omp_threads, sampler, flop_api, collect_detail=False):
    """
    Run a single GPU benchmark scenario and return a result dict.

    Parameters
    ----------
    size : str
        One of 'small', 'medium', 'large'.
    mode : int
        multiprocessor_mode (1=CPU-OpenMP, 2=GPU/C loop).
    omp_threads : int
        Value of OMP_NUM_THREADS (informational; caller sets env).
    sampler : _MemSampler
        Running memory sampler.
    flop_api : dict or None
        FLOP counter callables from _import_flop_api(), or None.
    collect_detail : bool
        If True, include per-kernel FLOP breakdown in result.

    Returns
    -------
    dict
    """
    cfg = SCENARIOS[size]
    tmpdir = tempfile.mkdtemp()

    try:
        domain = _create_domain(cfg['nx'], cfg['ny'], mode, tmpdir)
        n_tris = domain.number_of_triangles

        setup_rss_mb = sampler.current_mb()
        sampler.reset_peak()

        # Arm FLOP counters before evolve (mode=2 only; mode=1 counters stay zero)
        gpu_dom = None
        if flop_api is not None and mode == 2:
            try:
                gpu_dom = domain.gpu_interface.gpu_dom
                flop_api['reset'](gpu_dom)
                flop_api['enable'](gpu_dom, True)
            except AttributeError:
                gpu_dom = None

        if gpu_dom is not None:
            flop_api['start'](gpu_dom)

        t0 = time.perf_counter()
        for _ in domain.evolve(yieldstep=cfg['yieldstep'],
                               finaltime=cfg['finaltime']):
            pass
        wall_time_s = time.perf_counter() - t0

        if gpu_dom is not None:
            flop_api['stop'](gpu_dom)

        n_steps = domain.number_of_steps
        peak_rss_mb = sampler.peak_mb
        cells_per_s = (n_tris * n_steps / wall_time_s) if wall_time_s > 0 else 0.0

        # Collect FLOP stats
        gflops_per_second = 0.0
        total_gflops = 0.0
        flop_detail = {}
        if gpu_dom is not None:
            try:
                stats = flop_api['stats'](gpu_dom)
                gflops_per_second = stats.get('gflops_per_second', 0.0)
                total_gflops = stats.get('total_gflops', 0.0)
                if collect_detail:
                    detail_keys = [
                        'extrapolate', 'compute_fluxes', 'update',
                        'protect', 'manning', 'backup', 'saxpy',
                        'rate_operator', 'ghost_exchange',
                    ]
                    flop_detail = {k: stats[k] for k in detail_keys if k in stats}
            except Exception:
                pass

    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)

    rec = {
        'name': f'gpu_{size}_mode{mode}_t{omp_threads}',
        'size': size,
        'n_triangles': n_tris,
        'mode': mode,
        'omp_threads': omp_threads,
        'finaltime': cfg['finaltime'],
        'n_steps': n_steps,
        'wall_time_s': round(wall_time_s, 3),
        'cells_per_s': round(cells_per_s, 0),
        'gflops_per_second': round(gflops_per_second, 3),
        'total_gflops': round(total_gflops, 3),
        'setup_rss_mb': round(setup_rss_mb, 1),
        'peak_rss_mb': round(peak_rss_mb, 1),
        'mb_per_ktri': round(peak_rss_mb / (n_tris / 1000), 3) if n_tris else 0.0,
    }
    if collect_detail and flop_detail:
        rec['flop_detail'] = flop_detail

    return rec


# ---------------------------------------------------------------------------
# Metadata helpers
# ---------------------------------------------------------------------------

def _git_info():
    def _run(args):
        try:
            return subprocess.check_output(
                args, cwd=os.path.dirname(os.path.abspath(__file__)),
                stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return 'unknown'

    return {
        'commit': _run(['git', 'rev-parse', '--short', 'HEAD']),
        'branch': _run(['git', 'rev-parse', '--abbrev-ref', 'HEAD']),
        'commit_long': _run(['git', 'rev-parse', 'HEAD']),
    }


def _env_info(gpu_ext_present, gpu_device_available):
    omp = os.environ.get('OMP_NUM_THREADS', 'unset')
    offload = os.environ.get('OMP_TARGET_OFFLOAD', 'unset')
    try:
        import anuga
        anuga_version = anuga.__version__
    except Exception:
        anuga_version = 'unknown'

    return {
        'python_version': sys.version.split()[0],
        'platform': platform.system(),
        'hostname': platform.node().split('.')[0],
        'omp_num_threads_env': omp,
        'omp_target_offload_env': offload,
        'anuga_version': anuga_version,
        'gpu_ext_present': gpu_ext_present,
        'gpu_device_available': gpu_device_available,
        'cpu_only_mode': gpu_ext_present and not gpu_device_available,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='ANUGA GPU performance benchmark suite.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--sizes', default='small,medium',
        help='Comma-separated sizes to run: small,medium,large  (default: small,medium)',
    )
    parser.add_argument(
        '--modes', default='1,2',
        help='Comma-separated multiprocessor_modes (default: 1,2)',
    )
    parser.add_argument(
        '--output', default=None,
        help='Output JSON path.  Default: benchmarks/results/gpu_<branch>_<commit>_<ts>.json',
    )
    parser.add_argument(
        '--flop-detail', action='store_true',
        help='Include per-kernel FLOP breakdown in JSON output and table',
    )
    args = parser.parse_args()

    sizes = [s.strip() for s in args.sizes.split(',')]
    modes = [int(m.strip()) for m in args.modes.split(',')]

    for s in sizes:
        if s not in SCENARIOS:
            parser.error(f'Unknown size {s!r}. Choose from: {list(SCENARIOS)}')

    # --- probe GPU availability ---
    gpu_ext_present, gpu_device_available = _check_gpu_available()
    flop_api = _import_flop_api()

    # --- output path ---
    git = _git_info()
    env = _env_info(gpu_ext_present, gpu_device_available)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    omp_threads = int(os.environ.get('OMP_NUM_THREADS', 1))

    if args.output:
        outpath = args.output
    else:
        outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
        os.makedirs(outdir, exist_ok=True)
        branch_slug = git['branch'].replace('/', '_')
        fname = f"gpu_{branch_slug}_{git['commit']}_{timestamp}.json"
        outpath = os.path.join(outdir, fname)

    # --- header ---
    print(f'ANUGA GPU benchmark  commit={git["commit"]}  branch={git["branch"]}')
    print(f'Python {env["python_version"]}  OMP_NUM_THREADS={omp_threads}  '
          f'OMP_TARGET_OFFLOAD={env["omp_target_offload_env"]}')

    if gpu_device_available:
        print('GPU: device offloading ACTIVE')
    elif gpu_ext_present:
        print('GPU: sw_domain_gpu_ext present — running in CPU_ONLY_MODE '
              '(no real GPU target; set OMP_TARGET_OFFLOAD=mandatory to verify)')
    else:
        print('GPU: sw_domain_gpu_ext NOT present — mode=2 falls back to C RK loop on CPU')

    if flop_api is None:
        print('FLOP counters: unavailable (sw_domain_gpu_ext missing or too old)')
    else:
        print('FLOP counters: available')

    print(f'Sizes: {sizes}  Modes: {modes}')
    print(f'Output: {outpath}')
    print()

    sampler = _MemSampler(interval_s=0.1)
    sampler.start()

    results = []

    header = (
        f"{'Scenario':<32} {'tris':>10}  {'mode':>4}  {'thrd':>4}  "
        f"{'steps':>6}  {'wall(s)':>8}  {'cells/s':>12}  "
        f"{'GFLOP/s':>8}  {'peak MB':>8}  {'MB/Ktri':>8}"
    )
    rule = '-' * len(header)
    print(header)
    print(rule)

    for size in sizes:
        for mode in modes:
            name = f'gpu_{size}_mode{mode}_t{omp_threads}'
            sys.stdout.write(f'  Running {name} ... ')
            sys.stdout.flush()
            try:
                rec = run_one(size, mode, omp_threads, sampler,
                              flop_api, collect_detail=args.flop_detail)
                results.append(rec)
                print(
                    f"\r  {rec['name']:<32} {rec['n_triangles']:>10,}  "
                    f"{rec['mode']:>4}  {rec['omp_threads']:>4}  "
                    f"{rec['n_steps']:>6}  {rec['wall_time_s']:>8.2f}  "
                    f"{rec['cells_per_s']:>12,.0f}  "
                    f"{rec['gflops_per_second']:>8.2f}  "
                    f"{rec['peak_rss_mb']:>8.1f}  "
                    f"{rec['mb_per_ktri']:>8.3f}"
                )
                if args.flop_detail and 'flop_detail' in rec:
                    for k, v in rec['flop_detail'].items():
                        print(f"    {k:<28} {v:>14,} FLOPs")
            except Exception as exc:
                print(f'\r  {name:<32} FAILED: {exc}')
                import traceback
                traceback.print_exc()

    sampler.stop()
    print(rule)
    print()

    # --- mode=1 vs mode=2 speedup summary ---
    by_size_mode = {}
    for rec in results:
        by_size_mode[(rec['size'], rec['mode'])] = rec

    speedup_lines = []
    for size in sizes:
        r1 = by_size_mode.get((size, 1))
        r2 = by_size_mode.get((size, 2))
        if r1 and r2 and r1['wall_time_s'] > 0:
            speedup = r1['wall_time_s'] / r2['wall_time_s']
            speedup_lines.append(
                f"  {size:<10} mode2 vs mode1 speedup: {speedup:>6.2f}x  "
                f"({r2['cells_per_s']:>12,.0f} vs {r1['cells_per_s']:>12,.0f} cells/s)"
            )
    if speedup_lines:
        print('Speedup summary (mode=2 vs mode=1):')
        for line in speedup_lines:
            print(line)
        print()

    payload = {
        'timestamp': timestamp,
        'git': git,
        'env': env,
        'omp_threads': omp_threads,
        'results': results,
    }
    with open(outpath, 'w') as fh:
        json.dump(payload, fh, indent=2)
    print(f'Saved {len(results)} results → {outpath}')

    return 0


if __name__ == '__main__':
    sys.exit(main())
