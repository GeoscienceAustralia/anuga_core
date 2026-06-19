#!/usr/bin/env python3
"""
ANUGA Performance Benchmark Suite
----------------------------------
Measures wall time and peak memory for simulations of increasing size and
for each multiprocessor_mode (0=Python Euler, 1=Python RK2, 2=C RK2/GPU).

Results are written to a JSON file for later comparison with
compare_benchmarks.py.

Usage
-----
    # Run all scenarios (small + medium), save to results/
    python benchmarks/run_benchmarks.py

    # Quick sanity check (small only)
    python benchmarks/run_benchmarks.py --sizes small

    # All sizes including large (~360K tris, slow)
    python benchmarks/run_benchmarks.py --sizes small,medium,large

    # Only mode 0 and 2
    python benchmarks/run_benchmarks.py --modes 0,2

    # Custom output file
    python benchmarks/run_benchmarks.py --output my_results.json

Metrics
-------
- wall_time_s      : wall-clock seconds for the evolve loop only
- n_steps          : number of internal (CFL) timesteps taken
- cells_per_s      : n_triangles * n_steps / wall_time  (the primary figure)
- setup_rss_mb     : RSS after domain creation (before any evolve)
- peak_rss_mb      : peak RSS sampled during evolve (100 ms polling)
- mb_per_ktri      : peak_rss_mb / (n_triangles / 1000) — memory efficiency
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
# Memory sampler (background thread)
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
        # fallback: /proc/self/status on Linux
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
    Create a rectangular dam-break domain with nx*ny cells per side.

    The domain has no file output (store=False) and uses reflective boundaries.
    Initial condition: left half stage=2m, right half stage=0.5m.
    """
    import numpy as np
    import anuga
    from anuga import rectangular_cross_domain, Reflective_boundary

    domain = rectangular_cross_domain(nx, ny, len1=1000.0, len2=1000.0)
    domain.set_flow_algorithm('DE0')
    domain.set_low_froude(0)
    domain.set_name('bench')
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
# ---------------------------------------------------------------------------

SCENARIOS = {
    'small':  dict(nx=50,  ny=50,  finaltime=200.0, yieldstep=50.0),
    'medium': dict(nx=150, ny=150, finaltime=100.0, yieldstep=25.0),
    'large':  dict(nx=300, ny=300, finaltime=50.0,  yieldstep=12.5),
}


# ---------------------------------------------------------------------------
# Run one scenario
# ---------------------------------------------------------------------------

def run_one(size, mode, omp_threads, sampler):
    """
    Run a single benchmark scenario and return a result dict.

    Parameters
    ----------
    size : str
        One of 'small', 'medium', 'large'.
    mode : int
        multiprocessor_mode (0, 1, or 2).
    omp_threads : int
        Value of OMP_NUM_THREADS (informational only — caller must set env).
    sampler : _MemSampler
        Running memory sampler; peak is reset just before evolve starts.

    Returns
    -------
    dict
        Benchmark result record.
    """
    cfg = SCENARIOS[size]
    tmpdir = tempfile.mkdtemp()

    try:
        domain = _create_domain(cfg['nx'], cfg['ny'], mode, tmpdir)
        n_tris = domain.number_of_triangles

        # Memory snapshot after full domain setup (before any evolve)
        setup_rss_mb = sampler.current_mb()
        sampler.reset_peak()

        t0 = time.perf_counter()
        for _ in domain.evolve(yieldstep=cfg['yieldstep'],
                               finaltime=cfg['finaltime']):
            pass
        wall_time_s = time.perf_counter() - t0

        n_steps = domain.number_of_steps
        peak_rss_mb = sampler.peak_mb
        cells_per_s = (n_tris * n_steps / wall_time_s) if wall_time_s > 0 else 0.0

    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)

    return {
        'name': f'{size}_mode{mode}_t{omp_threads}',
        'size': size,
        'n_triangles': n_tris,
        'mode': mode,
        'omp_threads': omp_threads,
        'finaltime': cfg['finaltime'],
        'n_steps': n_steps,
        'wall_time_s': round(wall_time_s, 3),
        'cells_per_s': round(cells_per_s, 0),
        'setup_rss_mb': round(setup_rss_mb, 1),
        'peak_rss_mb': round(peak_rss_mb, 1),
        'mb_per_ktri': round(peak_rss_mb / (n_tris / 1000), 3) if n_tris else 0.0,
    }


# ---------------------------------------------------------------------------
# Metadata helpers
# ---------------------------------------------------------------------------

def _git_info():
    def _run(args):
        try:
            return subprocess.check_output(
                args, cwd=os.path.dirname(__file__),
                stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return 'unknown'

    return {
        'commit': _run(['git', 'rev-parse', '--short', 'HEAD']),
        'branch': _run(['git', 'rev-parse', '--abbrev-ref', 'HEAD']),
        'commit_long': _run(['git', 'rev-parse', 'HEAD']),
    }


def _env_info():
    omp = os.environ.get('OMP_NUM_THREADS', 'unset')
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
        'anuga_version': anuga_version,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='ANUGA performance benchmark suite.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--sizes', default='small,medium',
        help='Comma-separated sizes to run: small,medium,large  (default: small,medium)',
    )
    parser.add_argument(
        '--modes', default='0,1,2',
        help='Comma-separated multiprocessor_modes to test  (default: 0,1,2)',
    )
    parser.add_argument(
        '--output', default=None,
        help='Output JSON path. Default: benchmarks/results/<branch>_<commit>_<timestamp>.json',
    )
    args = parser.parse_args()

    sizes = [s.strip() for s in args.sizes.split(',')]
    modes = [int(m.strip()) for m in args.modes.split(',')]

    for s in sizes:
        if s not in SCENARIOS:
            parser.error(f'Unknown size {s!r}. Choose from: {list(SCENARIOS)}')

    # --- output path ---
    git = _git_info()
    env = _env_info()
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    omp_threads = int(os.environ.get('OMP_NUM_THREADS', 1))

    if args.output:
        outpath = args.output
    else:
        outdir = os.path.join(os.path.dirname(__file__), 'results')
        os.makedirs(outdir, exist_ok=True)
        fname = f"{git['branch'].replace('/', '_')}_{git['commit']}_{timestamp}.json"
        outpath = os.path.join(outdir, fname)

    print(f'ANUGA benchmark  commit={git["commit"]}  branch={git["branch"]}')
    print(f'Python {env["python_version"]}  OMP_NUM_THREADS={omp_threads}')
    print(f'Sizes: {sizes}  Modes: {modes}')
    print(f'Output: {outpath}')
    print()

    sampler = _MemSampler(interval_s=0.1)
    sampler.start()

    results = []
    header = f"{'Scenario':<28} {'tris':>8}  {'mode':>4}  {'thrd':>4}  {'steps':>6}  {'wall(s)':>8}  {'cells/s':>10}  {'setup MB':>9}  {'peak MB':>8}  {'MB/Ktri':>8}"
    rule = '-' * len(header)
    print(header)
    print(rule)

    for size in sizes:
        for mode in modes:
            name = f'{size}_mode{mode}_t{omp_threads}'
            sys.stdout.write(f'  Running {name} ... ')
            sys.stdout.flush()
            try:
                rec = run_one(size, mode, omp_threads, sampler)
                results.append(rec)
                print(
                    f"\r  {rec['name']:<28} {rec['n_triangles']:>8}  {rec['mode']:>4}  "
                    f"{rec['omp_threads']:>4}  {rec['n_steps']:>6}  "
                    f"{rec['wall_time_s']:>8.2f}  {rec['cells_per_s']:>10,.0f}  "
                    f"{rec['setup_rss_mb']:>9.1f}  {rec['peak_rss_mb']:>8.1f}  "
                    f"{rec['mb_per_ktri']:>8.3f}"
                )
            except Exception as exc:
                print(f'\r  {name:<28} FAILED: {exc}')

    sampler.stop()
    print(rule)
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
