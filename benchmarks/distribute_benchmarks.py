#!/usr/bin/env python3
"""Benchmark ANUGA parallel mesh distribution methods.

Compares four approaches for distributing a mesh across MPI ranks:

  distribute()              -- traditional: full Domain on rank 0, then distribute
  distribute_collaborative() -- shared-memory cooperative version of distribute()
  distribute_basic_mesh()   -- mesh-first: only Basic_mesh on rank 0 (no quantities),
                               quantities set locally after distribution
  dump()+load()             -- rank 0 partitions and writes files; all ranks read

Run with:

    mpirun -np 8 python benchmarks/distribute_benchmarks.py --size 500
    mpirun -np 8 python benchmarks/distribute_benchmarks.py --size 500 --scheme morton
    mpirun -np 8 python benchmarks/distribute_benchmarks.py --size 500 --reps 3

Options
-------
--size M       Grid size M: produces 4*M*M triangles (default 500)
--reps R       Repetitions for median timing (default 1)
--scheme S     Partition scheme: metis | morton | hilbert (default metis)
--interval T   Memory ticker sample interval in seconds (default 1.0)
--no-evolve    Skip the correctness evolve check after timing
"""

import argparse
import gc
import shutil
import statistics
import tempfile
import threading
import time

from mpi4py import MPI

import anuga
from anuga.parallel.parallel_api import distribute, distribute_collaborative
from anuga.parallel.parallel_api import distribute_basic_mesh
from anuga.parallel.sequential_distribute import (
    sequential_distribute_dump, sequential_distribute_load)

comm  = MPI.COMM_WORLD
myid  = comm.Get_rank()
nproc = comm.Get_size()


# ---------------------------------------------------------------------------
# Memory helpers
# ---------------------------------------------------------------------------

def _rss_mb():
    """Current RSS in MiB from /proc/self/status."""
    try:
        with open('/proc/self/status') as fh:
            for line in fh:
                if line.startswith('VmRSS:'):
                    return int(line.split()[1]) / 1024.0
    except OSError:
        pass
    return 0.0


def _pss_mb():
    """Proportional Set Size in MiB (shared pages counted proportionally).

    Summing PSS across all ranks gives the true physical-memory footprint
    of the job, unlike VmRSS which double-counts pages shared via MPI.Win.
    """
    try:
        with open('/proc/self/smaps_rollup') as fh:
            for line in fh:
                if line.startswith('Pss:'):
                    return int(line.split()[1]) / 1024.0
    except OSError:
        pass
    try:
        total = 0
        with open('/proc/self/smaps') as fh:
            for line in fh:
                if line.startswith('Pss:'):
                    total += int(line.split()[1])
        return total / 1024.0
    except OSError:
        pass
    return _rss_mb()


class MemoryMonitor:
    """Background thread sampling this process's RSS/PSS every `interval` s.

    Use as a context manager::

        with MemoryMonitor(label='distribute()', interval=1.0) as mon:
            fn(domain)
        peak_pss_mib = mon.peak_pss_mb
    """

    def __init__(self, label='', interval=1.0, print_ticker=False):
        self.label        = label
        self.interval     = interval
        self.print_ticker = print_ticker
        self.peak_rss_mb  = 0.0
        self.peak_pss_mb  = 0.0
        self._stop        = threading.Event()
        self._thread      = threading.Thread(target=self._run, daemon=True)

    def __enter__(self):
        self.peak_rss_mb = _rss_mb()
        self.peak_pss_mb = _pss_mb()
        self._t0         = time.time()
        self._stop.clear()
        self._thread.start()
        return self

    def __exit__(self, *_):
        self._stop.set()
        self._thread.join()

    def _run(self):
        while not self._stop.wait(self.interval):
            rss = _rss_mb()
            pss = _pss_mb()
            if rss > self.peak_rss_mb:
                self.peak_rss_mb = rss
            if pss > self.peak_pss_mb:
                self.peak_pss_mb = pss
            if self.print_ticker:
                elapsed = time.time() - self._t0
                print(f'  [{self.label}]  t={elapsed:6.1f}s  '
                      f'rank-0 RSS={rss:,.0f} MiB  PSS={pss:,.0f} MiB',
                      flush=True)


# ---------------------------------------------------------------------------
# Shared-memory diagnostic
# ---------------------------------------------------------------------------

def check_shmem():
    """Test whether MPI.Win.Allocate_shared works. Returns (ok, node_size, reason)."""
    import numpy as np
    node_comm = comm.Split_type(MPI.COMM_TYPE_SHARED)
    node_rank = node_comm.Get_rank()
    node_size = node_comm.Get_size()
    try:
        win = MPI.Win.Allocate_shared(
            8 if node_rank == 0 else 0, 8, MPI.INFO_NULL, node_comm)
        buf, _ = win.Shared_query(0)
        arr = np.ndarray((1,), dtype=np.float64, buffer=buf)
        if node_rank == 0:
            arr[0] = 42.0
        node_comm.Barrier()
        ok = float(arr[0]) == 42.0
        win.Free()
        node_comm.Free()
        return ok, node_size, ('Win.Allocate_shared OK' if ok
                               else 'shared query returned wrong value')
    except Exception as e:
        try:
            node_comm.Free()
        except Exception:
            pass
        return False, node_size, f'exception: {e}'


# ---------------------------------------------------------------------------
# Domain / mesh factories
# ---------------------------------------------------------------------------

def make_domain(grid_size):
    """Build a full Domain on rank 0 with a simple initial condition."""
    pts, verts, bnd = anuga.rectangular_cross(grid_size, grid_size,
                                              len1=float(grid_size),
                                              len2=float(grid_size))
    domain = anuga.Domain(pts, verts, bnd)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage',     lambda x, y: 1.0 + 0.01 * x)
    domain.set_quantity('friction',  0.03)
    return domain


def make_basic_mesh(grid_size):
    """Build a Basic_mesh on rank 0 only (no quantities)."""
    if myid != 0:
        return None
    from anuga.abstract_2d_finite_volumes.basic_mesh import (
        rectangular_cross_basic_mesh)
    return rectangular_cross_basic_mesh(grid_size, grid_size,
                                        len1=float(grid_size),
                                        len2=float(grid_size))


# ---------------------------------------------------------------------------
# Ghost-triangle statistics
# ---------------------------------------------------------------------------

def collect_ghost_stats(pd):
    """Return (ghost_sum, ghost_max, ghost_min) reduced to rank 0."""
    import numpy as np
    n_ghost  = int(np.sum(pd.tri_full_flag == 0))
    g_sum = comm.reduce(n_ghost, op=MPI.SUM, root=0)
    g_max = comm.reduce(n_ghost, op=MPI.MAX, root=0)
    g_min = comm.reduce(n_ghost, op=MPI.MIN, root=0)
    return (g_sum or 0), (g_max or 0), (g_min or 0)


# ---------------------------------------------------------------------------
# Timed benchmarks
# ---------------------------------------------------------------------------

def time_distribute(fn, grid_size, reps, ticker_interval, parameters=None):
    """Time distribute() or distribute_collaborative().

    Returns (times_s, pss_sum_mbs, rss_max_mbs, ghost_stats) where values
    are only valid on rank 0.
    """
    times       = []
    pss_sum_mbs = []
    rss_max_mbs = []
    ghost_stats = (0, 0, 0)
    label       = fn.__name__

    for _ in range(reps):
        gc.collect()
        comm.Barrier()
        domain = make_domain(grid_size)
        comm.Barrier()

        mon = MemoryMonitor(label=label, interval=ticker_interval,
                            print_ticker=(myid == 0))
        with mon:
            t0 = MPI.Wtime()
            pd = fn(domain, parameters=parameters)
            t1 = MPI.Wtime()

        elapsed  = t1 - t0
        peak_rss = max(mon.peak_rss_mb, _rss_mb())
        peak_pss = max(mon.peak_pss_mb, _pss_mb())
        ghost_stats = collect_ghost_stats(pd)

        wall    = comm.reduce(elapsed,  op=MPI.MAX, root=0)
        pss_sum = comm.reduce(peak_pss, op=MPI.SUM, root=0)
        rss_max = comm.reduce(peak_rss, op=MPI.MAX, root=0)

        if myid == 0:
            times.append(wall)
            pss_sum_mbs.append(pss_sum)
            rss_max_mbs.append(rss_max)

        del pd, domain

    return times, pss_sum_mbs, rss_max_mbs, ghost_stats


def time_distribute_basic_mesh(grid_size, reps, ticker_interval, parameters=None):
    """Time distribute_basic_mesh() -- mesh-first workflow.

    Rank 0 builds only a Basic_mesh (no quantities).  All ranks call
    distribute_basic_mesh(), then set quantities locally on the result.

    Returns (times_s, pss_sum_mbs, rss_max_mbs, ghost_stats).
    """
    times       = []
    pss_sum_mbs = []
    rss_max_mbs = []
    ghost_stats = (0, 0, 0)

    for _ in range(reps):
        gc.collect()
        comm.Barrier()
        bm = make_basic_mesh(grid_size)
        comm.Barrier()

        mon = MemoryMonitor(label='distribute_basic_mesh', interval=ticker_interval,
                            print_ticker=(myid == 0))
        with mon:
            t0 = MPI.Wtime()
            pd = distribute_basic_mesh(bm, parameters=parameters)
            t1 = MPI.Wtime()

        # Set quantities post-distribution (part of the workflow)
        pd.set_quantity('elevation', 0.0)
        pd.set_quantity('stage',     lambda x, y: 1.0 + 0.01 * x)
        pd.set_quantity('friction',  0.03)

        elapsed  = t1 - t0
        peak_rss = max(mon.peak_rss_mb, _rss_mb())
        peak_pss = max(mon.peak_pss_mb, _pss_mb())
        ghost_stats = collect_ghost_stats(pd)

        wall    = comm.reduce(elapsed,  op=MPI.MAX, root=0)
        pss_sum = comm.reduce(peak_pss, op=MPI.SUM, root=0)
        rss_max = comm.reduce(peak_rss, op=MPI.MAX, root=0)

        if myid == 0:
            times.append(wall)
            pss_sum_mbs.append(pss_sum)
            rss_max_mbs.append(rss_max)

        del pd, bm
        gc.collect()
        comm.Barrier()

    return times, pss_sum_mbs, rss_max_mbs, ghost_stats


def time_dump_load(grid_size, reps, ticker_interval, parameters=None):
    """Time sequential_distribute_dump (rank 0) + load (all ranks).

    Returns (times_dump, times_load, pss_sum_mbs, rss_max_mbs, ghost_stats).
    """
    times_dump  = []
    times_load  = []
    pss_sum_mbs = []
    rss_max_mbs = []
    ghost_stats = (0, 0, 0)

    for _ in range(reps):
        gc.collect()
        comm.Barrier()

        if myid == 0:
            tmpdir = tempfile.mkdtemp(prefix='anuga_bench_')
        else:
            tmpdir = None
        tmpdir = comm.bcast(tmpdir, root=0)

        domain      = make_domain(grid_size)
        domain_name = domain.get_name()
        comm.Barrier()

        # Dump phase (rank 0 only)
        if myid == 0:
            mon_dump = MemoryMonitor(label='dump', interval=ticker_interval,
                                     print_ticker=True)
            with mon_dump:
                t0 = MPI.Wtime()
                sequential_distribute_dump(domain, nproc,
                                           partition_dir=tmpdir,
                                           parameters=parameters)
                t1 = MPI.Wtime()
            dump_elapsed = t1 - t0
        else:
            dump_elapsed = 0.0

        del domain
        comm.Barrier()

        # Load phase (all ranks)
        mon_load = MemoryMonitor(label='load', interval=ticker_interval,
                                 print_ticker=(myid == 0))
        with mon_load:
            t0 = MPI.Wtime()
            pd = sequential_distribute_load(filename=domain_name,
                                            partition_dir=tmpdir)
            t1 = MPI.Wtime()

        load_elapsed = t1 - t0
        peak_rss = max(mon_load.peak_rss_mb, _rss_mb())
        peak_pss = max(mon_load.peak_pss_mb, _pss_mb())
        ghost_stats = collect_ghost_stats(pd)

        dump_wall = comm.bcast(dump_elapsed, root=0)
        load_wall = comm.reduce(load_elapsed, op=MPI.MAX, root=0)
        pss_sum   = comm.reduce(peak_pss,     op=MPI.SUM, root=0)
        rss_max   = comm.reduce(peak_rss,     op=MPI.MAX, root=0)

        if myid == 0:
            times_dump.append(dump_wall)
            times_load.append(load_wall)
            pss_sum_mbs.append(pss_sum)
            rss_max_mbs.append(rss_max)

        del pd

        if myid == 0:
            shutil.rmtree(tmpdir, ignore_errors=True)
        comm.Barrier()

    return times_dump, times_load, pss_sum_mbs, rss_max_mbs, ghost_stats


def run_evolve_check(grid_size, scheme):
    """Distribute a small mesh and run one evolve step to verify correctness."""
    import numpy as np

    small = min(grid_size, 50)
    if myid == 0:
        print(f'\n  Evolve check ({small}x{small} mesh, distribute_basic_mesh) ...')

    bm = make_basic_mesh(small)
    pd = distribute_basic_mesh(bm, parameters={'partition_scheme': scheme})

    pd.set_quantity('elevation', 0.0)
    pd.set_quantity('stage',     lambda x, y: np.where(
        (x > small * 0.4) & (x < small * 0.6), 1.0, 0.0))
    pd.set_quantity('friction',  0.03)
    Br = anuga.Reflective_boundary(pd)
    pd.set_boundary({t: Br for t in pd.get_boundary_tags()})

    t0 = time.perf_counter()
    for _ in pd.evolve(yieldstep=0.01, finaltime=0.01):
        pass
    dt = time.perf_counter() - t0

    if myid == 0:
        print(f'  Evolve check passed in {dt:.2f}s')


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def fmt_time(times):
    if not times:
        return '--'
    med = statistics.median(times)
    if len(times) == 1:
        return f'{med:.3f}s'
    return f'{med:.3f}s  (min {min(times):.3f}  max {max(times):.3f})'


def fmt_mem(mbs):
    if not mbs:
        return '--'
    med = statistics.median(mbs)
    return f'{med/1024:.2f} GiB' if med >= 1024 else f'{med:.0f} MiB'


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Benchmark ANUGA parallel mesh distribution methods.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument('--size',      type=int,   default=500)
    p.add_argument('--reps',      type=int,   default=1)
    p.add_argument('--scheme',    type=str,   default='metis',
                   choices=['metis', 'morton', 'hilbert'])
    p.add_argument('--interval',  type=float, default=1.0)
    p.add_argument('--no-evolve', action='store_true')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args     = parse_args()
    M        = args.size
    ntri     = 4 * M * M
    params   = {'partition_scheme': args.scheme}

    shmem_ok, node_size, shmem_reason = check_shmem()

    if myid == 0:
        shmem_status = ('yes' if shmem_ok
                        else f'NO -- falling back to Bcast ({shmem_reason})')
        W = 66
        print('=' * W)
        print('  ANUGA distribute benchmark')
        print(f'  Mesh:       rectangular_cross {M}x{M}  ({ntri:,} triangles)')
        print(f'  MPI ranks:  {nproc}   Repetitions: {args.reps}')
        print(f'  Ranks/node: {node_size}')
        print(f'  Scheme:     {args.scheme}')
        print(f'  Shared mem: {shmem_status}')
        print('=' * W)

    # Run all four methods
    times_std,   pss_std,   rss_std,   ghost_std   = time_distribute(
        distribute, M, args.reps, args.interval, params)

    gc.collect()
    comm.Barrier()

    times_col,   pss_col,   rss_col,   ghost_col   = time_distribute(
        distribute_collaborative, M, args.reps, args.interval, params)

    gc.collect()
    comm.Barrier()

    times_bm,    pss_bm,    rss_bm,    ghost_bm    = time_distribute_basic_mesh(
        M, args.reps, args.interval, params)

    gc.collect()
    comm.Barrier()

    times_dump, times_load, pss_dl, rss_dl, ghost_dl = time_dump_load(
        M, args.reps, args.interval, params)

    if not args.no_evolve:
        run_evolve_check(M, args.scheme)

    if myid == 0:
        times_dl = [d + load_t for d, load_t in zip(times_dump, times_load)]

        W = 32
        C = 22

        def row(label, *vals):
            cells = '  '.join(f'{v:<{C}}' for v in vals)
            print(f'  {label:<{W}}  {cells}')

        sep = '-' * (W + 2 + 4 * (C + 2))
        print(f'\n{sep}')
        row('', 'distribute()', 'collaborative()', 'distribute_basic_mesh()', 'dump()+load()')
        print(sep)
        row('Wall time (median)',
            fmt_time(times_std), fmt_time(times_col),
            fmt_time(times_bm),  fmt_time(times_dl))
        row('Peak PSS sum (physical total)',
            fmt_mem(pss_std), fmt_mem(pss_col),
            fmt_mem(pss_bm),  fmt_mem(pss_dl))
        row('Peak RSS max single rank',
            fmt_mem(rss_std), fmt_mem(rss_col),
            fmt_mem(rss_bm),  fmt_mem(rss_dl))
        print(sep)
        print('  Note: PSS sums shared pages proportionally.')
        if times_dump:
            med_dump = statistics.median(times_dump)
            med_load = statistics.median(times_load)
            print(f'  dump()+load() breakdown: '
                  f'dump (serial) {med_dump:.3f}s  +  '
                  f'load (parallel) {med_load:.3f}s')

        # Speedup summary
        print()
        candidates = [
            ('distribute()',             statistics.median(times_std) if times_std else None),
            ('collaborative()',          statistics.median(times_col) if times_col else None),
            ('distribute_basic_mesh()',  statistics.median(times_bm)  if times_bm  else None),
            ('dump()+load()',            statistics.median(times_dl)  if times_dl  else None),
        ]
        valid = [(n, t) for n, t in candidates if t is not None]
        best_name, best_time = min(valid, key=lambda x: x[1])
        print(f'  Fastest: {best_name}  ({best_time:.3f}s)')
        for name, t in valid:
            if name != best_name and t > 0:
                print(f'    vs {name}: {t / best_time:.2f}x slower')

        # Ghost triangle stats
        print()
        print(f'  Partition quality  ({args.scheme}, {nproc} ranks, {ntri:,} triangles):')
        for label, gs in [('distribute()', ghost_std),
                          ('collaborative()', ghost_col),
                          ('distribute_basic_mesh()', ghost_bm),
                          ('dump()+load()', ghost_dl)]:
            g_sum, g_max, g_min = gs
            if g_sum > 0:
                pct = 100.0 * g_sum / ntri
                avg = g_sum / nproc
                print(f'    {label:<28}  ghost={g_sum:,}  ({pct:.1f}%)  '
                      f'avg={avg:,.0f}/rank  min={g_min:,}  max={g_max:,}')
            else:
                print(f'    {label:<28}  (no ghost data)')
        print()


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        import sys
        import traceback
        print(f'\n[rank {myid}] ERROR: {e}', flush=True)
        traceback.print_exc()
        comm.Abort(1)
        sys.exit(1)
    MPI.Finalize()
