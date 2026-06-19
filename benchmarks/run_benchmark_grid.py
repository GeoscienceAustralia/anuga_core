#!/usr/bin/env python3
"""Run distribute_benchmarks.py over a grid of (np, scheme) values and
print a consolidated summary table.

Each (np, scheme) combination is run as a subprocess and its stdout saved to
<outdir>/bench_np<N>_<scheme>.txt.  Existing files are reused unless
--force is given, so partial grids can be resumed cheaply.

Usage
-----
    python benchmarks/run_benchmark_grid.py [options]

Options
-------
--size M         Mesh grid size M (4*M*M triangles, default 1000)
--reps R         Repetitions per function (default 1)
--interval S     Ticker sample interval in seconds (default 5.0)
--np LIST        Comma-separated process counts (default 10,20,30)
--schemes LIST   Comma-separated schemes (default metis,morton,hilbert)
--outdir DIR     Directory for per-run output files (default benchmarks/results/dist)
--mpirun CMD     MPI launcher (default mpirun)
--script PATH    Path to distribute_benchmarks.py (default: auto-detect)
--force          Re-run even if output file already exists
--dry-run        Print commands without executing them
--no-summary     Skip summary; just run and save
"""

import argparse
import re
import subprocess
from datetime import datetime
from pathlib import Path


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(add_help=True)
    p.add_argument('--size',       type=int,   default=1000)
    p.add_argument('--reps',       type=int,   default=1)
    p.add_argument('--interval',   type=float, default=5.0)
    p.add_argument('--np',         type=str,   default='10,20,30')
    p.add_argument('--schemes',    type=str,   default='metis,morton,hilbert')
    p.add_argument('--outdir',     type=str,   default=None)
    p.add_argument('--mpirun',     type=str,   default='mpirun')
    p.add_argument('--script',     type=str,   default=None)
    p.add_argument('--force',      action='store_true')
    p.add_argument('--dry-run',    action='store_true')
    p.add_argument('--no-summary', action='store_true')
    a = p.parse_args()
    a.np_values   = [int(x) for x in a.np.split(',')]
    a.scheme_list = [s.strip() for s in a.schemes.split(',')]
    if a.outdir is None:
        here = Path(__file__).resolve().parent
        a.outdir = str(here / 'results' / 'dist')
    return a


def find_benchmark_script(explicit=None):
    if explicit:
        return explicit
    here = Path(__file__).resolve().parent
    candidate = here / 'distribute_benchmarks.py'
    if candidate.exists():
        return str(candidate)
    raise FileNotFoundError(
        'Cannot find distribute_benchmarks.py; use --script to specify it.')


# ---------------------------------------------------------------------------
# Run one (np, scheme) combination
# ---------------------------------------------------------------------------

def output_path(outdir, np_val, scheme):
    return Path(outdir) / f'bench_np{np_val}_{scheme}.txt'


def run_one(mpirun, script, np_val, scheme, size, reps, interval,
            outfile, dry_run=False):
    cmd = [mpirun, '-np', str(np_val), 'python', script,
           '--size', str(size),
           '--reps', str(reps),
           '--interval', str(interval),
           '--scheme', scheme,
           '--no-evolve']
    print(f'  $ {" ".join(cmd)}')
    print(f'    -> {outfile}')
    if dry_run:
        return True
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600)
        text = result.stdout
        if result.returncode != 0:
            text += f'\n[STDERR]\n{result.stderr}'
        Path(outfile).write_text(text)
        if result.returncode != 0:
            print(f'    [FAILED -- exit code {result.returncode}]')
            return False
        return True
    except subprocess.TimeoutExpired:
        print('    [TIMED OUT]')
        return False
    except Exception as e:
        print(f'    [ERROR: {e}]')
        return False


# ---------------------------------------------------------------------------
# Parse one run's stdout into a metrics dict
# ---------------------------------------------------------------------------

def _first_float(pattern, text):
    m = re.search(pattern, text)
    return float(m.group(1)) if m else None


def _first_int(pattern, text):
    m = re.search(pattern, text)
    return int(m.group(1).replace(',', '')) if m else None


def parse_output(text):
    """Extract benchmark metrics from one run's stdout."""
    r = {}

    r['triangles'] = _first_int(r'(\d[\d,]+) triangles', text)
    r['ranks']     = _first_int(r'MPI ranks:\s+(\d+)', text)
    r['scheme']    = (re.search(r'Scheme:\s+(\S+)', text) or
                      type('', (), {'group': lambda s, n: 'metis'})()).group(1)

    # Wall time row: 4 columns separated by whitespace
    wall = re.search(
        r'Wall time[^\d]+([\d.]+)s\s+([\d.]+)s\s+([\d.]+)s\s+([\d.]+)s', text)
    if wall:
        r['wall_dist']  = float(wall.group(1))
        r['wall_col']   = float(wall.group(2))
        r['wall_bm']    = float(wall.group(3))
        r['wall_dl']    = float(wall.group(4))

    # PSS row: 4 columns
    pss = re.search(
        r'Peak PSS[^\d]+([\d.]+) GiB\s+([\d.]+) GiB\s+([\d.]+) GiB\s+([\d.]+) GiB',
        text)
    if pss:
        r['pss_dist'] = float(pss.group(1))
        r['pss_col']  = float(pss.group(2))
        r['pss_bm']   = float(pss.group(3))
        r['pss_dl']   = float(pss.group(4))

    # RSS row: 4 columns
    rss = re.search(
        r'Peak RSS[^\d]+([\d.]+) GiB\s+([\d.]+) GiB\s+([\d.]+) GiB\s+([\d.]+) GiB',
        text)
    if rss:
        r['rss_dist'] = float(rss.group(1))
        r['rss_col']  = float(rss.group(2))
        r['rss_bm']   = float(rss.group(3))
        r['rss_dl']   = float(rss.group(4))

    # dump+load breakdown
    dl = re.search(r'dump \(serial\) ([\d.]+)s.*?load \(parallel\) ([\d.]+)s', text)
    if dl:
        r['dump_time'] = float(dl.group(1))
        r['load_time'] = float(dl.group(2))

    # Ghost % (from distribute() column)
    g_tot = re.search(r'distribute\(\)\s+ghost=[\d,]+\s+\(([\d.]+)%\)', text)
    if g_tot:
        r['ghost_pct'] = float(g_tot.group(1))

    return r


# ---------------------------------------------------------------------------
# Summary printer
# ---------------------------------------------------------------------------

def _t(val):
    return f'{val:6.2f}s' if val is not None else '  --   '


def _g(val):
    return f'{val:6.2f}' if val is not None else '   --  '


def _pct(val):
    return f'{val:5.1f}%' if val is not None else '  --  '


def print_summary(results, args, ntri, outdir):
    lines = []

    def out(*a):
        s = ' '.join(str(x) for x in a)
        lines.append(s)
        print(s)

    W = 76
    out('=' * W)
    out('  ANUGA distribute() benchmark grid summary')
    out(f'  Mesh:      synthetic {args.size}x{args.size}  ({ntri:,} triangles)')
    out(f'  np values: {args.np_values}')
    out(f'  Schemes:   {args.scheme_list}')
    out(f'  Reps:      {args.reps}')
    out(f'  Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    out('=' * W)

    # Wall time table
    out()
    out('  Wall time (seconds)')
    hdr = (f'  {"np":>4}  {"scheme":<8}  '
           f'{"distribute()":>13}  {"collaborative()":>15}  '
           f'{"dist_basic_mesh()":>18}  {"dump+load":>10}  '
           f'{"dump":>7}  {"load":>6}')
    out(hdr)
    out('  ' + '-' * (len(hdr) - 2))
    for np_val in args.np_values:
        for scheme in args.scheme_list:
            r = results.get((np_val, scheme), {})
            out(f'  {np_val:>4}  {scheme:<8}  '
                f'{_t(r.get("wall_dist")):>13}  '
                f'{_t(r.get("wall_col")):>15}  '
                f'{_t(r.get("wall_bm")):>18}  '
                f'{_t(r.get("wall_dl")):>10}  '
                f'{_t(r.get("dump_time")):>7}  '
                f'{_t(r.get("load_time")):>6}')

    # Memory table
    out()
    out('  Peak PSS -- physical memory total across all ranks (GiB)')
    hdr2 = (f'  {"np":>4}  {"scheme":<8}  '
            f'{"distribute()":>13}  {"collaborative()":>15}  '
            f'{"dist_basic_mesh()":>18}  {"dump+load":>10}')
    out(hdr2)
    out('  ' + '-' * (len(hdr2) - 2))
    for np_val in args.np_values:
        for scheme in args.scheme_list:
            r = results.get((np_val, scheme), {})
            out(f'  {np_val:>4}  {scheme:<8}  '
                f'{_g(r.get("pss_dist")):>13}  '
                f'{_g(r.get("pss_col")):>15}  '
                f'{_g(r.get("pss_bm")):>18}  '
                f'{_g(r.get("pss_dl")):>10}')

    # Ghost triangle table
    if any('ghost_pct' in r for r in results.values()):
        out()
        out('  Partition quality -- ghost triangles (distribute() column)')
        out(f'  {"np":>4}  {"scheme":<8}  {"ghost %":>8}')
        for np_val in args.np_values:
            for scheme in args.scheme_list:
                r = results.get((np_val, scheme), {})
                out(f'  {np_val:>4}  {scheme:<8}  {_pct(r.get("ghost_pct")):>8}')

    out()
    out('=' * W)

    summary_path = Path(outdir) / 'summary.txt'
    summary_path.write_text('\n'.join(lines) + '\n')
    print(f'\n  Summary saved to {summary_path}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args   = parse_args()
    script = find_benchmark_script(args.script)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f'\nBenchmark grid: np={args.np_values}  schemes={args.scheme_list}')
    print(f'Output directory: {outdir.resolve()}')
    print(f'Script: {script}\n')

    ntri   = 4 * args.size * args.size
    failed = []

    for np_val in args.np_values:
        for scheme in args.scheme_list:
            outfile = output_path(outdir, np_val, scheme)
            if outfile.exists() and not args.force:
                print(f'[SKIP] np={np_val} scheme={scheme}  '
                      f'(file exists; use --force to re-run)')
                continue
            print(f'\n[RUN] np={np_val}  scheme={scheme}')
            ok = run_one(args.mpirun, script, np_val, scheme,
                         args.size, args.reps, args.interval,
                         outfile, dry_run=args.dry_run)
            if not ok:
                failed.append((np_val, scheme))

    if args.dry_run or args.no_summary:
        return

    results = {}
    for np_val in args.np_values:
        for scheme in args.scheme_list:
            outfile = output_path(outdir, np_val, scheme)
            if not outfile.exists():
                continue
            text = outfile.read_text()
            r = parse_output(text)
            if r.get('wall_dist') is not None:
                results[(np_val, scheme)] = r
            else:
                print(f'  [WARN] Could not parse {outfile}')

    if not results:
        print('\nNo parseable results found.')
        return

    if failed:
        print(f'\nWARNING: {len(failed)} run(s) failed: {failed}')

    print()
    print_summary(results, args, ntri, outdir)


if __name__ == '__main__':
    main()
