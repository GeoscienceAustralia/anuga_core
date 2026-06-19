#!/usr/bin/env python3
"""
ANUGA Benchmark Comparison Tool
---------------------------------
Compares two JSON results files produced by run_benchmarks.py and prints
a side-by-side delta table.

Usage
-----
    python benchmarks/compare_benchmarks.py before.json after.json

    # Show only speed changes > 5%
    python benchmarks/compare_benchmarks.py before.json after.json --threshold 5

    # List available result files
    python benchmarks/compare_benchmarks.py --list
"""

import argparse
import json
import os
import sys


def _load(path):
    with open(path) as fh:
        return json.load(fh)


def _pct(old, new):
    """Return signed % change from old to new, or None if old==0."""
    if old == 0:
        return None
    return (new - old) / old * 100.0


def _fmt_pct(pct):
    """Format a signed % change."""
    if pct is None:
        return '   n/a'
    return f'{pct:+.1f}%'


def _arrow(pct, invert=False):
    """▲/▼/= based on direction and whether lower-is-better."""
    if pct is None or abs(pct) < 0.5:
        return ' '
    improved = (pct < 0) if invert else (pct > 0)
    return '↑' if improved else '↓'


def compare(before, after, threshold=0.0):
    b_results = {r['name']: r for r in before['results']}
    a_results = {r['name']: r for r in after['results']}

    all_names = sorted(set(b_results) | set(a_results))

    print(f"\nBenchmark comparison")
    print(f"  Before : {before['git']['commit']}  branch={before['git']['branch']}  "
          f"t={before['timestamp']}")
    print(f"  After  : {after['git']['commit']}  branch={after['git']['branch']}  "
          f"t={after['timestamp']}")
    print()

    # Header
    col = '{:<28}  {:>7}  {:>4}  {:>10}  {:>8}  {:>10}  {:>8}  {:>10}  {:>8}'
    hdr = col.format(
        'Scenario', 'tris', 'mode',
        'wall(s)', 'Δwall',
        'cells/s', 'Δspeed',
        'peak MB', 'Δmem',
    )
    rule = '─' * len(hdr)
    print(hdr)
    print(rule)

    any_printed = False
    for name in all_names:
        b = b_results.get(name)
        a = a_results.get(name)

        if b is None:
            print(f'  {name:<28}  (only in after)')
            any_printed = True
            continue
        if a is None:
            print(f'  {name:<28}  (only in before)')
            any_printed = True
            continue

        d_wall = _pct(b['wall_time_s'], a['wall_time_s'])
        d_speed = _pct(b['cells_per_s'], a['cells_per_s'])
        d_mem = _pct(b['peak_rss_mb'], a['peak_rss_mb'])

        # Skip if all changes are below threshold
        def _mag(v):
            return abs(v) if v is not None else 0.0

        if threshold > 0 and max(_mag(d_wall), _mag(d_speed), _mag(d_mem)) < threshold:
            continue

        w_arrow = _arrow(d_wall, invert=True)
        s_arrow = _arrow(d_speed, invert=False)
        m_arrow = _arrow(d_mem, invert=True)

        print(col.format(
            name,
            f"{b['n_triangles']:,}", b['mode'],
            f"{a['wall_time_s']:.2f}",
            f"{_fmt_pct(d_wall)}{w_arrow}",
            f"{a['cells_per_s']:,.0f}",
            f"{_fmt_pct(d_speed)}{s_arrow}",
            f"{a['peak_rss_mb']:.1f}",
            f"{_fmt_pct(d_mem)}{m_arrow}",
        ))
        any_printed = True

    if not any_printed:
        print(f'  (no changes exceed {threshold}% threshold)')

    print(rule)
    print()

    # Summary statistics
    common = [(b_results[n], a_results[n])
              for n in all_names if n in b_results and n in a_results]
    if common:
        speed_pcts = [_pct(b['cells_per_s'], a['cells_per_s'])
                      for b, a in common if b['cells_per_s'] > 0]
        mem_pcts = [_pct(b['peak_rss_mb'], a['peak_rss_mb'])
                    for b, a in common if b['peak_rss_mb'] > 0]

        def _mean(lst):
            lst = [v for v in lst if v is not None]
            return sum(lst) / len(lst) if lst else None

        ms = _mean(speed_pcts)
        mm = _mean(mem_pcts)

        if ms is not None:
            dir_s = 'faster' if ms > 0 else 'slower'
            print(f'  Average speed change : {ms:+.1f}% ({dir_s})')
        if mm is not None:
            dir_m = 'less memory' if mm < 0 else 'more memory'
            print(f'  Average memory change: {mm:+.1f}% ({dir_m})')
        print()


def list_results():
    results_dir = os.path.join(os.path.dirname(__file__), 'results')
    if not os.path.isdir(results_dir):
        print('No results/ directory found. Run run_benchmarks.py first.')
        return

    files = sorted(
        f for f in os.listdir(results_dir) if f.endswith('.json')
    )
    if not files:
        print('No result files in benchmarks/results/')
        return

    print(f'Result files in {results_dir}:')
    print()
    col = '  {:<50}  {:>10}  {:>10}  {:>6}'
    print(col.format('File', 'Commit', 'Branch', 'Runs'))
    print('  ' + '─' * 80)
    for fname in files:
        path = os.path.join(results_dir, fname)
        try:
            data = _load(path)
            commit = data.get('git', {}).get('commit', '?')
            branch = data.get('git', {}).get('branch', '?')
            n = len(data.get('results', []))
            print(col.format(fname, commit, branch, n))
        except Exception as exc:
            print(f'  {fname:<50}  (error: {exc})')
    print()


def main():
    parser = argparse.ArgumentParser(
        description='Compare two ANUGA benchmark result files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('before', nargs='?', help='Before JSON file')
    parser.add_argument('after', nargs='?', help='After JSON file')
    parser.add_argument(
        '--threshold', type=float, default=0.0,
        help='Only show rows where any metric changes by at least this %% (default: 0)',
    )
    parser.add_argument(
        '--list', action='store_true',
        help='List available result files in benchmarks/results/',
    )
    args = parser.parse_args()

    if args.list:
        list_results()
        return 0

    if not args.before or not args.after:
        parser.error('Provide both before and after JSON files, or use --list')

    before = _load(args.before)
    after = _load(args.after)
    compare(before, after, threshold=args.threshold)
    return 0


if __name__ == '__main__':
    sys.exit(main())
