#!/usr/bin/env python3
"""Compare baseline vs prototype Ns=0 benchmark records."""
import json, sys

def load(p):
    with open(p) as f:
        return json.load(f)

def main():
    if len(sys.argv) < 3:
        sys.exit('usage: compare.py baseline.json prototype.json [gate_pct]')
    b, p = load(sys.argv[1]), load(sys.argv[2])
    gate = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0

    print(f"  {'':22} {'baseline':>12} {'prototype':>12} {'delta':>10}")
    print('  ' + '-' * 58)
    rows = [('wall min (s)', 'wall_min', 3),
            ('wall median (s)', 'wall_median', 3),
            ('wall stdev (s)', 'wall_stdev', 4),
            ('cells/s (max)', 'cells_per_s_max', 0)]
    for label, key, dp in rows:
        bv, pv = b[key], p[key]
        d = 100 * (pv - bv) / bv if bv else 0.0
        print(f"  {label:22} {bv:12.{dp}f} {pv:12.{dp}f} {d:+9.2f}%")

    print()
    # correctness fingerprints must be identical at Ns=0
    same = True
    for key in ('n_triangles', 'n_steps', 'volume', 'stage_min', 'stage_max'):
        bv, pv = b.get(key), p.get(key)
        ok = (bv == pv)
        same &= ok
        flag = 'OK ' if ok else 'DIFF'
        print(f"  {flag} {key:16} {bv!r} vs {pv!r}")

    # the headline number: use min, which is the least noise-prone estimator
    slow = 100 * (p['wall_min'] - b['wall_min']) / b['wall_min']
    noise = 100 * max(b['wall_stdev'] / b['wall_median'],
                      p['wall_stdev'] / p['wall_median'])
    print()
    print(f"  Ns=0 slowdown (min wall): {slow:+.2f}%   [run-to-run noise ~{noise:.2f}%]")
    print(f"  Gate: <= {gate:.1f}%  ->  {'PASS' if slow <= gate else 'FAIL'}")
    if not same:
        print("  WARNING: numerical fingerprints differ - Ns=0 must be bit-identical")
    if abs(slow) < noise:
        print(f"  (|slowdown| < noise: effect not resolvable at this sample size)")

if __name__ == '__main__':
    main()
