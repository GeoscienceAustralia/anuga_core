#!/usr/bin/env python3
"""
Focused Ns=0 regression benchmark for the tracer flux-kernel prototype.

Answers exactly one question: does adding tracer support to the flux kernel
slow down an ordinary (no-tracer) run?

Matches benchmarks/run_benchmarks.py domain setup so numbers are comparable,
but adds repeats + warmup because the effect we are looking for (~1%) is the
same order as run-to-run noise.

Usage:
    OMP_NUM_THREADS=1 python bench_tracer.py --size large --mode 1 --repeats 5 \
        --label baseline --out baseline.json
"""
import argparse, json, os, statistics, sys, tempfile, time, shutil

SCENARIOS = {
    'small':  dict(nx=50,  ny=50,  finaltime=200.0, yieldstep=50.0),
    'medium': dict(nx=150, ny=150, finaltime=100.0, yieldstep=25.0),
    'large':  dict(nx=300, ny=300, finaltime=50.0,  yieldstep=12.5),
}


def create_domain(nx, ny, mode, tmpdir):
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


def one_run(size, mode):
    cfg = SCENARIOS[size]
    tmpdir = tempfile.mkdtemp()
    try:
        domain = create_domain(cfg['nx'], cfg['ny'], mode, tmpdir)
        n_tris = domain.number_of_triangles
        t0 = time.perf_counter()
        for _ in domain.evolve(yieldstep=cfg['yieldstep'],
                               finaltime=cfg['finaltime']):
            pass
        wall = time.perf_counter() - t0
        steps = domain.number_of_steps
        # a cheap correctness fingerprint: total volume + stage extremes
        stage = domain.quantities['stage'].centroid_values
        vol = float((stage * domain.areas).sum())
        return dict(n_triangles=n_tris, n_steps=steps, wall_s=wall,
                    cells_per_s=n_tris * steps / wall,
                    volume=vol,
                    stage_min=float(stage.min()), stage_max=float(stage.max()))
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--size', default='large', choices=list(SCENARIOS))
    ap.add_argument('--mode', type=int, default=1)
    ap.add_argument('--repeats', type=int, default=5)
    ap.add_argument('--label', default='run')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()

    import anuga
    print(f"# {a.label}: anuga {anuga.__version__}", flush=True)
    print(f"# size={a.size} mode={a.mode} OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS','unset')}",
          flush=True)

    print("  warmup...", flush=True, end=' ')
    w = one_run(a.size, a.mode)
    print(f"{w['wall_s']:.2f}s  {w['n_triangles']} tris  {w['n_steps']} steps", flush=True)

    runs = []
    for i in range(a.repeats):
        r = one_run(a.size, a.mode)
        runs.append(r)
        print(f"  run {i+1}/{a.repeats}: {r['wall_s']:.3f}s  "
              f"{r['cells_per_s']:,.0f} cells/s", flush=True)

    walls = [r['wall_s'] for r in runs]
    cps = [r['cells_per_s'] for r in runs]
    res = dict(
        label=a.label, size=a.size, mode=a.mode,
        omp_threads=os.environ.get('OMP_NUM_THREADS', 'unset'),
        anuga_version=anuga.__version__,
        n_triangles=runs[0]['n_triangles'], n_steps=runs[0]['n_steps'],
        wall_min=min(walls), wall_median=statistics.median(walls),
        wall_stdev=statistics.stdev(walls) if len(walls) > 1 else 0.0,
        cells_per_s_max=max(cps), cells_per_s_median=statistics.median(cps),
        volume=runs[0]['volume'],
        stage_min=runs[0]['stage_min'], stage_max=runs[0]['stage_max'],
        raw_walls=walls,
    )
    spread = 100 * res['wall_stdev'] / res['wall_median'] if res['wall_median'] else 0
    print(f"\n  min {res['wall_min']:.3f}s  median {res['wall_median']:.3f}s  "
          f"stdev {res['wall_stdev']:.3f}s ({spread:.2f}%)")
    print(f"  volume fingerprint {res['volume']:.10e}")

    if a.out:
        with open(a.out, 'w') as f:
            json.dump(res, f, indent=2)
        print(f"  -> {a.out}")


if __name__ == '__main__':
    main()
