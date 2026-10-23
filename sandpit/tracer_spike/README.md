# Tracer spike — test scripts

Verification for the generic-tracer work on `spike/tracer-flux-ns0-benchmark`.
Run from a directory OUTSIDE the anuga source tree (these import `anuga`).

    OMP_NUM_THREADS=1 python test_recon.py       # 6/6  edge reconstruction
    OMP_NUM_THREADS=1 python test_tracer_ns1.py  # 5/5  flux kernel, Ns=1
    OMP_NUM_THREADS=1 python test_ns2.py         # 3/3  multi-class striding
    OMP_NUM_THREADS=1 python test_time_int.py    # 6/6  time integration via evolve()
    OMP_NUM_THREADS=1 python test_add_tracer.py  # 20/20 add_tracer() registration API

Benchmark (Ns=0 regression gate):

    OMP_NUM_THREADS=1 python bench_tracer.py --size large --mode 2 --repeats 5 \
        --label baseline --out baseline.json
    python compare.py baseline.json prototype.json 1.0

`test_add_tracer.py` covers the registration API rather than the kernels: name
-> index mapping, array shapes/dtype/contiguity, growth from Ns=1 to Ns=2
without losing the existing tracer, `c`/`m = h*c` consistency, and the guards.
Its section B is the regression test for the cached-C-struct trap — a tracer
added after `evolve()` has built the struct must still move.

The four kernel scripts above still wire the six arrays by hand; that is
deliberate, since they are the verification record for the kernels themselves.
New work should use `domain.add_tracer()`.

`--size large` is 360k triangles. Mode 1 = legacy CPU, mode 2 = unified/GPU.
Run baseline and prototype on the SAME build configuration — see the rebuild
warning in ../../../Projects/Sediment_Transport/HANDOVER.md.

## Where the tests went

Everything in this directory that was a `test_*.py` is now a pytest module in
`anuga/shallow_water/tests/`, so CI runs it:

    test_tracers.py          test_tracers_gpu.py       (tracers, via develop #270)
    test_sediment_*.py                                 (sediment, 15 modules)

Run them with

    pytest anuga/shallow_water/tests/test_tracers*.py
    pytest anuga/shallow_water/tests/test_sediment_*.py --run-fast

The `*_gpu.py` modules collect every mode 1 / mode 2 comparison. They skip
themselves on a GPU-offload build unless `ANUGA_GPU_TESTS_ISOLATED=1` is set --
the NVHPC runtime aborts a process that builds many mode-2 domains.

## What is left here

The Ns=0 benchmark harness, which is a TIMING gate rather than a correctness
test and so does not belong in the pytest suite:

    OMP_NUM_THREADS=1 python bench_tracer.py --size large --mode 1
    python compare.py baseline.json prototype.json

It exists to answer one question: does carrying the tracer machinery cost
anything when no tracers are registered? It has earned its keep -- it caught a
+2.88% regression from hoisting pointer loads out of a guard, at a point when
every correctness test was green.

`run_parallel_tracer.py` drives the MPI halo exchange, which pytest does not
cover either.
