# C/Cython Extension Audit — TODOs

Audit date: 2026-06-09 (branch `develop`). Re-verified 2026-09-17 and filed as
GitHub issues #323–#343 (one per open item; the issue number follows each
title below). Two items were found already fixed and are ticked.

**Status (2026-09-17): complete.** Every item is closed, in four PRs by
priority — P1 correctness #345, P2 #346, P3 performance #347, P4 hygiene
#349 (plus #348, a benchmark follow-up). The one thing left outside the code
is #333's self-hosted GPU runner: the workflow is on develop, the runner has
to be registered. `scripts/c_static_analysis.sh` (gcc -fanalyzer, run by
CI) now guards against the P1 class of finding recurring. Line numbers
below are as of the audit and have drifted; follow the PR for each item.

Scope: all 20 compiled
extensions (meson targets), their C sources, the Python↔C boundary, and test
coverage.

## Extension inventory (context)

20 extensions built via meson; 19 actively imported. Two execution modes share
`gpu/core_kernels.c` (mode 1 compiles it with `-DCPU_ONLY_MODE`). 9 extensions
use OpenMP. `quad_tree.c` is compiled into both `quad_tree_ext` and
`fitsmooth_ext`.

---

## P1 — Correctness / safety (do first)

- [x] **Unchecked division by zero in `_gradient2`** — #323 — DONE in #345 —
  `anuga/utilities/util_ext.h:146`: `det = xx*xx + yy*yy; //FIXME catch det == 0`
  then divides by `det` unconditionally. Degenerate (zero-length) edges produce
  NaN that propagates into friction/gravity terms. Guard with `D->epsilon` (or
  return an error). The FIXME has been in the source for years.

- [x] **CG solver divide-by-zero** — #324 — DONE in #345 — `anuga/utilities/cg.c:147`
  (`z[i]=1.0/D[i]*x[i]` — Jacobi preconditioner with no `D[i]==0` check) and
  `anuga/utilities/cg.c:256` (`alpha = rTr/cg_ddot(M,d,q)` — singular/indefinite
  matrix gives 0 denominator). Both silently produce NaN/Inf instead of a
  convergence failure. Affects `kinematic_viscosity_operator`.

- [x] **No shape validation when wiring domain pointers** — #325 — DONE in #345 —
  `anuga/shallow_water/sw_domain_openmp_ext.pyx` (`get_python_domain_pointers`,
  ~lines 145–515) extracts 80+ raw pointers with no checks that
  `centroid_values.shape[0] == number_of_elements`, `edge_values.shape == (N,3)`,
  etc. A malformed quantity array means out-of-bounds reads/writes in the flux
  kernels. Same pattern in the `manning_friction_*` wrappers
  (sw_domain_openmp_ext.pyx:1032–1077): only `w.shape[0]` is read; sibling
  arrays are never length-checked. Add cheap assertions at the boundary — this
  is once-per-call, not per-element, so cost is negligible.

- [x] **Static substep counters are shared mutable state** — #326 — DONE in #345 —
  `anuga/shallow_water/sw_domain_openmp.c:104–106`: `static anuga_int call,
  timestep_fluxcalls, base_call` inside `_openmp_compute_fluxes_central`. Breaks
  with two domains evolving in one process (counters interleave) and races if
  called from multiple threads with the GIL released. Move into `struct domain`.
  (Check whether the standalone/CMake `libanuga_sw` build inherited this.)

- [x] **Thread-unsafe globals in URS reader** — #327 — DONE in #345 — `anuga/file/urs.c:24–28`:
  `static int32_t *fros, *lros; static struct tgsrwg* mytgs0; static anuga_int
  numDataMax`. Concurrent `read_mux2` calls corrupt each other. Refactor to a
  context struct passed through the call chain (or document single-threaded-only).

- [x] **GPU error paths leak mapped memory** — #328 — DONE in #345 —
  `anuga/shallow_water/gpu/gpu_halo.c:37–97`: sequential mallocs with no cleanup
  of earlier allocations on failure. `gpu_inlet_operator.c:116–177`: early
  `fprintf(...); return` paths between `omp target enter data` and the matching
  `exit data` leave device buffers mapped forever. Use a goto-cleanup pattern.

- [x] **Error codes swallowed at the Cython layer** — #329 — DONE in #345 —
  e.g. `sw_domain_openmp_ext.pyx:1204–1218`: `gravity()`/`gravity_wb()` turn a
  `-1` return into a silent `return None`. Sweep all wrappers: a C error code
  should raise a Python exception with context, never return None.

- [x] **`exit()` in library code** — #330 — DONE in #345 — `anuga/utilities/quad_tree.c:9–12`,
  `anuga/utilities/sparse_dok.c:9–12` (`emalloc` calls `exit(EXIT_FAILURE)` on
  OOM), killing the whole Python interpreter (and every MPI rank). Also raw
  `malloc` without NULL checks at `sparse_dok.c:54`, `quad_tree.c:516,528`.
  Return error codes up to Cython and raise `MemoryError`.

## P2 — Silent-wrong-answer and divergence risks

- [x] **Silent GPU→CPU fallback** — #331 — DONE in #346 —
  `shallow_water_domain.py:5107–5150`: `set_gpu_interface()` auto-downgrades
  `multiprocessor_mode` 2→1 when no GPU is found, with only a printed warning.
  On HPC batch runs this means a "GPU job" silently burns CPU hours. Add a
  strict mode (e.g. `set_multiprocessor_mode(2, strict=True)` or an env var)
  that raises instead. Same for the per-boundary runtime fallbacks
  (shallow_water_domain.py:3102, 3384, 3883, …) — warn once at init, not
  mid-evolve.

- [x] **Finish the CPU/GPU kernel unification** — DONE 2026-06-09.
  `core_gravity_wb` was a stub silently falling back to `core_gravity`; real
  well-balanced body ported. `core_gravity` read stale `height_centroid_values`;
  now uses `stage - bed`. `core_fix_negative_cells` had diverged (broader
  `minimum_allowed_height` threshold, ignored `tri_full_flag`); rewritten to
  match the CPU reference. New `core_manning_friction_sloped_semi_implicit_edge_based`
  ported with dual-mode pragmas. `_openmp_*` versions are now thin wrappers.
  Explicit (non-semi-implicit) Manning variants kept as legacy CPU-only — not
  reachable from production code. `_openmp_evaluate_reflective_segment` left
  alone (GPU has its own boundary architecture). NOTE: the fix_negative_cells
  and gravity changes alter GPU-mode behaviour (now matches CPU reference) —
  re-validate a GPU run.

- [x] **GPU mode applies WRONG friction for sloped-Manning domains** — DONE: `gpu_kernels.c` now branches on `GD->use_sloped_mannings` (lines 220, 244, 419) and `friction.py` dispatches on it; verified 2026-09-17. (found
  during unification): `gpu_kernels.c:181` `gpu_manning_friction` always calls
  the *flat* semi-implicit variant regardless of `use_sloped_mannings`, and
  `friction.py:39–48` falls back to the CPU edge-based function against
  possibly-stale host arrays while data is device-resident. Fix (4 steps):
  add `use_sloped_mannings` to `struct domain` in `sw_domain.h`; populate it in
  `sw_domain_gpu_ext.pyx` init; branch in `gpu_manning_friction` to the new
  `core_manning_friction_sloped_semi_implicit_edge_based`; remove the CPU
  fallback in friction.py. (`bed_edge_values` is already mapped to device.)

- [x] **Epsilon drift between C and config.py** — #332 — DONE in #346 (the `eta > 1.0e-16` literal is gone; the duplicated constants and hfactor magic numbers remain) —
  `sw_domain_openmp.c:231` tests `eta > 1.0e-16` (a literal) while the domain
  carries `D->epsilon` (1.0e-12 from config). Decide which threshold is intended
  and name it. Also hoist the duplicated `seven_thirds`/`one_third` constants
  (4 copies across the Manning functions) and the wet-dry hfactor magic numbers
  in `gpu/core_kernels.c:30–33` (`a_tmp=0.3, b_tmp=0.1, …`) into one header /
  the domain struct so they're tunable and single-sourced.

- [x] **GPU test coverage is conditional and thin** — #333 — DONE in #346 (workflow merged; a self-hosted GPU runner still has to be registered) ((a) is done — the CPU build runs the mode-2 path in-process; (b) the hardware runner remains) — the only CPU-vs-GPU
  consistency tests (`anuga/shallow_water/tests/test_DE_gpu_omp.py`, parallel
  `test_parallel_sw_flow_gpu_de*.py`) skip entirely when no GPU is present, so
  CI never exercises mode 2. Two TODOs: (a) run the mode-2 code path with
  `OMP_TARGET_OFFLOAD=DISABLED`/CPU fallback in CI so the dispatch + kernels are
  at least compiled and executed; (b) add a nightly/manual GPU runner that runs
  the consistency suite on real hardware.

- [x] **No direct unit tests for `quantity_openmp_ext`** — #334 — DONE in #346 — gradient
  computation, extrapolation, and the limiters (`limit_edges_by_all_neighbours`
  etc.) are only tested through whole-domain evolution. Add small fixed-mesh
  unit tests with hand-checkable values; these kernels are where limiter bugs
  hide.

## P3 — Performance

- [x] **False sharing / shared(D) in older OpenMP loops** — DONE: no `shared(D)` pragmas remain in `sw_domain_openmp.c`; verified 2026-09-17. — newer
  `core_kernels.c` correctly hoists `restrict` pointers and uses
  `firstprivate`; older functions in `sw_domain_openmp.c` (e.g. the friction
  loops, `#pragma omp parallel for simd default(none) shared(D)`) dereference
  the domain struct inside the loop. Hoist hot pointers to locals — measurable
  on multi-socket nodes. Benchmark with `benchmarks/run_benchmarks.py`
  before/after.

- [x] **malloc in the fitting hot loop** — #335 — DONE in #347 — `anuga/utilities/quad_tree.c:82–94`
  mallocs a 3-double sigma array per query point inside the OpenMP fit loop
  (`fitsmooth.c:199–235`), which also contains a critical section. Pass a
  stack/caller buffer instead (the FIXME in the source says exactly this).
  Affects `fit_to_mesh` / set_quantity-from-points performance.

- [x] **Unparallelized polygon loop** — #336 — DONE in #347 — `anuga/geometry/polygon.c:699–700`
  main separate-points loop has a `// TODO, JLGV: Use OpenMP` comment. Matters
  for large point clouds in `inside_polygon` during region setup.

- [x] **Per-kernel microbenchmarks** — #337 — DONE in #347 — `benchmarks/` measures whole-evolve
  throughput (cells/s) only. Add timed harnesses for compute_fluxes /
  extrapolate / distribute / protect individually so kernel-level regressions
  are attributable. (The gpu_flop.c instrumentation could feed this.)

- [x] **GIL held in quantity/mesh extensions** — #338 — DONE in #347 — extern C calls in
  `quantity_openmp_ext.pyx` aren't declared `nogil`. Low impact today
  (single-threaded driver), but blocks any future threaded ensemble use and is
  cheap to fix while touching the file.

## P4 — Cleanup / hygiene

- [x] **Delete or quarantine dead code** — #339 — DONE in #349 — `mannings_operator_ext` is built but
  imported nowhere (friction lives in the sw extensions); remove the meson
  target + sources or document why it stays. Uncompiled orphans:
  `anuga/fit_interpolate/p_test.c`, `ptinpoly.c/.h`, `rand48.c`.

- [x] **Modernize Python-2-isms in .pyx** — #340 — DONE in #349 — `xrange` in
  `fitsmooth_ext.pyx:80,93,149,166,239,241` and
  `sparse_matrix_ext.pyx:53,77`. Works only because Cython optimizes typed
  loops; will break the day a loop variable becomes untyped.

- [x] **Replace fprintf/printf with proper error propagation** — #341 — DONE in #349 — 15+
  `fprintf(stderr, ...)` sites under `gpu/`, plus `print_*_array` helpers in
  `util_ext.h:189–200`. On HPC, per-rank stderr is often discarded; errors
  should surface as return codes → Python exceptions.

- [x] **Add static analysis to CI** — #342 — DONE in #349 — a `clang-tidy` / `gcc -fanalyzer` pass
  over `anuga/**/ *.c` would have caught most of the P1 items mechanically.
  One-time setup, ongoing payoff.

- [x] **Harvest in-source FIXMEs** — #343 — DONE in #349 — besides those above:
  `operators/kinematic_viscosity_operator.c:7` ("replace with library call"),
  `operators/mannings_operator.c:44,92` (Taylor expansion),
  `geometry/polygon.c:671` ("pass rtol/atol from Python").

## Verified non-issues (don't re-flag)

- `polygon.c:747` point-in-polygon division: the crossing test guarantees
  `py_j != py_i`; no divide-by-zero despite appearances.
- NumPy C API: builds with `-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION`; no
  deprecated macros or `->data` access found.
- Contiguity at the Cython boundary: `mode="c"` / `[::1]` is consistently
  enforced, so non-contiguous input fails loudly rather than corrupting.
- `_openmp_fix_negative_cells` uses a proper `reduction(+:...)` clause, not a
  critical section.
