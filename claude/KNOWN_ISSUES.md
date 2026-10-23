# Known Issues and Gotchas

Things discovered during development sessions that are surprising, non-obvious,
or require caution when working in specific areas.

---

## Sediment

### `'depth_slope'` shear closure is self-amplifying on an evolving bed (2026-09-18)

`[T-7]` `tau_b = rho g h S` takes `S` from the bed itself. Once the bed evolves,
any non-uniform erosion steepens local slopes, which raises `tau_b`, which erodes
faster: a tank of still water on a 1e-3 slope scoured 0.8 m in 60 s and a
fraction with `tau_c* = 10` (30x above the stress the imposed slope gives) was
entrained. anugaSed keeps this in check with the undocumented domain-global
clamp `S <- min(S, mean(S)/2)` that PHYSICS_SPEC D1a describes; ANUGA does not
apply that. Since 2026-09-20 `set_shear_closure` offers two explicit bounds:
`max_slope` (a per-cell cap) and `freeze_slope=True` (S from the bed at setup,
kept for the run), the latter being what the sediment_erosion validation case
uses to evolve its bed. Without either, use `'quadratic_drag'` (default) or
`'energy_slope'` when the bed moves.

## Structures

### A culvert in still water amplifies roundoff (not a mode-1 vs mode-2 problem)

A Boyd structure whose inlets sit in near-still water is an unstable configuration: the
operator reads the enquiry cells and writes the inlet cells, and the two are close enough
to close a feedback loop. With no head across the structure the discharge is decided by
roundoff, and the perturbation grows exponentially — from 1 ULP to ~1e-2 m of stage in a
few seconds on a 200 m x 50 m test domain.

**It is not a GPU port problem.** Perturbing a single enquiry cell by 1e-15 in mode 1
diverges from unperturbed mode 1 exactly as fast as mode 2 does (measured: 5.4e-06 at
t=0.25 s, 5.3e-03 at t=1 s — the same as |mode1 - mode2| over the same run). A CPU/GPU
comparison over such a configuration measures the instability, not the port.

Consequences when writing tests:

* A GPU-vs-CPU culvert test needs a **genuine driving head**, so the discharge is
  deterministic and dominant. `Test_GPU_LargeInlet` in `test_DE_gpu_omp.py` does this
  (ponded upstream, low downstream) and then agrees to ~1e-14. An earlier version of it
  used still water and was hostage to this effect.
* A well-balancedness test (issue #229) must be kept **short** (~1 s) or use enquiry
  points well clear of the inlets, or the instability grows into the measurement.

With the #229 fix (stage leveling) the old uniform-depth perturbation is gone, so this
effect is no longer masked — it is now the dominant behaviour of a structure in still
water, and the constraints above apply to any test touching one.

---

## Build

### Building with GPU offloading (NVIDIA HPC SDK / nvc)

ANUGA's GPU extension (`sw_domain_gpu_ext`, `multiprocessor_mode=2`) requires
`nvc` from the NVIDIA HPC SDK — GCC 15's nvptx backend ICEs on `core_kernels.c`
(ompdevlow GIMPLE pass segfault in `core_extrapolate_second_order_edge`).

**One-time setup (Ubuntu, requires sudo):**
```bash
# Add NVIDIA HPC SDK apt repo
curl -fsSL https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK \
  | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' \
  | sudo tee /etc/apt/sources.list.d/nvhpc.list
sudo apt-get update -y && sudo apt-get install -y nvhpc   # ~5 GB
```

**GPU build (RTX 5070 = Blackwell cc120; adjust gpu_arch for other GPUs):**
```bash
NVC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3/compilers/bin/nvc
conda run -n anuga_env_3.14 bash -c "CC=$NVC pip install --no-build-isolation -v -e . \
  -Csetup-args=-Dgpu_offload=true \
  -Csetup-args=-Dgpu_arch=cc120"
```

Meson auto-detects nvc as `nvidia_hpc`; the build uses `-mp=gpu,multicore -gpu=cc120`.
The build dir must be clean if switching from a prior GCC build (`rm -rf build/cp314`).

**Verify GPU works:**
```bash
conda run -n anuga_env_3.14 pytest anuga/shallow_water/tests/test_DE_gpu_omp.py -v
```
All 56 tests pass on the RTX 5070.

**Switching back to CPU-only build:**
```bash
rm -rf build/cp314
conda run -n anuga_env_3.14 pip install --no-build-isolation -e .
```

### `--no-build-isolation` is recommended

`pip install --no-build-isolation -e .` is the recommended build approach.
It is not strictly required in all environments, but is preferred because
it ensures meson-python uses the Cython/numpy already installed in the conda
environment rather than fetching isolated build dependencies.

### Generated C files appear as untracked in `git status`

`sw_domain_openmp_ext.c` and other generated `.c` files are listed in
`.gitignore` but still show up as untracked. This is expected — they are
build artifacts.

### A reused meson build dir does NOT re-detect MPI (2026-06-25)

`sw_domain_gpu_ext` is built with real C MPI (`HAVE_MPI4PY=True`, multi-rank
GPU halo exchange) only when meson finds MPI **at configure time**. Detection is
two-stage: `dependency('mpi', language: 'c')` first, then — because meson's
`mpi` dependency does **not** match conda's `mpich` pkg-config name — an
**mpi4py fallback** that parses `mpicc -show` to locate `mpi.h`
(`anuga/shallow_water/meson.build`). The fallback runs via `run_command`, which
is evaluated only on a **fresh configure**.

Gotcha: if you `pip install -e .` *before* MPI/mpi4py is in the env, the gpu
extension is compiled against the single-process stubs (`gpu_mpi_stubs.h`),
`gpu_has_mpi()` returns False, and the four
`anuga/parallel/tests/test_parallel_sw_flow_gpu_*` tests **skip**
("GPU extension built without C MPI"). Installing MPI afterwards and re-running
plain `pip install -e .` does **not** fix it — meson-python reuses the cached
build dir and just relinks the no-MPI `.so` (it never re-runs the fallback).
This is independent of `gpu_offload`: it bites the standard gcc CPU-only build,
where these tests otherwise run mode-2 on the host.

Fix — force a fresh configure by pointing at a new build dir (or deleting the
cached one):

```bash
CC=gcc pip install --no-build-isolation -e . \
  -Csetup-args=-Dgpu_offload=false \
  -Cbuild-dir=build/cp314-mpi -v
```

Verify (any one is sufficient):
- meson logs `GPU extension will be built WITH MPI support (multi-GPU enabled)`
  (visible with `-v`);
- `readelf -d <sw_domain_gpu_ext...so> | grep NEEDED` lists `libmpi.so.*`;
- `python -c "from anuga.shallow_water import sw_domain_gpu_ext as e; print(e.gpu_has_mpi())"` → `True`.

Then the `test_parallel_sw_flow_gpu_*` files run (`real_gpu_available()` stays
False on a CPU build, so the per-test "needs N GPUs" guards do not fire either).

### A reused meson build dir keeps the old compiler — switching gcc↔nvc needs `rm -rf build/cp*` (2026-07-01)

Same root cause as the MPI note above, but for the **compiler**. meson-python
reuses `build/cp<ver>` and only reads `CC` on the **first** configure of a dir;
a later build just runs `meson setup --reconfigure`, which keeps the originally
detected compiler. So building for GPU (`CC=nvc pip install -e . -Dgpu_offload=true`)
in a tree that already has a gcc-configured `build/cp314` **stays on gcc**, and
`anuga/shallow_water/meson.build`'s guard aborts:

```
C compiler for the host machine: cc (gcc 15.2.0)
ERROR: gpu_offload=true is not supported with gcc ... rm -rf build/cp314 required when switching compiler
```

The reverse bites too (an nvc-configured dir stays on nvc for a later gcc CPU
build). **Fix: remove the build dir before switching compiler** so `CC` is read
on a clean configure:

```bash
rm -rf build/cp*                     # force a fresh meson configure
CC=$(which nvc) pip install --no-build-isolation -e . \
    -Csetup-args=-Dgpu_offload=true -Csetup-args=-Dgpu_arch=cc120
```

`tools/install_anuga_nvc.sh` now does this `rm -rf build/cp*` automatically
before the nvc build. Verified end to end: fresh nvc build succeeds and the
isolated GPU runner reports 65/65 passed. See also `SESSION_GUIDE.md` → "CPU and
GPU are separate builds".

---

## Testing

### `str.find()` returns 0 for first-position match (2026-03-26)

In `anuga/pmesh/mesh.py::_generateMesh_impl`, the old check `not self.mode.find('Q')`
was buggy: `str.find()` returns 0 when 'Q' is at position 0, and `not 0` is `True`,
so the check was treating a 'Q' at position 0 as "not found". Fixed by using
`'Q' not in self.mode`.

### Triangle library prints to stdout during pytest

The triangle C library writes to stdout when in verbose mode ('V' flag). Since
pytest `-s` does not suppress stdout, these appear as noise during test runs.
Fixed by ensuring `_generateMesh_impl` adds 'Q' (quiet) when `verbose=False`.

### `test_verbose_does_not_raise` triggers logging output

`anuga/abstract_2d_finite_volumes/tests/test_pmesh_to_mesh.py::test_verbose_does_not_raise`
intentionally calls with `verbose=True`. This triggers `General_mesh:` log output.
Fixed by wrapping with `logging.disable(logging.CRITICAL)` / `logging.disable(logging.NOTSET)`.

### Parallel tests run as subprocesses

Tests in `anuga/parallel/tests/` spawn `mpiexec` subprocesses. They cannot be
parallelised with `pytest-xdist` and must run serially. They are marked slow
and skipped by `--run-fast`.

### GPU build: `test_DE_gpu_omp.py` aborts mid-file (NVHPC target present-table) (2026-06-15)

On a GPU build (`-Dgpu_offload=true`, nvc), running the whole
`anuga/shallow_water/tests/test_DE_gpu_omp.py` file in one process **aborts
silently** (exit 1, no traceback, not a SIGSEGV, GPU idle/no-OOM) partway
through — around the 9th–11th test (`Test_GPU_InletOperator::test_inlet_operator_basic`).
The NVHPC OpenMP-target runtime calls `exit()`.

**Not caused by the mode-2 session changes** — reproduces identically on the
pre-change commit (`d96ae357`) rebuilt with nvc.

**It is NOT a simple cumulative-resource leak.** Diagnostics:
- Each test class *alone*, and `test_inlet_operator_basic` alone, pass.
- Creating 16–20 GPU domains in a loop — both dropped each iteration *and* kept
  simultaneously live — works fine.
- The crash only appears with the file's specific mix of low-level kernel tests,
  `set_multiprocessor_mode(1)`↔`(2)` switching, `sync_to/from_device`, and
  operator setup, accumulated across ≥9 tests.
- Forcing finalization between tests (`gc.collect()`, or nulling
  `domain.gpu_interface`) makes it **worse**: introduces assertion *failures*
  before the abort.

**Root cause (diagnosed, not yet fixed):** device arrays are bound with
`#pragma omp target enter data map(to: host_ptr...)` keyed on the *host* pointer,
and released with `map(delete:)` (`gpu_domain_unmap_arrays` /
`gpu_domain_finalize`). The OpenMP present table is reference-counted and
host-pointer-keyed, so repeated map/unmap of arrays whose host addresses numpy
recycles across domains corrupts the table (stale entry reused, or a live
entry deleted) — hence both the "leak then abort" and the "eager-unmap then
assertion failure" signatures. Two reference cycles
(`domain ↔ gpu_interface`, and `gpu_dom → python_domain → domain`) also defer
`GPUDomain.__dealloc__`/`gpu_domain_finalize` to the cyclic GC, so finalization
timing is non-deterministic — but breaking either cycle does not fix the
underlying present-table issue (the other cycle still pins the domain, and eager
finalize corrupts).

It also crashes the *whole* `pytest --pyargs anuga.shallow_water` run (in either
`ANUGA_DEFAULT_COMPUTE_MODE`) because `test_DE_gpu_omp.py` collects early and
these GPU tests set mode 2 explicitly regardless of the env default — the abort
at ~3% kills the run before the rest of the suite executes.

**On a GPU build, `ANUGA_DEFAULT_COMPUTE_MODE=unified` over the full suite is not
viable even with the GPU file excluded:** every default domain then offloads to
the GPU, so the whole suite churns hundreds of mode-2 GPU domains in one process
and aborts early (~20%). This is the same root cause. Validate the unified
default over the full suite on the **gcc CPU build** (`-Dgpu_offload=false`),
where it is the documented 2657-passed run; on a GPU build, drive the bulk suite
with `ANUGA_DEFAULT_COMPUTE_MODE=legacy` and cover GPU paths via the per-class
runner below.

**Impact:** production use (a single, or a few sequential, mode-2 GPU domains)
is unaffected — single-domain evolve and each test class pass. Only running many
GPU-domain tests in *one* process trips it.

**Auto-skip:** on a GPU-offload build, `test_DE_gpu_omp.py` skips itself at
collection (module-level `pytest.skip`, gated on `anuga.gpu_offload_supported()`
and `not ANUGA_GPU_TESTS_ISOLATED`), so a normal `pytest --pyargs anuga` no longer
crashes — the file is reported as skipped with a message pointing here. On a CPU
build the guard is inert and the file runs in-process as usual.

**Workaround — run the GPU tests in isolated processes:**
```bash
# one fresh process per CLASS (fast):
bash anuga/shallow_water/tests/run_gpu_tests_isolated.sh
# one fresh process per TEST FUNCTION, with a per-test timeout (most robust;
# turns a genuine hang into a reported TIMEOUT). Works on any pytest target.
# Installed as `anuga_run_isolated_tests` (scripts/, via meson); in a source
# checkout run scripts/anuga_run_isolated_tests.py directly:
anuga_run_isolated_tests [TARGET] [--timeout S] [-k EXPR]
```
Both set `ANUGA_GPU_TESTS_ISOLATED=1` to bypass the auto-skip; all tests pass
this way (verified 65/65 per-function on the nvc build). Then run the rest of the
suite normally (it does not trip the issue) under the **legacy** default:
```bash
ANUGA_DEFAULT_COMPUTE_MODE=legacy \
  pytest anuga/shallow_water/tests/ --ignore=anuga/shallow_water/tests/test_DE_gpu_omp.py
```
**Do NOT use `pytest --forked`** for these tests: CUDA contexts are fork-unsafe,
so forking from a GPU-initialised parent poisons every child (it turns the abort
into ~53 spurious failures). Isolation must be *fresh* processes (separate
`python -m pytest` invocations), not `os.fork()`.

A real fix (FUTURE_WORK P1.10) needs either per-test fresh-process isolation
baked into the GPU test file, strict 1:1 map/unmap reference-count discipline per
domain, or device-pointer allocation (`omp_target_alloc` + `is_device_ptr`)
instead of host-pointer-keyed `map(to:)`. (`omp target enter/exit data` cleanup
is reference-counted and host-pointer-keyed; forcing finalization between tests
removes the abort but still yields ~7 aliasing failures, so clean teardown alone
is not sufficient.)

### GPU build: `anuga.shallow_water` is green under `unified` via the isolated runner (2026-06-17)

The per-function isolated runner now passes the **entire** `anuga.shallow_water`
set under the unified default on a GPU-offload build:

```bash
anuga_run_isolated_tests --pyargs anuga.shallow_water -cm unified
# 410 collected -> pass=408 skip=2 (2 skips are pre-existing legacy-default guards)
# (-cm/--compute-mode sets ANUGA_DEFAULT_COMPUTE_MODE for every child; omit to
#  inherit the environment.)
```

This works because each test runs in its own fresh process (no mode-2 domain
accumulation -> no NVHPC abort), **and** because 11 tests that probed mode-1-only
host state are now pinned to `legacy` (`domain.set_compute_mode('legacy')`):

- 9 white-box tests call `compute_forcing_terms()` / `compute_fluxes()` and assert
  on the host `semi_implicit_update` / `explicit_update` arrays, which mode-2 GPU
  computes on-device and never syncs back (so the host arrays read stale zeros) —
  in `test_forcing.py`, `test_friction.py`, `test_physics_sw.py` (Manning friction
  cases) and `test_data_manager.py::test_sww_extrema` (extrema monitoring).
- 2 numerical tests compare against legacy-recorded references and diverge at the
  ~1e-6 level under mode-2's different reduction/eval order
  (`test_regression_snapshots.py::test_dam_break_DE1_stage_snapshot` and
  `test_sww_interrogate.py::test_get_maximum_inundation_de0`). The two
  regression-snapshot domain helpers are pinned so that whole file stays
  deterministic under any `ANUGA_DEFAULT_COMPUTE_MODE`.
- (Session 52) `test_negative_cells_warning.py` — its `make_domain()` helper calls
  `domain.update_conserved_quantities()` **directly, outside `evolve()`** and reads host
  `centroid_values`; pinned to `legacy`. Under the unified default this had failed with
  `AttributeError: 'NoneType' object has no attribute 'update_conserved_quantities_kernel'`
  because `gpu_interface` is only built during `evolve` (commit `641db3bd`). The method
  itself was also hardened: it now calls `_ensure_gpu_interface()` first, like every other
  mode-2 entry point, so any direct call builds the interface (or falls back to legacy)
  rather than crashing on a `None` (commit `c80bc457`).

These are test-harness artifacts, not solver bugs; the pins are no-ops for the
distribution-default legacy path. Mode-2 numerical fidelity remains covered by the
mode1-vs-mode2 comparison tests in `test_DE_gpu_omp.py`. Note this complements —
does not replace — the guidance above: the *full* `pytest --pyargs anuga.shallow_water`
(non-isolated) under `unified` on a GPU build still aborts; use the isolated runner.
Commit `0c50947d`.

### Targeted `--cov=anuga.submodule` runs corrupt numpy's `_NoValue` sentinel

Running `pytest --cov=anuga.structures.structure_operator` (or any sub-package
path) causes test failures with:
```
TypeError: float() argument must be a string or a real number, not '_NoValueType'
```

**Root cause:** coverage.py calls `importlib.util.find_spec('anuga.structures.structure_operator')`
inside a `sys_modules_saved()` context (in `inorout.py`). This auto-imports parent
packages (including `anuga/__init__.py` → `shallow_water_domain.py` → numpy),
then purges all newly-imported modules from `sys.modules`. The subsequent real
import re-executes `numpy/__init__.py`. Since numpy's C extension (`_multiarray_umath`)
was already initialized, the reload guard fires and a new `_NoValue` singleton is
created — but C extensions hold references to the old one, breaking identity checks.

**Workaround:** Always use `--cov=anuga` (not a sub-path). For per-module numbers:
```bash
pytest --run-fast --cov=anuga anuga/structures/tests/ -q 2>&1 | grep structure_operator
```

**Not fixable** from conftest.py: pytest-cov creates `CovPlugin` (which starts
coverage) inside `pytest_load_initial_conftests(tryfirst=True)` — before conftest.py
is even loaded.

---

## Numerical

### `== None` vs `is None` with numpy arrays

Using `== None` on a numpy array raises `ValueError: The truth value of an array
is ambiguous`. Always use `is None` / `is not None` throughout the codebase.

### `epsilon = 1.0e-6` wet/dry threshold

`anuga/config.py` defines `epsilon` as the wet/dry threshold. Many conditional
checks use `depth > epsilon` rather than `depth > 0`. Be aware of this when
writing new flux/operator code.

### `minimum_allowed_height = 1.0e-05`

Cells below this height are treated as dry. Negative depths are clipped.

### mode-1 vs mode-2 differ by ~1 ULP at the wet/dry margin (2026-07-20)

Both modes are **stage-primary** (`height == stage - bed` exactly in each) — there
is *no* stage-vs-height representation flip. The residual is a **1-ULP difference
in the core's near-dry stage value itself**: for a cell sitting essentially *at*
the bed that gets lifted a hair off it (e.g. rain on dry ground, terrain ~343 m +
~1e-5 m of rainfall), mode-2's stage going into the fractional step is already ~1
ULP off from mode-1, and it is **masked by the dry-cell bed-clamp** (`protect` sets
`stage = bed` exactly in both) until the lift exposes it. Only in the just-wetted
cells, zero in momentum; chaos-amplifies to mm over hours. Ruled out along the way:
the rate operator's arithmetic (rate inputs are bit-identical; rain on a *fully
wet* domain is bit-identical) and the device sync (a plain memcpy). This is why
**towradgi's mode-1/mode-2 divergence is rain-on-dry-cells** — remove rain and the
whole run is bit-identical. Same family as the #200 dry-cell gap; benign (both
modes valid). The exact operation was not isolated (it needs instrumenting the C
RK loop's device state mid-step). The session-51 in-process double-precision
localization harness (build two domains from the same setup, restartable lockstep
evolve, diff `centroid_values` in double precision) is the tool if it's ever chased.

### DE_ader2's steady state depends on dt, so it never quite settles (2026-10-01)

DE1 and DE2 evaluate every flux on a state reconstructed from Q alone, and dt
only multiplies R(Q): their discrete steady state is R(Q) = 0 whatever the
timestep (dt 0.04 vs 0.08 s: steady depths agree to 1e-15). DE_ader2 evaluates
the fluxes on the predicted state Q + (dt/2) Q_t, where Q_t is the local
Cauchy-Kovalewski derivative from each cell's own gradients, not the discrete
residual. At a discrete steady state R = 0 but Q_t != 0 wherever the flow is not
uniform, so the steady state solves R(Q + (dt/2) Q_t) = 0 and moves with dt. This
is the familiar Lax-Wendroff / MUSCL-Hancock property. Measured on the
rundown_mild_slope stable case (100 m, n = 0.06), dt 0.04 vs 0.08 s moves the
depth by 4e-3 h0 in the cells next to the inflow and outflow, 7e-4 at 10-20 m
and 2e-5 in the interior. Uniform flow itself is dt-independent since #404 put
the friction in the predictor.

So every change in dt sends a small pulse downstream from the boundary cells.
The predictor also uses the *previous* step's dt (`prev_dt * 0.5` in
`evolve_one_ader2_step`), which adds a mismatch whenever dt changes. There are
two sources of dt changes:

- **Clipping to output times.** The step before each yield is cut short. On a
  stable flow this alone keeps the inlet at about 1e-4 h0.
- **Feedback in unstable flow.** Above the roll-wave threshold (Froude > 1.5 for
  Manning), the start-up transient's dt changes seed roll waves. The waves then
  swing the global CFL dt by about 17%, which kicks the inlet again. The loop
  sustains itself: about 1e-3 h0 near the inlet even with no output clipping.
  DE1 and DE2, under the same dt swings, stay steady to 1e-14 and grow no roll
  waves unless seeded.

The effect is far below any physical disturbance, and it only shows in flows
that are physically unstable. That is why `rundown_mild_slope` gives every solver
the same 1% inflow seed (PR #407). A fixed dt that divides yieldstep exactly
gives a DE_ader2 steady state to 1e-13.

Prototypes tried (scratch only, not adopted):

- **Predictor with the current step's dt** (CFL pass on Q^n first). This removes
  the prev_dt mismatch but not the dt-dependence: noise falls about 3x with
  output clipping (1.6e-3 to 4.8e-4 h0) and is unchanged without it (9e-4 to
  8.4e-4).
- **Held dt.** Keep dt fixed until the CFL limit forces it down, reset to
  MARGIN x the CFL dt, and meet output times with equal steps r/ceil(r/held).
  Between held changes the inlet relaxes exponentially (3e-4 to 1e-10 h0 in
  25 s). With MARGIN 0.7-0.85 and no raising, the unseeded roll-wave case goes
  steady (1e-13). It is fragile, though: MARGIN 0.9, or 0.75 with a rule to raise
  dt, lets the start-up kicks establish the waves, and the loop then persists.
  It costs about 1/MARGIN more steps, which gives back much of ADER-2's
  one-flux-call advantage.

---

## API

### `numpy` imported as `num` (not `np`)

This is a project-wide convention — do not change it to `np` in existing files.

### `anuga/__init__.py` is the single public API surface

All public names must be both imported and listed in `__all__` in `anuga/__init__.py`.
The file is ~1000 lines; search carefully before adding to avoid duplicates.

### camelCase methods in `pmesh/mesh.py` are deprecated

As of 2026-03-24, camelCase public methods have snake_case equivalents.
The camelCase versions emit `DeprecationWarning`. Prefer snake_case in new code.

### `get_CFL` / `set_CFL` are deprecated in `generic_domain.py`

Use `get_cfl()` / `set_cfl()` instead.

---

## Memory and Performance

### `psutil` is optional

`anuga/utilities/system_tools.py::memory_stats()` tries `psutil` first and falls
back to parsing `/proc/self/status` via `_VmB('VmRSS:')`. If neither works it
returns `'mem=?'`. The `psutil` package is not in the conda environment files
by default.

### Kinematic viscosity operator is slow

`test_kinematic_viscosity_operator.py` runs 4 tests that take 2–5 seconds each.
These are marked `@pytest.mark.slow` at module level.

---

## Meshes

### A couple of sliver triangles can pin the global timestep (towradgi, 2026-08-11)

`run_small_towradgi.py` at the default `scale=1` (256688 triangles) reports a
`delta t` that looks suspiciously frozen — `delta t in [0.10136885, 0.10136885]`,
identical to 8 decimal places across 1184 steps, drifting only in the 6th decimal
between yieldsteps. It is **not** a fixed timestep, a stale reduction, or a GPU
artifact (the same numbers come out on cc75/cc89/cc120 and on the CPU build).

Reconstructing the per-cell CFL limit from the SWW centroid values identifies the
culprit: **cell 72200, the smallest triangle in the mesh** (area 1.457 m²,
inradius 0.425 m vs a 1.99 m median), sitting in a **static pond** — h = 0.923 m,
|v| = 0.000 m/s, unchanged from t=0 to t=600. It binds `min(r/(|v|+√(gh)))` at
every step, so the global minimum barely moves while the flood evolves elsewhere.
The tightest dozen cells are all slivers clustered in a ~50 m patch around easting
307664–307713, northing 6193697–6193735, where `Model/Creeks/creeks.csv` and
`Model/Bdy/CreekBanks.csv` nearly coincide.

The cost is real: relaxing the worst 5 cells would raise `dt` from 0.141 to 0.231
in the reconstruction (~1.6× fewer steps for the whole run).

Confirmed by refinement — at `-sc 0.1` (1636238 triangles) those regions are
regenerated without the pathological slivers and `dt` **rises** to ~0.1176 (the
opposite of what resolution alone does) and becomes genuinely adaptive again
(`delta t in [0.11755643, 0.11783941]`, min ≠ max).

**Lesson:** a `delta t` that is constant to many decimal places is a signal to go
looking for a mesh artifact, not evidence of a well-behaved solver. A handful of
slivers in still water can tax every timestep of a run.

---

## Structures

### `RiverWall` tests require full mesh with breaklines

`anuga/structures/riverwall.py` — tests were deferred because `RiverWall`
requires a domain with a mesh that has breaklines (specific mesh construction).
Simple rectangular domains don't suffice.

### RESOLVED (2026-06-15): "riverwall flux divergence" was really a DE0 boundary bug

**Symptom (now fixed):** a riverwall simulation under `multiprocessor_mode=2`
diverged from legacy — on `run_parallel_riverwall.py` (sequential), stage drifted
from 0 at t=0 to ~0.095 m by t≈100 s.

**Misdiagnosis → real cause.** It was *not* the riverwall flux. The riverwall
kernel (`core_compute_fluxes_central` elevation override + Villemonte weir) is
correct: with a GPU-supported boundary (e.g. `Dirichlet`), mode-1 vs mode-2
riverwall results are **bit-identical (0.0)**. The actual bug was the **boundary**:
`run_parallel_riverwall.py` uses `Transmissive_momentum_set_stage_boundary`
(*not* in `GPU_BOUNDARY_TYPES`), and it was **euler-specific** —
`evolve_one_euler_step()` dispatched straight to `_evolve_one_euler_step_c`, which
handles only GPU boundary types in C and **skips `update_boundary()` entirely**,
so the Transmissive boundary was silently never evaluated (stale edge values →
drift). rk2/rk3/ader2 already fell back to a Python-orchestrated `_gpu` loop for
non-GPU boundaries; **euler had no such fallback**. Confirmed by DE1/DE2/DE_ader2
matching (0.0) while DE0 diverged.

**Fix:** added `_evolve_one_euler_step_gpu()` (host evaluation of non-GPU
boundaries via `evaluate_segment` + `sync_boundary_values`, mirroring rk2/rk3),
and `_evolve_one_euler_step_c()` now delegates to it when
`not self._gpu_all_on_gpu`. DE0 + `Transmissive_momentum_set_stage` + riverwall is
now bit-identical to legacy; `run_parallel_riverwall.py` is **un-pinned** (passes
under `ANUGA_DEFAULT_COMPUTE_MODE=unified` again). Regression test:
`test_DE_gpu_omp.py::Test_GPU_NonGPUBoundaryFallback` (DE0/DE1/DE2/DE_ader2 with a
Transmissive_momentum_set_stage boundary).

**General lesson:** in mode 2, a boundary type not in `GPU_BOUNDARY_TYPES` is only
correct if the active step path falls back to host evaluation. All four DE
algorithms now do. If you add a new evolve path, replicate the
`if not self._gpu_all_on_gpu: return self._evolve_one_*_step_gpu(...)` fallback.

### RESOLVED (2026-07-26): mode-2 shared one global `Time_boundary` value across all time-boundary edges

**Symptom (now fixed):** `validation_tests/analytical_exact/avalanche_wet`
diverged catastrophically under `ANUGA_DEFAULT_COMPUTE_MODE=unified` (xvelocity
L¹ error 0.93 vs legacy 0.006; momentum ran to ~230 vs ~62). `avalanche_dry`,
with the *same* physics, passed — the tell was that dry uses a **single**
`Time_boundary` (the other end is Transmissive) while wet uses **two** (left and
right), and on the sloped bed their absolute stages differ by ~10 m.

**Root cause.** The GPU time boundary stored a **single global** `(stage, xmom,
ymom)` (`struct time_boundary` in `gpu/gpu_domain.h`) applied to *every*
time-boundary edge. `init_time_boundary` (`sw_domain_gpu_ext.pyx`) lumps the
edges of all `Time_boundary` tags into one list, and the evolve loop did
`for B in self._gpu_time_boundaries: set_time_boundary_values(gpu_dom, q0, q1,
q2)` — each call **overwrote** the global, so the last boundary won and
`gpu_evaluate_time_boundary` wrote that one value to all edges. With two
differing boundaries the other one was corrupted. Ablation confirmed it:
`Reflective + slope` was bit-identical (1e-13), only `Time_boundary + slope`
diverged.

**Fix.** Per-edge value arrays, mirroring `file_boundary`: `time_boundary` now
holds `stage_values/xmom_values/ymom_values[num_edges]` (mapped to device, pushed
each step via `omp target update`). New helper
`Domain._push_gpu_time_boundary_values()` builds the per-edge array by evaluating
**each** `Time_boundary` over its own edges in `boundary_map` order (matching
`init_time_boundary`), replacing the 10 clobbering loops. After the fix,
avalanche_wet mode-1 vs mode-2 is bit-identical. Files: `gpu/gpu_domain.h`,
`gpu/gpu_boundaries.c`, `gpu/gpu_domain_core.c`, `sw_domain_gpu_ext.pyx`,
`shallow_water_domain.py`. **Requires a C/Cython rebuild.** Regression:
`test_DE_gpu_omp.py::Test_GPU_TimeBoundary` (two differing Time_boundaries on a
slope, mode 1 == mode 2). **General lesson:** any GPU boundary/operator that
holds a per-edge quantity must store it per-edge, not as one scalar shared across
tags — the existing single-substep `Test_GPU_TimeBoundarySubstep` used a single
boundary and so missed this.

### RESOLVED (2026-07-26): mode-2 GPU fractional operators clobbered by the CPU-sync bracket

**Symptom (now fixed):** `validation_tests/behaviour_only/bridge_hecras2` drained
to the bed under `unified` (peak_max_stage −0.0067 vs baseline 1.19; HEC-RAS
correlation −0.39 vs 0.997), while `bridge_hecras` passed. Ablation localized it:
removing the bridge made the modes match, and the largest stage divergence was at
the **inflow** (y≈11), not the bridge (y≈480–520) or the outflow.

**Root cause.** GPU-accelerated fractional operators (`Inlet_operator` /
`Parallel_Inlet_operator`, `Rate_operator`) apply their update straight to the
**device** arrays in mode 2. When a **CPU-only** fractional operator is also
present (here the `Internal_boundary_operator` bridge),
`apply_fractional_steps()` brackets the operator loop with
`sync_from_device()` … `sync_to_device()`. The trailing host→device sync then
**overwrote the GPU operator's device write** with host data that never received
it — silently dropping the inflow (~600 m³/step here). On a CPU-only build
(host == device) the sync is a no-op so this stays hidden; it bites only on a
real GPU-offload build (device ≠ host).

**Fix.** In `__call__`, skip the device fast-path while
`domain._gpu_host_writes_suppressed` is set (that flag marks the sync-bracketed
region), falling through to the host path so the batch `sync_to_device()` carries
the change. `bridge_hecras2` uses the *parallel* factory, so the guard is needed
in `Parallel_Inlet_operator.__call__` (`anuga/parallel/parallel_inlet_operator.py`)
as well as the base `Inlet_operator` (`anuga/structures/inlet_operator.py`) and
`Rate_operator` (`anuga/operators/rate_operators.py`). Python-only, no rebuild.
Regression: `test_DE_gpu_omp.py::Test_GPU_InletWithCpuOnlyOperator` (inlet +
a no-op CPU-only operator, mode 1 == mode 2, inflow retained). **General
lesson:** a GPU-path operator that writes conserved quantities on-device must
route through the host path whenever `_gpu_host_writes_suppressed` is set, or the
batch host→device sync will discard its work. `collect_max` is exempt — it writes
only the separate `max_*` arrays, not synced conserved quantities.

---

## SWW GUI / animate.py

### `replace_all=True` in Edit tool can change more than intended

When reverting a colormap from `terrain` → `Greys_r` with `replace_all=True`, the `_elev_frame` and `save_elev_frame` default arguments (which must stay `terrain`) were also reverted — requiring a second manual fix. Always check every occurrence of the target string in the file before using `replace_all`.

### Worker must accept all params even when a save method doesn't use them

`worker_frame` in `_animate_worker.py` calls `save_fn(frame=..., show_elev=..., elev_levels=..., show_mesh=...)` for every quantity. If a `save_*` method (e.g. `save_elev_frame`) doesn't declare those params, it raises `TypeError`. All `save_*` methods must accept `show_elev`, `elev_levels`, `show_mesh` even if they ignore the values.

### Double overlays when baked + canvas overlay both active

If Show Elev or Show Mesh is ticked during generation (baked into PNGs) and the canvas overlay is also active, contours/mesh appear twice. The canvas overlay methods check `self._last_gen_show_elev` / `self._last_gen_show_mesh` and return early when already baked. This guard must be maintained if either system is extended.

### Live mesh viewer redraw requires `ax.cla()` + full re-draw

When toggling the Basemap checkbox in `_show_mesh`, a simple `ax.set_visible()` or artist removal is not sufficient — the basemap tiles are added by `contextily` as Axes-level patches. The only reliable approach is `ax.cla()` (clear axis), re-draw the triplot, conditionally call `_add_basemap`, call `mesh_fig.tight_layout()`, then `mesh_canvas.draw()`.

---

## Scenario Module

### `anuga/scenario/` depends on `spatialInputUtil`

The scenario module (`prepare_data.py`, `setup_boundary_conditions.py`, etc.)
imports `spatialInputUtil`, a compiled C extension not included in the main repo.
Meaningful unit tests require this extension plus real shapefile/Excel test data.
Tests for this module are deferred.

---

## Hydrata Current-State Assessment (2026-02-28)

These are known issues identified in the Hydrata fork analysis that also apply to anuga-community.

### `pyproject.toml` declares only `numpy` as a dependency

Despite the codebase importing scipy, netCDF4, matplotlib, meshpy, dill, pymetis, pyproj,
and affine, `pyproject.toml` only lists `numpy>=2.0.0`. This means `pip install anuga`
on a clean venv will produce a package that fails at runtime.

**Fix:** Add the missing dependencies to `[project].dependencies`.

### Phantom dependencies: `cartopy` and `openpyxl`

These appear in code paths but are never actually imported at runtime. Their presence in
any install documentation is misleading.

### GDAL remnants on `remove-gdal` branch

GDAL was partially removed but remnants remain. The `remove-gdal` branch has the work
in progress. Merge not yet complete in anuga-community.

### `setup.py` still present alongside meson-python

Both `setup.py` and `pyproject.toml` (meson-python) exist. The `setup.py` is a
legacy artifact and should be removed once meson-only builds are confirmed in CI.

### Test isolation problems

- **47 `set_datadir('.')` calls** — many tests write files relative to CWD rather than
  a temp directory. Running tests from a non-repo directory can fail or pollute the tree.
- **198 `tempfile.mktemp()` uses** — `mktemp()` is a security risk (TOCTOU) and deprecated.
  Should be replaced with `tmp_path` fixture (pytest) or `tempfile.mkdtemp()`.
- **7+ tests write `domain.sww` to CWD** — parallel test runs step on each other.

### Code duplication (~7,700 redundant lines)

- Three quantity kernels share ~90% code: `quantity_ext.pyx`, `quantity_ext_openmp.pyx`, `quantity_ext2.pyx`
- Five parallel operator wrappers are near-identical to their `structures/` counterparts
- `Culvert_operator` and `Culvert_operator_Parallel` have near-identical logic
- `system_tools.py` is 750 lines with overlap against `numerical_tools.py`

### No linting or type annotations

Zero pre-commit hooks, no ruff/flake8 config, 4,189 functions with no type annotations.
Current approach is manual `pyflakes` / `autopep8` before commits.

### GPU build forced to CPU (`set_gpu_offload(False)` / `-ngo`) is slow — nvc limitation

A `gpu_offload=true` (nvc `-mp=gpu,multicore`) build forced onto the host runs the
`omp target teams distribute` regions through nvc's host fallback, which **does not
scale with threads** (it gets *slower* with more threads). Microbenchmark (40M-element
memory-bound loop, 60 iters, RTX 5070 box, HPC SDK 26.3):

| config | 1t | 8t | 16t |
|--------|----|----|-----|
| nvc `-mp=gpu,multicore` + `OMP_TARGET_OFFLOAD=disabled` | 0.91s | 4.48s | 3.00s (pathological) |
| nvc `-mp=multicore` (multicore-only build) | 0.92s | 0.73s | 0.64s |
| gcc `-fopenmp` (`#pragma omp parallel for`) | 0.92s | 0.76s | 0.61s |
| nvc GPU offload | — | — | 0.12s |

Neither `OMP_TARGET_OFFLOAD=disabled` nor `CUDA_VISIBLE_DEVICES=` engages the good
multicore variant of the dual build — the host always gets the GPU variant's serial-ish
fallback. towradgi small (256k tri, -ft 200 -ys 50) confirms it: GPU 6.35s, `-ngo`
60–100s (1.7× scaling 1→16 threads), vs a gcc `gpu_offload=false` build at ~17s.

**Implication:** a GPU build is not a substitute for a CPU build. `set_gpu_offload(False)` /
`-ngo` is for **correctness A/B** (verify GPU and CPU give identical results — they are
bit-identical) only, NOT timing. For CPU-multicore performance, build with
`-Dgpu_offload=false` (gcc → host-optimised `omp parallel for` via the `CPU_ONLY_MODE`
macros). `set_gpu_offload(False)` warns about this on a GPU build.

**This is a confirmed, documented NVHPC limitation, not an ANUGA bug** (investigated
2026-06-13). The NVHPC Reference Guide defines `-mp=gpu` as "compiled for GPU execution
*as well as host fallback to the CPU*" — and that host fallback runs **single-threaded**.
A peer-reviewed compiler comparison (IPDPSW 2023, Iowa State) measured NVHPC host fallback
at "OMP 1" (1 CPU thread) and found "the GPU code version on CPU in host fallback mode
performs worse than the CPU version with 1 thread". NVHPC is also documented to handle
nested/inner parallel regions poorly (NVIDIA recommends the `loop` directive over
`teams distribute parallel for` for this reason). The fast multicore variant only exists
when `-mp=multicore` is the *sole* mode; `-mp=gpu,multicore` does not let the runtime pick
it on the host (verified across `OMP_TARGET_OFFLOAD=disabled`, `CUDA_VISIBLE_DEVICES=`,
argument order, and `ACC_DEVICE_TYPE` — the last hangs). So a single nvc binary cannot be
fast on both GPU and CPU; the two-build split (gcc CPU / nvc GPU) is required.

References:
- NVHPC Compilers Reference Guide 26.3 — https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-ref-guide/index.html
- "OpenMP Offload Features and Strategies for High Performance across Architectures and
  Compilers", IPDPSW 2023 — https://swapp.cs.iastate.edu/files/inline-files/OpenMP_Offload_Features_and_Strategies_for_High_Performance_across_Architectures_and_Compilers-ipdpsw-may-2023.pdf
- OMP_TARGET_OFFLOAD, OpenMP 5.0 spec — https://www.openmp.org/spec-html/5.0/openmpse65.html

## `struct domain` fields must be explicitly initialised in sw_domain_openmp_ext

(2026-08-24) `Domain_C_struct.__cinit__` allocates the C struct with
`PyMem_Malloc` — **uninitialized memory** — and fills every field explicitly in
`get_python_domain_parameters` / `get_python_domain_pointers`. Adding a field
to `struct domain` (sw_domain.h) without also setting it there leaves garbage
in the legacy path. Symptom when it happened (`reconstruct_edge_bed`): the
flux kernel read a garbage flag, well-balance broke, and
`test_ader2::test_still_water_flat_bed` spun forever in its evolve loop at
100% CPU (timestep collapse — not a hang in any one kernel). The GPU-path
`GPUDomain` is a cdef-class member and IS zero-initialised, which makes the
asymmetry easy to miss. Both pyx files also redeclare the struct member list,
so a new field must be added there before Cython can assign it.

## ninja does not reliably re-track `sw_domain.h` — clean-rebuild after struct changes

(2026-08-24) After editing `struct domain` (sw_domain.h), incremental `ninja`
in `build/cp313` rebuilt only a subset of objects (and once reported "no work
to do" straight after a header revert). The resulting mixed-layout objects
produce baffling failures: OpenMP "partially present on the device" fatals at
map time, silent `exit(1)` during mode-2 init, and test failures that migrate
as you touch unrelated files. If a change touches any shared struct header,
run `ninja -t clean && ninja` (74 targets, ~1 min) before trusting any test
result — an afternoon was lost bisecting phantom bugs that were stale objects.

## OMP_TEAMS_THREAD_LIMIT silently corrupts nvc `target teams loop` kernels

(2026-08-25, GPU build, nvhpc 25.9) Setting `OMP_TEAMS_THREAD_LIMIT` (32/64/96)
on the unified GPU kernels produced spectacular "speedups" (up to 1.6x) with
SILENTLY WRONG results — the lake-at-rest well-balance check exploded from
1e-14 to O(10) momentum, and at 48 the kernels appear not to run at all
(exact 0.0 outputs). nvc's `omp target teams loop` codegen evidently bakes in
blocking assumptions the runtime override violates. Never benchmark this knob
without a physics gate; `OMP_NUM_TEAMS` / `OMP_THREAD_LIMIT` were separately
measured as safe-but-useless (nvc defaults are already right).

## Bedload is first order in the bed

The bedload edge flux is centred with Rusanov dissipation scaled by the
bed-wave speed (`core_apply_bedload`, pass 2). That makes it monotone and
exactly conservative, but first order: a migrating bed form is smoothed at
O(dx). The `sediment_bed_hump` validation case measures 9.9% L1 error of the
hump volume at 10 m cells and 5.5% at 5 m. The bare centred flux it replaced
was unstable on that case (scour hole at the upstream toe, overshoot at the
crest, noise fed back into the flow). The upgrade, if sharper bed forms are
ever needed, is a limited (MUSCL) reconstruction of z at the edge feeding the
same Rusanov flux.

The bed update runs in a pass of its own (pass 3): the Rusanov term reads the
neighbour's bed elevation, so writing the bed inside the divergence loop made
the two sides of an edge see different states and broke conservation. Any
future term that reads a neighbour's bed must keep that separation.
