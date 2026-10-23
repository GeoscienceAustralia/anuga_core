# Plan — make `multiprocessor_mode=2` + `gpu_offload=false` the standard distribution

Created: 2026-06-12. Owner: Stephen Roberts. Status: **in progress** (step 1 in review as PR #144).

## Goal

Ship ANUGA so that a standard `pip install` builds the unified `sw_domain_gpu_ext`
kernels for **CPU multicore** (`gpu_offload=false`) and runs the solver **and operators**
through them by default (`multiprocessor_mode=2`). This routes the rainfall / culvert /
inlet operators through C OpenMP kernels instead of serial Python, closing nearly the
entire OpenMP→MPI gap for single-node users with no MPI setup.

## Why (evidence)

Towradgi small (~256k tri, OMP_NUM_THREADS=16), gpu_offload=false build:

| Config | Time | Note |
|--------|------|------|
| mode=1 (Python ops) + metis_rcm | 18.64 s | current default path |
| mode=2 (C ops) + no reorder | 17.03 s | C operators alone |
| **mode=2 (C ops) + metis_rcm** | **12.27 s** | both — within ~11% of MPI-16 (11.08 s) |

Profiling showed operators cost ~40% of a mode-1 OpenMP run as serial Python that does
not scale with threads — the structural reason MPI beats OpenMP. The C kernels already
exist in `anuga/shallow_water/gpu/` (`gpu_rate_operator.c`, `gpu_culvert_operator.c`,
`gpu_inlet_operator.c`) and compile to `#pragma omp parallel for` under
`-DCPU_ONLY_MODE` (see `gpu_omp_macros.h`). Validated numerically: `test_DE_gpu_omp.py`
56/56 mode-1-vs-mode-2 equivalence tests pass on a `gpu_offload=false` build.

## Naming note

`multiprocessor_mode=2` is labelled "GPU" but really means "use the unified gpu_ext
kernels" — which target a GPU when `gpu_offload=true` and CPU multicore when
`gpu_offload=false`. The terminology should change (see step 4).

---

## Steps

### Step 1 — Deferred interface build  ✅ done, in review (PR #144)
Branch `feat/defer-gpu-interface-build`. Decouples *choosing* mode 2 from *building*
the device interface: mode recorded immediately, interface built eagerly if boundaries
are ready, else lazily at first `evolve()`. Removes the "boundaries before mode" ordering
constraint so mode 2 can be selected at `__init__`. Awaiting Jorge's review.

### Step 2 — Audit fall-back for kernel-less operators  ✅ DONE 2026-06-12 (audit + Inlet fix)
The equivalence tests only cover operators that HAVE C kernels (rate, inlet, culvert,
weir). Before any default switch, confirm every other operator behaves correctly in
mode 2 — **graceful fall-through to Python, never a silent no-op**.

**How dispatch works (verified).** `Domain.apply_fractional_steps()`
(`shallow_water_domain.py:4391`) is overridden in mode 2. It runs **every** operator's
`__call__` in the loop — nothing is skipped or no-op'd by mode. Operators with no
mode-2 branch simply execute their normal Python body. `_has_cpu_only_fractional_operators()`
(`:4316`) classifies each op; any operator not in the known-GPU-safe set
(Rate / boundary_flux / Inlet / Collect_max / Boyd-via-manager) is flagged "CPU-only"
and the loop is wrapped in `sync_from_device()` … `sync_to_device()`.

**CPU-multicore (`gpu_offload=false`) is safe.** In that build `set_gpu_interface()`
still creates a `GPU_OMP_interface` (`:5109`), so `gpu_interface is not None` and the
sync path *does* run — but the C sync (`gpu_domain_sync_to/from_device`,
`gpu_domain_core.c:1098/1111`) is `#pragma omp target update` on shared host memory,
i.e. a no-op under `CPU_ONLY_MODE`, guarded by `if (!GD->gpu_initialized) return;`.
Host == device, so every Python operator reads/writes the same arrays the gpu_ext
kernels use. **No silent no-ops, no corruption for the CPU default.**

**Audit table** (operator → status in mode 2):

| Operator | Mode-2 branch? | Status |
|----------|---------------|--------|
| `Rate_operator` (scalar/t/quantity/centroid_array) | yes | **C kernel**; graceful `_init_gpu` returns → CPU fallback |
| `Rate_operator` (`rate_spatial` / `rate_xarray`) | guarded out | **Falls back to Python** (good pattern) |
| `Inlet_operator` | yes | **C kernel**; ⚠ `_init_gpu` *re-raises* on failure (no graceful fallback) — see fix below |
| `Boyd_box` / `Boyd_pipe` / `Weir_orifice_trapezoid` | via `GPUCulvertManager` | **C kernel** (`is_boyd_operator` covers all three) |
| `Internal_boundary_operator` | no | Python `__call__` runs; classified CPU-only → sync wrap ✓ |
| `Collect_max_quantities_operator` | yes | **C kernel**; falls through to NumPy when `_gpu_initialized` False ✓ |
| `Collect_max_stage_operator` | no | Python ✓ |
| `Kinematic_viscosity_operator` | no | Python / C-CG (+MPI) on host arrays ✓ CPU; ⚠ GPU edge sync (real-GPU only) |
| `Erosion_operator` + subclasses (`Bed_shear`, `Circular`, `Polygonal`, `Flat_slice`, `Flat_fill_slice`) | no | Python ✓ CPU; ⚠ modifies **elevation**, not covered by centroid sync (real-GPU only) |
| `Sanddune_erosion_operator` | no | Python ✓ CPU; ⚠ elevation sync (real-GPU only) |
| `Mannings_operator`, `Wind_stress_operator` | no | Python ✓ |
| `set_elevation` / `set_friction` / `set_quantity` / `set_stage` / `set_w_uh_vh` | no | Python ✓ CPU; `set_elevation` ⚠ elevation/edge sync (real-GPU only) |
| `boundary_flux_integral_operator` | classified skip | GPU-safe (reads `boundary_flux_sum` only) ✓ |
| `elliptic_operator` | no | Python ✓ |

**Conclusion for the CPU-multicore default (this plan's target): nothing is BROKEN.**
Every operator either has a CPU-OpenMP C kernel or runs its Python body correctly on the
shared host arrays.

**Fix applied (robustness):** `Inlet_operator._init_gpu` (`inlet_operator.py`) used to
**re-raise** on any failure, so an inlet whose GPU init failed — missing/partial
`gpu_interface`, or `MAX_INLET_OPERATORS=32` slot limit exceeded — **crashed** instead of
falling back to Python. Rewritten to mirror `Rate_operator._init_gpu`: graceful
precondition guard (`gpu_interface`/`gpu_dom` missing → silent `return`), and
`warnings.warn(...) + return` on init exception or `op_id < 0`, leaving
`_gpu_initialized=False` so `__call__`'s existing `if self._gpu_initialized:` guard falls
through to the Python path. This also makes `_has_cpu_only_fractional_operators` correctly
classify a failed-init inlet as CPU-only → it gets the sync wrap the Python path needs.
Verified: mode-2 inlet with `gpu_interface=None` returns silently and runs the Python
`__call__` without crashing; `test_inlet_operator.py` 14/14 pass.

**Rate fix applied (same class of issue):** `Rate_operator._init_gpu`
(`rate_operators.py`) previously hard-`raise`d on the `MAX_RATE_OPERATORS=64` slot limit
(and an unguarded kernel import/init could propagate too). Its precondition guards were
already graceful; the init call + slot-limit check are now wrapped to `warnings.warn(...)
+ return` (leaving `_gpu_initialized=False`, `_gpu_op_id=None`), matching Inlet. Verified:
`test_rate_operators.py` 37/37 pass; mode-2 rate with `gpu_interface=None` falls through
to the Python path without crashing. Both kernel-backed fractional operators now share one
graceful-fallback contract.

**Real-GPU follow-up (out of scope here, gate for step 5 under `gpu_offload=true`):**
the host/device sync only covers centroid values (stage/xmom/ymom/height). Operators
that modify **elevation** (erosion, sanddune, set_elevation) or read **edge** values
(kinematic viscosity) would not have those changes reflected on the device. Harmless in
CPU-multicore (host == device); must be addressed before erosion/KV operators are used
in real-GPU runs. Track separately from the CPU-default switch.

### Step 3 — Build/packaging: default `gpu_offload=false`  🔄 IN PROGRESS 2026-06-12
- **`gpu_offload` already defaults to `false`** (`meson_options.txt:8`) — a stock
  `pip install` already builds the CPU-multicore extension. GPU users opt in with
  `-Dgpu_offload=true -Dgpu_arch=...`. ✅ (no change needed)
- ⬜ Confirm `sw_domain_gpu_ext` builds/imports on a minimal **no-MPI** conda env (meson
  builds it with MPI stubs — verify on a clean env, since a mode-2 default hard-fails on
  import if the extension is missing). **Still to do — needs a clean no-MPI env / CI run.**
- ⬜ Verify CI builds and conda recipes pass with the CPU-multicore extension.

### Step 4 — Compute model: two orthogonal knobs  ✅ DONE 2026-06-13 (API shipped; CLI/docs follow)
**Design correction (user decision).** The first cut exposed `legacy`/`cpu`/`gpu` as three
peer *per-domain* modes — but `cpu` and `gpu` are the **same** mode-2 path differing only
by GPU offload, which is a **process-global** OpenMP setting, not a per-domain property.
That scope-mismatch is what produced the cross-domain offload bug. The model is now split
into two orthogonal knobs:

**Knob 1 — per-domain compute mode** (`Domain`): **`legacy` / `unified`** (mode-2 named
`unified` — user choice).
- `set_compute_mode(mode, verbose=False)` — `'legacy'`→mode 1 (`sw_domain_openmp_ext` +
  serial-Python ops); `'unified'`→mode 2 (unified `gpu_ext` C kernels, solver + operators).
  Per-domain: different domains may use different modes. No offload side-effects.
- `get_compute_mode()` → `'legacy'|'unified'`. `set_multiprocessor_mode(1|2)` → `'legacy'`
  / `'unified'`.
- Parallel guard retained: under MPI, `'unified'` on a no-MPI gpu_ext build falls back to
  `'legacy'` (rank-0 warning) — see below.

**Knob 2 — process-global GPU offload** (module functions, exported from `anuga`):
- `set_gpu_offload(enable=True, verbose=True)` — process-wide; decides whether `'unified'`
  domains offload to a GPU or run CPU-multicore. Enabling on a CPU-only build (or no
  device) warns and stays off; never hard-fails. Must be called **before the first
  `evolve()`** (OpenMP reads `OMP_TARGET_OFFLOAD` at init). Returns the resolved state.
- `gpu_offload_enabled()` — resolved process state (explicit request ∧ build support).
- `gpu_offload_supported()` — wraps `sw_domain_gpu_ext.gpu_available()` (False in
  `CPU_ONLY_MODE`, no device, or launch-time `OMP_TARGET_OFFLOAD=disabled`).
- Module-global `_GPU_OFFLOAD_REQUESTED` (None=build default) replaces the old per-instance
  `_auto_disabled_offload` flag — fixing the bug where domain B asking to offload couldn't
  clear a disable set by domain A.

So `cpu` = `unified` + offload off; `gpu` = `unified` + offload on — *compositions of the
two knobs*, not primitives. `set_gpu_interface()` now sets `gpu_offload_active =
gpu_offload_enabled()`; banner reads "'unified' CPU multicore … no offload" when off, and
the misleading "GPU<->CPU sync every RK2 step" warning stays suppressed on CPU (host==device).

**Why offload is process-global, not per-domain:** OpenMP target offload is a process-level
runtime ICV — one process cannot run the same unified kernels on GPU for domain A and CPU
for domain B. Hence knob 2 is a module function.

**Offload mechanism (what actually works — learned the hard way).** Two parts are needed;
neither alone is enough:
1. **`OMP_TARGET_OFFLOAD=disabled` is the real lever.** It is what keeps the `omp target`
   regions (solver AND operators) on the host. Confirmed on towradgi (256k tri, -ft 200):
   GPU = 6.35 s; `set_gpu_offload(False)` *without* the env var still ran on the GPU
   (6.35 s — the default-device flag alone does NOT stop nvc offloading); with
   `OMP_TARGET_OFFLOAD=disabled` it ran on CPU (≈60 s). So `set_gpu_offload()` sets/clears
   this env var. The OpenMP runtime reads it at the first target region, so call
   `set_gpu_offload()` **before the first evolve()/domain build**.
2. **Device-id consistency** still matters so data mapping and operator device selection
   agree on the host. The kernels carry no `device()` clause and map arrays via
   `omp target enter data map` to the **default device**, chosen in `gpu_domain_init`.
   - `gpu_domain_init` routes to the host (`device_id=-1`) when offload is off (the
     `g_gpu_offload_enabled` flag, set via `set_offload_enabled`).
   - **Critical bug fixed:** the inlet/culvert operators re-set the default device with
     `omp_set_default_device(GD->device_id)`; with `device_id=-1` that is an invalid device
     number that silently re-routed work back to the GPU. Added `gpu_compute_device(GD)`
     (`device_id>=0 ? device_id : omp_get_initial_device()`) and used it at all 7 sites.
   - New C exposure: `set_default_device`/`get_default_device`/`get_initial_device` and
     `set_offload_enabled`/`get_offload_enabled`.
- Contract: call `set_gpu_offload()` **before building the first 'unified' domain**;
  toggling an already-mapped domain is unsupported. Runtime re-enable after a disable is
  best-effort (the OpenMP runtime may not re-read `OMP_TARGET_OFFLOAD`).

**`set_omp_num_threads` is now process-level too.** OpenMP thread count is a process ICV,
so it moved to a module function `anuga.set_omp_num_threads()` (one `omp_set_num_threads`
covers both the legacy and unified kernels). `Domain.set_omp_num_threads` stays as a thin
backward-compatible wrapper that delegates and records `self.omp_num_threads`.

**Validated on real GPU hardware (2026-06-13).** GPU build via
`CC=nvc … -Dgpu_offload=true -Dgpu_arch=cc120` (RTX 5070 Laptop, 1 device, HPC SDK 26.3).
`gpu_available()`=True, `gpu_has_mpi()`=True. Verified end-to-end:
- `test_DE_gpu_omp.py` 56/56 (GPU path), `test_compute_mode.py` 12 pass + 1 skip
  (CPU-only-warning test correctly skipped), 172 broader domain/operator tests pass.
- Three-way run of the same scenario — `legacy` / `unified`+GPU / `unified`+offload-off —
  agree to machine precision: GPU vs legacy = 1.5e-15, CPU vs legacy = 1.5e-15, **GPU vs
  CPU = 0.0 (bit-identical)**.
- The robust toggle works on the GPU build: `set_gpu_offload(False)` routed the unified
  domain to the host (banner "CPU multicore", `device_id=-1`, `gpu_offload_active=False`);
  `set_gpu_offload(True)` back to the GPU (`device_id=0`) — re-enableable, no env mutation.
- Multi-rank parallel GPU tests skip cleanly on a 1-GPU box (need 1 GPU/rank).

**Parallel (MPI) safety guard.** `'unified'` exchanges ghosts at the C level
(`exchange_ghosts`), a *silent no-op* unless `sw_domain_gpu_ext` was compiled with MPI
(`gpu_has_mpi()` — build-time `HAVE_MPI4PY`). `set_compute_mode('unified')` with
`numprocs > 1` and no MPI build warns (rank 0) and falls back to `'legacy'` (correct
Python/mpi4py exchange). Supported matrix:

| build | serial | MPI (numprocs>1) |
|-------|--------|------------------|
| CPU-only, **no** MPI | legacy, unified | legacy only (unified→legacy with warning) |
| CPU-only, **with** MPI | legacy, unified | legacy, **unified** (hybrid MPI+OpenMP) |
| GPU build, with MPI | legacy, unified (+offload) | legacy, unified (+offload) |

Introspection: `Domain.compute_capabilities()` → `{gpu_offload, num_gpu_devices, mpi,
modes}` (`modes` = `legacy` always, `unified` when gpu_ext imports).

Tests: `anuga/shallow_water/tests/test_compute_mode.py` (13 tests — legacy/unified, int
API, cpu/gpu rejected as modes, capabilities, parallel-no-MPI fallback, parallel+MPI keeps
unified, and a `Test_gpu_offload` class for the process toggle), registered in
`tests/meson.build`, robust to both build types. Regression green: `test_DE_gpu_omp.py`
56/56, `test_shallow_water_domain.py` + `test_sw_domain_openmp.py` (106), parallel mode-2.

⬜ **Follow-up:** wire `--compute-mode {legacy,unified}` and a `--gpu-offload` flag into the
standard arg parser (`-mpm` stays); update help and user docs (rolls into Step 6).

### Step 5 — Flip the default  ⬜ (the actual switch)
- In `Domain.__init__` (shallow_water), default to mode 2 **with an auto-fallback**:
  if `sw_domain_gpu_ext` failed to import OR the build lacks the needed kernels, fall
  back to mode 1 and warn once. Never hard-fail for basic users.
- Gate on import success:
  ```python
  try:
      from anuga.shallow_water import sw_domain_gpu_ext   # noqa: F401
      _default_mode = MULTIPROCESSOR_GPU
  except Exception:
      _default_mode = MULTIPROCESSOR_OPENMP
  ```
- Run the FULL suite (`pytest --pyargs anuga`, not just --run-fast) with the new default
  on a `gpu_offload=false` build. Many unit tests build domains without driving a full
  evolve / setting boundaries — verify the deferred build (step 1) means they don't
  trip the mode-2 setup. Fix any that assume mode 1 internals.
- Re-run validation_tests/ (analytical + experimental) under the new default.

### Step 6 — Docs & comms  ⬜
- README / install docs: explain CPU-multicore is the default, GPU is opt-in.
- Note `OMP_NUM_THREADS` controls parallelism; recommend `-ro metis_rcm` (or wire a
  sensible default reorder — see open question).
- Migration note for users who relied on mode-1-specific behaviour.

---

## Open questions / risks
- **Determinism**: mode 2 reductions (OpenMP `reduction(+:)`) may differ from mode 1 at
  the ULP level → bit-for-bit reproducibility across thread counts is not guaranteed.
  Confirm regression-snapshot tests tolerate this (they passed at 56/56, but the full
  snapshot suite under a forced default should be checked in step 5).
- **Should reorder be automatic?** metis_rcm gives a big chunk of the mode-2 win. Consider
  a default `reorder='metis_rcm'` at domain build (or first evolve) rather than requiring
  `-ro`. Separate decision; don't couple to this plan unless cheap.
- **Windows**: mode-2 CPU build needs an OpenMP-capable compiler (mingw) — same constraint
  as the existing `sw_domain_openmp_ext`. Verify the gpu_ext builds under the Windows CI.
- **mode 1 retirement**: if mode 2/CPU is strictly better and fully covering, `sw_domain_
  openmp_ext` (mode 1) eventually becomes redundant. Out of scope here; revisit later.

## Key references
- PR #144 (step 1): https://github.com/anuga-community/anuga_core/pull/144
- Benchmark + validation detail: `claude/SESSION_GUIDE.md` → "CPU multicore via the unified
  gpu_ext C kernels" section.
- C kernel macros: `anuga/shallow_water/gpu/gpu_omp_macros.h`
- Mode dispatch: `anuga/shallow_water/shallow_water_domain.py` (`set_multiprocessor_mode`,
  `set_gpu_interface`, `_boundaries_ready`, evolve lazy hook); `anuga/config.py`
  (`MULTIPROCESSOR_OPENMP=1`, `MULTIPROCESSOR_GPU=2`).
- Equivalence tests: `anuga/shallow_water/tests/test_DE_gpu_omp.py`.
