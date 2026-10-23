# GPU/OpenMP Offloading — Development Plan (v4.0.0 / SC26)

Written 2026-04-01 after full architecture review of `anuga/shallow_water/gpu/`.

---

## Architecture Summary

The GPU pipeline is:

```
Python (sw_domain_gpu_omp.py)
  → Cython bridge (sw_domain_gpu_ext.pyx, 1 797 lines)
      → C kernels (gpu/*.c, ~12 files)
          → #pragma omp target teams loop  [GPU mode]
          → #pragma omp parallel for simd  [CPU_ONLY_MODE]
```

The macro header `gpu_omp_macros.h` switches the backend at compile time. All
kernels are shared between GPU and CPU-multicore modes — no duplicate code paths.

**Per-element FLOP costs (Gordon Bell counting):**

| Kernel             | FLOPs/elem | Notes |
|--------------------|------------|-------|
| extrapolate         | 220        | Dominant cold-start cost |
| compute_fluxes      | 400        | Dominant hot-loop cost |
| update              | 21         |  |
| protect             | 5          |  |
| Manning friction    | 18         |  |
| SAXPY (RK2 stages) | 9          |  |
| rate_operator       | 8          |  |

Full RK2 is executed in C (`gpu_evolve_one_rk2_step`) — no Python round-trip
per step. MPI Allreduce picks global minimum timestep.

**Operator status:**

| Operator | GPU-accelerated? | Notes |
|----------|-----------------|-------|
| Rate_operator | Yes (kernel) | Up to MAX_RATE_OPERATORS=64 |
| Inlet_operator | Yes (small D2H/H2D) | ~6 KB gather/scatter, up to 32 |
| Boyd_box / Boyd_pipe culverts | Physics CPU, gather/scatter GPU | ~2 KB/timestep |
| Gate operators (Weir_orifice, etc.) | **No — CPU only** | Gap |
| File_boundary | **No** | Gap |
| Flather / Compound boundaries | **No** | Gap |

---

## Phase 1 — Correctness and Test Coverage (weeks 1–4)

These are blockers for any production use or paper submission.

### 1.1 Boundary condition completeness

**Currently supported on GPU:** Reflective, Dirichlet, Transmissive,
Transmissive_n_zero_t, Time_boundary.

**Not yet supported:** File_boundary, Flather_boundary, compound boundaries.

Priority: **File_boundary** — it is the standard open-ocean forced boundary
in all real ANUGA tsunami simulations. Without it, GPU mode cannot run most
real models.

Approach:
- Add a `struct file_boundary` to `gpu_domain.h` with a time-interpolated
  stage/xmom/ymom triplet (same pattern as `time_boundary` but driven by an
  array rather than a scalar — Python computes the interpolated value and
  pushes it to the C struct each timestep, identical to `time_boundary`).
- The boundary kernel is trivial (same as Dirichlet, one-value-per-edge).
- Time interpolation stays in Python (no GPU work needed).

### 1.2 Device memory check before data mapping

`gpu_domain_core.c` maps all arrays to device without checking whether enough
device memory is available. On large meshes (>5 M elements) this silently
corrupts results or crashes.

Approach:
- Add `gpu_check_device_memory()` before the first `omp target enter data` call.
- Use `omp_target_alloc` with a canary size to probe free memory.
- Print a clear error and fall back to CPU_ONLY_MODE (or abort) if insufficient.

### 1.3 Operator slot limits

Hard limits (`MAX_RATE_OPERATORS=64`, `MAX_INLET_OPERATORS=32`,
`MAX_CULVERTS=64`) are static arrays in the C struct. Real models often exceed
these. Current behaviour: silent truncation.

Approach:
- Add an assertion + meaningful error in `gpu_rate_operator.c`,
  `gpu_inlet_operator.c`, and `gpu_culvert_operator.c` when the limit is
  exceeded.
- Medium-term: replace static arrays with heap-allocated dynamic lists.

### 1.4 Test coverage gaps

Current `test_DE_gpu_omp.py` (899 lines) covers:
- Kernel-by-kernel checks (flux, extrapolate, update, Manning, rate, inlet)
- Single-process evolve tests

Missing:
- Multi-rank (2- and 4-process) GPU tests — need to verify halo exchange
  correctness on non-trivial partitions.
- Culvert operator test in GPU mode.
- Test that `multiprocessor_mode=2` gives the same final state as `mode=1`
  after N full timesteps (end-to-end regression, not just one RK step).
- CPU_ONLY_MODE test (compile flag) to catch macro regressions.

Approach:
- Add `@pytest.mark.slow` multi-rank tests using `mpirun -n 2` subprocess
  (same pattern as existing parallel tests in `anuga/parallel/tests/`).
- Add an evolve regression test: run 10 s of simulation in mode 1 and mode 2,
  compare stage/momentum to tolerance 1e-10.

---

## Phase 2 — Performance Validation (weeks 5–10)

Required for SC26 submission. Target: demonstrate strong/weak scaling.

### 2.1 Benchmark suite

Create `examples/gpu_benchmark/` with:
- **Small** (100 K triangles) — fits on one GPU, tests baseline throughput.
- **Medium** (2 M triangles) — typical real model, multi-rank (4–8 GPUs).
- **Large** (20 M triangles) — SC26 size, designed for ~64 GPUs.

Each benchmark should:
1. Run N timesteps.
2. Print Gordon Bell FLOP/s (already instrumented via `gpu_flop.c`).
3. Report wall-clock time, mesh size, MPI rank count.

### 2.2 GPU-aware MPI

Currently `-DGPU_AWARE_MPI` is a compile-time option that enables direct
device-pointer MPI calls (eliminates the D2H/H2D halo copy). This is a
significant win on NVLink / InfiniBand systems.

Actions:
- Verify correctness on NV systems with GPU-aware MPI (CUDA-aware OpenMPI or
  MPICH).
- Add CMake/meson option `gpu_aware_mpi=true/false` with a runtime check that
  the MPI library supports it (detect via `MPIX_CUDA_AWARE_SUPPORT` or
  `MPIX_GPU_QUERY_CUDA_SUPPORT`).
- Document in `docs/source/parallel/` that this requires a CUDA-aware MPI.

### 2.3 NVTX / OMPT profiling hooks

Add profiling region markers around key kernels to enable `nsys` / `ncu`
profiling without recompilation:

```c
// In gpu_kernels.c — wrap major kernels
#ifdef ENABLE_NVTX
nvtxRangePush("compute_fluxes");
#endif
// ... kernel ...
#ifdef ENABLE_NVTX
nvtxRangePop();
#endif
```

These should be gated on a compile-time flag so they add zero overhead in
production.

### 2.4 Weak scaling experiment

The most important SC26 metric is weak scaling: keep elements-per-GPU constant
as GPU count increases. Target: >80% parallel efficiency up to 64 GPUs.

Key potential bottlenecks to measure:
- MPI Allreduce for global timestep (called every RK2 step).
- Halo exchange volume (grows with partition boundary length ≈ √N).
- Culvert MPI communication (per-culvert point-to-point messages).

---

## Phase 3 — Feature Completeness for v4.0.0 (weeks 11–20)

### 3.1 Gate operators on GPU

`Weir_orifice_trapezoid_operator` and `Broad_crested_weir_operator` are
currently CPU-only. Note: `gpu_adjust_edgeflux_with_weir()` **already exists**
in `gpu_device_helpers.h` as a device function implementing Villemonte (1947)
weir discharge — the hard physics is already ported.

What's missing is the operator registration path (analogous to culvert batch
gather/scatter). Approach:
- Add `struct weir_operator_info` to `gpu_domain.h` (same pattern as
  `culvert_indices`).
- Register in Cython during `init_gpu_domain`.
- Apply in a GPU kernel per timestep (the physics call is already in device code).

### 3.2 Riverwall support

Riverwalls define weir-like structures along mesh edges. The physics helper
`gpu_adjust_edgeflux_with_weir` is already in `gpu_device_helpers.h`.
The flux kernel `gpu_kernels.c` needs to check a per-edge riverwall flag and
call the adjustment for flagged edges. This is a medium-complexity addition.

### 3.3 Dynamic operator slot limits

Replace static `MAX_*` arrays with heap allocation:

```c
// gpu_domain.h
struct rate_operators {
    struct rate_operator_info *ops;  // heap-allocated
    int capacity;
    int num_operators;
    int initialized;
};
```

Required for large models with many rainfall patches.

### 3.4 GPU documentation

Add to `docs/source/`:
- `gpu_mode.rst` — how to enable multiprocessor_mode=2, build flags, hardware
  requirements.
- Benchmark results page (once 2.1 is done).
- GPU-specific boundary and operator limitations.

---

## Phase 4 — SC26 Paper Preparation (months 4–6)

### 4.1 Gordon Bell metrics

The FLOP profiling infrastructure (`gpu_flop.c`) already measures:
- Total FLOPs per timestep (extrapolate + fluxes + update + Manning + SAXPY + …)
- Wall-clock time → FLOP/s
- MPI-reduced global totals across all ranks

Need to add:
- Peak theoretical FLOP/s for the target GPU (compile-time constant or queried
  at runtime via `omp_get_device_num()` + device-specific value).
- Roofline model comparison (memory-bound vs compute-bound analysis).
- Time-breakdown per kernel (currently total only — add per-kernel timing).

### 4.2 Validation against physical benchmarks

GPU mode must produce results matching known analytical solutions:
- Thacker paraboloid (rotating flow, exact solution exists).
- Dam break (Ritter solution).
- Tide benchmark (compare to measured tide gauge data).

These already exist in `validation_tests/` — add a GPU-mode run path to the
validation scripts.

### 4.3 Multi-node strong scaling

Demonstrate that a fixed-size problem (e.g. 20 M triangles) scales from 1 to
64 GPUs. Expected runtime reduction: ~50× (accounting for communication
overhead). This is the key SC26 performance claim.

---

## Implementation Order (Recommended)

| Priority | Item | Effort | Impact |
|----------|------|--------|--------|
| P0 | 1.4 End-to-end regression test (mode1 vs mode2) | 1 day | Confidence in correctness |
| P0 | 1.3 Slot limit assertions + errors | half day | Prevent silent failures |
| P1 | 1.1 File_boundary GPU support | 2 days | Enables real models in GPU mode |
| P1 | 1.2 Device memory check | 1 day | Prevents silent crashes on large runs |
| P1 | 1.4 Multi-rank halo exchange test | 2 days | Validates MPI correctness |
| P2 | 2.1 Benchmark suite | 3 days | SC26 results |
| P2 | 2.2 GPU-aware MPI validation | 2 days | Scaling performance |
| P2 | 2.4 Weak scaling experiment | 1 week | Core SC26 result |
| P3 | 3.1 Weir/gate operators on GPU | 3 days | Feature parity with CPU mode |
| P3 | 3.2 Riverwall support | 2 days | Real-world model support |
| P3 | 3.4 GPU documentation | 2 days | v4.0.0 release requirement |
| P4 | 4.1 Gordon Bell metrics | 1 week | Paper submission |
| P4 | 4.2 Validation against physical benchmarks | 1 week | Paper rigour |
| P4 | 4.3 Multi-node strong scaling | ongoing | Core SC26 claim |

---

## Quick Wins Available Right Now

These can be done without hardware — they run in CPU_ONLY_MODE:

1. **End-to-end regression test** (P0 above) — compare mode=1 and mode=2 over
   10 s of simulation. Catches any future kernel regressions immediately.

2. **Slot limit assertions** (P0 above) — 30-minute change, prevents silent
   bugs in large models.

3. **File_boundary skeleton** (P1 above) — the struct + Python push path can be
   written and tested in CPU_ONLY_MODE before GPU hardware is available.

4. **Culvert test in GPU mode** — the infrastructure exists; test doesn't.
