# Session Guide

How to orient a new Claude session for ANUGA development work.

---

## Quick orientation

```bash
git branch          # see all branches
git log --oneline -10   # recent commits
git status          # current state
```

Key files to read first:
- `CLAUDE.md` — build system, test commands, architecture overview
- `claude/PROGRESS.md` — what has been done and what remains
- `claude/DECISIONS.md` — why things are the way they are
- `claude/KNOWN_ISSUES.md` — surprises and gotchas

---

## Release roadmap

| Milestone | Branch | Status |
|-----------|--------|--------|
| **v3.3.2** | `develop` → `main` | **SHIPPED 2026-04-05** — tagged, PyPI + conda-forge published; propagated to GA remote |
| **v3.3.8** | cherry-picks → `main` | **SHIPPED 2026-07-11** — latest 3.3.x patch; tagged + PyPI (see session 48) |
| **v4.0.0** | `develop` → `main` | **SHIPPED 2026-08-30** — tagged, PyPI (see `RELEASE_PLAN_4.0.0.md`) |
| **v4.0.1** | `main` | **Tagged 2026-09-17** |
| **v4.1.0** | `develop` → `main` | Next. Users test tracers, sediment and docs first. The first release with validation reports + a Zenodo draft (`RELEASE_VALIDATION_REPORTS.md`) |

**v3.3.2:** Shipped. Includes EPSG/CRS support, utm→pyproj replacement, sww_merge fixes,
sww2vtu converter, pyproj DeprecationWarning fixes, ruff linting, riverwall throughflow,
NPY002 fixes, GDAL removal, regression snapshot tests.

**v3.3.8:** Shipped (2026-07-11). Patch release on the 3.3.x line: parallel
structure-operator logging (accumulated_flow column, unique operator numbering,
sequential-matching filenames) and the weir critical-depth gravity fix
(hardcoded 9.81 → `domain.g`) — both cherry-picked to `main` via PR #182 — plus
the SWW large-t store crash (#149), osx-64 wheels (#142) and checkout 6→7 (#146).
**Release procedure:** the publish workflow uploads to PyPI on a *published
GitHub Release*, **not** on a bare tag push — so create the annotated tag
(`git tag -a 3.3.8 -m "ANUGA 3.3.8"`, bare version, no `v`), push it, then
`gh release create 3.3.8 --verify-tag`. conda-forge follows via its feedstock bot.

**v4.0.0:** `feat/sc26` has been merged into `develop` (2026-04-01). `develop` is now
the active working branch. feat/sc26 contains GPU/OpenMP-offloading work
(`multiprocessor_mode=2`) forming the basis of a **Supercomputing 2026 (SC26)** paper.

## Active branches

| Branch | Purpose |
|--------|---------|
| `main` | Stable — v3.3.1 release |
| `develop` | Active development for v4.0.0 — contains GPU work + ADER-2 |
| `develop_sc26` | Working branch for GPU/SC26 incremental improvements |
| `develop_gpu` / `develop_cupy` | Earlier GPU experiments (CuPy-based) |
| `experiment/claude_culvert_refactor` | Culvert structure refactoring experiment |

Target PR branch is `develop` for all new work going into v4.0.0.

`develop_ader` merged into `develop` 2026-04-29.

---

## Common tasks

### Run tests
```bash
pytest --pyargs anuga                    # full suite (~163s)
pytest --pyargs anuga --run-fast         # skip slow tests (~41s)
pytest --pyargs anuga -m slow            # only slow tests
pytest anuga/shallow_water/tests/test_shallow_water_domain.py  # single file

# Per-test process isolation (required on a GPU build; works on any build).
# -cm legacy|unified sets the default compute mode for every child.
anuga_run_isolated_tests --pyargs anuga.shallow_water -cm unified   # 437 pass, 2 skip
anuga_run_isolated_tests                                            # the GPU file
```

⚠️ **The fast suite SKIPS `test_DE_gpu_omp.py` (86 tests) on a GPU-offload build.** A green
`--run-fast` on the nvc build is *not* evidence for anything in that file — it never ran it.
Tell the two apart by the counts:

| build | `pytest --pyargs anuga --run-fast` |
|---|---|
| GPU offload (nvc) | 2610 pass / 218 skip ← the 86 are skipped |
| CPU-only (gcc) | **2696 pass** / 219 skip ← the 86 actually run |

So **validate any change to `test_DE_gpu_omp.py` on a CPU-only build** — which is what CI
builds, and which is how a green local run still turned CI red on #192:

```bash
pip install --no-build-isolation -e . -Csetup-args=-Dgpu_offload=false   # CI's build
```

On a GPU build, use `anuga_run_isolated_tests` (it opts in via `ANUGA_GPU_TESTS_ISOLATED`)
to exercise the file — but that is a *complement* to the CPU-build run, not a substitute:
some assertions come out differently on the two builds (see the order-sensitivity trap in
session 49).
See `CLAUDE.md` → "Testing a GPU-offload (nvc) build" for the full GPU recipe.

### Build
```bash
conda activate anuga_env_3.14
pip install --no-build-isolation -e .

# CPU multicore (distribution default — fast on CPU). gcc, host omp parallel for:
CC=gcc pip install --no-build-isolation -e . -Csetup-args=-Dgpu_offload=false

# GPU offload (fast on GPU). nvc from the NVIDIA HPC SDK:
CC=nvc pip install --no-build-isolation -e . \
    -Csetup-args=-Dgpu_offload=true -Csetup-args=-Dgpu_arch=cc120
```
**CPU and GPU are separate builds** — one nvc binary can't be fast on both (NVHPC's OpenMP
host fallback is single-threaded; see KNOWN_ISSUES). For CPU performance use the gcc
`gpu_offload=false` build; the GPU build is for GPU runs (and correctness A/B on CPU only).

### Compute model (v4.0) — two orthogonal knobs
```python
# Per-domain: which compute path (legacy mode 1 vs unified C kernels mode 2)
domain.set_compute_mode('legacy')    # = set_multiprocessor_mode(1): openmp_ext + Python ops
domain.set_compute_mode('unified')   # = set_multiprocessor_mode(2): unified gpu_ext C kernels

# Process-global (call before the first evolve()):
anuga.set_gpu_offload(True/False)    # offload unified kernels to GPU (GPU build only)
anuga.set_omp_num_threads(16)        # OpenMP thread count for the whole process
anuga.gpu_offload_enabled()          # resolved offload state
domain.compute_capabilities()        # {gpu_offload, num_gpu_devices, mpi, modes}
```
`cpu` = unified + offload off; `gpu` = unified + offload on (compositions, not modes).
CLI: `-mpm 1|2`, `-nt N`, `-go`/`-ngo` (gpu offload on/off), `-ro metis_rcm` (reorder).
Rationale + the nvc finding: `claude/DECISIONS.md` and `claude/KNOWN_ISSUES.md`.

### Check code quality
```bash
pyflakes anuga/path/to/module.py
autopep8 anuga/path/to/module.py
```

### Transmissive riverwalls (`Cd_through`) — issue #32

Riverwalls can leak *through* the wall body (below the crest), not just overtop —
via a per-riverwall `Cd_through` discharge coefficient in the hydraulic-parameter
dict (alongside `Qfactor, s1, s2, h1, h2`):

```python
riverWall_Par = {'fence': {'Cd_through': 0.5}}   # 0.0 (default) = impermeable
domain.riverwallData.create_riverwalls(riverWall, riverWall_Par)
```

- Physics: submerged orifice `Q = Cd_through · h_eff · √(2g·|Δstage|) · sign(Δstage)`,
  `h_eff` = upstream (driving-side) submerged depth below the crest (so it flows
  even when the downstream side is dry). Momentum contribution is zero (conservative).
- Applied **on top of** the weir/overtopping discharge, so a transmissive wall does
  seepage *and* overtopping. `Cd_through=0` reproduces the old impermeable behaviour.
- Lives in the **shared** flux kernel `anuga/shallow_water/gpu/core_kernels.c`
  (`gpu_adjust_edgeflux_with_throughflow`, hydraulic-properties **column 5**), which
  the legacy `_openmp_compute_fluxes_central` (mode 1) *and* the unified path (mode 2)
  both call — so it works in **both compute modes with bit-identical results**. Do
  NOT look for it in a `sw_domain_openmp_ext.c` (that's a build artifact and doesn't
  exist). Added in commit `6ebb4453`; documented on issue #32.

---

## Benchmark timings — Towradgi small (MSI laptop, RTX 5070, AMD Ryzen 9, 2026-06-11)

Case: `run_small_towradgi.py -ft 200 -ys 50`, ~256k triangles, DE1 algorithm.

### Baseline and reference

| Mode | Config | Time (s) | Speedup |
|------|--------|----------|---------|
| Serial | 1 rank / 1 thread | 96.27 | 1× |
| OpenMP (no reorder) | `OMP_NUM_THREADS=16`, mode=1 | 22.73 | 4.2× |
| MPI | `mpiexec -np 16`, mode=1 | 12.13 | 7.9× |
| MPI + RCM reorder | `mpiexec -np 16 -ro rcm` | 11.08 | **8.7×** |
| GPU | mode=2 (RTX 5070, cc120) | 6.25 | **15.4×** |
| GPU + Hilbert reorder | mode=2, `-ro hilbert` | 5.62 | **17.1×** |
| GPU + metis_rcm reorder | mode=2, `-ro metis_rcm` | 5.79 | 16.6× |

### Reordering comparison (all `OMP_NUM_THREADS=16`, mode=1)

| Reorder method | CLI | Time (s) | vs no-reorder |
|----------------|-----|----------|---------------|
| None | — | 22.73 | — |
| Metis-16 | `-ro metis` | 19.74 | −13% |
| Metis-Hilbert-16 | `-ro metis_hilbert` | 18.73 | −18% |
| Hilbert | `-ro hilbert` | 18.53 | −18% |
| Morton | `-ro morton` | 18.36 | −19% |
| RCM (global) | `-ro rcm` | 18.24 | −20% |
| **Metis-RCM-16** | **`-ro metis_rcm`** | **17.43** | **−23%** |
| Metis-RCM-24 | `-ro metis_rcm -rn 24` | 18.28 | −20% |
| Metis-RCM-32 | `-ro metis_rcm -rn 32` | 18.50 | −19% |

**Best reorder by mode: OpenMP → `metis_rcm` (17.43 s, 5.5×); MPI → `rcm` (11.63 s); GPU → `hilbert` (5.62 s, 17.1×).**

### Multi-GPU strong scaling — gadi (NVIDIA V100, GPU-aware MPI), 2026-06-17

First multi-GPU run via the GPU-aware-MPI build (`-Dgpu_aware_mpi=true`, nvc,
cc70). Times are the **evolve-loop** wall time for `run_small_towradgi.py -mpm 2
-ro hilbert -ft 3600` (18× the laptop table's `-ft 200`). The absolute seconds
are therefore NOT comparable to the 5.62 s reference above — only the V100
scaling below is meaningful.

| Config | Evolve (s) | Speedup vs 1×V100 | Parallel eff. |
|--------|-----------:|-------------------|---------------|
| RTX 5070 (1 GPU, local) | 200 | — | — |
| V100 ×1 | 151 | 1.00× | — |
| V100 ×2 | 105 | 1.44× | 72% |
| V100 ×4 | 83 | 1.82× | 45% |

- A single **V100 (151 s) beats the RTX 5070 (200 s)** — this solver is
  memory-bandwidth-bound and V100 HBM2 (~900 GB/s) > 5070 GDDR7.
- Strong scaling tails off (72% → 45%): "small" Towradgi (~257k tris) gives only
  ~64k cells/rank at 4 GPUs, so halo-exchange/compute ratio rises and the GPUs
  underutilise. A larger mesh should scale substantially better.

### CPU multicore via the unified gpu_ext C kernels (mode=2, gpu_offload=false)

The operators (rainfall/culverts/inlet) run as **serial Python** in mode=1 — profiling showed
they cost ~8.5 s of a 21 s run (~40%), and do NOT scale with OMP threads. This is the main
reason MPI beats OpenMP: MPI splits operators across ranks.

The `shallow_water/gpu/` C kernels (`gpu_rate_operator.c`, `gpu_culvert_operator.c`,
`gpu_inlet_operator.c`) implement these operators with `OMP_PARALLEL_LOOP_*` macros
(`gpu_omp_macros.h`) that compile to `#pragma omp parallel for` when built with
`gpu_offload=false` (`-DCPU_ONLY_MODE`). They are dispatched only when `multiprocessor_mode=2`.
So a `gpu_offload=false` build + `-mpm 2` runs the *entire* step — solver AND operators —
as multicore CPU OpenMP.

Build: `CC=gcc pip install --no-build-isolation -e . -Csetup-args=-Dgpu_offload=false`

| Config | Time (s) | Notes |
|--------|----------|-------|
| mode=1 (Python ops) + metis_rcm | 18.64 | previous best OpenMP |
| mode=2 (C ops) + no reorder | 17.03 | C operators alone |
| **mode=2 (C ops) + metis_rcm** | **12.27** | **both — nearly matches MPI** |
| MPI-16 + rcm (reference) | 11.08 | — |

The two optimisations are **super-additive**: C operators alone save ~1.6 s, reorder alone
~4 s, but together ~6.4 s — because in mode=2 the gpu_ext *solver* kernels are also strongly
cache-sensitive (17.03→12.27 s with reorder) while the operators become parallel C.
This closes the OpenMP→MPI gap to within ~11%.

**Re-confirmed 2026-06-13** on the RTX 5070 laptop with a fresh gcc `gpu_offload=false`
build (this session's compute-model work):

| `-mpm 2` config | no reorder | `-ro metis_rcm` | reorder win |
|-----------------|-----------:|----------------:|------------:|
| `-nt 1`  (serial)   | 95.12 s | 46.22 s | **2.06×** |
| `-nt 16` (16 cores) | 19.18 s | 11.27 s | **1.7×**  |

11.27 s now *matches* MPI-16+rcm (11.08 s) with no MPI setup. Reorder helps more
serially (cache locality dominates) than at 16 threads (closer to memory-bandwidth
bound on this single-NUMA box). Thread scaling 1→16: ~5.0× (no reorder), ~4.1×
(reorder). Validates a stock single-node `pip install` + `-mpm 2 -ro metis_rcm` as
the CPU path the migration targets.

**Note:** `gpu_offload=false` overwrites the GPU build. Rebuild with
`-Dgpu_offload=true -Dgpu_arch=cc120` (and `CC=nvc`) to restore GPU mode.

**Migration plan:** making `mode=2 + gpu_offload=false` the standard distribution
default is tracked in `claude/PLAN_default_mode2_cpu.md`. Step 1 (deferred interface
build) is in review as PR #144; step 2 (audit operator fall-back) is next.

Optimal reorder differs by execution model: CPU sequential traversal benefits from RCM
graph-bandwidth minimisation; GPU warp-parallel execution benefits from Hilbert's tight
spatial clustering for coalesced memory access. metis_rcm (5.79s) is slightly worse than
hilbert (5.62s) for GPU but better than no reorder (6.25s).

### MPI partition scheme comparison + reorder

| Config | Time (s) | vs default |
|--------|----------|------------|
| MPI-16, default (`-ps metis`) | 12.13 | — |
| MPI-16, `-ps morton` | 13.29 | +10% |
| MPI-16, `-ps hilbert` | 13.68 | +13% |
| MPI-16, `-ps rcm` | 15.59 | +29% |
| **MPI-16, `-ps metis -ro rcm`** | **11.08** | **−9%** |

**Metis must remain as the partition scheme** — it explicitly minimises inter-rank edge
cuts, which determines communication cost. RCM/Hilbert/Morton as partition schemes hurt
because they optimise spatial/topological locality but ignore communication.

The winning combination is orthogonal: Metis to partition (minimise communication) +
RCM to reorder within each rank (minimise cache misses). Benefits stack: 11.08 s vs
12.13 s baseline.

### OMP binding flags and numactl — do not use

Thread pinning (`OMP_PROC_BIND`, `OMP_PLACES`) slows culvert/polygon setup, which uses
scipy/numpy internally and benefits from free core migration. `numactl --interleave=all`
caused only 2 active threads → 57 s (worse than serial); likely a GCC OpenMP
affinity-detection quirk on single-NUMA systems.

Hardware: **single NUMA node**, 32 CPUs, 30 GB — no NUMA effects to exploit.

### Why MPI-16 still beats OpenMP-16+metis_rcm (~40% gap)

With a single NUMA node the gap is NOT memory topology. Remaining causes:
- **False sharing**: OpenMP thread chunk boundaries share 64-byte cache lines in the flux
  arrays. MPI processes have separate address spaces so this cannot happen.
- **Single-threaded operators**: culverts, rainfall, boundary conditions run serially.
- **Reduction overhead**: global stage max, wet count, water volume require barriers.

The false sharing fix requires padding flux array chunks to cache-line boundaries in the
C kernel — beyond Python-side mesh reordering. **17.43 s with metis_rcm is the ceiling
for Python-side OpenMP optimisation on this hardware.**

Key findings:
- RCM beats Hilbert/Morton because it minimises graph bandwidth (neighbour-index distance),
  directly targeting the flux kernel's `domain.neighbours[i]` access pattern.
- Metis-RCM beats global RCM because Metis produces compact sub-graphs; local RCM wavefronts
  stay tight rather than snaking across the whole mesh.
- Sweet spot for OpenMP is `n_procs = OMP_NUM_THREADS`: more partitions introduce seams
  within each thread's working set and hurt more than they help.
- For MPI, reordering applies after distribute() to each rank's local ~16k-triangle mesh;
  benefit is smaller but positive (5.5%).

---

## GPU hardware comparison — RTX 5070 vs V100, per-kernel (nsys, 2026-07-16)

Case: `run_small_towradgi.py -mpm 2 -go -ft 240`, ~256k tri, **OpenMP-target on both**
(`nvc -mp=gpu`), single rank/GPU. Identical 4736-instance counts on both profiles ⇒
same work, directly comparable. RTX 5070 = laptop (8GB GDDR7, ~450 GB/s); V100 = gadi
SXM2 (HBM2, ~900 GB/s), profiled by Steve on gadi.

**Headline: V100 only 1.16× faster in TOTAL GPU kernel time (4.64 s vs 5.40 s) — the 5070
delivers 86% of V100 kernel throughput, nowhere near the ~2:1 raw-DRAM-bandwidth ratio.**

| kernel | 5070 ms | V100 ms | faster |
|---|--:|--:|:--|
| **compute_fluxes_central** (dominant) | 2303.9 | 2690.2 | **5070 1.17×** |
| extrapolate_second_order_edge | 1746.1 | 938.1 | V100 1.86× |
| update_conserved_quantities | 585.8 | 182.7 | V100 3.21× |
| manning_friction | 279.7 | 120.4 | V100 2.32× |
| rate_operator_apply | 239.2 | 353.0 | 5070 1.48× |
| **TOTAL GPU kernel** | **5399.5** | **4635.4** | **V100 1.16×** |

Two camps. **Memory-bound** kernels (extrapolate, update_conserved, manning — gather /
scatter / elementwise) go to the V100 by 1.9–3.2× — that IS the HBM2 bandwidth advantage.
But the **compute-bound** `fluxes_central` (Riemann flux math, the single largest cost —
43% of 5070 time, 58% of V100 time) goes to the **5070**: it is not bandwidth-limited, so
Blackwell's newer/higher-clocked SMs win the arithmetic. Because that kernel dominates, the
5070's win there nearly cancels the V100's bandwidth wins elsewhere.

**Lesson (and a correction to record).** When Steve estimated the 5070 at ~75% of the V100
from wall-clock (200 s vs 150 s at ft=3600), I pushed back that the real DRAM ratio is ~50%
(V100 2×). I was right about the *spec sheet* — the memory-bound kernels confirm ~2× — but
wrong about what matters for *this* code: towradgi is not bandwidth-bound end to end (its
heaviest kernel is compute-bound), so the achieved ratio is 86%, not 50%. Steve's effective
estimate was closer than my bandwidth correction. **Measure the kernels; do not infer GPU
throughput from either the spec sheet or wall-clock.** Also note the kernel ratio (V100
1.16×) is *smaller* than the wall-clock ratio (1.33×) — the ~1.15× remainder is host-side
(gadi's server CPUs beat the laptop on the 30–70% non-kernel setup / Python / I/O portion).

**Method, reusable:** for identical work, `time_A / time_B` per kernel *is* the inverse
achieved-throughput ratio (no need for absolute GB/s). Match kernels by base name after
stripping the `nvkernel_` prefix and `_F<n>L<n>_<n>` source-location suffix (the two builds'
line numbers differ). Extract with `nsys stats --report cuda_gpu_kern_sum --format csv`.

### FP64 is the hidden bottleneck on GeForce — ncu, 2026-07-18 (issue #199)

Followed up the "memory-bound" kernels with `ncu`, and the bandwidth framing was **wrong**.
On the RTX 5070 these DP kernels are **FP64-COMPUTE-bound, not bandwidth-bound**:

| kernel (5070) | DRAM % | FP64 pipe % | verdict |
|---|--:|--:|---|
| update_conserved_quantities | **12%** | **89%** | FP64-pipe-bound |
| extrapolate_second_order_edge | 22% | 81% | FP64-pipe-bound |

GeForce FP64 is **1/64 rate** (vs datacenter V100 **1/2**), so ANUGA's double-precision
kernels — full of divisions — become **FP64-throughput-limited** on consumer cards even
though they *look* memory-bound and are memory-bound on a V100. This is **why the 5070/V100
kernel ratios varied** (update_conserved 3.2× — worse than the ~2× DRAM ratio — is the most
FP64-bound one): the "achieved 55% of peak bandwidth" number was a red herring; the kernel
is slow because it is doing FP64 division, not because it fails to move memory.

**Actionable lever: cut FP64 divisions.** `update_conserved_quantities` did 6 FP64 divisions
per element; reformulating `cv/(1-dt*siu/c) == cv*c/(c-dt*siu)` (halving to 3) gave **2.85×**
on the 5070 (123.7 -> 43.4 µs/call, closing the gap to the V100's 38.6 µs), and helps the
**CPU unified path** too (shared kernel). Committed `2d7893ee`. `extrapolate` (still 81% FP64)
is the obvious next candidate — same treatment (count and reduce its DP divisions) should pay
off similarly on GeForce. **General rule:** on a GeForce/consumer GPU, profile the DP kernels'
`sm__pipe_fp64_cycles_active` before assuming they are bandwidth-bound; division count, not
bytes, is often the lever. (Datacenter GPUs with fast FP64 are unaffected — they stay
bandwidth-bound, so this is a consumer-hardware optimisation, not a portable win.)

**Trap that cost time here, recorded for reuse:** a roundoff-level change to a *shared* kernel
perturbs the chaotic trajectory of every mode-1-vs-mode-2 comparison test. It collapsed the
#192 order-sensitivity test (whose signal was pure chaotic amplification — the flat-bed
culvert was a no-op, so operator order never mattered *deterministically*). Fix was to rebuild
that test on a setup where order matters by a **large deterministic margin** (sloped bed ->
active culvert; rate confined to one inlet), giving a chaos-free 1e-15-vs-7e-3 fix/bug signal.
**Lesson (again): never let a guard depend on chaotic-amplification magnitude** — assert the
deterministic property (mode-2 [rate,culvert] matches mode-1 [rate,culvert], not [culvert,rate]).

---

## Weak-scaling benchmark — rectangular, constant timestep (MSI laptop, RTX 5070, AMD Ryzen 9 32C/30GB, 2026-06-26)

Case: `examples/parallel/run_parallel_rectangular_weak_scaling.py` — the extent grows with
the grid (`length = cell_size * sqrtN`, `cell_size=0.002`) so the triangle size, and hence
the explicit CFL timestep, stays ~constant as the triangle count grows (true weak scaling).
Each `-sn` doubling is ×4 triangles. DE0, `finaltime=0.015` (~276–318 steps), no reorder.
`dt` drifts only `4.72→5.44e-5` across the range (the fixed central bump is a smaller
fraction of larger domains). `-mp 2` offloads to the GPU; `-sn 1900` (14.44M, ~6.7 GB) is
the largest that fits the 8 GB RTX 5070 — `-sn 2000` (16M) OOMs the GPU.

16 cores where applicable: MPI = `mpiexec -np 16` × 1 thread/rank; OpenMP = serial `-mp 1`,
`OMP_NUM_THREADS=16`; GPU = `-mp 2`, np=1.

### Evolve time (s)

| triangles | MPI-16 | OpenMP-16 | GPU |
|-----------|-------:|----------:|----:|
| 0.90M (sn 475)  | 5.01  | 5.75  | 2.25 |
| 3.61M (sn 950)  | 22.03 | 22.56 | 8.03 |
| 14.44M (sn 1900)| 92.52 | 96.84 | 31.21 |

### Cost per (million-triangles · step) — flat = ideal weak scaling

| triangles | MPI-16 | OpenMP-16 | GPU |
|-----------|-------:|----------:|----:|
| 0.90M  | 17.5 ms | 20.1 ms | 7.86 ms |
| 3.61M  | 21.5 ms | 22.0 ms | 7.83 ms |
| 14.44M | 23.2 ms | 24.3 ms | 7.84 ms |
| **weak-scaling eff.** | 100→75% | 100→83% | 100→**~100%** |

### GPU evolve speedup

| triangles | vs MPI-16 | vs OpenMP-16 |
|-----------|----------:|-------------:|
| 0.90M  | 2.23× | 2.55× |
| 3.61M  | 2.74× | 2.81× |
| 14.44M | **2.96×** | **3.10×** |

Findings:
- **GPU weak-scales near-perfectly** (~7.84 ms/(Mtri·step), <0.5% spread over 16× size),
  2.2–3.1× faster than 16 CPU cores, lead widening with size.
- **MPI-16 evolve is only ~2–13% faster than OpenMP-16** here — this is a *pure-solver* case;
  the larger MPI>OpenMP gap on Towradgi came from serial operators (culverts/rainfall) that
  this rectangle has none of. Both are memory-bandwidth bound and degrade with size.
- **MPI setup overhead is large and grows**: `distribute` 1.3 → 5.8 → **29.9 s** plus rank-0
  `creation` 0.4 → 1.9 → 9.6 s. End-to-end at 14.44M: MPI 132.1 s vs OpenMP 104.4 s vs GPU
  38.7 s — for these short (~276-step) runs MPI's faster evolve does NOT repay its ~30 s
  distribute; on production-length runs the one-time setup amortizes and MPI's evolve edge wins.

---

## Recent session summaries (sessions 21–58)

**Session 58 (2026-10-01/03):** Roll waves in the validation suite, why
DE_ader2 never quite settles, and validation reports as part of every
release. (Work between sessions 57 and 58, including 4.0.0 shipping and the
sediment cases, is not summarised here; see the git log and `ROADMAP.md`.)

- **rundown_mild_slope shows roll waves (PR #407, merged).** The oscillation
  DE_ader2 showed down the slope is not the outflow boundary: q = 0.2,
  n = 0.03 gives Froude 2.26, above the Manning roll-wave threshold of 1.5.
  In 100/200/400 m channels the onset stays the same distance from the
  inflow. The case now runs a stable n = 0.06 case (Froude 1.21, analytic)
  and a 200 m roll-wave case with a 1% inflow seed (period 3.3 s, chosen to
  avoid aliasing with the 2 s yieldstep). `compare_solvers.py` runs the
  roll-wave case for DE0/DE1/DE2/DE_ader2: with the same seed the
  second-order schemes grow it at about 0.055/m, and DE0 about 8x slower.
- **DE_ader2's leftover unsteadiness (KNOWN_ISSUES).** Unseeded, DE0, DE1 and
  DE2 reach a bit-exact steady state, while DE_ader2 does not, because its
  fluxes are taken on Q + (dt/2) Q_t. Its steady state therefore moves with
  dt (4e-3 h0 at the boundary cells for dt 0.04 vs 0.08), and every dt
  change (output clipping; the roll waves swinging the CFL dt) sends a pulse.
  Two prototypes were tried and not adopted: using the current dt in the
  predictor (cuts the noise only about 3x), and holding dt (works, but fragile,
  and costs 15-40% more steps).
- **Validation reports per release (PR #408, merged).**
  `validation-report.yml` builds DE0/DE1/DE_ader2 reports from the tag on a
  published release, attaches them, and for final releases makes a Zenodo
  DRAFT (`.github/scripts/zenodo_validation_reports.py`, never
  auto-publishes). It is a gate, which needed `produce_report()` and
  `validations_produce_results.py` to stop swallowing failures. Zenodo:
  personal accounts plus the `anuga` community. Tokens are set, and the
  sandbox dry run (run locally) passed both the create and new-version paths.
  Setup and steps: `claude/RELEASE_VALIDATION_REPORTS.md`.
- **Combined report (PR #409, merged).** One report with every case once and
  each differing figure as DE0 | DE1 | DE_ader2 panels, plus a run-time
  summary. `validations_report_body.tex` is now shared by both reports
  (add new cases there). `validations_combine_algs.py collect/combine`.
  Report runs use `reports/matplotlibrc` (13 pt, 200 dpi) so legends survive
  the third-width panels. `subcritical_depth_expansion` was being run but had
  never been in the report; it is now.
- **Citation and Zenodo software DOI (#410, #412, merged).** Zenodo's GitHub
  integration is on, so each published release is archived as software from
  `CITATION.cff`. The file was invalid (the licence must be SPDX `Apache-2.0`)
  and is fixed, with the title "ANUGA: hydrodynamic modelling of floods, storm
  surges and tsunamis". The lint workflow now validates it (`cffconvert`).
- **4.1.0 preparation.** The call for testing is pinned issue #411. The
  release notes are completed (#413): an opening, highlights, DE_ader2's
  changed results, the audit fixes, Requirements, Known issues, and Thanks
  with a placeholder for the testers. #406 (the corrected DE0 wave-damping
  figures) is merged. The merged branches are cleaned up.
- **Zaragoza pier wake (#414, merged).** Right behind the pier the flume
  scours where the model deposits. Suspension changes the bed by under 1 mm,
  and a resolved wake (103k triangles, St 0.20, recirculation) still fills
  the lee, so the lee scour is 3D. The case now adds a wake tail on the
  shear amplification (A_w = 2 over 5 D, a calibration): ring 7.9/11.8 mm
  against 10.4/12.6, RMS 3.2/3.6.
- **Lessons.** `workflow_dispatch` only works for workflows on the default
  branch (`main`), so test a new release workflow's script locally. The GPU
  build in unified mode is no faster on the validation suite (17-24 min per
  algorithm against 19 on 8 CPU threads), because the meshes are small. For small
  meshes (the Zaragoza case, 18.8k triangles: 30 s on the CPU) say so before
  choosing the CPU; the user expects the GPU for runs.
  `/tmp` is a 16 GB tmpfs with a quota, so full report runs belong in an on-disk worktree.
  `gh pr merge` is blocked by the permission classifier, so the user merges.

**Session 57 (2026-08-21/23):** **4.0.0 release engineering** — Phases 1 and 2
complete, RC published, and three separate upstream breakages absorbed along the
way. None of the four defects found were in ANUGA's solver.

- **Review + release plan.** An engineering review of the whole repo (published
  as an Artifact) produced `claude/RELEASE_PLAN_4.0.0.md`: go straight to 4.0.0,
  no 3.4 bridge, because `main` and `develop` already share the packaging
  foundation and develop still defaults to `legacy` compute mode — so there is
  no behavioural cliff. The scary parts of "v4" (mode-2 default flip, forcing
  removal) are explicitly deferred to 4.1.
- **Phase 1: nine verification gates, all green.** Full suite both builds
  (2937 / 2948), unified-mode one-process (2947), GPU isolated runners on TWO
  architectures (local RTX 5070 cc120 and cloud g6.xlarge Ada L4 — identical
  counts), validation suite (120 passed), Towradgi mode-1 vs mode-2, MPI,
  wheels, fresh conda 3.10/3.14, docs.
- **Phase 2: notes, guide, citation, RC.** `RELEASE_NOTES.md` and
  `UPGRADING_TO_4.0.md` written and approved; `CITATION.cff` was stale at
  version 3.1.9 / 2022-06-26 with a 2022 commit pin. `4.0.0rc1` published to
  **TestPyPI** — 21 artifacts, cp310–cp314, Linux/macOS(arm64+x86_64)/Windows —
  built via the `ANUGA_VERSION` override so **no rc tag** was needed on develop.
  Verified installable in a clean venv by the exact command the announcement
  gives people.

**Four defects found, none in the solver:**

1. **`run_gpu_tests_isolated.sh` swallowed failures** — `rc=$?` after a pipeline
   read `tail`'s status, so every failing class reported "ok" and the script
   exited 0. This is why several "all GPU classes green" claims earlier were
   worthless.
2. **`anuga_run_isolated_tests` crashed against an installed anuga** — all 106
   tests CRASH in a container. pytest emits node ids relative to ITS rootdir
   (the common ancestor of the targets); the runner resolved them against its
   own cwd-derived ROOTDIR. In a checkout the two coincide, which is why it was
   invisible locally. Now resolves against target-derived base dirs and **fails
   loudly** when an id resolves to no file.
3. **gcc/nvc disagree on line-region triangle selection** (#231) — a culvert
   exchange line lying exactly on mesh edges selected 32 vs 54 triangles for the
   same script. Test made geometry-robust.
4. **Concurrent writes to one SWW silently interleave** (#232) — `domain_name`
   is hardcoded in the case-study scripts, so two runs corrupt each other's
   output with a non-monotonic time axis and plausible-looking numbers. Cost a
   full invalid gate-5 result. Found via the frame-index pattern (writer A at
   0,2,4-8; B at 1,3).

**Three upstream breakages, none our doing:**

- **conda-forge mingw build 11** broke every Windows job for two days: binaries
  compile and then abort in the stack protector (`__stack_chk_fail` on an empty
  `main`). Took five CI cycles to converge because the ground kept moving — the
  compilers were rebuilt mid-investigation, so a sysroot-only pin worked for a
  day and then failed. **Only one self-consistent set exists**: gcc/gxx_impl
  build 2 + sysroot chain build 10 + `libwinpthread` build 10 + `libstdcxx`
  build 2. Pinning one half looks like it works (wheels build) and then every
  test process dies silently at exit 127, because libwinpthread bites at RUN
  time. Removal tracked in **#235**; upstream is
  conda-forge/m2w64-sysroot-feedstock#23, where our reproduction is cited.
- **Cython 3.3.0** broke all 20 wheel builds: `dict.items()` can no longer be
  assigned to a `cdef list`. `sparse_matrix_ext.pyx` had one such assignment
  that was **never used**, plus a second on the next line that would have failed
  next. Local builds stayed green on Cython 3.2.5.
- **Branch protection** turned out to be enforced by *rulesets*, not classic
  protection — removing the classic setting did nothing, which only surfaced
  when a merge was refused.

**Also:** #237 (bare `ModuleNotFoundError` from an unbuilt tree) fixed;
#238 filed (make the meshpy import lazy — Triangle's licence makes it a hard
non-free dependency at import time); `claude/PPA_FEASIBILITY.md` written for
#25 (Ubuntu has every dependency except meshpy and pymetis; recommend against a
PPA since pip wheels, conda-forge and Docker now cover the original ask).

**RESUME STATE (as of 2026-08-23):**
- `develop` == `origin/develop` == `stoiver/develop` at `743b6369`. All CI green.
- **Waiting on**: RC soak feedback from Ole, Rudy, Petar, David and Jorge; and
  conda-forge/admin-requests#2297 landing so the Windows pin (#235) can go.
- **Petar's verdict on Towradgi can reopen the #229 physics decision**, not just
  the wording — it is his case study and it is the evidence behind the release
  notes' numbers. Reproduction recipe kept in
  `claude/RELEASE_4.0.0rc1_REPLY_PETAR.md`.
- **RC drift** (commits after 4.0.0rc1, tracked in the release plan): #237's
  import message and the Cython 3.3.0 fix. Both build/import level; neither can
  change results, so the soak stands. Anything that COULD change results means
  rebuilding the RC, not adding a row.
- Phase 3 is short and scripted in `claude/RELEASE_PLAN_4.0.0.md`: release PR to
  `main`, set `date-released` in CITATION.cff, annotated bare-version tag, then
  the GitHub Release triggers PyPI / Docker / the conda-forge bot.
- Open from this session: #231, #232, #235, #238. Closed: #189, #214, #224,
  #225, #229, #237.

**Session 56 (2026-08-20/21):** **Issue triage sweep → five issues closed**, one
of them the long-standing well-balancedness defect in every structure operator.

- **Triage of all 24 open issues** against the code. #214 (TOML interior holes +
  erosion) was already implemented on develop — closed citing the commits. #25
  got a comment pointing at the docker/ images + GHCR workflow (PPA half stays
  open). #143 awaits samcom12's re-run. The rest confirmed still open.
- **#224 riverwall crest → device (PR #227, merged).** `set_elevation()` wrote
  only the host array; mode-2 kernels kept the crest they were mapped with — a
  gate operated mid-run was silently ignored. New `gpu_sync_riverwall_to_device()`
  pushed by `evolve()` **at yieldstep boundaries only** (the user asked for
  that timing): once before the loop, once on resuming from each yield, gated by
  a `device_data_dirty` flag set by the three RiverWall setters. Also warns if
  `create_riverwalls()` runs after the interface exists (fresh arrays — the
  device still points at the old ones; not fixable by a push).
- **#189 degenerate-timestep protection (PR #228, merged).** Mode 2 never runs
  it: the C loops return before `update_timestep()`, and the Python-orchestrated
  GPU loops reach it but the stale host `max_speed` (device never syncs it back)
  makes its `<10` guard return — absent on every path, one of them by accident.
  Option 1 from the issue: warn once per domain, and the routine itself now
  warns-and-returns in mode 2 so the absence is a decision.
- **#225 MAX_INLET_TRIANGLES (PR #228, merged).** The 64-triangle inlet cap was
  fixed arrays in `culvert_indices` — host-side staging only; the device side
  was already dynamic. Now heap arrays sized by the real count, cap deleted.
  `Test_GPU_LargeInlet`: 79 triangles/inlet, modes agree ~1e-14.
- **#229 structure well-balancedness (PR #230, merged) — took two attempts.**
  Filed first: `Inlet.set_depths()` writes a uniform DEPTH, which on a sloping
  bed tilts a lake at rest by half the inlet's bed range per step at zero Q.
  **Attempt 1 (surface shift, e899973f) reverted** (f4e2ab11): broke the pipe
  culvert on all 15 CI jobs — flattening is the smoothing the culvert feedback
  loop needs; a dry receiving inlet grew a ridge and Q stalled to 0.
  **Attempt 2 (stage leveling): water finds its level** — fills raise the
  lowest stages, drawdowns lower the highest, clamped at the bed, bisection (no
  sort, same loop in numpy and one GPU team). Flat bed reproduces the old write
  **exactly**, so the pipe test passes by construction; zero transfer is
  bit-exact (lake at rest 6.8e-02 → 3.7e-09). Second trap: one momentum value
  written uniformly over leveled depths gave a near-dry cell 1.7e6 m/s and
  collapsed the timestep — `set_average_momenta()` now writes it depth-weighted
  (uniform VELOCITY). Full details in DECISIONS.md.
- **`run_gpu_tests_isolated.sh` lied (8662239c).** `rc=$?` after a pipeline read
  `tail`'s status — every failing class reported "ok", which is how attempt 1's
  failure reached CI instead of this desk. Fixed + verified against a deliberate
  failure. Lesson applied: attempt 2 was validated with the full suite exactly
  as CI runs it (2937 passed), not the fast subset.
- **New gotcha (KNOWN_ISSUES):** a structure in near-still water amplifies
  roundoff exponentially (1 ULP → ~1e-2 m stage in seconds) via the
  enquiry↔inlet feedback. NOT mode-related — a 1-ULP mode-1 perturbation
  diverges as fast as mode 2. Un-masked now #229 is fixed: CPU/GPU culvert
  comparisons need a driven head; well-balancedness tests stay short.
- **Process:** merges to `develop` don't fire `Fixes #N` (default branch is
  `main`) — #189/#224/#225/#229 all closed by hand with commit citations.

*(superseded by session 57's resume state above)*

**RESUME STATE (as of 2026-08-21):**
- `develop` == `origin/develop` == `stoiver/develop` at `7f98999a` (merge of
  PR #230). Topic branch `fix/229-well-balanced-inlets` deleted after merge.
- Merged this session: PRs #227, #228, #230. Closed: issues #189, #214, #224,
  #225, #229.
- **Still open from the triage:** #223 (gpu-aware MPI is host-staged), #170
  (per-substep boundaries in the C RK loop), #190 (blocked on draft PR #188),
  #143 (awaiting samcom12 re-run), #77, #141, and the older wishlist items.
- Uncommitted, deliberately: `towradgi.toml` arithmetic demo. Untracked:
  `sandpit/culvert_issue/` (48 MB reporter model — do not commit),
  `docs/source/reference/generated/`.

**Session 55 (2026-08-13):** A regression of our own making, the culvert crash,
and hardening the installers so neither can happen quietly again.

- **gpu culvert stack overflow (#217, PR merged via #218).** 101 culverts wedged
  at `Time = 0` then segfaulted on a NaN-bit-pattern pointer.
  `gpu_culverts_apply_all()` sized its per-step buffers as fixed stack arrays of
  `MAX_CULVERTS` (64) but looped over `num_culverts`: 24240 bytes into 15360, on
  every timestep. `MAX_CULVERTS` is only the *initial* capacity — registration
  grows by doubling. Buffers moved to the heap. Regression test
  `Test_GPU_ManyCulverts` (80 culverts) hangs on the old code, passes in 1.3 s.
- **`-Dgpu_arch` never selected anything (PR #220).** nvc device-links at the
  LINK step and `gpu_link_args` was empty, so it built for the build host's GPU,
  or — with no GPU visible, as in a container — every arch the CUDA supports.
  The "multi-arch" images were portable *by accident*. Fixed by repeating
  `-gpu=<arch>` at link.
- **…which then broke Petar's RTX 30xx (PR #221).** With the flag working,
  `install_anuga_nvc.sh`'s hardcoded `GPU_ARCH=cc120` produced sm_120-only
  builds that compile cleanly and crash at every kernel launch. Now detected
  from `nvidia-smi`.
- **Compatibility is ONE-DIRECTIONAL** (measured, not assumed): nvc embeds PTX
  and the driver JIT-compiles it *forward*, so a cc86 build runs on sm_120;
  a cc120 build can never run on sm_86. JIT costs +0.38 s at first launch, no
  throughput difference — matters where the cache is cold (isolated runner,
  containers, MPI ranks). `cc70`/V100 needs a CUDA ≤ 12.x in the SDK; CUDA 13
  dropped Volta.
- **Installers (PR #222, open).** New `tools/anuga_build_report.py` (also for bug
  reports) prints version/install type/GPU build/MPI/SM archs/GPU present, and
  `--check` fails a build that cannot run here. Both scripts log a transcript.
  `GPU_ARCH=auto` probes what the toolchain can build rather than assuming (nvc
  on Gadi is a module wrapper with no SDK layout beneath it). Explicit arch lists
  are never silently narrowed. `install_miniforge.sh` is re-runnable, detects a
  GPU and *advises* (does not switch compilers — the nvc host fallback is
  several times slower than gcc), `GPU=1` opts in. GPU tests are skipped when no
  GPU is visible (they were running on a Gadi login node).
- **`anuga:gpu-mpi` image** (2.69 GB, in ECR): HPC-X OpenMPI 4.1.7a1, CUDA-aware,
  mpi4py + pymetis, `-Dgpu_aware_mpi=true`. Four container gotchas solved and
  documented: `OMPI_CC=gcc`, `OPAL_PREFIX`, `OMPI_MCA_ess_singleton_isolated=1`
  (singleton `MPI_Init` hangs — bites *serial* runs, since `import anuga` pulls
  in mpi4py), and container transports (`pml=ob1, btl=self,smcuda,tcp`).
- **`.dockerignore` excludes `docker/Dockerfile*`** — they landed in the
  `COPY . .` layer, so a one-line Dockerfile edit invalidated the ~45 min ANUGA
  compile. Now 2.6 s with 14 layers cached.
- **Issue #223 (open):** `-Dgpu_aware_mpi=true` is not GPU-aware — `gpu_halo.c`
  stages every halo through host memory in *both* paths (UCX `uct_mm` SIGSEGVs
  on device pointers). Measure halo cost before implementing real GPUDirect.

*(Session 56's resume state below is superseded by session 57's above; and
session 55's by session 56's — PR #222
merged, and the issue landscape has moved.)*

**RESUME STATE (as of 2026-08-13):**
- `develop` = `stoiver/develop`, **8 commits ahead of origin/develop**, all in
  **PR #222** (installer hardening). CI was green on `5d53d602`; later commits
  only touch shell scripts, which CI does not exercise.
- Merged this session: #218, #220, #221. Closed: #219 (superseded — its fix rode
  in with #218), issue #217.
- **Open: PR #222, issue #223.**
- ECR (ap-southeast-2, `067589750115`): `anuga:gpu-slim` 460 MB (cc70–cc120,
  portable, AMD-verified on g6) and `anuga:gpu-mpi` 599 MB.
- Uncommitted, deliberately: `towradgi.toml` arithmetic demo. Untracked:
  `sandpit/culvert_issue/` (48 MB reporter model — do not commit).
- Gadi: build works from a login node; the open question is whether that
  `nvhpc` module can build `cc70` — run `tools/anuga_build_report.py --check`
  on a `gpuvolta` node to find out.

**Session 54 (2026-08-11):** Portable GPU image proven on AMD, C-level output
captured in the logfile, and the first GPU scaling numbers at 1.6M triangles.
- **`-tp=haswell` fix never actually built.** 394b61b8 set `CXXFLAGS=-tp=haswell`
  too, but only `CC` is nvc — `CXX` stays g++, which rejects `-tp` and kills
  meson's compiler detection (`Unable to detect linker for compiler c++ …`).
  That is why the ECR image was still the non-portable AVX-512 build. Fixed in
  `54b18f69`: **CFLAGS only** (`meson.build:166` already avoids `-march=native`,
  so the g++-built C++ extensions were always portable).
- **Portability PROVEN on AMD** — the thing that was outstanding all of session
  53. Run landed on **g6.xlarge** (AMD EPYC Zen 4 + L4): `rc=0`, no SIGILL,
  real GPU offload. g5 had zero capacity in Sydney all session; **g6 is the type
  that actually launches** there.
- **Logging: the logfile never contained C output.** `TeeStream` tees the
  `sys.stdout` *object*, so everything the C extensions printf to fd 1 (the GPU
  domain banner, Triangle, ~99 printf sites in the GPU kernels) reached the
  terminal but not the file — batch logs had no record of the GPU config used.
  `set_logfile()` now installs an **fd-level tee** (fd 1 → pipe → terminal +
  file). Needed `setvbuf`/`reconfigure(line_buffering=True)` (fd 1 is no longer
  a tty, so libc *and* Python block-buffer), and the logging file handler now
  writes **through the pipe** so the file has a single writer and exact write
  order. `file_only()` consequently mutes C output too. Adds `close_logfile()`;
  `set_logfile()` now **truncates** (logs used to accumulate every run ever).
  `anuga/utilities/tests/test_log.py` — subprocess-based, because pytest's
  capture replaces `sys.stdout` with something that never touches fd 1.
- **`-sc/--scale`** added to `run_small_towradgi.py` (10 = 135237 tri, 1 =
  256688, 0.1 = 1636238). Deliberately **in the script, not anuga's standard
  parser**: containerised runs take the script from the staged S3 inputs but
  anuga from the image, so a parser-level flag needs an image rebuild. First
  attempt (`c109e65d`) made exactly that mistake and failed on AWS with
  `unrecognized arguments: -sc 0.1`; reverted in `39e3e429`.
- **GPU scaling — small meshes flatter slow GPUs.** Same 600 s sim, mode 2:

  | GPU | 257k mesh | 1.64M mesh |
  |---|---|---|
  | RTX 5070 laptop (cc120) | 24.77 s / 61.4 M cell-updates/s | 111.52 s / **73.4 M/s** |
  | L4, g6.2xlarge (cc89) | 31.89 s / 47.7 M/s | 120.32 s / **68.1 M/s** |
  | T4, g4dn (cc75) | 53.06 s / 28.6 M/s | — |

  At 257k the 5070 looks 1.29× the L4; at 1.64M it is only **1.08×** — both
  cards are under-saturated at 257k and the L4 recovers more (+43% vs +20%
  throughput). Expect **g6 ≈ 93% of the laptop** on production-size problems.
  Memory at 1.64M is only **3.4 GB host / 762 MB GPU** (per-cell cost is small;
  the 1.6 GB at 257k is mostly fixed overhead — do not extrapolate linearly).
  Results **bit-identical** across cc120/cc89, native vs container.
- **towradgi's constant `delta t` is a mesh artifact, not physics** — see
  KNOWN_ISSUES.

**RESUME STATE (as of 2026-08-11, end of session 54):**
- `develop` pushed to **`stoiver`**; **PR #218** open to `anuga-community/develop`
  (14 commits at open; 2 more since — `-sc` flag add + revert/move).
- ECR `anuga:gpu-slim` is current: portable multi-arch build **including** the
  logging fixes (`sha256:ed7f18dcde0b…`).
- Uncommitted working-tree edit: `towradgi.toml` arithmetic *demo*
  (`finaltime="3600*24"`) — still not for committing.
- Open follow-ups: PR #218 CI is a **gcc `gpu_offload=false`** build (differs
  from the local nvc one); the terminal — not the logfile — is only
  approximately ordered under the fd tee; REDIST trimming in the slim image.

**Session 53 (2026-08-09 – 08-11):** Parallel bug fixes + Docker slim image +
first real AWS GPU run of towradgi.
- **MPI deadlock in `update_conserved_quantities`** (negative-cells "loss of
  conservation" warning called the collective `get_water_volume()` under a
  per-rank `if num_negative_ids>0` branch → `-np 1` worked, `-np 2` hung;
  faulthandler pinned it). Fixed: warn **serial-only** (no collective in the hot
  per-substep path; two earlier attempts — local volume, then unconditional
  Allreduce — each caused a *new* hang, incl. `test_parallel_boyd_box_operator`).
  → **PR #216, merged to develop.** Regression test
  `test_parallel_negative_cells_deadlock.py` added.
- **anuga_toml_run**: per-rank MPI logs; ranks≠GPUs unified guard (fail fast, not
  hang); end-of-run water-balance V0 fix (capture at first yieldstep, not before
  evolve — pre-evolve `get_water_volume` gives a bogus unclamped V0); `--emit-script`
  (eject an editable standalone run script). *(on branch add-toml-scenario-features / PR #213)*
- **Erosion**: default elevation to time-varying storage when erosion present —
  TOML parser (warn only on explicit static) **and** core operators
  (`_mark_domain_erosion`; `initialise_storage` warns if reset to static).
- **Riverwall verbose output** condensed to one line per wall.
- **TOML numeric fields** accept quoted arithmetic (`finaltime = "5*60"`, safe ast eval).
- **Docker**: slim multi-stage GPU image (CUDA **-base** + NVHPC REDIST + venv) =
  **2.15 GB / 459 MB compressed, ~37x smaller** to pull; validated locally. `-tp=haswell`
  portability fix (native nvc build SIGILLs on AWS g5/Zen2 — no AVX-512). `aws_run_gpu.sh`
  gained **`--ecr`** (private in-region image + role ECR-read + instance login).
- **AWS**: validated `run_small_towradgi.py` end-to-end via ECR slim image
  (`aws_run_gpu.sh --ecr --instance g4dn.xlarge`, rc=0, ~2 min). Coords in memory
  `reference_aws_towradgi_run.md`.

*(Session 53's resume state is superseded by session 54's above: the commits are
pushed, PR #213 merged, and the ECR image rebuilt portable.)*


**Session 52 (2026-07-25):** **Fix the CI unified-compute-mode failure + harden
`update_conserved_quantities`.** The `github CI` "Test package (unified compute mode,
full — Linux / Python 3.14)" step (CPU build, `ANUGA_DEFAULT_COMPUTE_MODE=unified`)
went red with 5 failures in `Test_negative_cells_warning`, all
`AttributeError: 'NoneType' object has no attribute 'update_conserved_quantities_kernel'`
at `shallow_water_domain.py:2948`.
- **Root cause.** Those white-box tests call `domain.update_conserved_quantities()`
  **directly, outside `evolve()`**, then inspect host `centroid_values`. Under the unified
  default the domain is mode 2 (`MULTIPROCESSOR_GPU`) but `gpu_interface` is only built during
  `evolve`, so the direct call hit `self.gpu_interface.…` with `gpu_interface` still `None`.
- **Fix 1 — pin the tests (commit `641db3bd`).** Added `domain.set_compute_mode('legacy')`
  in the test's `make_domain()` helper — the documented white-box pattern
  (`test_forcing.py`/`test_friction.py`). Build-agnostic: legacy just selects the CPU openmp
  solver per-domain, never touches GPU offload, so it's correct on the GPU-offload build too
  (and sidesteps the NVHPC mode-2 domain-count abort). Verified 8/8 pass in both default and
  unified. **CI confirmed fully green** (run 30145027044).
- **Fix 2 — harden the method (commit `c80bc457`).** `update_conserved_quantities` was the lone
  mode-2 entry point that read `self.gpu_interface` **without** first calling
  `_ensure_gpu_interface()` (its siblings `compute_fluxes`/`compute_forcing_terms`/extrapolate/
  distribute all guard). Added the guard: it builds the interface when boundaries are ready,
  else falls back to legacy. No-op during normal evolve. Verified a direct unified call now runs
  mode 2 instead of crashing; sw-domain evolve/conservation tests 14/14 green under unified.
- Both commits on `develop`, pushed to `origin` (anuga-community).

**Session 51 (2026-07-19/20):** **Unify the culvert kernel (mode-1 == mode-2), diagnose
Towradgi's residual as rain-on-dry, tighten comparison tolerances.** Rebuilt CPU-only
(gcc, `-Dgpu_offload=false`), confirmed the full suite green in both modes (legacy 2698 /
unified 2697), then closed out the mode-1-vs-mode-2 culvert difference and chased what's left.
- **Culvert unification (commit `b7a010ae`).** Extracted the per-culvert discharge +
  semi-implicit transfer into one C routine `culvert_compute_one()` plus a shared host inlet
  gather `culvert_gather_inlet_host()` (`gpu_culvert_operator.c`). Mode-2's batch calls it
  (behaviour-preserving); **mode-1's Python operators now route their per-step update through it
  via a Cython bridge** (`culvert_apply_one_host` / `culvert_gather_inlet_host_py`), gated to
  fully-local culverts (cross-boundary MPI keeps the Python path). Wired into **both**
  `Structure_operator` and — the one that matters — `Parallel_Structure_operator`. Result:
  **mode-1 == mode-2 bit-for-bit** for every culvert config (box/pipe, velocity head, blockage),
  1 and 16 threads. De-dups the Python/C physics on the runtime path.
  - **The seed was the inlet-average gather**, not momentum/ordering/FMA: mode-1 summed with
    numpy, mode-2 with a C loop, 1 ULP apart for multi-cell inlets.
  - **Two wrong-turn traps worth remembering:** (1) `anuga.Boyd_box_operator` is a **factory**
    returning `Parallel_Boyd_box_operator` — a *separate* class hierarchy from `Structure_operator`
    — even in serial, so a monkeypatch of `Structure_operator.__call__` silently hit an unused
    class and produced two invalid "proofs". Instrument the class that actually runs. (2) An FMA
    hypothesis was moot: the CPU build is generic x86-64 (no `-march`/`-mfma`).
- **Towradgi still diverges → it's rain-on-dry, NOT culverts.** Built an **in-process
  double-precision localization harness** (two domains from the same setup, restartable lockstep
  evolve, diff `centroid_values` in double). It pinned the seed to `Rate_operator` on **dry**
  cells: remove rain → bit-identical; rain on a **fully-wet** domain → bit-identical; rain on
  **dry** → 1 ULP. **Corrected framing:** both modes are *stage-primary* (`height == stage - bed`
  in each); the residual is a 1-ULP difference in the core's near-dry **stage itself**, masked by
  the dry-cell bed-clamp (`protect` sets `stage=bed`) until rain lifts the cell off the bed. Rate
  inputs and the device sync (a plain memcpy) are ruled out; exact operation not isolated (needs
  C-RK-loop device-state instrumentation). Same family as #200; **benign, accepted**. Red
  herring: the #200 init dry-cell clamp asymmetry — forcing mode-1 to clamp too left the
  divergence unchanged, so the initial state is not the cause. Documented in FUTURE_WORK /
  KNOWN_ISSUES (commit `aadbd8e7`).
- **GPU build validated.** Rebuilt nvc (`gpu_offload=true`, cc120, RTX 5070). Towradgi
  mode-1(CPU) vs mode-2(GPU): same-order divergence dominated by rain-on-dry (**wet-only ΔStage
  ≈ 0 at t=120** — only dry cells differ; wet cells diverge later via chaos). Isolated GPU test
  runner: **all classes pass** on the device, including the culvert/weir/dry-cell/fractional-step
  tests that touch the changed code.
- **Tightened mode-1/mode-2 comparison tolerances (commit `7a42bb2d`).** The culvert/weir
  CPU-vs-GPU tests used **atol=0.02/0.05** with comments blaming Python-vs-C ~3% drift — obsolete
  now that both run the same C kernel. Re-measured and tightened: culvert 0.02/0.05→1e-10/1e-9,
  weir 0.02→1e-10, inlet 1e-8→1e-10, forcing helper (wind/pressure/**rate**) 1e-8→1e-10 on GPU.
  Measured ≤1.2e-15 everywhere (rain bit-identical); kept 5–6 orders of headroom. These now catch
  a mode-1/mode-2 regression instead of silently passing multi-percent drift. All 17 pass under
  the isolated GPU runner. **Lesson:** when a comparison test's loose tolerance is justified by a
  "Python vs C FP drift" note, re-measure after any unification — the rationale may be stale by
  orders of magnitude. All commits on `develop`, pushed to `stoiver`.

**Session 50 (2026-07-16/19):** **mode-1 vs mode-2 difference notebook → two real GPU bugs.**
Built `validation_tests/case_studies/towradgi/compare_mode1_mode2.ipynb` (uses
`anuga.SWW_plotter`; white-background `YlOrRd` log difference maps, 3×1 stacked panels)
to show *where/when* legacy and unified diverge on Towradgi small. Chasing the plot
features found two genuine bugs, both fixed on `develop`.
- **BUG (#199, fixed `2d7893ee`): `update_conserved_quantities` FP64-bound on GeForce.**
  ncu (not nsys) is decisive: DRAM 12% / **FP64 pipe 89%** on the RTX 5070 — the kernel is
  *compute*-bound on double-precision division, not bandwidth-bound (GeForce FP64 = 1/64
  rate vs V100 1/2). Reformulating `cv/(1-dt·siu/c) == cv·c/(c-dt·siu)` halves the divisions
  (6→3) → **2.85× faster** on the 5070 (123.7→43.4 µs/call), and helps the CPU unified path
  too (shared kernel). See the "GPU hardware comparison" and "FP64 is the hidden bottleneck"
  sections above.
- **BUG (#200, fixed `b3a7da9a`): mode-2 loses startup forcing on deeply-dry cells.**
  This is the one the notebook's t≈120 divergence was really showing. Towradgi sets initial
  `stage=0`, but the creek-inlet bed is 215–320 m, so those cells start *deeply* dry
  (stage ≪ bed). Mode 1 reconciles them to `stage=bed` in one step (its per-step protect);
  the mode-2 device path converges only gradually — **halving the deficit each step, ~12
  steps** — and the inlet's inflow during that window raises the sub-bed stage instead of
  making depth, and is **permanently lost** (~24 m³ on Towradgi). Fix: clamp `stage` up to
  `bed` for dry cells in `set_gpu_interface()` before the initial device sync (mass-neutral,
  no-op for wet cells; mode 1 never calls it so is untouched). After the fix,
  mode-1-vs-mode-2 max|dStage| at t=120 drops **2.8 cm → 7.6e-6** — the difference is finally
  the roundoff-seeded chaos the notebook always *claimed* it was.
  - **The chase, 4 wrong layers before the truth** (each ruled out by measurement, worth
    reusing): "roundoff chaos" → no, systematic & localized; "inlet under-injects" → no, its
    accounting says applied=Q·t; "whole device state frozen" → no, only inlet-cell *depth*;
    "deep-dry alone" → no, needs deep-dry **at a bed gradient** on the real unstructured
    mesh. Instrument the *device* (sync_from_device), not the host — the host-side protect
    reconciles the host copy regardless, so a host read hides the bug.
  - **Fine output cracked it open.** The default 120 s yieldstep made it look like a t=120
    event; `-ys 0.1` showed the freeze is the **first ~1 s / ~12 internal steps** and the
    coarse output was hiding it. Lesson: when a divergence sits suspiciously on an output
    time, re-run with fine output before theorising.
  - **Regression guard:** `Test_GPU_DryCellStartupReconciliation` checks the *device* stage
    is reconciled after interface build (fails with fix reverted, GPU build only — CPU builds
    share host/device arrays so cannot express it). The full dynamic mass loss is
    **mesh-dependent** (a regular grid with the same bed does not reproduce it — only the real
    unstructured mesh), so it is validated on Towradgi, not re-created as a unit test.
- **Ruled out an operator bug at t=360** (user asked): toggled rain and culverts off — the
  mode-1/mode-2 divergence is *bit-identical* with/without them, distributed domain-wide,
  smooth (no step at t=360), sub-mm. It is chaotic amplification; the rain ramp only makes
  chaos *visible earlier* (lift-off t~90 with rain, t~300 without), not a systematic jump.
- **RTX 5070 vs V100 per-kernel** (nsys, both OpenMP-target): V100 only **1.16× faster in
  total kernel time** (5070 = 86% of V100), because the dominant `fluxes_central` is
  compute-bound and Blackwell wins it, while only the memory-bound kernels show the ~2× DRAM
  ratio. Recorded in the "GPU hardware comparison" benchmark section; corrected the earlier
  spec-sheet "50%" claim — measure kernels, don't infer from spec or wall-clock.

**Session 49 (2026-07-13/14):** **Operator parity audit, mode-1 vs mode-2 vs parallel** —
which found three real bugs, two of them in the *legacy* path, not the GPU one. Started
from a question about the CPU/GPU divergence in the PR-188 validation: the differences grow
sharply around t=360 s, so was an operator implemented differently on the two paths?
- **The t=360 divergence was NOT a bug** — worth recording, because the negative result is
  the point. The rainfall time series **steps up 3.5x at exactly t=360** (0.000375 ->
  0.001314 m/s); more rain ⇒ more newly wetted cells ⇒ more wet/dry fronts ⇒ faster chaotic
  amplification. Both paths see the identical rate. The clincher: the growth table showing
  the t=360 jump was **OpenACC vs OpenMP-target — both mode 2, sharing identical operator
  code**. Two runs with the *same* operators still show the jump, so it cannot be a CPU/GPU
  operator mismatch. **Method:** when a divergence correlates with a time, check the
  *forcing* before you go hunting in the code.
- **Operator physics does match.** `Rate_operator`: both paths evaluate `t`, `timestep`,
  `rate` and `factor` identically and apply the same update. `Inlet_operator`: both use the
  time-averaged `Q = 0.5*(Q(t) + Q(t+dt))`. Culverts: harmonised back in session 48.
- **BUG (#191, fixed `52705240`): mode-2 `Rate_operator` counted ghost cells in its mass
  tracking.** The CPU masks its influx sum with `full_indices`; the GPU kernel summed over
  every index. Under MPI a rainfall polygon straddling a partition boundary appears on
  several ranks, so the reported influx (and `fractional_step_volume_integral`) was inflated
  by the halo: **+4.8% at np=2, +9.1% at np=4** — which is exactly the discrepancy Steve had
  seen on gadi with 4 GPUs/4 ranks. The insidious part: `full_indices`/`num_full` were
  plumbed from Python into the C struct, malloc'd, memcpy'd, freed — and **never read by any
  kernel**. The code *looked* like it handled ghosts. Fix bakes the mask into the area array
  at init (`areas` -> `mass_areas`, zero for unowned cells; free, since `areas` was used for
  nothing else) and **deletes the dead fields** so the trap cannot reappear.
- **BUG (#193, fixed `c843ee78`): parallel inlet over-counted the mass balance by
  x(number of ranks).** `fractional_step_volume_integral` is a per-rank LOCAL accumulator
  (allreduce-summed on read), but the parallel inlet added the master's GLOBAL volume on
  *every* participating rank: **np=1 -> 40.0, np=2 -> 80.0, np=4 -> 160.0** against a true
  `Q*dt = 40.0`. **Not a GPU bug** — mode 1 and mode 2 gave identical wrong numbers, so this
  had been corrupting `Water_volume_statistics` in *legacy* production MPI runs. Root defect
  was a **convention clash**: `Rate_operator` contributes its local share, the inlet
  contributed the global total, into the same accumulator. Fixed by routing all nine inlet
  sites through `_add_fractional_step_volume()` (master-only) and documenting the
  convention. Audited the other contributors: parallel **structure operators never touch the
  accumulator at all** (internal transfers, mass-neutral), so they are clean.
- **BUG (#194, PR #195): the mode-2 banner printed the RANK count labelled "GPU(s)".** A
  4-rank run on a 1-GPU box reported "4 GPU(s)". It concealed the one thing the banner was
  best placed to catch — ranks map to devices round-robin (`rank % num_devices`), so
  oversubscription silently shares a device.
- **Can you run several ranks on one GPU? Measured, and the answer is useful** (160k-tri
  mode-2 evolve, one RTX 5070):

    | ranks | no MPS | with MPS |
    |---|---|---|
    | 1 | 3.80 s | 3.51 s |
    | 2 | 7.13 s | **3.48 s** |
    | 4 | 11.21 s | **3.72 s** |

  Without MPS the ranks **time-slice** the device: ~3x slower at np=4. **NVIDIA MPS erases
  the penalty but never beats one rank per GPU** — a single rank with the whole mesh already
  saturates the device, so splitting it adds no parallelism; MPS just lets the fragments
  overlap again. *MPS is a way to stop losing, not a way to go faster.* Note it did **not**
  hang and every rank count gave a **bit-identical checksum** — so the known hang/garble
  failure is real (Steve has hit it on gadi) but **not universal**, which is why the banner
  now leads with the measured slowdown and reports hangs as a possibility. A warning that
  predicts a crash that does not come is one people learn to ignore — which is how the
  original banner failed.
- **BUG (#192, fixed — PR to `develop`): mode 2 IGNORED fractional-step operator
  registration order.** Mode 1 runs operators in registration order; mode 2 hoisted every
  Boyd culvert to the **front** (batched via `GPUCulvertManager`) before the loop. Mode 2
  therefore returned **bit-identical** results whether a `Rate_operator` was registered
  before or after a culvert, while mode 1 gave different answers:

    | order sensitivity, max abs stage diff | mode 1 | mode 2 |
    |---|---|---|
    | pre-fix | responds (9.75e-05) | **0.0 — ignores order** |
    | post-fix | responds (9.75e-05) | responds (9.78e-05) |

  Towradgi registers its 22 culverts first, so the orders coincided *by luck*.
  - **Steve asked the right design question: should the CPU adopt culverts-first instead?**
    **No** — and the investigation shows why it would not even have worked. The batched GPU
    path is **gather-all -> compute-all -> scatter-all**, so there are *two* divergences:
    the hoisting, AND mode 2 applying culverts **simultaneously from a snapshot** where mode
    1 applies them **sequentially** (culvert 2 sees culvert 1's stage change). Reordering the
    CPU fixes the first and leaves the second — you would move every calibrated production
    model's answers and *still* not have parity. Also: **mode 1 is the oracle** we validate
    GPU work against; changing the oracle to match the thing under test is backwards.
  - **Fix:** fire the batched cycle **at the position of the first culvert** rather than
    before the loop. Contiguous culverts (what every real script does) then reproduce
    registration order *exactly*, batching untouched, **zero performance cost**. Interleaved
    culverts still get pulled forward — so it now **warns** instead of diverging silently.
  - **The measurement trap — READ THIS BEFORE TOUCHING THE GUARD.** My first check said the
    fix made things *worse* (656 -> 1076 cells differing). It had not. Comparing
    mode-1-vs-mode-2 *divergence* is the wrong instrument: a separate, pre-existing
    mode-1/mode-2 culvert discrepancy (~1e-4) dominates, and in this setup it is nearly
    **equal and opposite** to the order effect — so the *wrong* order coincidentally looked
    closer (1.4e-06) than the right one (9.77e-05, the floor). The correct instrument is
    order **sensitivity** (does swapping the registration order change the answer?), which is
    unambiguous: exactly 0.0 with the bug. The test docstring carries this warning so nobody
    later "improves" the guard back into the misleading form.
  - **Still open, separate:** the snapshot-vs-sequential culvert difference above, and the
    ~1e-4 mode-1/mode-2 Boyd_box discrepancy that forms the floor. Neither is caused by #192
    or by its fix (the culvert-first control is byte-for-byte unchanged across it).
  - **NEAR-MISS + a real lint gap.** Inserting `_warn_if_culverts_interleaved()` after the
    operator loop silently swallowed the trailing
    `if gpu_mode and needs_cpu_sync: sync_to_device()` into the *new method*, where those
    locals do not exist. `apply_fractional_steps()` then stopped pushing host work back to
    the device, so CPU-only operators (wind stress) never reached the GPU — the exact bug
    class fixed earlier this session. **The unified suite caught it (3 wind-stress
    failures); `ruff check` did NOT.** The project's ruff config does not enable **F821
    (undefined name)** — `ruff check` reports "All checks passed" on code with two undefined
    names, while `ruff check --select F821` flags them. **112 F821s exist repo-wide**, so
    enabling it is its own piece of work, but *lint here cannot catch a typo'd or orphaned
    variable reference.* Do not rely on it to; run the suite.
- **Testing lessons worth keeping.**
  - ⚠️ **On a GPU-offload build the fast suite SKIPS `test_DE_gpu_omp.py` ENTIRELY — all 86
    tests.** This bit hard: the #192 guard was written, "verified" with a green fast suite on
    the nvc build, pushed — and CI went red on 8 jobs, because the fast suite had never
    *run* the tests being added. Counts make it obvious once you know:

        GPU-offload build (nvc):  2610 pass / 218 skip     <- test_DE_gpu_omp skipped
        CPU-only build (gcc):     2696 pass / 219 skip     <- the extra 86 ARE those tests

    **Any change to `test_DE_gpu_omp.py` must be validated on a CPU-only build**
    (`pip install --no-build-isolation -e . -Csetup-args=-Dgpu_offload=false`), which is
    what CI builds. The isolated runner (`anuga_run_isolated_tests`) *does* exercise the file
    on a GPU build, but the *fast suite alone is not evidence* for it.
  - **Do not assert on a mode-1-vs-mode-2 divergence MAGNITUDE — it is build-dependent.**
    The first version of the #192 guard required mode-2's order-sensitivity to match mode-1's
    within 5%. It passed on the GPU build (ratio 1.00) and failed on CI's CPU build (ratio
    2.0), because the magnitude is confounded by the pre-existing mode-1/mode-2 culvert
    discharge discrepancy, which differs per build. Assert on the *property* (mode 2 obeys
    registration order — exactly 0.0 sensitivity when it does not), not on a number that a
    separate discrepancy contaminates. **This is the same trap already documented two bullets
    up, and I walked straight into it anyway** — which is precisely why it is written down
    twice.
  - A parallel bug can often be made **serially catchable**: #191's guard fakes the partition
    by clearing `tri_full_flag` (the only thing `set_full_indices()` reads) and asserts
    against the analytic influx — so it runs in the ordinary non-MPI suite. #193's could
    *not* be (at np=1 there is only the master, so the bug is structurally invisible), hence
    a real `mpicmd` np=3 test.
  - **Extract the untestable into a pure function.** `gpu_startup_banner()` only ever fires
    on hardware CI lacks, so it was pulled out as a pure function of
    `(numprocs, num_devices, device_id, offload_active)`; the whole matrix is now covered
    serially.
  - **`ruff check` does not catch undefined names here (issue #197).** F821 is not enabled,
    so an orphaned/typo'd variable reference lints clean — it took the test suite to catch
    one. 112 F821s exist repo-wide, 74 of them in the single broken module
    `anuga/tsunami_source/eqf_v2.py`. Do not trust lint for this class of error.
  - Every fix this session was **proven by re-injecting the bug** and watching the new test
    fail. Do that; a guard you have not seen fail is not a guard.

**Session 48 (2026-07-11/12):** Culvert/weir mode-1 vs mode-2 reporting parity,
two real physics bugs, the **3.3.8 patch release**, then **partition-save
performance** (parallel dump, single-file layout) and docs. A user reported the parallel
(MPI) culvert/trap logs differed from sequential. Root-caused and fixed on
`develop` (branch `fix/culvert-parallel-mode2-reporting`, PRs pushed direct to
both develops, four commits `d55e1984`→`467255f0`):
- **Parallel structure logging** — parallel `.log` files were missing the
  `accumulated_flow` column (computed on master, never written) and mis-numbered
  operators (counter advanced only on the master proc, so a culvert and trap with
  different masters collided). Fixed: write the column + header, advance the
  counter on every rank in lockstep (participants in the ctor, non-participants
  in the factory `return None` branch), and drop the `_P<rank>` filename suffix so
  names match sequential. Verified np=4 and np=16 bit-identical to serial.
- **Mode-2 (GPU) stats write-back** — in mode 2 the batched `GPUCulvertManager`
  runs the Boyd/weir physics in C and skips the Python `__call__`, so the logs
  were **all zeros**. Added per-culvert report fields to `culvert_state`
  (C `gpu_culvert_operator.c`), a `gpu_culverts_get_report` getter, and a
  write-back in `apply_all` mirroring mode-1's accumulation.
- **boyd_pipe critical-depth bug** (real): the GPU C translation **divided** by
  `(bf·d)^2.5` where the Python reference **multiplies** → flow_area/velocity ~8%
  off (discharge matched only because it was inlet-controlled; outlet-controlled
  flow would corrupt discharge too). Fixed both dcrit pairs; added a boyd-pipe
  cross-mode regression test (the existing one only used boyd_box).
- **weir critical-depth gravity** (real, issue #181, closed): the Python
  `weir_orifice_trapezoid_function` hardcoded **9.81** in the trapezoid
  critical-depth Newton solve, vs the C using `domain.g` (9.8) → ~0.034% velocity
  offset in partly-full flow. Fixed to derive `g` from `domain.g` on both sides
  (Python grows a `g` arg, discharge routines pass `self.domain.g`); note this
  shifts mode-1 weir results ~0.034%. Guard test added.
- **3.3.8 release.** Cherry-picked the two **legacy** (non-GPU) fixes — parallel
  logging + weir gravity — to `main` via **PR #182** (GPU-only pipe/mode-2 fixes
  omitted; that machinery isn't on 3.3.x). Also in 3.3.8 since 3.3.7: SWW
  large-t store crash (#149), osx-64 wheels (#142), checkout 6→7 (#146). Tagged
  **annotated `3.3.8`** on the merge commit (bare version, no `v`, per the
  convention), published the GitHub Release → PyPI (21 files, trusted OIDC);
  conda-forge follows via feedstock bot.
- **Gotcha (memory):** MPI parallel + mode-2 needs `#ranks == #GPUs` or it hangs/
  garbles — validate GPU mode-2 in serial on a 1-GPU box.

Then, **partition-save performance** (a 173M-triangle / 18,400-partition run
measured ~1,000 s load + ~4,000 s partition + **~37,000 s to save the
partitions** — ~88% of the job; the dump routines wrote files in a serial loop):
- **Four stacking speedups (PR #183).**
  - **Opt-in parallel dump** (`num_workers`) for `sequential_distribute_dump`
    (domain) and `sequential_mesh_dump` (mesh). `extract_submesh` only *reads*
    the shared submesh, so a **fork** `ProcessPoolExecutor` shares the
    partitioned mesh **copy-on-write** (never pickled), and `gc.freeze()` keeps
    the big inherited structure out of the workers' GC scans (avoids both the
    rescan cost and COW page-dirtying). Default 1 = today's serial, memory-frugal
    path (it releases each rank as it writes; the parallel path keeps the whole
    submesh live — the memory-for-speed trade-off).
  - **Single-file domain dump** (`single_file=True`, new default): inline
    points/triangles/quantities into the one per-partition pickle instead of
    separate `.npy` files → **8 files/partition → 1** (~147k → 18.4k files at
    18,400 ranks), the dominant cost on metadata-bound Lustre/GPFS. The loader
    reads **both** layouts (filename-string vs inline array), so old dumps still
    load; `single_file=False` restores the legacy layout.
  - **Batched serial gc** — was collecting every rank (18,400 full graph
    traversals); now every 64 (big arrays are freed by refcounting anyway).
  - **Vectorized mesh `boundary_tag` write** — `_write_mesh_partition` wrote the
    NetCDF string variable **one row per boundary edge** (one HDF5 write each).
    cProfile showed this loop dominating; assembling the char array and writing
    it once took the mesh write **9.45 → 2.57 ms/file (~3.7×)**, byte-identical.
- **Partition file format — pickle vs NetCDF (decision).** Benchmarked on real
  partitions. **Domain → keep pickle**: writes ~6–8× faster than NetCDF at ~equal
  size, and it carries arbitrary Python state (boundary map, geo_reference) that
  NetCDF cannot hold without a blob. **Mesh → keep NetCDF**: its content is 100%
  arrays, it is **1.6× smaller** than pickle *and* portable/inspectable
  (`ncdump -h`), safe to load, and write-once/read-many. Prototyped a **hybrid**
  (bulk arrays as NetCDF variables + a small pickled config blob, one file) — it
  round-trips exactly but writes ~6× slower, so **not adopted**. Two formats for
  two genuinely different jobs; the mesh's apparent 34× gap was mostly the
  boundary_tag bug above, not NetCDF itself.
- **Docs (PRs #184, #185).** Documented `num_workers` / `single_file` and added
  "Performance for large partition counts" sections. Then **harmonized** the
  domain/mesh offline-partitioning pages: `sequential_distribute_dump` had a
  **stub docstring** (no `Parameters` block) and `sequential_distribute_load` had
  **none**, so the domain page's autodoc API rendered sparse next to the mesh
  page — both now have full numpydoc + `See Also`, and each dump states its
  on-disk format (domain → pickle, mesh → NetCDF). Also reordered the
  **Appendices** (builds → advanced usage → GPU internals → theory → Contributing)
  and shortened two long nav labels.
  **Gotcha:** #183's final commit landed on the branch *after* the merge captured
  it, so it silently missed `develop` — recovered via #184. Check the merge picked
  up the head you expect.
- **`pymetis` → core dependency** (was a `parallel` extra). `partition_mesh` needs
  it, and the *sequential* offline-partitioning preprocessing uses that without
  mpi4py — so gating it behind `parallel` meant a plain install could not
  partition a mesh.
- **GPU install docs.** Stated the real requirement is **a compiler with working
  OpenMP offloading** — the kernels are standard OpenMP `target`, not CUDA, so
  nvc is simply the best option *at the moment* (GCC's nvptx backend ICEs on the
  ANUGA kernels; LLVM/Clang offload, AMD AOMP, Intel icx untested but feasible).
- **`tools/install_anuga_nvc.sh` preflight (user report).** A user ran the script
  in a conda env **not** created from `environments/environment_<PY>.yml` and it
  died with `BackendUnavailable: Cannot import 'mesonpy'` — *after* nvc was found,
  so it read like a GPU/compiler failure when it was just a missing build dep.
  Cause: the build uses `pip install --no-build-isolation`, so pip does **not**
  create a temp build env — the **meson-python** backend and the rest of
  pyproject's `build-system.requires` must already be in the env. Added a
  **preflight** that checks `mesonpy`/`cython`/`pybind11`/`numpy` are importable
  and `meson`/`ninja` are on PATH, and on failure names the missing packages and
  gives the fix (`conda install -c conda-forge …`, or recreate from the matching
  `environments/environment_<ver>.yml`). It also prints the env's **actual**
  Python version — the script prefers an *already-activated* env and ignores `$PY`,
  so the banner could say `PY=3.14` while building into a 3.12 env (confusing in
  the user's log).
  **Triage tip:** `Cannot import 'mesonpy'` from this script is never an nvc
  problem — it is a missing build backend in the target env.
- **Test isolation: `Test_File_Conversion` (issue #186, closed).**
  `test_ferret2sww3` passed in isolation but ERRORed in a full run. `setUp` and
  several tests in the class create files with **bare relative names**
  (`most_small_*.nc`, `test_*.nc`, `test.sww`, `*.tms`, `*.asc`/`.dem`) in the
  **CWD**, and clean them up *inline at the end of each test body* — so a test
  that fails or is interrupted leaves stale files behind and the **next** run dies
  overwriting them with `NetCDFFile(..., 'w')` (triggered by a stale `tmp/` left by
  a previously failed `install_anuga.py`). Fixed by giving each test its **own
  temp working directory** (`setUp`: `mkdtemp` + `chdir`; `tearDown`: `chdir` back
  + `rmtree`), so cleanup is **unconditional** even on failure. Chose the temp-CWD
  approach over patching ferret3 alone because it fixes *every* CWD-writer in the
  class at once and stops the suite polluting the directory it runs from; safe
  because the class reads no fixture files by relative path (it creates everything
  it reads). Verified it still passes with **corrupt + read-only** leftovers
  planted in the CWD (the exact failure mode), leaving them untouched.
  **Pattern for similar failures:** "passes alone, fails in a full run" + files
  written to the CWD ⇒ isolate the whole TestCase in a temp CWD rather than
  chasing individual `os.remove` calls. (An earlier fix, `55d91479`, had used
  per-file `mkstemp` for `test_sww_extent` but never touched the ferret tests.)
- **`set_quantity()` did not reach the device in mode 2 — a REAL silent
  correctness bug** (`d22bb53f`). Triaging unified-mode test failures reported on
  PR #187 (which we first suspected were just the known unified-in-one-process
  aborts — they were **not**: no aborts, and they reproduced *in isolation*).
  `test_sww2dem_verbose_True` stored an SWW whose stage was `[-0.6, -0.1]` after
  being set to `1.0`. Root cause: once the GPU interface exists the **device**
  holds the authoritative centroid state, but `Domain.set_quantity()` updated only
  the **host** arrays — so the device kept evolving stale/default values. With the
  initial condition removed, mode 2 evolved **stage = 0** (the device default),
  proving `set_quantity` never reached it. **Reachable from ordinary scripts**,
  because `set_boundary()` *and* `distribute_to_vertices_and_edges()` both call
  `_ensure_gpu_interface()`:

      domain.set_boundary(...)           # builds the device interface
      domain.set_quantity('stage', ...)  # host only -> device never sees it
      domain.evolve(...)                 # evolves stale values -> WRONG output

  i.e. GPU runs could be silently wrong, no error or warning. Session 39 had fixed
  exactly this for `set_boundary()`; `set_quantity()` was missed. Fix mirrors it:
  sync host **from** device, apply the change, sync **to** device (no-op in legacy;
  off the hot path — no operator/structure calls `set_quantity()` per timestep).
  Guard: `Test_GPU_SetQuantityReachesDevice` (fails with the fix reverted).
  **Lesson:** when a host array and the device can both hold state, *every* public
  mutator needs the sync, not just the one that bit you. Worth auditing other
  host-mutating APIs (e.g. direct `quantity.set_values()`) for the same gap.
- **`test_urs2sts` unified failures were NOT a bug** (`b20bfca0`). Three tests
  assert on `quantity.boundary_values`, a *host* array mode 2 never populates (it
  evaluates boundaries on-device and only syncs centroids back). Verified the
  **physics is mode-2-correct**: skipping only those host-array asserts, the tests'
  own fbound-vs-Dirichlet comparisons all pass under unified. So they are white-box
  mode-1-only-state checks ⇒ pinned to legacy per `CONVENTIONS.md`. `test_2pts`
  already passed under unified.
  **Triage rule that paid off:** before "fixing" a mode-2 test failure by pinning it
  to legacy, check whether the *physics* still agrees — if it doesn't, you are
  hiding a real bug. PR #187 proposed pinning these tests, which would have masked
  the `set_quantity` corruption above.
- **Audit of the other host-mutating paths, and the choke-point fix** (`c6a6f503`).
  Followed the lesson above. `set_quantity_vertices_dict()` routes through
  `set_quantity()` ⇒ covered; `set_boundary()` ⇒ session 39; fractional-step
  operators ⇒ covered (non-GPU-safe ops are classified CPU-only and
  `apply_fractional_steps()` brackets them with a sync pair); no direct
  `centroid_values[:] =` on conserved quantities exists; checkpoint restore drops
  the interface and rebuilds lazily. **Two residuals**, both actioned:
  - Direct `Quantity.set_values()` still bypassed the sync. Fixed by moving the
    round-trip **down from `Domain.set_quantity()` into `Quantity.set_values()`** —
    the single choke point for host-side quantity writes, so both routes are
    covered by one mechanism and neither double-syncs. Two guards keep it off the
    hot paths: `GPU_SYNCED_QUANTITIES` limits the sync to the four centroid arrays
    the C actually round-trips (`stage`, `xmomentum`, `ymomentum`, `height` — see
    `gpu_domain_sync_to_device()`), so elevation/friction/user tracers pay nothing;
    and `apply_fractional_steps()` sets `_gpu_host_writes_suppressed` inside its
    already-bracketed region, without which an operator that sets a quantity would
    cost a full host↔device transfer pair *per operator per timestep*.
    **Rejected:** the deferred-push ("dirty flag, flush at next step") variant — it
    races with the yieldstep's `sync_from_device()`, which would silently drop the
    pending host write and reintroduce this exact bug class.
  - **Issue #189** — `apply_protection_against_isolated_degenerate_timesteps()`
    never runs in mode 2. It hangs off `update_timestep()` (line 2718), but all
    three mode-2 step functions return early into their C counterparts *before*
    reaching it. Default-off (`config.py:151`) so nothing is silently wrong today;
    the sharp edge is that enabling it under GPU gives no protection *and no
    warning*. Suggested fix: warn, as mode 2 already does for unsupported forcing.
  - **Verification gotcha worth remembering:** these sync guards only mean anything
    on a **GPU-offload build** (`gpu_offload=True`, nvc). On a CPU build the
    `omp target update` pragmas are no-ops and the host arrays *are* the "device"
    arrays, so the tests pass whether or not the fix exists. Always confirm the
    build (`build/*/meson-info/intro-buildoptions.json`) before trusting a green
    run as proof. With the sync disabled on the GPU build, the direct-Quantity test
    and the `set_quantity` guard both fail while the "interface not yet built"
    control still passes — the right signature.
- **PR #188 (JorgeG94) — OpenACC offload back end: reviewed and numerically
  validated.** Adds `gpu_backend=openmp|openacc` (meson), routing every kernel
  through `gpu/gpu_omp_macros.h` so the same source compiles for OpenMP-target,
  OpenACC (`nvc -acc=gpu`) or CPU multicore. Verified the load-bearing claim rather
  than trusting it: **no raw `#pragma omp target`/`omp parallel` and no `reduction(`
  survives** in `gpu/*.c` — important, because under `-acc=gpu -mp=multicore` a
  leftover `omp target` would silently run on the **host** while the data sat on the
  GPU. Every host↔device transfer also genuinely routes through a draining macro.
  - **Numerical validation (the gate the PR left open).** Ran `run_small_towradgi.py`
    (`-mpm 2 -go`, 256k tri, ft=3600) on both back ends **built from the same branch**,
    so the back end was the only variable. Result: **OpenACC validates.**

    | comparison | max ΔStage | RMS | >1mm |
    |---|---|---|---|
    | OpenACC run1 vs run2 (determinism) | **0.000e+00** | 0 | 0 |
    | OpenACC vs OpenMP-target | 7.29e-02 | 1.26e-04 | 500 |
    | OpenMP-target vs legacy CPU (*already shipped*) | 8.15e-02 | 3.70e-04 | 1650 |

    **The yardstick row is the third** — that is the divergence we already ship and
    accept between the current GPU back end and the CPU solver. ACC-vs-OMP is *smaller*
    on every measure. Divergence starts at **exactly zero at t=0** (so the
    `acc enter data copyin` mapping is right), is **5.96e-08 = 2⁻²⁴ (one ULP of the
    `.sww`'s float32 storage) at t=120**, and only then amplifies — textbook chaotic
    sensitivity in a wetting/drying model seeded by a different FP reduction order.
    **Reusable method:** a *bitwise* run-to-run repeat is what separates "chaotic but
    correct" from "data race" — both produce the same growth curve, but only a race
    breaks reproducibility. Worth reaching for on any future GPU-backend change.
  - **Perf: the kernel win is real but does not reach wall clock.** Steve spotted that
    the PR's "fluxes-central 2.0s → 1.62s" had to be a short run; re-measured with
    **nsys at `-ft 240`** and it reproduces — 2.304 s → 1.843 s (**−20%**, same 4736
    instances); total GPU kernel time 5.40 → 4.92 s. Yet wall clock is a dead heat at
    ft=240 (17.59 s vs 17.58 s) and ~1% *slower* at ft=3600 (206.2 s vs 208.3 s). The
    CUDA API breakdown locates it: **the async queue issues 2.5× MORE
    `cuStreamSynchronize` than OpenMP-target** (554,260 calls / 5.46 s vs 225,038 /
    3.00 s) — the opposite of its stated motivation — masked only because OpenACC also
    does far less D2H (3.61 s → 0.93 s), leaving total host API time a wash (7.13 vs
    7.15 s). That wash is what eats the 20%. **Headroom, not a defect.**
    **Lesson:** "kernel is faster" and "run is faster" are different claims; always ask
    which was measured before agreeing *or* disagreeing.
  - **Real bug found: `set_gpu_offload(False)` / `-ngo` silently no-ops on OpenACC.**
    The host-fallback routes target regions to the initial device
    (`omp_set_default_device(omp_get_initial_device())`); the ACC shim turns that into
    `acc_set_device_num(-1, ...)`, and a negative devicenum in OpenACC is not "run on
    the host" — the runtime just reverts to implementation-defined default. OpenACC has
    no equivalent of OpenMP's initial-device semantics. So the run stays on the GPU with
    no error. **Issue #190**, scoped as a merge blocker for #188 rather than a live bug
    (the back end is only on the feature branch) so it survives a rebase/split of that PR.
    Fix is likely to *hard-error* rather than implement: the OpenMP idiom has no OpenACC
    counterpart at all (host execution is a device *type*, not a device *number*), and the
    separate `g_gpu_offload_enabled` flag may already carry most of the weight.
    **Note the validation above does not cover this path** — both runs used the GPU, so a
    green `.sww` diff gives false comfort here.
  - Also flagged: the PR silently changes the **default** build's flags
    (`-mp=gpu,multicore` → `-mp=gpu`), which bears on the known nvc host-fallback
    pathology; and `acc_free` at `gpu_halo.c:121` frees device buffers without draining
    the queue, inconsistent with the header's own `exit data` reasoning.

**Session 47 (2026-07-07/08):** Documentation overhaul — restructure, API
cross-linking, meta-pages, and a warning-free Read the Docs build.
- **Issue #32 "Make riverwalls transmissive."** Verified the `Cd_through`
  submerged-orifice throughflow is already active in **both** compute modes:
  legacy mode-1's `_openmp_compute_fluxes_central` delegates to the shared
  `core_compute_fluxes_central` in `core_kernels.c` (no mode gating), so the
  earlier "not in legacy" reading was wrong. Confirmed empirically (mode-1 and
  mode-2 give bit-identical throughflow). Documented on the issue and here.
- **Docs overhaul (PR #157).** Split into a standard-user **Contents** vs
  advanced **Appendices** structure; added a landing quick-start, a
  **Conventions & units** primer, and an evolve **Stability/blow-ups** section;
  moved Parallelisation into Contents and TOML/ANUGA-Viewer/QGIS into the
  standard sections; led the appendices with the developer + new **GPU install**
  pages. Converted the narrative "Reference" blocks to compact autosummary tables
  linking into the **API Reference** (expanded API consolidated there), and made
  every class page's method summary link to per-method signatures
  (`autodoc_default_options={'members':True}` + `sphinx.ext.napoleon`). Reframed
  the **ANUGA Viewer** as the recommended fast viewer for large `.sww` (dropped
  "legacy"). Content-review fixes (typos, heading levels, missing pip path,
  smoke-test, quantity units). New **Citing**, **Contributing**, **Glossary**
  meta-pages (surfacing `CITATION.cff`, `CONTRIBUTING.rst`, Apache-2.0).
- **Docstring + tooling fixes.** Reformatted malformed-RST docstrings surfaced by
  autodoc `members` (`Quantity`/`Domain` ×6 in #158, `internal_boundary_operator`
  ×2 in #160). `install_anuga_nvc.sh` now builds into an already-activated conda
  env if present (#158).
- **Warning-free build (PRs #159, #161).** Fixed the local warnings —
  `html_static_path` (`_static` dir) and the `Geo_reference.epsg` duplicate
  (napoleon renders a class-docstring *Attribute* **and** the real property →
  `napoleon_use_ivar = True`). Then **watched Read the Docs** (it builds
  `develop`): it had surfaced 62 warnings the local build hid — 56 `ipython3`
  Pygments-lexer (RTD lacks IPython → add `ipython` to `docs/requirements.txt`)
  plus the operator docstrings. RTD `develop` now builds **clean** (only the
  harmless MPI-less `Could not import mpi4py`). Added a Contributing "Building
  the documentation" note on reproducing RTD locally (#162). **Lesson:** a
  locally-installed IPython and `-D nbsphinx_execute=never` hid warnings RTD
  shows — build docs from a clean `docs/requirements.txt` env.
- **Branch policy (PR #163 + memory).** Recorded in `ROADMAP.md`: **do not merge
  `develop` → `main` until the team cuts v4.0.0**; all work/syncs stay on
  `develop`.
- **PDF structure + Conventions polish.** Noted the Unix-epoch time origin
  (seconds from 00:00:00 UTC 1 Jan 1970; `set_starttime` via `datetime`/
  `zoneinfo`) in the Conventions page (#166). Fixed the PDF putting every
  section under a single "Quick start" chapter — the `.. toctree::` sat *inside*
  the Quick-start section on the landing page, so moved Quick start to its own
  `quickstart.rst` (first in the toctree); verified against the live RTD PDF that
  Quick start and the other sections are now sibling chapters (#167).
- **Versioning diagnosis + wheel-repair cherry-pick (2026-07-08).**
  `git describe` on `develop` showed `1.3.1-4093-g…`, but the package version is
  correct — `anuga.__version__ = 3.3.6.dev767+g…` (`_git_version.py` uses
  `git describe --tags`). Two causes: the `3.x` tags are **lightweight**, so a
  plain `git describe` skips them and falls back to the old **annotated** `1.3.1`
  (`git describe --tags` → `3.3.6-…`); and `3.3.7` is a release tagged on `main`
  **not reachable from `develop`**. That 3.3.7 release had exactly one commit
  `develop` lacked — `e4700a7b` *"Repair wheels on all platforms so they are
  self-contained"* (delocate/delvewheel/auditwheel so mac/win/linux wheels are
  self-contained) — **cherry-picked into `develop`** (`014451da`; pushed **direct
  to `anuga-community/develop`** as an authorised one-off, no PR). Recorded a
  **tagging convention** in `ROADMAP.md` (PR #168): future release tags must be
  **annotated**, bare version, no `v` prefix (`git tag -a 3.3.8 -m "ANUGA 3.3.8"`)
  so plain `git describe` is correct.
- **nvc GPU rebuild + GPU sanity check.** Rebuilt via `tools/install_anuga_nvc.sh`
  (active-env path from #158) → `3.3.6.dev772+g08d4c10b`, 65/65 isolated GPU
  tests, and a dam-break sanity sim (GPU offload engaged, mass conserved to ~1e-13,
  GPU vs CPU-legacy ~4e-7).
- **GPU `Time_boundary` mode-1/mode-2 divergence — root-caused and fixed
  (PR #171, issue #170).** A rising-tide flood diverged ~4e-3 between mode-1 and
  mode-2. Isolated it: **single-substep** algorithms (DE0, DE_ader2) agree to
  machine precision, **multi-substep** (DE1/DE2) diverge — and only with
  **Python-evaluated boundaries** (Time/File/wave/Flather); reflective/steady are
  exact. Cause: the single-call **C RK loop** (`_evolve_one_rk*_step_c`) sets
  time-varying boundaries on the device **once per step**, reusing that value for
  every RK substep, whereas mode-1 calls `update_boundary()` **per substep**
  (an O(dt) boundary-forcing error). Fix (**option B**): route domains with any
  Python-evaluated boundary to the Python-orchestrated GPU loop (refreshes per
  substep → bit-matches mode-1; benchmarked GPU cost ≤~4%, within noise). Added
  `Test_GPU_TimeBoundarySubstep`; 69/69 GPU tests green. Filed **issue #170** for
  the proper per-substep C-loop fix (option A).
- **Fractional-step operator/structure timing (PRs #174, #175, #177; issue #176).**
  Checked operators/structures for the yieldstep-vs-inner-step concern — they are
  fine (applied **every inner timestep** in both modes: 13/13 and 46/46 evals ≫
  yieldsteps). But found a related **operator-evaluation-time** bug: fractional-step
  operators (applied by the evolve loop *before* it advances time to t+dt) should
  see the pre-step time **t**, and DE0/DE_ader2 do — but **DE1 (rk2)** and **DE2
  (rk3)** left `relative_time` advanced, so operators evaluated forcing at t+dt
  ("one step too far"), ~4e-4 for a time-varying rate. For **DE1** only mode-1 was
  affected (diverged from mode-2), fixed by restoring the pre-step time in the
  mode-1 rk2 body (**PR #174**). For **DE2** *both* modes advanced (mode-1 body +
  mode-2 C **and** GPU loops), so it stayed self-consistent and the cross-mode
  check missed it — fixed all three rk3 paths (**PR #177**, closes #176). Now all
  four algorithms evaluate operators at t in both modes. Guards: `test_operator_timing.py`
  (CPU, mode-1 — runs on any build) + `Test_GPU_OperatorTimeAlignment` (cross-mode);
  the GPU-only guard was split out to a CPU test in **PR #175**. No prior test used
  a time-varying operator (coverage gap). CPU suite 2610 pass, GPU 73/73.
- All merged to `anuga-community/develop` (PRs **#157–#177**, admin-merged by
  number, except the one-off direct cherry-pick `014451da`); RTD `develop` HTML
  **and PDF** confirmed clean.

**Session 46 (2026-07-06):** Issue #33 (memory) documentation + measurement,
plus follow-ups.
- **Issue #33 "Reduce Memory usage".** Documented on the issue everything already
  implemented for v4.0.0: the **Quantity per-type allocation** (QM1–QM7 — each
  quantity allocates only the arrays its `qty_type` needs instead of the blanket
  9; lazy `vertex_values`; centroid-primary elevation; lazy gradients/`phi`;
  shared gradient workspace `22559a5b`; ~54–58% off quantity memory) and the
  **domain C work-array reduction** (DM1: 9 dead arrays removed + deferred to
  first evolve; DM2: riverwall arrays lazy; ~740 MB at 2.25M tris), plus the
  exported `memory_stats`/`quantity_memory_stats`/`domain_memory_stats`/
  `domain_struct_stats` instrumentation. Then **re-ran the issue's exact
  benchmark** (`run_parallel_rectangular.py`, `mpiexec -np 2`, proc-0 Max RSS):
  **710→511 MB (−28%), 2.5 GB→1.37 GB (−45%), 5.1 GB→2.74 GB (−46%)** — RSS
  roughly **halved at 2.25M triangles**, the saving growing with N. Recorded the
  numbers in `PROGRESS_ARCHIVE.md` (PR #155) and posted them to the issue. The
  `print_domain_memory_stats` breakdown shows `river wall 0.00 MB` (DM2 lazy) and
  the trimmed work arrays (DM1). Remaining lever: rank-0's peak building the full
  domain before `distribute`.
- **`tools/install_ubuntu.sh`.** Combined the ~90%-identical `install_ubuntu_2X_04.sh`
  scripts into one version-aware `install_ubuntu.sh` (22.04/24.04/26.04 `case`,
  auto-derived python version); fixed the earlier cleanup's broken 20.04 dispatch
  and README mislabel (PR #154).
- **PR triage.** Reviewed **#148** (multi-compiler flag tuning, GCC 15/NVHPC/ICX)
  — recommended gating the FP-semantics flags (`-ffinite-math-only`,
  `-fassociative-math`, `-Mfprelaxed`) behind an opt-in `-Dfast_math` option so
  the default `pip install` keeps strict IEEE, while keeping the safe flags +
  `pow→cbrt` rewrite as default; asked @samcom12 to make that change.
- **Synced `develop` → `anuga-community`** (PRs #153/#154/#155); upstream CI green.

**Session 45 (2026-07-02 – 07-05):** Geodata CI breakage (libjxl), NumPy-2.5
warnings, PR triage, and a `tools/` cleanup.
- **Geodata CI failure.** A red CI on PR #150 (SeanWong's timestep clamp) turned
  out to be unrelated: conda-forge omits `libjxl` on Linux/macOS for py3.11–3.14,
  so `libgdal-core 3.12.3` can't load (`libjxl.so.0.11: cannot open shared object
  file`) and the whole geodata stack (fiona/rasterio/shapely) is down. Fixes,
  merged to `develop`: (1) skip-guard the one un-guarded geodata test
  (`test_rain_with_polygon_csv`) — **PR #151**; (2) a `spatialInputUtil`
  diagnostic that surfaces the real import error instead of the blind
  `except ImportError` (this is what revealed `libjxl`); (3) pin `libjxl <0.12`
  in `environment_3.10..3.14.yml` to force the `0.11` soname — **PR #152**.
  Restored ~20 geodata tests (upstream CI: 2667 passed, **0** `requires rasterio`
  skips). Confirmed transient-vs-persistent by re-running CI and by a fresh
  `conda env create` locally (imports clean — the same versions work; CI's solve
  was the outlier).
- **NumPy 2.5 `arr.shape=` deprecations.** Surfaced once the geodata tests
  un-skipped. Converted ANUGA's own 14 sites to `arr.reshape(...)`; suppressed
  the remaining one (rasterio ≤1.5.0 doing it internally in `raster.read()`) via
  a `filterwarnings` ignore, removable after a rasterio update.
- **PR triage.** Reviewed **#147** (BFS locality partitioning — recommended
  *not* merging: functionally equivalent to existing RCM but 3–4× slower, per the
  author's own benchmarks). Reviewed **#150** (timestep clamp at yield/final
  boundaries — a good dedup + real negative-timestep fix; LGTM with minor notes;
  merged).
- **`tools/` cleanup.** Renamed `install_gpu_anuga.sh` → `install_anuga_nvc.sh`
  (it's the nvc GPU build); removed 12 obsolete scripts (Travis/AppVeyor, Python
  3.8 era, Ubuntu 20.04, Travis-era conda, dead `old_div` helper); added
  `tools/README.md` naming the canonical path (`install_miniforge*` +
  `pip install -e .`) and `environments/*.yml`/`CLAUDE.md` as authoritative.
  `tools/` went 22 → 11 files (commits `dac3d8a5`, `2cf16e91`, `8bce344a`,
  `a02fe190`).

**Session 44 (2026-06-29 – 07-01):** SWW writer crash fix, laptop-guide
reconciliation, and GPU install-script fixes.
- **SWW crash on `main`.** `Write_sww.store_quantities()` raised `IndexError:
  index 0 is out of bounds for axis 0 with size 0` deep into a long run
  (t≈6960 s). The checkpoint/overwrite path locates the existing time slot with
  an *absolute* `1e-14` tolerance — below the float64 ULP (~9e-13) for any t
  beyond a few tens of seconds — so `numpy.where` returns empty and `check[0][0]`
  throws, killing the run at the write. Scaled the tolerance to the time
  magnitude (`1e-9*max(1,|t|)`) and append-with-warning when nothing matches;
  reproduced at t≈6960 and verified (commit `d04fa3ef`). Cherry-picked to `main`
  as **PR #149**.
- **Laptop guide.** Reconciled `cdac_script/ANUGA_Laptop_Guide.docx`
  benchmarking claims against this file and fixed the docx: the "dual-CCD NUMA
  limits OpenMP" claim contradicted the measured **single-NUMA-node** finding
  (the OpenMP-vs-MPI gap is false sharing in the flux arrays + serial operators);
  corrected "~3x"→"~4–5x" on 16 cores; added GPU/MPI speedup figures; repaired a
  garbled `<12 GB VRAM` bullet (commit `4fc2f134`).
- **GPU install script (`tools/install_anuga_nvc.sh`).** (1) Run
  `scripts/anuga_run_isolated_tests.py` instead of `pytest test_DE_gpu_omp.py`,
  which auto-skips on a GPU build so the test step ran nothing (commit
  `99bbc29f`). (2) `rm -rf build/cp*` before the nvc build — meson-python reads
  `CC` only on the *first* configure of a build dir, so a leftover
  gcc-configured `build/cp314` kept gcc and failed the `gpu_offload=true` guard
  (commit `ec079cd7`). Verified end to end: fresh nvc build succeeds, isolated
  runner reports **65/65 GPU tests pass**. Documented the compiler-switch gotcha
  in `KNOWN_ISSUES.md` (commit `89624513`).
- Gitignored the generated `validation_tests/case_studies/towradgi/MODEL_OUTPUTS/`
  artifacts (commit `6cfb4e62`).

**Session 43 (2026-06-26):** GPU test skip visibility + weak-scaling benchmark.
`test_DE_gpu_omp.py`'s module-level skip (on a GPU-offload build, to dodge the
NVHPC mode-2 abort) now also emits a `UserWarning` with the reason, so a plain
`pytest` shows *why* it skipped in the warnings summary without `-rs` (commit
`e37c0e76`). New `examples/parallel/run_parallel_rectangular_weak_scaling.py`
(based on `run_parallel_rectangular.py`): grows the extent with the grid
(`length = cell_size*sqrtN`) **and** normalises the sloped bed by the extent, so
triangle size, water depth, and hence the explicit CFL timestep all stay constant
as the triangle count grows — true weak scaling (commit `490cc00a`). Ran a
constant-`dt` sweep (sn 475/950/1900 → 0.9M/3.6M/14.4M triangles) comparing
**MPI-16 vs OpenMP-16 vs GPU** (full table above, "Weak-scaling benchmark"): GPU
weak-scales near-perfectly (~7.84 ms/(Mtri·step), 2.2–3.1× faster than 16 CPU
cores); MPI≈OpenMP for this pure-solver case (the Towradgi MPI>OpenMP gap was
serial operators, absent here); MPI's growing `distribute` overhead (→30 s at
14.4M) isn't repaid on these short runs.

**Session 42 (2026-06-25):** Test robustness, PR #144 conflict resolution, and an
MPI build gotcha. Made the run_toml end-to-end test directory-independent: locate
the checkout via `cwd` (not just `__file__`, which lives in site-packages for an
installed pkg), then self-contained — prefer the installed `anuga_toml_run`
console command + inline-generated inputs, so it runs from any directory and
skips only when no runner exists (commits `316620e5`, `e9400c6d`). Resolved the
conflicts on **PR #144** ("Defer GPU/offload interface build to first evolve"):
develop's compute-mode refactor (`set_compute_mode`/`_ensure_gpu_interface`/
`_boundaries_ready`, commit `289bd5c7`) had independently *superseded* the PR's
mechanism, so both conflicted files were resolved to develop's version; the PR's
only remaining net contribution is the operator-audit (graceful GPU-init fallback
in `rate_operators.py`/`inlet_operator.py`); merge pushed to the PR branch
(`e73adbfa`). Diagnosed why `gpu_has_mpi()` was False on the gcc CPU build (the
four `test_parallel_sw_flow_gpu_*` skipping): `sw_domain_gpu_ext` was built
against the single-process MPI stubs because a **reused meson build dir** never
re-ran MPI detection; fixed by rebuilding into a fresh `-Cbuild-dir` (mpi4py
fallback then finds `mpicc`/`mpi.h`), documented in `KNOWN_ISSUES.md` (commit
`a07216dd`). Those four tests then pass (8/8).

**Session 41 (2026-06-23):** `anuga_toml_run` TOML scenario runner — examples,
georeferencing, and an end-to-end test. Added `examples/run_toml/` with
self-contained `simple/` (dam break) and `complex/` (floodplain) scenarios + per-
scenario READMEs (commits `0c3c9546`, `d2f85ad4`). Added `"EPSG:<code>"`
projection support to the TOML parser and made the runner resolve the scenario
CRS to an EPSG code and call `domain.set_epsg()`, so zone/hemisphere/EPSG
propagate into the SWW (commits `ca0eedd4`, `43010ab9`, `c4ef89f5`; cairns uses
`EPSG:32755`, simple `-56`→zone 56/EPSG:32756). Consolidated the real-DEM Cairns
scenario as a third TOML example `examples/run_toml/cairns/` (TOML-only; the
~9 MB DEM moved to the shared `examples/data/cairns/`; `cairns_toml_excel/` kept
as the legacy Excel front-end) and added the runner smoke test
`anuga/scenario/tests/test_run_toml_end_to_end.py` (commit `036cb484`). Fixes: a
parallel race opening the run log before its dir exists (`acc826f5`), `Make_Geotif`
misreading a MaskedArray as a filename (`7fe595b0`); plus `-ro/-rn` reorder args
on `run_parallel_rectangular` (`39f4ccc6`) and GUI font scaling (`c0b599b3`,
`f0c7418d`).

**Session 40 (2026-06-17):** Mode-2 ('unified') triage on the **GPU build** + the
isolated runner became a first-class installed tool. Running
`anuga_run_isolated_tests --pyargs anuga.shallow_water` under the unified default on
a GPU-offload build surfaced 11 failures — all GPU-offload artifacts, not solver
regressions (Session 39's all-green was the CPU build, where device memory == host
memory). Two groups: (1) **9 white-box tests** call `compute_forcing_terms()` /
`compute_fluxes()` and assert on the host `semi_implicit_update` / `explicit_update`
arrays, which mode-2 GPU computes on-device and never syncs back (`test_forcing.py`,
`test_friction.py`, `test_physics_sw.py` Manning cases, `test_data_manager.py::
test_sww_extrema`); (2) **2 numerical tests** compare against legacy-recorded
references and diverge at ~1e-6 from mode-2's reduction/eval order
(`test_regression_snapshots.py::test_dam_break_DE1_stage_snapshot`,
`test_sww_interrogate.py::test_get_maximum_inundation_de0`). Fix: pin each to legacy
with `domain.set_compute_mode('legacy')` (snapshot helpers pinned so the whole file
is deterministic); no-op for the legacy default. The shallow_water set is now
**408 pass / 2 skip** under `-cm unified` on the GPU build (commit `0c50947d`).
Then **moved** the harness `anuga/shallow_water/tests/run_isolated_tests.py` →
`scripts/anuga_run_isolated_tests.py`, installed it to bindir via `meson.build`
(`configure_file`, matching the other `anuga_*` scripts), and made it install-safe
(importlib-resolved default target; cwd-seeded rootdir; `_abs_nodeid` passes through
absolute/`--pyargs` ids) — commit `37eccc6d`. Added **`-cm`/`--compute-mode
{legacy,unified}`** to set `ANUGA_DEFAULT_COMPUTE_MODE` for every child (omit to
inherit; banner prints the resolved mode) — commit `34401cde`. Docs: `KNOWN_ISSUES.md`
(green-run note + new command/flag), `CLAUDE.md` (testing section), new
`CONVENTIONS.md` → "Compute mode in tests" (commit `f57d0532` + this session's docs).

**Session 39 (2026-06-15):** Mode-2 ('unified') unit-suite triage — drove the fast suite
to **zero failures** under `ANUGA_DEFAULT_COMPUTE_MODE=unified` (2657 passed; 2658 in
legacy). Four genuine code fixes (commit `6bec9f8b`): (1) `recorded_min/max_timestep`
now records the CFL step before the yield cap (new `GD.recorded_flux_timestep` exposed by
the C evolve kernels, read by the Python step wrappers); (2) `gpu_manning_friction`
dispatches sloped-vs-flat on `domain.use_sloped_mannings` (new GD flag); (3) Domain
pickling restored — `__getstate__`/`__setstate__` drop the non-picklable cdef GPUDomain
and rebuild lazily; (4) `set_boundary()` invalidates+rebuilds the mode-2 device interface
so mid-run boundary changes (Reflective→Dirichlet) reach the device — fixed a progressive
runup-inundation divergence. Test/doc adaptations (commit `ff9083b4`): pinned legacy on the
deprecated forcing-function tests (Rainfall/Inflow are mode-1-only, skipped in mode 2),
skipped `test_default_is_legacy` under the env override, and pinned
`run_parallel_riverwall.py` to legacy. **Known gap:** mode-2 riverwall flux diverges from
legacy (~0.095 m on the riverwall case) — root-cause/fix tracked as FUTURE_WORK P1.9.

**Session 38 (2026-06-11):** Mesh reordering suite for OpenMP/MPI/GPU cache locality.
Added `rcm_partition()` (scipy `reverse_cuthill_mckee` on triangle adjacency graph),
`_rcm_within_partition()` helper, and `metis_rcm`/`metis_hilbert` hybrid methods to
`anuga/parallel/partitioning.py`. `reorder_domain()` gained `n_procs` parameter
(defaults to `OMP_NUM_THREADS`); for metis variants this controls partition count.
Added `-rn`/`--reorder_nprocs` CLI arg to standard parser so partition count can be
set independently of `OMP_NUM_THREADS`. Added `metis_rcm`, `metis_hilbert`, `rcm` to
`-ro` CLI choices. `hilbert_order_from_points()` refactored to use new
`hilbert_codes_from_points()` helper (raw codes without argsort, needed for
per-partition Hilbert sorting). Best OpenMP result: **`-ro metis_rcm` at 17.43 s (5.5×)**
vs 22.73 s baseline — 23% improvement. Sweet spot is `n_procs=OMP_NUM_THREADS=16`.
MPI + RCM: 11.63 s vs 12.31 s baseline (5.5% improvement). GPU + Hilbert: 5.62 s vs
6.25 s baseline (10% improvement, 17.1×). GPU + metis_rcm: 5.79 s (better than no
reorder but worse than hilbert — GPU warps prefer Hilbert's spatial clustering for
coalesced access over RCM's sequential bandwidth minimisation). Hardware is single NUMA
node so NUMA effects are absent; remaining OpenMP vs MPI gap is false sharing in flux
arrays. OMP binding flags and numactl both harmful. Best reorder per mode:
OpenMP→metis_rcm, MPI→rcm, GPU→hilbert. Full benchmark table in SESSION_GUIDE.
Commits `b5ed2b6d`, `93670929`.

**Session 37 (2026-06-11):** CLI improvements to `run_small_towradgi.py` and standard
arg parser. Added `--multiprocessor_mode`/`-mpm` (choices 1/2, default 1) to standard
parser; `run_small_towradgi.py` reads it from `args`. Added `-os`/`--outputstep` to
standard parser (SUPPRESS default, scripts default to `yieldstep`). Wired `-ys`/`-ft`
overrides in `run_small_towradgi.py` via `getattr(args, key, default)` — removed
`verbose = False` hardcode that was silently swallowing `--verbose`. Fixed
`Rate_operator` crash on empty `local_rates[fid]` array (NumPy 2.x raises `ValueError`
on `.max()` of zero-size array; MPI ranks with no polygon triangles hit this). Towradgi
benchmark results recorded (see table above). Commits `a51eec3d`–`7e6b8f40`.

**Session 36 (2026-06-11):** GCC CPU build fix after `gcc-15-offload-nvptx` install:
`-foffload=disable` added to GCC `openmp_c_args` and link args in `meson.build` so
`sw_domain_openmp_ext` (CPU-only) doesn't trigger the nvptx mkoffload pass. `--threads`
argument added to `scripts/anuga_benchmark_omp.py` (comma-separated OMP_NUM_THREADS
override, e.g. `--threads 1,2,4,8,16`). Parallel GPU test hang fix: added
`gpu_get_num_devices()` to `gpu_domain_core.c`/`gpu_domain.h`, exposed as
`get_num_gpu_devices()` in `sw_domain_gpu_ext.pyx`; parallel GPU tests now skip per-rank
when `real_gpu_available() and get_num_gpu_devices() < N` — CPU_ONLY_MODE builds are
unaffected and continue to run. Towradgi benchmark data downloaded (~86 MB). Full test
suite: 2852 passed. Commits `a96c3ca6`, `c5eb80eb`.

**Session 35 (2026-06-10):** GPU build on RTX 5070 laptop (Blackwell GB206M, 8 GB,
cc120). GCC 15 nvptx backend ICEs on `core_kernels.c` (`ompdevlow` GIMPLE pass
segfault in `core_extrapolate_second_order_edge`) — unfixable at source level.
Solution: NVIDIA HPC SDK 26.3 (`nvc`) installed via apt; `meson` auto-detects it
as `nvidia_hpc`; build command: `CC=nvc pip install --no-build-isolation -e .
-Csetup-args=-Dgpu_offload=true -Csetup-args=-Dgpu_arch=cc120`. All 56 GPU tests
pass. New `tools/install_anuga_nvc.sh`: auto-detects nvc under
`/opt/nvidia/hpc_sdk/Linux_x86_64/`, configurable via `PY`/`GPU_ARCH`/`NVHPC_ROOT`.
Fixed: `pytest-regressions` missing from all 5 intel conda environment YMLs (already
present in non-intel variants; pip-installed in existing env). `KNOWN_ISSUES.md`
updated with full GPU build recipe. Full test suite: 2852 passed. Commits
`6a10f5d0`, `f47967d2`.

**Session 21 (2026-04-15):** Domain work-array memory reduction — ~740 MB saved
at N=2.25M. DM1: 9 dead C work arrays removed from `_ensure_work_arrays()`;
only 3 live arrays remain. DM2: `edge_flux_type`/`edge_river_wall_counter` lazy
for non-riverwall simulations; `sw_domain.h` NULL guard added. DM3:
`domain_memory_stats`/`print_domain_memory_stats`/`domain_struct_stats`/`print_domain_struct_stats`
added to `system_tools.py` and exported from `anuga`.

**Session 22 (2026-04-21):** `anuga_animate_sww_gui` major feature release.
Parallel frame generation (ProcessPoolExecutor, fork on Linux, up to 4 workers)
via new `_animate_worker.py`. Zoom region (Set Zoom / Reset Zoom). `elev`
quantity: static or time-varying (erosion); terrain colormap; timeseries panel.
Fix View Mesh multiple windows, Cancel button, app-close hang. Sphinx docs with
automated screenshot capture. Commit `ebc68c37`.

**Session 23 (2026-04-24):** `anuga_sww_gui` GUI overhaul (renamed from
`anuga_animate_sww_gui`). Baked overlays (elevation contours + mesh at correct
z-order). Multi-point timeseries picking (tab10 palette, legend, CSV export,
Clear button). Save Frame time-selection dialog. 3-tab ttk.Notebook UI. Basemap
checkbox for mesh viewer and Save Mesh dialog. Sphinx RST, help, screenshots
updated. Commit `49c5b7d8`.

**Session 24 (2026-04-24):** P2.5 `Rate_operator` factories — `rainfall()` and
`inflow()` classmethods, input validation (TypeError/ValueError), 13 new tests.
P1.4/P1.8: `gauge.py` print hygiene (all `log.info()`), `file_function.py`
FIXMEs resolved. P3.6 erosion delta-bed view (`elev_delta` quantity, RdBu_r
colormap, symmetric auto-limits, 6 new tests). All P1 FUTURE_WORK items done.

**Session 25 (2026-04-25):** P2.3 `create_riverwalls` refactor — extracted 3
helpers, orchestrator ~50 lines. P2.2 `Generic_Domain.__init__` refactor —
extracted 4 helpers, `__init__` ~25 lines. Split `test_shallow_water_domain.py`
into 5 files; cleanup −101 lines. Fix 383 pytest warnings (animate.py,
rate_operators.py, pyproject.toml). P2.8/P2.9 scenario validation and TOML
docs done. claude/ rationalisation.

**Session 26 (2026-04-26):** P3.3 `fit_interpolate` L-curve alpha auto-selection.
`Fit.select_alpha()`: 20 log-spaced candidates (1e-6 … 100), numerically stable
RSS via normal equations, max-curvature corner detection, fallback to DEFAULT_ALPHA.
`dok_to_csr` added to `fitsmooth_ext.pyx` (non-destructive DOK→CSR). `alpha='auto'`
wired in `Fit.fit()`. Removed dead `_RawCSR`/`_SumRawCSR`. 13 new tests covering
row_ptr extension, multi-attribute, degenerate/interior paths. fit.py 85→92%.
P2.7 gauge modernisation fully done (session continuation). Commit `12864187`.

**Session 27 (2026-04-27):** `Kinematic_viscosity_operator` MPI-parallel (Option B
distributed CG). Phase 1: removed Apple OpenMP guards from 4 C files. Phase 2:
`parabolic_solve` serial path routed through C CG (`cg_solve_c_precon`) with Jacobi
preconditioner. Phase 3: full distributed CG — `_exchange_ghost_vector` (MPI tag 198
non-blocking), `_distributed_dot` (Allreduce), `_parabolic_matvec_distributed` (ghost
exchange before SpMV), `_parabolic_solve_distributed` (n_full-length CG loop). `parallel_safe()`
returns True. New `run_parallel_kv_operator.py` + `test_parallel_kv_operator.py`
(serial-vs-parallel xvel, max diff 8.6×10⁻⁶). New `run_parallel_kv_unit_tests.py`
+ `test_parallel_kv_unit_tests.py` (4 in-process MPI assertions for each primitive).
Bug fix: `test_select_alpha_degenerate_falls_back_to_default` platform-dependent on
Windows py3.10/3.11/3.13 — now uses `return_curve=True` to branch on actual kappa.
Commits `61418742`, `5498f98d`. All CI passed.

**Session 32 (2026-05-03):** CI fixes. macOS `sw_domain_gpu_ext` import failure
resolved: `_omp_target_is_present` missing from conda-forge libomp — fixed by adding
`static inline` stubs + macro redirects for all OpenMP 4.5 target-alloc API functions
(`omp_target_alloc/free/memcpy`, `omp_target_is_present`, `omp_get_initial_device`) in
`gpu_omp_macros.h` under `CPU_ONLY_MODE`. Also: MPI ABI fix for Linux conda envs (use
`mpicc -show` compile flags, not hardcoded sys.prefix), Windows MPI detection (Library/
include path, msmpi), macOS rpath for libmpi.dylib (`@executable_path/../lib` +
`DYLD_LIBRARY_PATH` in CI), macOS simd-reduction warnings (dropped `simd` from reduction
macros), CI Node.js 24 upgrades (setup-miniconda v4, setup-python v6), skip-reason
diagnostics in test_DE_gpu_omp.py. macOS GPU tests now run (not skipped).
Commits `94c4d74f`–`7cc69fac`.

**Session 31 (2026-04-29):** ADER-2 GPU wiring + optimisation. `gpu_ader_ck_predictor` /
`gpu_evolve_one_ader2_step` added to `gpu_kernels.c` + declared in `gpu_domain.h` (fix for
missing header causing Windows/3.12 CI failure). `evolve_one_ader2_step` dispatches to
`_evolve_one_ader2_step_c` / `_evolve_one_ader2_step_gpu` in GPU mode. `DE_ader2` flow
algorithm added to `set_flow_algorithm()` (DE1 settings + ader2 timestepping). Fused
predict-extrapolate C-K loop (P3.8): merged extrapolation and C-K predictor into a single
kernel pass — eliminates the second extrapolation entirely → **1.75× faster than DE1**.
`FLOPS_ADER_PREDICTOR=105` constant. Fix: `create_sts_boundary`
in `sts.py` now calls `fid.close()` in a `try/finally` block (Windows WinError 32).
`develop_ader` merged into `develop`. Commits `e9d15803`–`3b00dc79`.

**Session 30 (2026-04-29):** ADER-2 timestepping via Cauchy-Kovalewski predictor.
`core_ader_ck_predictor()` in `core_kernels.c`: recovers cell slopes from the 2×2
edge system (edges 0 and 1), evaluates well-balanced SWE time derivatives locally
(dz/dx = dw/dx − dh/dx from reconstruction), advances centroids by dt in-place.
`evolve_one_ader2_step()` in `shallow_water_domain.py`: backup → CFL step →
C-K predictor(dt/2) → midpoint flux → restore Q^n → update. `boundary_flux_integral_operator`
updated for 'ader2'. 10 tests: well-balance (flat/sloped), mass conservation, dam-break
consistency, non-negative depths. All 2656 fast-suite tests pass. Commit `825f1e5f`.

**Session 29 (2026-04-28):** Investigated numpy `_NoValue` reload issue triggered
by `--cov=anuga.submodule` targeted coverage runs. Root cause: coverage.py's
`sys_modules_saved()` context in `inorout.py` calls `importlib.util.find_spec()`
on the subpackage, auto-importing the parent chain including numpy, then purging
all new imports from sys.modules. The second real numpy import fires the reload
guard and corrupts the C-level `_NoValue` singleton. Fix: removed `numpy._pytesttester.PytestTester`
from `anuga/__init__.py` (replaced with plain `def test()`) and all 23 subpackage
`__init__.py` files (was unused boilerplate). Targeted submodule `--cov=` runs are
unfixable from conftest.py (pytest-cov starts before conftest.py loads);
use `--cov=anuga` always. Documented in KNOWN_ISSUES.md. Commit `af71f10b`.

**Session 28 (2026-04-27):** P2.6 fast-suite coverage continued. `anuga/file/`:
`test_netcdf_nc.py` (10 tests, netcdf.py 34%→100%), `test_sts.py` (11 tests,
sts.py 47%→89%). `anuga/structures/`: 9 new tests in `test_inlet_operator.py`
(inlet_operator.py 45%→64%); 16 new tests in `test_structure_operator.py`
(structure_operator.py 65%→96% — enquiry getters, setters, skew 4-point, error
paths, print/timestepping stats, non-constant elevation warning). Overall fast
suite: 58.13% → 58.68%.

---

## File locations for common operations

| Task | Files |
|------|-------|
| Add public API export | `anuga/__init__.py` (import + `__all__`) |
| Add slow test marker | `@pytest.mark.slow` decorator or module-level `pytestmark` |
| Configure pytest options | `conftest.py` (repo root), `pyproject.toml` `[tool.pytest.ini_options]` |
| Per-test process isolation runner | `scripts/anuga_run_isolated_tests.py` (installed `anuga_run_isolated_tests`; `-cm legacy\|unified`) |
| Default compute mode | `ANUGA_DEFAULT_COMPUTE_MODE` env var; `domain.set_compute_mode('legacy'\|'unified')` |
| Memory reporting | `anuga/utilities/system_tools.py::memory_stats()` |
| Timestepping output | `anuga/abstract_2d_finite_volumes/generic_domain.py::timestepping_statistics()` |
| Triangle quiet/verbose | `anuga/pmesh/mesh.py::_generateMesh_impl()` |
| TOML scenario config | `anuga/scenario/`, `scripts/anuga_toml_run.py`, `examples/run_toml/` (simple/complex/cairns; shared DEM in `examples/data/cairns/`); legacy Excel front-end in `examples/cairns_toml_excel/` |
| Single-process benchmark | `benchmarks/run_benchmarks.py` + `benchmarks/compare_benchmarks.py` |
| MPI distribution benchmark | `benchmarks/distribute_benchmarks.py` + `benchmarks/run_benchmark_grid.py` |
| Validation report (one algorithm) | `validation_tests/reports/validations_produce_results.py -alg X [-l] [-nr]` |
| Combined validation report | `validation_tests/reports/validations_combine_algs.py collect -alg X` then `combine` |
| Add a case to the reports | `validation_tests/reports/validations_report_body.tex` (`\inputresults{../dir/case}`) |
| Release reports + Zenodo | `.github/workflows/validation-report.yml`, `.github/scripts/zenodo_validation_reports.py`, `claude/RELEASE_VALIDATION_REPORTS.md` |

---

## Key reference documents

| Document | URL |
|----------|-----|
| Hydrata refactor plan | https://github.com/Hydrata/anuga_core/blob/anuga-4.0-refactor-plan/REFACTOR_PLAN.md |
| anuga-community GitHub | https://github.com/anuga-community/anuga_core |
| Hydrata fork | https://github.com/Hydrata/anuga_core |

## Next priorities

See `claude/FUTURE_WORK.md` for the full prioritised list (P1–P3).

**SC26 (needs GPU hardware):** G4.1 Gordon Bell metrics, G4.2 physical benchmark validation, G4.3 multi-node strong scaling (scripts in `benchmarks/` and `scripts/hpc/` are ready).

**Best standalone value:** P2.6 fast-suite coverage, P2.7 gauge module modernisation, P2.4 culvert compute_rates deduplication.
