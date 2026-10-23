# 4.0.0 Phase-1 verification gates — results log

Freeze commit: `55ef856d` (was b90d2e40; +1 test-robustness fix) (develop, 2026-08-21). Plan: `claude/RELEASE_PLAN_4.0.0.md`.
Machine: local dev box (Ubuntu 26.04, RTX 5070, nvc GPU build) unless noted.

| # | Gate | Status | Result |
|---|------|--------|--------|
| 1 | Full suite (as CI runs it) | **GREEN** both builds | nvc: 2937 passed/10 skipped, 74 s. gcc venv: 2948 passed/104 skipped (no mpi4py in venv; MPI covered by the nvc run) |
| 2 | Unified-mode suite, CPU build, one process | **GREEN** | FULL suite, unified default, one process, gcc build: 2947 passed/105 skipped. GPU test file 106/106 after 55ef856d |
| 3 | GPU build, isolated runners | **GREEN** (local + cloud) | local RTX 5070 (cc120): 25 classes + unified sweep 465 pass/2 skip. Cloud **g6.xlarge (Ada L4)**: 106/106 GPU file, then 465 pass/2 skip unified — identical counts on a second GPU architecture. rc=0, container 723 s |
| 4 | Validation suite | **GREEN** | 120 passed / 0 failed (1895.9 s). NOTE: the structure regression baselines passed *unchanged* — the existing validation set is blind to #229 (flat beds at every inlet) |
| 5 | Towradgi mode 1 vs mode 2 | **GREEN** (rerun) | isolated dirs, monotonic axes asserted: final max\|Δ\|=3.52e-02 m, localised (9/128539 cells >1 cm), mass agrees 3.4e-06 rel. See "Gate 5 rerun" |
| 6 | MPI smoke | **GREEN** | parallel MPI tests (mpirun subprocess spawns) ran inside gate 1's 2937 |
| 7 | Wheel smoke (3 OSes) | **GREEN** (all 3) | Windows was red 08-20 → 08-22 on a conda-forge mingw toolchain regression (stack-protector ABI vs the build-11 CRT); fixed by pinning the sysroot chain to build 10, PR #233 / `6d5e7c4b`. 20 wheel builds + sdist green |
| 8 | Fresh conda installs (3.10, 3.14) | **GREEN** | env create + `pip install --no-build-isolation` + evolve smoke: py3.10.20/numpy 2.2.6 and py3.14.7/numpy 2.5.2, both at `55ef856d` |
| 9 | Docs build clean | **GREEN** | `make html` build succeeded, 0 warnings/errors; `generated/` now gitignored |

**Found by the gates so far:** gcc/nvc disagree on line-region triangle selection when an
exchange line lies on mesh edges (32 vs 54 triangles) — filed as #231, test made
geometry-robust in `55ef856d`. This is the process working: a red gate, a diagnosis, a
fix, an issue.


---

## Deltas from #229 (for the Phase-2 release notes)

Method: the pre-#229 behaviour is reproduced **in-process** by patching
`Inlet.set_average_depth` / `set_average_momenta` (and the `Parallel_Inlet`
overrides) back to the uniform-depth write — same build, same commit, only the
#229 change toggled. Cleaner than rebuilding an old checkout and isolates
exactly the write-back effect. Harness: `scratchpad/delta_runner.py`,
`delta_sweep.sh`, `wb3.py`.

**1. Controlled case — lake at rest, idle culvert on a sloping bed** (the defect
in isolation; `wb3.py`):

| bed slope | t | old (uniform depth) | new (leveling) |
|---|---|---|---|
| 1/50 | 1 s | 6.805e-02 m | 3.664e-09 m |
| 1/50 | 5 s | 1.306e-01 m | 7.353e-04 m |
| 1/200 | 1 s | 1.797e-02 m | 2.423e-07 m |
| 1/200 | 5 s | 1.284e-01 m | 7.278e-03 m |

The "new" growth over seconds is the still-water roundoff amplification
documented in KNOWN_ISSUES, not the tilt.

**2. Validation structure cases** (capped at 3600 s simulated, both arms):

| case | final max\|Δ\| | peak max\|Δ\| | final RMS |
|---|---|---|---|
| bridge_hecras | 1.3547e-02 m | 1.3545e-02 m | 3.4143e-04 m |
| bridge_hecras2 | 0 | 0 | 0 |
| lateral_weir_hecras | 0 | 0 | 0 |
| tides_hecras | 0 | 0 | 0 |

Three of four are **bit-identical**, confirming the design claim: on a flat bed
the leveling write reproduces the old uniform-depth write exactly. Only
`bridge_hecras` moves (1.4 cm peak, 0.34 mm RMS) — worth checking *why* it
differs from `bridge_hecras2` before writing the notes (likely a sloping or
non-uniform bed under one inlet; **not yet confirmed**).

**3. Towradgi** (22 Boyd culverts on real terrain, 1 h, 128 539 points):

| metric | value |
|---|---|
| final max\|Δ\| | **3.1817e-01 m** |
| 99th percentile | 1.5450e-03 m |
| cells >1 mm | 2340 / 128539 |

**Why so much larger than the validation cases:** the towradgi inlets sit on
real terrain. Bed-elevation spread across the 44 inlets of its 22 structures:
**median 0.774 m, max 4.954 m, and 41 of 44 inlets exceed 1 cm** (worst:
Branch_5_Collins_St_Culverts inlet 1, 4.95 m across 12 triangles). Since the
old write-back tilted a level surface by ~half the bed range across an inlet,
the predicted worst-case tilt is ~2.48 m; the observed 0.32 m is well inside
that, as expected when only some of those steep inlets carry water in a 1 h run.

This is the release-notes headline: the synthetic flat-bed validation cases were
**blind** to #229 (3 of 4 bit-identical), while a real-terrain model moves by
decimetres. Analysis script: `scratchpad/inlet_slopes.py`.

---

## Cloud gate blocked — the isolated runner crashes inside the container

Two AWS runs, both g4dn.xlarge (no g6 capacity in any Sydney AZ today):

1. `gate3-g6/` — rc=2 in 15 s: **the slim image ships no pytest**
   (`No module named pytest`). Worth a Dockerfile note or a `[test]` extra.
2. `gate3-g6b/` — with `pip install pytest` prepended: rc=1, **all 106 tests
   report CRASH at ~0.2 s each**.

**Reproduced locally in the same image**, so it is NOT the T4 and not the cloud:

```
docker run --rm --gpus all anuga:gpu-slim-4rc \
  sh -c "pip install -q pytest && anuga_run_isolated_tests"   # 106x CRASH
```

But the same test run *directly* in the same container passes:

```
docker run --rm --gpus all -e ANUGA_GPU_TESTS_ISOLATED=1 anuga:gpu-slim-4rc \
  sh -c "pip install -q pytest && python -m pytest --pyargs \
         anuga.shallow_water.tests.test_DE_gpu_omp -k test_flux_kernel -q"   # 1 passed
```

So the runner's *child spawning* is what breaks in a container, not the tests.
State of the investigation at handoff:

* `run_one()` (scripts/anuga_run_isolated_tests.py:232) spawns
  `[sys.executable, -m, pytest, nodeid, -p, no:cacheprovider, -q, --no-header,
  -o, addopts=]` with `cwd=str(ROOTDIR)` and `env=_base_env()`.
* Running that exact command by hand in the container **works** (rc=0).
* `classify()` maps "non-zero exit, no pytest summary" to CRASH, so CRASH here
  may be a *misclassification* of a child that failed for an environmental
  reason (e.g. `ROOTDIR` resolution in a site-packages install — `_find_rootdir`
  walks up from `Path.cwd()`; in the container there is no source checkout).
* **Next step:** print a failing child's captured stdout/stderr (the runner
  keeps it in `run_one`'s return) rather than inferring — e.g. run the runner
  with a single `-k` selection in the container and dump the output.

**RESOLVED — fixed in `69c586f6`.** Not the T4, not CUDA, not the NVHPC
runtime. `-v` showed the children were handed unrunnable node ids
(`file or directory not found: test_DE_gpu_omp.py::Test::test_flux_kernel`):
pytest emits ids relative to *its own* rootdir (the common ancestor of the
targets), while the runner resolved them against its own cwd-derived `ROOTDIR`.
In a checkout the two coincide — which is exactly why this never showed up
locally; against an installed anuga they do not, so every child failed to find
its test and `classify()` mapped the result to CRASH. (`--rootdir` is not the
fix: it makes pytest emit ids with an *empty* path component.)

The runner now resolves ids against target-derived base directories (a file
target's parent, a directory target itself, a `--pyargs` module/package via
importlib), trying `ROOTDIR` first so in-checkout behaviour is bit-identical,
and accepting only candidates that name a real file. It also now **errors
loudly** — naming the id and the directories searched — when an id resolves to
no file, instead of spawning N children that each report the same thing as a
test CRASH.

**Second, separate gap: the slim image ships no test dependencies.** Run 1 died
on `No module named pytest`; run 3 also needed `pytest-regressions` or the
snapshot tests error at fixture setup (`fixture 'num_regression' not found`).
Worth a `[test]` extra in the image or a documented line in `docker/README.md`,
so the image can verify itself.

Cost: four instances, ~15 min of GPU time total.


---

## Gate 5 is invalid (ran, but do not use the number)

The run completed (mode 1 13.6 min at 16 threads, mode 2 3.5 min on the 5070)
and printed `final-stage max|d| = 3.52e-02 m`. **Do not quote that number.**
The two SWW files do not have comparable time axes:

```
mode1  ntimes 9  times [0, 60, 600, 120, 1200, 1800, 2400, 3000, 3600]   <-- non-monotonic
mode2  ntimes 7  times [0, 600, 1200, 1800, 2400, 3000, 3600]
```

**Root cause identified — filed as #232.** The frame indices show two writers
alternating on the shared unlimited dimension: my run occupies indices
0,2,4,5,6,7,8 and a second writer (yieldstep=60) occupies 1,3. The towradgi
script hardcodes `domain_name` (line 129) and `MODEL_OUTPUTS` (line 382), so
*every* run of that case study writes the same path; gate 5's mode-1 arm
(18:19-18:33) overlapped the validation suite, still writing at 18:28.
`sww_merge` is NOT involved (it only runs when `numproc > 1`; these were
serial) — an early guess of mine that the index pattern disproved.

**Rerun requirement:** point `--datadir` at a clean, empty directory per arm,
and assert monotonic `time` before differencing.

## Towradgi #229 delta — not obtained (harness bug, one-line fix)

The chained old-write arm died immediately:

```
File "run_small_towradgi.py", line 12: from project import *
ModuleNotFoundError: No module named 'project'
```

`runpy.run_path()` does not put the script's own directory on `sys.path`. The
four HEC-RAS cases were unaffected because they import no local module. Fix:
`sys.path.insert(0, case_dir)` in `scratchpad/towradgi_runner.py` (and
`delta_runner.py`, which has the same latent bug). So the Towradgi row of the
#229 delta table is **still missing** — the four validation cases above stand.


---

## Gate 5 rerun (clean) — 2026-08-21 19:18-19:44

Three arms, each forced into its own datadir (`Domain.set_datadir` patched),
with monotonic-time and matching-axes assertions before differencing.

| arm | wall time |
|---|---|
| mode 1 (CPU, 16 threads) | 11.9 min |
| mode 2 (GPU, RTX 5070) | 3.6 min — **3.3x faster** |
| old write-back (mode 1) | 11.1 min |

`GATE5 mode1-vs-mode2 final max|Δ| = 3.5246e-02 m` — **identical to the number
the corrupted file produced.** The original figure was in fact correct; the
corrupt file simply was not valid evidence for it. Rejecting it was still right.

Character of that 3.5 cm: 9 of 128 539 cells exceed 1 cm (0.007%), 99th pct
0.54 mm, total water agrees to 3.4e-06 relative, worst cell in 16 cm of water
(shallow flow, not a dry cell). Consistent with the accepted mode-1/mode-2
divergence documented in FUTURE_WORK (chaos amplification of a 1-ULP wet/dry
seed). **Caveat kept for the record:** FUTURE_WORK records an earlier towradgi
figure of 2.5e-3; this max is ~14x larger, but the configurations differ
(duration, mesh scale, yieldstep) so it is not a like-for-like comparison and
should not be read as a regression without one.
