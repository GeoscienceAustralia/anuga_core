# Future Work Recommendations

Generated: 2026-04-24 (session 23). Last updated: 2026-10-02 (session 58).
Based on codebase investigation cross-referenced against 50 sessions of completed work.

Items marked ~~strikethrough~~ have been invalidated (see notes).

> **See also:** `claude/C_EXTENSION_AUDIT_TODOS.md` (2026-06-09) — prioritised
> TODOs from an audit of all 20 C/Cython extensions (correctness, GPU/CPU
> kernel unification, performance, cleanup).
>
> **Active multi-step plan:** `claude/PLAN_default_mode2_cpu.md` (2026-06-12) —
> migrate the standard distribution to `multiprocessor_mode=2` + `gpu_offload=false`
> (CPU-multicore C operators by default). Step 1 in review as PR #144; **step 2
> (audit operator fall-back) is the next action.**

---

## Session 50 follow-ups (2026-07-19) — GPU mode-2 / OpenACC

Open items from the mode-1-vs-mode-2 investigation (session 50). Seven bugs from that
thread are already fixed on `develop` (#191–#194, #197, #199, #200); these remain.

**P1 — #190 OpenACC `set_gpu_offload(False)` / `-ngo` silently stays on the GPU.** MERGE
BLOCKER for PR #188. The host-fallback idiom `omp_set_default_device(omp_get_initial_device())`
maps to `acc_set_device_num(-1, ...)`, and a negative devicenum is NOT "run on host" in
OpenACC — the run silently stays on the GPU with no error. OpenACC has no device-*number*
for host execution (it is a device-*type* concept). Likely fix: **hard-error** on the OpenACC
build rather than implement, and lean on the separate `g_gpu_offload_enabled` flag. Only on
the `#188` branch, so not live on `develop`. (Full analysis in the issue and session guide.)

**P2 — #189 mode-2 never runs `apply_protection_against_isolated_degenerate_timesteps()`.**
It hangs off `update_timestep()`, which the mode-2 step path returns early past (all three
`evolve_one_*_step` functions dispatch to `_evolve_one_*_step_c` before it). Default-OFF
(`config.py:151`), so low priority — but a user who enables it under GPU gets no protection
*and no warning*. Cheapest fix: warn in mode 2 when it is enabled (mode 2 already warns for
other unsupported forcing).

**P2 — PR #188 (OpenACC backend) is WIP and needs the author.** Reviewed + numerically
validated in session 50: bit-reproducible, agrees with OpenMP-target at the mode-1-vs-mode-2
tolerance (see `validation_tests/case_studies/towradgi/compare_openmp_openacc.ipynb`), and the
per-kernel `fluxes-central` ~20% win is real. Before it leaves WIP: fix #190; explain the
unmentioned `-mp=gpu,multicore` → `-mp=gpu` change to the *default* build; drain the queue at
the `acc_free` teardown in `gpu_halo.c`; add tests; and rebase onto `develop` (it lacks #199,
#200). The kernel win does **not** reach wall-clock (the async queue issues 2.5× more
`cuStreamSynchronize` than OpenMP-target) — headroom, not a blocker, but worth noting.

**~~P2 — Unify the culvert implementation so mode-1 and mode-2 are bit-identical.~~ DONE
(session 51, uncommitted in working tree as of write-up).**

**What was done.** There is now **one** implementation of the Boyd/weir per-culvert update —
`culvert_compute_one()` in `shallow_water/gpu/gpu_culvert_operator.c` — plus a shared host inlet
gather `culvert_gather_inlet_host()`. Mode-2's batch was refactored to call it (behaviour
preserving); mode-1's Python operators (both `Structure_operator` and the *default*
`Parallel_Structure_operator`) now route their per-step update through it via a Cython bridge
(`culvert_apply_one_host` / `culvert_gather_inlet_host_py` in `sw_domain_gpu_ext.pyx`), gated to
fully-local culverts (cross-boundary MPI culverts keep the Python+MPI path). Files touched:
`gpu_culvert_operator.c/.h`, `sw_domain_gpu_ext.pyx`, `structures/structure_operator.py`,
`parallel/parallel_structure_operator.py`. Result: **mode-1 == mode-2 bit-for-bit** for every
culvert config (box/pipe, velocity head, blockage), at 1 and 16 threads. Full suite green in
legacy (2698) and unified (2697); 230 structure/GPU tests pass. De-dups the Python/C physics on
the runtime path (single source of truth).

**Correction to the earlier diagnosis (don't repeat the wrong turns):** the seed was **not**
momentum, ordering, or FMA. It was the **inlet-average gather** — mode-1 summed with numpy
(`num.sum(v*a)/area`), mode-2 with a C loop; for a multi-cell inlet these round 1 ULP apart in a
value-dependent way (plus a dry-cell depth clamp the C gather does and the Python global-average
path didn't). Two earlier "proofs" were invalid because a monkeypatch hit the **unused serial
`Structure_operator`** while the default operator is `Parallel_Structure_operator` (a *separate*
class hierarchy). Lesson: `anuga.Boyd_box_operator` is a **factory** returning the parallel
class even in serial — instrument the class that actually runs.

---

**Towradgi still diverges — and it is NOT culverts (diagnosed, recommend ACCEPT).**
After the culvert fix, towradgi mode-1 vs mode-2 is unchanged (7.6e-6 → 2.5e-3). An in-process
double-precision localization harness (two domains from the same setup, restartable lockstep
evolve, double diff — reconstructable, see session 51) pinned the real seed:

- The seed is **rainfall (`Rate_operator`) falling on DRY cells.** Remove the rate operators from
  both domains → **zero divergence** (bit-identical); rain on a **fully wet** domain →
  bit-identical; rain on **dry** cells → 1 ULP.
- It is **not** the rate operator's arithmetic (numpy `arr+scalar` == a C loop; wet cells prove
  it, and the rate inputs `local_rate/timestep/rate/factor` are bit-identical between modes) and
  **not** the sync (a plain `omp target update` memcpy of stage). Both modes are in fact
  **stage-primary** (`height == stage - bed` exactly in each) — there is no `stage = bed + height`
  reconstruction to flip. The residual is a **1-ULP difference in the core's near-dry stage
  itself**: mode-2's stage going *into* the rate op is already ~1 ULP off from mode-1, produced by
  the unified C RK loop's dry-cell handling and **masked by the bed-clamp** (`protect` sets
  `stage = bed` exactly in both) right up until rain lifts the cell off the bed and exposes it.
  That's why bare runs are bit-identical and only *rain-on-dry* diverges. Same family as the
  **#200** dry-cell issue (a wet/dry-margin roundoff gap between the legacy and unified paths).
  The exact operation was not isolated — it needs instrumenting the C RK loop's device state
  mid-step.
- A red herring ruled out along the way: mode-2 also clamps dry cells to bed at
  `set_multiprocessor_mode(2)` (#200) while mode-1 defers to first protect — but forcing mode-1
  to clamp too left the divergence unchanged, so the *initial* state is not the cause.

**Recommendation: accept it.** One ULP in the stage of dry ground under a hair of rain —
physically meaningless, chaos-amplified to mm over hours, both modes equally valid. A real fix
means reconciling stage-primary vs height-primary at the wet/dry margin (core-level, #200
family, large blast radius, no physical payoff). Unifying the `Rate_operator` would **not** help
(it's already identical on wet cells). If ever pursued, the localization harness is the tool.

**P3 — process, not code:**
- **`develop` is ~833 commits ahead of `main`.** *Every* session-50 fix (including the #200
  startup mass-loss and #193 parallel-inlet mass-balance correctness fixes) is unreleased,
  gated behind the "no develop→main until v4.0.0" rule. Worth a conscious call on when v4.0.0
  is cut — real correctness fixes are sitting unreleased.
- **Branch-protection bypasses.** Merges/pushes to `develop` this session went through with
  `--admin` ("Bypassed rule violations — changes must be made through a pull request"). Decide
  whether the PR-required rule on `develop` should be enforced or relaxed.

---

## Priority 1 — High value, low effort (1–3 days each)

~~**P1.1 Delete or absorb `boyd_box_operator_Amended3.py`**~~ — Done (session 25).

~~**P1.2 Add tests for the `rain` module**~~ — Done (session 25).

~~**P1.3 Add tests for `simulation/` and `validation_utilities/`**~~ — Done. `test_simulation.py` and `test_validation_utilities.py` exist.

~~**P1.4 Fix `gauge.py` verbose/print hygiene**~~ — Done (P2.7 session 24). No bare `print` calls remain; all logging via `log.info()`/`log.warning()`.

~~**P1.5 Add deprecation warnings to legacy forcing classes**~~ — Done (session 25). `DeprecationWarning` added to `Inflow`, `Rainfall`, `Wind_stress`, `Barometric_pressure` in `shallow_water/forcing.py`; `filterwarnings` in `pyproject.toml` suppresses them in the test suite.

~~**P1.6 Remove local-timestepping dead infrastructure**~~ — Done (session 23). Removed
`max_flux_update_frequency`, `flux_update_frequency`, `update_next_flux`,
`update_extrapolation`, `edge_timestep`, and `allow_timestep_increase` from Python domain,
C header, Cython wrapper, scenario system, and tests. Deleted
`test_local_extrapolation_and_flux_updating.py`. See P3.1 for the future implementation plan.

~~**P1.7 Write tests for `anuga/utilities/animate.py`**~~ — Done. `test_animate.py` exists under `anuga/utilities/tests/`.

~~**P1.8 Clean up `file_function.py` FIXMEs**~~ — Done (session 24). FIXMEs already resolved; deleted dead commented-out blocks, replaced raw `fid.xllcorner`/`fid.yllcorner`/`fid.zone` reads with `Geo_reference(NetCDFObject=fid)`, cleaned up redundant `.csv` branch and trailing NOTE comment.

**P1.10 Fix GPU-build `test_DE_gpu_omp.py` mid-file abort (NVHPC present-table)** — Running
the whole file in one process aborts (~9th–11th test) on a GPU build; the NVHPC OpenMP-target
runtime calls `exit()`. Pre-existing (reproduces on `d96ae357`). Not a simple leak: each class
alone passes and 16–20 looped/live domains are fine; forcing finalization between tests makes it
worse (assertion failures). Root cause is host-pointer-keyed, reference-counted OpenMP present-table
management (`map(to:)`/`map(delete:)` in `gpu_domain_core.c`) corrupted by numpy host-address reuse
across domains, plus two reference cycles deferring finalization. Fix options: per-test process
isolation for the GPU test file (e.g. `pytest-forked`/subprocess), strict 1:1 map/unmap reference
discipline per domain, or switch to `omp_target_alloc` + `is_device_ptr` device-pointer allocation.
Production (single/few sequential GPU domains) is unaffected. **Interim workaround in place:**
`anuga/shallow_water/tests/run_gpu_tests_isolated.sh` runs the file one fresh process per class
(all classes pass); `--forked` does NOT work (CUDA is fork-unsafe). See `claude/KNOWN_ISSUES.md`
("GPU build: `test_DE_gpu_omp.py` aborts mid-file").

~~**P1.9 Root-cause and fix mode-2 ('unified') riverwall flux divergence**~~ — DONE
(2026-06-15). Misdiagnosed: the riverwall flux is correct (bit-identical mode-1 vs mode-2
with a GPU-supported boundary). The real cause was a DE0/Euler-specific boundary bug —
`evolve_one_euler_step()` went straight to the C Euler loop and silently ignored non-GPU
boundary types (`run_parallel_riverwall.py` uses `Transmissive_momentum_set_stage_boundary`),
while rk2/rk3/ader2 already fell back to host evaluation. Fixed by adding
`_evolve_one_euler_step_gpu()` and delegating to it from `_evolve_one_euler_step_c()` when
`not self._gpu_all_on_gpu`. `run_parallel_riverwall.py` un-pinned; regression test
`Test_GPU_NonGPUBoundaryFallback`. See `claude/KNOWN_ISSUES.md` ("RESOLVED … DE0 boundary bug").

---

## Priority 2 — Medium effort (1–2 weeks each)

~~**P2.10 Remove the deprecated forcing-function classes**~~ — Done (2026-09-18).
`Wind_stress`, `Rainfall`, `Inflow`, `Barometric_pressure`, the `_fast` variants and
`General_forcing` deleted from `shallow_water/forcing.py` (only `Cross_section` remains);
`anuga.Inflow/Rainfall/Wind_stress` exports and the `pyproject.toml` filterwarnings
removed; tests ported to the operators (`test_forcing.py` rewritten). File-driven fields
(`use_coordinates=False` + `file_function`) moved into `Wind_stress_operator` /
`Barometric_pressure_operator`. `manning_friction_semi_implicit` kept as a forcing term.

~~**P2.1 Type hints on the public API**~~ — Done (session 33). Annotated ~130 public methods
across four files using `from __future__ import annotations` (PEP 563, Python 3.10+ compatible).
`base_operator.Operator` (11 methods, full coverage), `quantity.Quantity` (~28 methods),
`structure_operator.Structure_operator` (~28 methods including all `get_enquiry_*` getters),
`shallow_water_domain.Domain` (~60 methods including `__init__`, `evolve`, `set_flow_algorithm`,
`get_wet_elements`, `timestepping_statistics`, all inundation/volume queries). Uses
`numpy.typing.ArrayLike`, `collections.abc.Callable`, `TYPE_CHECKING` guards for heavy imports.
Commits `6c16986e`, `8f39e645`.

~~**P2.2 Refactor `Generic_Domain.__init__` (367 lines)**~~ — Done (session 25). Extracted
`_init_mesh()`, `_init_quantities()`, `_init_parallel()`, `_init_timestepping()`.
`__init__` is now ~25 lines. 743 domain/shallow-water tests pass.

~~**P2.3 Refactor `create_riverwalls` (300 lines)**~~ — Done (session 25). Extracted
`_validate_riverwall_inputs()`, `_match_edges_to_segments()`, `_build_hydraulic_properties()`
from the 300-line monolith. `create_riverwalls` is now a ~50-line orchestrator. All 43
riverwall tests pass.

~~**P2.4 Delete the `anuga/culvert_flows/` package**~~ — Done (session 34). Deleted the
entire package (5 412 lines removed, commit `b151fa66`): `culvert_class.py`,
`culvert_routines.py`, `culvert_polygons.py`, all tests and test data. Updated
`run_open_slot_wide_bridge.py` to drop legacy imports and show a `Boyd_box_operator`
equivalent. Removed dead culvert_flows references from `test_failure.py` and five
parallel test files. `subdir('culvert_flows')` removed from `anuga/meson.build`.
All 2 633 fast tests pass.

~~**P2.5 Improve `Rate_operator` usability**~~ — Done (session 24). Added `Rate_operator.rainfall(domain, rate_mm_hr)` and `Rate_operator.inflow(domain, rate_m3s)` factory classmethods; input validation (bad rate type → TypeError, region+polygon conflict → ValueError); updated `__init__` docstring pointing to factories. 13 new tests.

**P2.6 Raise fast-suite coverage threshold**
Fast suite at 54.66% against `fail_under=57` set for the full suite. Either set separate
thresholds in `.coveragerc` for fast vs full, or add targeted tests in `anuga/file/`,
`anuga/fit_interpolate/`, and `anuga/structures/` to lift the fast-suite baseline. Session
20 added ~90 tests as a model.

~~**P2.7 Modernise `sww2timeseries` / gauge module**~~ — Done (sessions 27–28).
- `gauge_get_from_file` rewritten with `csv.DictReader` (case-insensitive, whitespace-tolerant)
- `open().close()` file-existence checks replaced with `os.path.isfile()`
- `_generate_figures` marked `# pragma: no cover` (matplotlib/LaTeX display dependency)
- `plot_polygons` in `geometry/polygon.py` fixed: replaced `matplotlib.use('Agg')` with
  `plt.switch_backend('Agg')` (safe post-import); added defensive try-except around import
  block and plot body — resolves the matplotlib 3.10 / numpy 2.x `_NoValueType` crash
- `test_gauge.py` at 41 tests covering gauge.py at 99% (only lines 177-178, read-permission
  error path, remain uncovered)

Speculative future work: add EPSG/`Geo_reference` coordinate support to
`gauge_get_from_file` (accept optional EPSG code, convert to domain projection).

~~**P2.8 Scenario system input validation**~~ — Done (session 25). Schema validation added to TOML inputs; detailed error messages naming bad fields and expected types; range checks for physical parameters.

~~**P2.9 Document the scenario/TOML system**~~ — Done (session 25). Sphinx reference page added listing all supported TOML keys with types, defaults, and examples.

---

~~**P2.11 Overlaid algorithm comparisons for the key validation cases**~~ (session 58). **Not pursued** (decided 2026-10-02): the combined report's DE0 | DE1 | DE_ader2 panels are sufficient.
The combined report (PR #409) shows each figure as DE0 | DE1 | DE_ader2 panels,
which covers every case but compares less well than one plot with the three
algorithms overlaid, as `rundown_mild_slope/compare_solvers.py` does. Convert
the cases where the algorithms differ most (dam breaks, deep_wave, run-up,
sediment) one at a time. The open design point: the CI combine job sees only
the collected figures, so either each case's plot step writes a small data
file (centreline profile, gauge series) that `collect` also copies, or the
overlay is drawn inside the case by running all algorithms there, as
compare_solvers does.

## Priority 3 — Larger initiatives (weeks to months)

**P3.1 Implement local timestepping (GPU-compatible redesign)**
The 3.1.9 implementation (tag `anuga_core_3.1.9`, `swDE1_domain.c`) used per-edge power-of-2
update frequencies computed by a 3-pass algorithm (lines 482–641). The skip logic in
`_compute_fluxes_central` checked `update_next_flux[ki] != 1` before computing each edge flux.
**Why this design is not GPU-compatible**: uses `already_computed_flux[k,i]` as a per-edge
mutex (race condition under parallel), accumulates a static `local_timestep` across skipped
steps, and processes edge pairs sequentially. The current GPU kernel
(`core_kernels.c:_compute_fluxes_central`) uses `#pragma omp target teams distribute parallel
for` with reductions — incompatible with per-edge skip flags.

**GPU-compatible redesign**: per-triangle activity mask (grouped sub-cycling). Slow triangles
sit out for multiple steps, but the flux loop remains fully data-parallel. Requires: (1) CFL
criterion mapping triangle velocity+size → activity level, (2) grouped sub-cycle scheduler,
(3) conservation validation against analytical solutions. Estimated 2–5× speedup on domains
with large dry/slow areas. The 3.1.9 source in `/home/steve/anuga_core_3.1.9` is the
algorithmic reference.

**P3.0 Debris-flow mode from the Iverson & George (2014) dilatancy model (issue #391, added 2026-09-24)**
The D-Claw physics: five unknowns (h, hu, hv, hm, basal pore pressure p_b), dilatancy
`tan(psi) = m - m_crit/sqrt(1+N)`, Darcy dilation, pore-pressure relaxation on
`alpha mu h^2/(3k)`, and Coulomb-Terzaghi friction `(rho g_z h - p_b) tan(phi + psi)` with a
static threshold. With kappa = 1 the hyperbolic part is the SWE plus two contact fields, so
the flux solver is reused as is; hm and h p_b go on the tracer machinery, the closures in a
fractional-step operator like the sediment one. **The risk is static friction**: holding a
pile at rest on a 30 degree slope until the driving force exceeds the cap, which the
semi-implicit Manning factor cannot do (D-Claw needed a purpose-built Riemann treatment).
Stage 1 (Coulomb friction + static pile + Savage-Hutter dam break) is the go/no-go; stages 2
(m, p_b, consolidation solution, USGS gate-release flume) and 3 (mobilisation, loose/dense
bifurcation) follow. Papers in `~/Projects/Sediment_Transport/References/`; design document
linked from the issue. Does not bear on the suspended-sediment trench discrepancy (#389).

**P3.2 Higher-order spatial reconstruction**
Current extrapolation is linear (second-order smooth, first-order near gradients). Limited
third-order reconstruction (MUSCL-Hancock or ADER) would improve accuracy for long-distance
tsunami propagation. The consolidated `quantity_openmp_ext.pyx` (session 14, H3.1) is the
right place. Requires careful monotonicity limiting.

~~**P3.8 ADER-2 / MUSCL-Hancock fused predict-extrapolate kernel (single extrapolation)**~~ — Done (session 31). Fused C-K predictor + extrapolation into a single kernel pass, eliminating the second extrapolation. **1.75× faster than DE1.**

**P3.3 Improve `fit_interpolate` accuracy and performance** *(partial — sessions 25–26)*
Session 25–26: `Fit.select_alpha()` added with L-curve criterion (20 log-spaced candidates
1e-6 … 100, scipy sparse solves, max-curvature corner detection, fallback to DEFAULT_ALPHA).
`dok_to_csr` added to `fitsmooth_ext.pyx` for non-destructive DOK→CSR conversion.
`alpha='auto'` wired in `Fit.fit()`.  fit.py coverage 78→92%.
Remaining: (2) validation suite against known surfaces, (3) profile and vectorise inner loops
in `interpolate.py` (1200 lines, multiple "DESIGN ISSUES" comments, 82% coverage).

**P3.4 Parallel load-balancing monitoring**
Static METIS decomposition doesn't adapt as the wet front advances in inundation simulations.
Add runtime imbalance reporting (which rank is the bottleneck, imbalance ratio) and explore
dynamic repartitioning. The weak-scaling scripts from session 14 provide the benchmarking
framework.

**P3.5 GPU memory ceiling for large domains**
Current GPU offloading caps at ~2.25M triangles on typical 16 GB VRAM. The quantity memory
reduction work (session 13, ~58% saving) helps but is not sufficient for continental-scale
runs. Options: CUDA Unified Memory (`cudaMallocManaged`) or selective quantity transfer (only
GPU what's needed per sub-step).

~~**P3.6 `anuga_sww_gui` erosion delta-bed view**~~ — Done (session 24). Added `elev_delta` quantity: `elev_delta` property on `SWW_plotter`, `_elev_delta_frame`/`save_elev_delta_frame` methods (RdBu_r colormap, symmetric auto-limits), `_elev_delta_frame_count`, worker entry in `_animate_worker.py`, full GUI wiring in `anuga_sww_gui.py`. 6 new tests.

**P3.7 Streaming SWW reads for very long simulations**
`SWW_plotter` currently reads the full time dimension into memory on load. For very long
simulations (thousands of timesteps at high resolution) this can become the memory bottleneck.
A lazy/chunked reader using NetCDF4 variable slicing would allow the GUI and animate.py to
work without loading all data upfront.

---

## Speculative / Long-term

**S1 ML-fitted friction coefficients** — Replace fixed Manning's n with spatially varying
values trained on observed water levels. Potentially high accuracy gain for urban flood
applications. Requires calibration data and adjoint or ensemble methods.

**S2 Adaptive mesh refinement** — Dynamically refine triangles near the wet/dry front or
structures during simulation. Significant algorithmic complexity (remeshing, quantity
projection, parallel redistribution) but would reduce element count for long-range runs.

**S3 Real-time web visualisation** — Replace the desktop GUI with a WebGL viewer that can
stream frames from a running simulation. Lower barrier for classroom and stakeholder use.
Existing frame-generation pipeline as backend.

**S4 Operator composition / scenario DSL** — Higher-level description language for scenarios
(e.g., "rainstorm at 50 mm/hr for 2 hours, then tidal forcing") that auto-composes
`Rate_operator`, `File_boundary`, etc. Would reduce boilerplate for common operational setups.

---

**Held timestep for DE_ader2** (session 58; see KNOWN_ISSUES "DE_ader2's
steady state depends on dt"). Holding dt until the CFL limit forces it down,
and meeting output times with equal steps, stops DE_ader2's self-sustained
unsteadiness in unstable flow (unseeded roll waves went steady to 1e-13 at
MARGIN 0.7-0.85). Not adopted: it is fragile (MARGIN 0.9, or a raise rule,
re-establishes the loop), costs about 1/MARGIN more steps, and would have to go
into the legacy, unified C and GPU step paths. The prototype is described in
KNOWN_ISSUES. Revisit only if the effect matters in a real application.

## Invalidated suggestions

~~**HDF5/Zarr output format**~~ — ANUGA uses NetCDF4 (HDF5-backed), which has no 2 GB
per-variable size limit. The NetCDF3 classic restriction does not apply. (Invalidated
2026-04-24.)

---

## Summary

| Priority | Total | Remaining | Effort | Biggest payoff |
|----------|-------|-----------|--------|----------------|
| P1 — Quick wins | 8 | 0 ✅ | 1–3 days | All done |
| P2 — Medium | 9 | 1 | 1–2 weeks | Test coverage |
| P3 — Initiatives | 8 | 6 | 1–3 months | Performance, scalability, accuracy |
| Speculative | 4 | 4 | 6+ months | Strategic differentiation |

**Top 3 near-term recommendations:**
1. **P2.6** — Raise fast-suite coverage threshold (currently ~58–59%; next targets in `fit_interpolate/` and `structures/`)
2. **P3.3** — `fit_interpolate` accuracy: validation suite against known surfaces, profile `interpolate.py` inner loops
3. **P3.1** — Local timestepping GPU-compatible redesign (per-triangle activity mask sub-cycling)
