# Key Technical Decisions

Decisions from two sources: **anuga-community** sessions (2026-03-23 to 2026-03-26)
and the **Hydrata** refactor plan (2026-02-28, [REFACTOR_PLAN.md](https://github.com/Hydrata/anuga_core/blob/anuga-4.0-refactor-plan/REFACTOR_PLAN.md)).

Decisions made during Claude sessions that are not obvious from reading the code.
Each entry includes what was decided, why, and when.

---

## EPSG / CRS in Geo_reference (2026-03-26)

### `epsg` is a property, auto-computed for WGS84 UTM

**Decision:** `epsg` is a computed property — for WGS84 UTM it returns `32600+zone`
(north) or `32700+zone` (south) without storing anything. Explicit storage (`_epsg`)
only kicks in when the user supplies a non-auto-computable EPSG (e.g. EPSG:28992).

**Why:** Avoids redundant storage for the overwhelmingly common UTM case, while still
allowing any EPSG to be stored explicitly.

### Non-UTM EPSG: zone stays DEFAULT_ZONE (-1), use `is_located()` to distinguish from wavetank

**Decision:** For national grid EPSG codes (28992, 27700, etc.), `zone` remains -1.
The new `is_located()` method returns `True` when `zone != -1 OR _epsg is not None`,
cleanly separating "real CRS, no UTM zone" from "wavetank / no CRS".

**Why:** `zone = -1` previously meant both "no UTM zone" and "wavetank" — adding a
non-UTM EPSG conflated a located simulation with an unlocated one.

### Use pyproj to populate datum, projection, false_easting, false_northing for non-UTM EPSG

**Decision:** `_populate_from_pyproj()` is called for non-UTM EPSG codes. Uses
`to_cf()` first (works for BNG but not RD New), falls back to `to_dict()` (works for
RD New). The PROJ-string conversion warning from pyproj is suppressed via
`warnings.catch_warnings()`.

**Why:** Without pyproj, `datum='wgs84'` and `projection='UTM'` would be written to
SWW files for RD New / BNG — plainly wrong metadata.
`pyproj` is in the `[data]` optional extras, so it may not always be present; the
code degrades gracefully to leaving the defaults.

### Fixed pre-existing bug in `read_NetCDF`: `self.zone == 'southern'`

**Decision:** Lines 234 and 250 of the original code compared `self.zone` (an `int`)
to the strings `'southern'` and `'northern'` — always `False`. Corrected to
`self.hemisphere == 'southern'` / `'northern'`.

**Why:** The warnings about non-standard false easting/northing were silently never
firing. Fixed as part of the EPSG work.

### `pyproject.toml` optional extras structure (2026-03-26)

**Decision:**
- Core deps: `numpy`, `dill`, `matplotlib`, `netCDF4`, `scipy`, `meshpy`
- `[parallel]`: `mpi4py`, `pymetis`
- `[data]`: `pyproj`, `affine`, `xarray`, `pandas`, `openpyxl`, `tomli` (Py<3.11)
- `[dev]`: `pytest`, `psutil`
- `cartopy` excluded entirely — never actually imported anywhere in the codebase.

**Why:** `pyproject.toml` previously only declared `numpy`, meaning `pip install anuga`
produced a broken package on a clean venv. Confirmed hard deps by tracing top-level
imports triggered by `import anuga`.

---

## Testing

### `--run-fast` skips slow, not the other way around (2026-03-26)

**Decision:** Default `pytest` runs ALL tests (including parallel MPI tests). `--run-fast`
skips tests marked `@pytest.mark.slow` and tests under `anuga/parallel/tests/`.

**Why:** User wanted parallel tests included in the standard CI run. The flag was originally
`--run-slow` (to opt-in to slow tests) but then inverted so the common default is the complete
run and `--run-fast` is the quick-feedback mode.

**Result:** Default ≈163s (~1600 tests); `--run-fast` ≈41s (~1500 tests, 101 skipped).

---

### `pytest.mark.slow` criteria (2026-03-26)

Tests are marked slow if they:
- Live under `anuga/parallel/tests/` — auto-marked by `conftest.py` path check
- Take more than ≈5 seconds individually — manually decorated with `@pytest.mark.slow`

Current slow tests (manually marked):
- `anuga/operators/tests/test_kinematic_viscosity_operator.py` — all 4 tests (module-level mark)
- `anuga/tsunami_source/tests/test_tsunami_okada.py` — all tests (module-level mark)
- `test_culvert_class.py`: `test_that_culvert_runs_rating`, `test_momentum_jet`
- `test_new_culvert_class.py`: `test_that_culvert_runs_rating`, `test_that_culvert_flows_conserves_volume`, `test_momentum_jet`
- `test_shallow_water_domain.py`: `test_inflow_using_flowline` (13s)
- `test_model_tools.py`: `test_create_culvert_bridge_operator`
- `test_csv2sts.py`: `test_run_via_commandline`

---

## Docstring Style

### Use NumPy docstring style, not Sphinx `:param:` style (2026-03-26)

**Decision:** All new docstrings in ANUGA use NumPy style.

**Example:**
```python
Parameters
----------
domain : Domain
    The ANUGA domain to operate on.
rate : float or callable
    Flow rate in m³/s. Negative values extract water.
```

**Not:**
```python
:param domain: The ANUGA domain to operate on.
:type domain: Domain
```

**Why:** User preference; NumPy style is more readable for scientific code and consistent
with NumPy/SciPy ecosystem conventions.

---

## Verbose/Quiet Output

### Triangle library quiet flag (2026-03-26)

**Decision:** `_generateMesh_impl` in `anuga/pmesh/mesh.py` now honours the `verbose`
parameter to choose between 'Q' (quiet) and 'V' (verbose) triangle flags when neither is
already in `self.mode`.

**Why:** `str.find()` returns 0 for a first-position match, making `not self.mode.find('Q')`
True when 'Q' is at position 0. The old conditional `if 'Q' not in self.mode` → add 'V'
was running even with `verbose=False`, causing triangle library output to appear in `pytest -s`.

**Fix location:** `anuga/pmesh/mesh.py` in `_generateMesh_impl`.

---

### Suppress logging in tests that intentionally call verbose=True (2026-03-26)

**Decision:** `test_verbose_does_not_raise` in `test_pmesh_to_mesh.py` wraps the call with
`logging.disable(logging.CRITICAL)` / `logging.disable(logging.NOTSET)` to prevent expected
verbose output from polluting test output.

**Why:** The test's purpose is to confirm no exception is raised — not to check the log output.

---

## Memory Reporting

### Memory display format (2026-03-26)

**Decision:** Memory stats show as `mem=NNNmb` below 1024 MB, and `mem=N.NNgb` above.
- Below 1024 MB: `f'mem={rss_mb:.0f}MB'`
- At or above 1024 MB: `f'mem={rss_mb / 1024:.2f}GB'`

Uses `psutil` if available, falls back to `_VmB('VmRSS:')` on Linux.
Added to `timestepping_statistics()` output in `generic_domain.py`.

---

## API Design

### `basic_mesh_from_mesh_file()` as a function, not a class method (2026-03-26)

**Decision:** Implemented as a standalone function in `basic_mesh.py`, not as `Basic_mesh.from_file()`.

**Why:** User explicitly requested a function, not a class method. Functions are easier to
import and use without knowing the class.

**Returns:** A `Basic_mesh` instance with extra instance attributes:
- `vertex_attributes` — ndarray or None
- `vertex_attribute_titles` — list of strings
- `triangle_tags` — list

---

### `is None` vs `== None` (2026-03-26)

**Decision:** Always use `is None` / `is not None` in this codebase.

**Why:** `== None` can raise `ValueError` for numpy arrays (ambiguous truth value).
`is None` tests identity and always works correctly.

---

## Deprecation Strategy

### camelCase → snake_case in `pmesh/mesh.py` (2026-03-24)

**Decision:** 39 public methods renamed to snake_case. camelCase names kept as wrappers
that emit `DeprecationWarning`. All internal callers updated to snake_case immediately.

**Why:** The module had 49 camelCase alongside 31 snake_case methods — inconsistent and
confusing. Keeping camelCase as deprecated wrappers avoids breaking external user code.

---

### `get_CFL` / `set_CFL` in `generic_domain.py` (2026-03-24)

**Decision:** `set_cfl` promoted to the real method. `set_CFL` / `get_CFL` are deprecated
wrappers with `DeprecationWarning`.

---

## SWW GUI design decisions (2026-04-24)

### Baked overlays vs canvas overlays — two coexisting systems

**Decision:** Show Mesh and Show Elev use two independent rendering paths that coexist:
- **Canvas overlay** (quick preview): drawn in Tkinter pixel-space on top of the imshow `PhotoImage`. Fast, no re-generation. Tracks `_ax` coordinates.
- **Baked overlay** (permanent output): rendered into the matplotlib `Figure` during `save_*_frame` calls in `animate.py` using data-space triangulations at the correct z-order (elev contours z=3, mesh z=4). Permanent in the PNG.

When a quantity is generated with Show Elev/Mesh ticked, `_last_gen_show_elev` / `_last_gen_show_mesh` are set and the canvas overlay methods return early, preventing a double-render.

**Why:** Canvas overlays can only go on top of the pre-rendered opaque imshow; baked overlays integrate correctly at any z-order. Both are useful: canvas for cheap live preview while browsing, baked for output files.

### `_nice_contour_levels` rounds steps to 1/2/5/10/20/50/100 × magnitude

**Decision:** Contour levels are not equally-spaced floats — they are computed by finding the nearest round step value (factor of 1, 2, 5, 10 … × the order of magnitude) such that at most `n × 1.5` levels span the data range.

**Why:** Arbitrary equal-spacing (e.g. 7.3 m intervals) produces ugly tick labels. Round steps like 5 m or 50 m are immediately readable.

### Multi-point timeseries — list, not single triangle; tab10 palette for colour consistency

**Decision:** `_ts_triangles` is a list; picks are appended (deduplicated). Both frame-star overlays (`_update_pick_overlay`) and timeseries lines (`_update_timeseries`) call `self._ts_color(i)` for the same index, guaranteeing identical colours in both panels without any cross-communication.

**Why:** A single picked index is a common case, but users asked for the ability to compare multiple points. The tab10 palette is categorical and perceptually distinct up to 10 points — sufficient for practical use.

### Tabbed UI layout — 3 tabs + always-visible bars

**Decision:** The 9 flat rows of the original UI were reorganised into:
- **Plot tab**: quantity, vmin/vmax, colormap, overlays
- **Generate tab**: output path, DPI, stride, EPSG, basemap
- **Output tab**: save/export/animation/mesh buttons
- Always-visible: file row, Generate Frames bar, playback bar, frame slider, status bar

**Why:** The always-visible bars give access to the most-used controls (generate, play, scrub) without tab switching. Less-used settings (DPI, EPSG, basemap) move to the Generate tab; save/export buttons move to Output. The flat layout became too wide to fit without scrolling as features accumulated.

### `use_basemap=None` sentinel in `_render_and_save_mesh`

**Decision:** The parameter defaults to `None`, which means "inherit from last generation state" (`self._gen_used_basemap and sp.epsg is not None`). An explicit `True` or `False` overrides this.

**Why:** The `_save_mesh` dialog needs to pre-populate the checkbox from the last gen state, but the user can then change it. Passing the checkbox value as an explicit bool ensures the override is honoured. `None` as sentinel keeps the old default behaviour for any caller that doesn't pass the argument.

---

## Hydrata Refactor Plan Decisions

From [REFACTOR_PLAN.md](https://github.com/Hydrata/anuga_core/blob/anuga-4.0-refactor-plan/REFACTOR_PLAN.md) (2026-02-28). These are adopted as guidance for future work.

### Quantity kernels: single source with compile-time OpenMP toggle

**Decision:** Unify `quantity_ext.pyx`, `quantity_ext_openmp.pyx`, `quantity_ext2.pyx` into
a single source file. OpenMP pragmas become no-ops when compiled without `-fopenmp`.

**Why:** The three files share ~90% code, creating a maintenance burden and risk of divergence.

### Parallel operators: MPI awareness in base class, not wrapper files

**Decision:** Eliminate the 5 thin parallel wrapper files. Move `self.parallel` flag into base
operator classes so serial/parallel behaviour is a runtime switch, not a separate class.

**Why:** The wrappers add complexity without meaningful behaviour difference.

### Linting: incremental ruff enforcement, never bulk-format

**Decision:** Use `ruff` (not pyflakes/pep8) for linting. Format only files being modified —
never run a bulk formatter across the whole repo.

**Why:** Bulk formatting creates massive diffs that conflict with upstream PRs and obscure
meaningful changes in git history.

### Testing: pytest over unittest for new tests, convert selectively

**Decision:** New tests written as plain pytest functions. Convert existing `unittest.TestCase`
classes selectively when touching a file, not as a bulk migration.

**Why:** pytest fixtures (especially `tmp_path`) solve the test isolation problem more cleanly
than the current `set_datadir('.')` / `tempfile.mktemp()` pattern.

### Coverage: 65% target enforced in CI, 80% on changed lines via diff-cover

**Decision:** Configure `.coveragerc` with `branch=true, fail_under=65`. Install `diff-cover`
to enforce 80% coverage on lines changed in each PR.

**Why:** A blanket 65% target is achievable; the 80% diff-cover rule ensures new code is
always well-tested without requiring a wholesale coverage lift immediately.

### Fork hygiene: rebase monthly, squash-merge phases, never force-push main

**Decision (Hydrata-specific, relevant for upstream sync):**
- Rebase feature branches against `anuga-community/main` at least monthly
- Squash-merge completed phases rather than accumulating merge commits
- Never bulk-format to keep diffs clean for upstream PRs

**Why:** Hydrata's fork must stay mergeable with upstream. Coordination on
`pyproject.toml` dependency changes is especially important.

## Compute model: per-domain mode + process-level offload/threads (2026-06-13)

Reached while building the `multiprocessor_mode=2` migration (see
`claude/PLAN_default_mode2_cpu.md`). The execution model has **two orthogonal knobs** —
one per-domain, one process-global — because they have genuinely different scopes.

### Per-domain: `set_compute_mode('legacy' | 'unified')`

**Decision:** A domain picks its compute path. `'legacy'` = `multiprocessor_mode=1`
(`sw_domain_openmp_ext` solver + serial-Python operators); `'unified'` =
`multiprocessor_mode=2` (the unified `gpu_ext` C kernels: solver *and* operators).
`set_multiprocessor_mode(1|2)` is kept as a thin int alias (1→legacy, 2→unified).

**Why:** Which solver/extension a domain uses is genuinely per-domain — two domains can
differ with no interaction. `'unified'` was chosen over `'cpu'`/`'gpu'` because those are
*not* per-domain choices (see next).

### Process-global: `set_gpu_offload(bool)`, `set_omp_num_threads(int)`

**Decision:** GPU offload on/off and the OpenMP thread count are **module-level functions**
(`anuga.set_gpu_offload`, `anuga.set_omp_num_threads`), not domain methods.
`Domain.set_omp_num_threads` stays as a backward-compatible wrapper.

**Why:** Both are process-level OpenMP ICVs. One process cannot run the same unified
kernels on a GPU for domain A and on the CPU for domain B; and `omp_set_num_threads`
covers the whole process. Encoding `cpu`/`gpu` as per-domain "modes" mixed the two scopes
and produced a real cross-domain bug (a second domain couldn't re-enable offload). So
`cpu` = `unified` + offload off and `gpu` = `unified` + offload on are *compositions* of
the two knobs, not primitives.

### Offload toggle mechanism

**Decision:** `set_gpu_offload(False)` sets `OMP_TARGET_OFFLOAD=disabled` (the real lever
that forces `omp target` regions onto the host) **and** a process-global C flag
(`g_gpu_offload_enabled`) that makes `gpu_domain_init` choose the host device and the
inlet/culvert operators route there (`gpu_compute_device(GD)`, never `device_id=-1`).
Must be called **before the first `evolve()`/domain build** (OpenMP reads the env at its
first target region; data is mapped at init).

**Why:** The C flag / `omp_set_default_device(host)` alone does NOT stop nvc offloading —
the solver kept running on the GPU (towradgi `-ngo` was 6.35s = GPU speed) until
`OMP_TARGET_OFFLOAD=disabled` was added. The flag is still needed for device-id/data-map
consistency and to fix operators calling `omp_set_default_device(-1)`.

### Two builds, not one binary, for CPU vs GPU

**Decision:** Ship CPU and GPU as **separate builds** selected by the `gpu_offload` meson
option: gcc `-Dgpu_offload=false` (fast CPU multicore via `CPU_ONLY_MODE` → `omp parallel
for`) and nvc `-Dgpu_offload=true` (fast GPU). The distribution default is
`gpu_offload=false`. `set_gpu_offload(False)` / `-ngo` on a GPU build is for **correctness
A/B only** (results are bit-identical to GPU), not performance.

**Why:** A single nvc `-mp=gpu,multicore` binary cannot be fast on both — NVHPC's OpenMP
target *host fallback runs single-threaded* (a confirmed, documented NVHPC limitation, not
an ANUGA bug — see `claude/KNOWN_ISSUES.md` for the microbenchmark and citations). No
runtime knob re-routes the host execution to the fast multicore variant.

### CLI

**Decision:** Standard arg parser (`anuga.get_args`) gains `-nt/--omp_num_threads`,
`-go/--gpu_offload`, and `-ngo/--no-gpu_offload` (tri-state; unset = follow the build).
Scripts apply them right after parsing, before building the domain.

### Forcing-function classes → operators: confirm removal, add mode-2 safety (2026-06-14)

**Context:** `forcing.py`'s `Wind_stress`, `Rainfall`, `Inflow`, `Barometric_pressure`
were already **deprecated** in session 25 (`DeprecationWarning`, suppressed in tests via
`pyproject.toml`). Operator equivalents exist: `Rate_operator.rainfall()`/`inflow()`,
`Wind_stress_operator`, `Barometric_pressure_operator`. Manning friction
(`manning_friction_semi_implicit`) stays a forcing term and is NOT deprecated.

**Decision:** Mode 2 confirms the direction toward **removal** (the next phase after
deprecation). Interim safety warning added now; remove after a release (see FUTURE_WORK).
Also verify the `_fast` variants (`Wind_stress_fast`, `Barometric_pressure_fast`) carry
the same deprecation.

**Why:** Forcing terms are evaluated *inside* the timestep (`compute_forcing_terms`).
The mode-2 C step loop does forcing in C and only handles Manning, so Python forcing
terms are **silently skipped** in mode 2 (surfaced by
`test_rainfall_forcing_with_evolve_1`). Honouring them would require a
`sync_from_device → f() → sync_to_device` every step — the serial-Python-per-step cost
the mode-2 work exists to eliminate. Operators run via `apply_fractional_steps`
(fractional-step splitting), which mode 2 already supports — C kernels where they exist
(rate/inlet/culvert) and a clean one-sync Python fallback where they don't
(wind/barometric). So operators are the architecturally correct path for mode 2, and the
forcing classes are redundant.

### Mode-2 testing: pin white-box host-state tests to legacy (2026-06-17)

**Context:** Running `anuga.shallow_water` under the unified default *on a GPU build*
(via `anuga_run_isolated_tests -cm unified`) failed 11 tests that pass on the CPU build.
On a CPU build mode-2 device memory *is* host memory (`CPU_ONLY_MODE` stubs), so the
divergence only appears with real offload.

**Decision:** Pin such tests to legacy with `domain.set_compute_mode('legacy')` rather
than teaching mode-2 to sync or loosening tolerances. Two kinds qualify: (1) **white-box**
tests that call `compute_forcing_terms()` / `compute_fluxes()` and assert on the host
`semi_implicit_update` / `explicit_update` arrays (mode-2 computes those on-device and
never syncs back → stale zeros); (2) **numerical reference/snapshot** tests recorded under
legacy that diverge at ~1e-6 from mode-2's reduction/eval order. End-to-end evolve tests
are NOT pinned — they pass in both modes. Documented as a convention in `CONVENTIONS.md`.

**Why:** These probe mode-1 internals or a legacy-recorded baseline, not behaviour that
mode 2 must reproduce bit-for-bit; mode-1-vs-mode-2 equivalence is covered separately by
`test_DE_gpu_omp.py`. Adding per-step device syncs purely to satisfy a white-box assertion
would reintroduce exactly the cost mode 2 exists to remove. The pin is a no-op for the
distribution-default legacy path.

**Caveats:** (1) Manning must stay in-step (semi-implicit, stability) — the
`forcing_terms` mechanism is kept for it; only the optional classes are deprecated.
(2) Forcing-as-forcing-term (in-step) vs operator (end-of-step fractional split) are NOT
bit-identical; for slowly-varying rain/wind/barometric the difference is negligible and
operators are the modern standard, but it is a documented numerics change.

**Interim safety (done 2026-06-14):** until the classes are removed, mode 2 warns once
(`Domain._warn_unsupported_mode2_forcing`) when `forcing_terms` contains anything other
than Manning, turning the silent skip into a loud, actionable message.

## Partition file I/O

### Domain partitions stay pickle; mesh partitions stay NetCDF (2026-07-12)

**Context:** Saving the partition files dominates offline partitioning at scale — a
173M-triangle / 18,400-partition run measured ~1,000 s to load, ~4,000 s to partition,
and **~37,000 s to save** (~88% of the job). While parallelising the write
(`num_workers`) and collapsing the domain dump to one file per partition, the obvious
question was whether both dumps should share one format: should the domain dump move to
NetCDF (like the mesh dump), or the mesh dump to pickle (like the domain dump)?

**Decision:** Keep **two formats**. `sequential_distribute_dump` (full domain: mesh +
all quantities) writes **pickle**; `sequential_mesh_dump` (mesh topology + halo only)
writes **NetCDF**. A **hybrid** was prototyped — bulk arrays as NetCDF variables plus a
small pickled config blob, still one file per partition — and **rejected**.

**Why:** Measured on real ~10k-triangle partitions (local SSD):

| Path | write | read | size |
|------|------:|-----:|-----:|
| domain: single pickle (chosen) | **0.36 ms** | 0.11 ms | 868 KiB |
| domain: hybrid NetCDF + blob | 2.09 ms | 0.98 ms | 881 KiB |
| mesh: NetCDF (chosen) | 2.57 ms* | 1.20 ms | **294 KiB** |
| mesh: pickle | 0.28 ms | 0.08 ms | 465 KiB |

The two dumps do genuinely different jobs:

- **Domain = fast internal cache.** It carries arbitrary Python state (boundary map with
  condition objects, `Geo_reference`, domain flags) that NetCDF cannot hold without an
  opaque pickled blob — so NetCDF buys little while costing ~6× on write, which *is* the
  bottleneck. Size is ~equal. The hybrid round-trips exactly but is still ~6× slower to
  write, because the bulk arrays still go through the NetCDF/HDF5 container.
- **Mesh = portable, reusable preprocessing artifact** (partition once, run many
  scenarios). Its content is **100% arrays**, so NetCDF is fully self-describing
  (`ncdump -h`) with *no* opaque blob, is **1.6× smaller** than pickle, and is safe to
  load — pickle is an arbitrary-code-execution risk for shared/archived files. Its slower
  write is amortised over many reuses and is now parallelisable via `num_workers`.

\* The mesh NetCDF write was **9.45 ms/file** until profiling showed
`_write_mesh_partition` was writing the `boundary_tag` string variable **one row per
boundary edge** (one HDF5 call each). Vectorising to a single write took it to 2.57 ms
(~3.7×, byte-identical output). The apparent "NetCDF is 34× slower than pickle" gap was
mostly that bug, not NetCDF — which is exactly why the fix was to **optimise the writer,
not change the format**. Lesson: benchmark the *implementation* before blaming the format.

**Caveats:** (1) The domain dump's `single_file=True` (new default) inlines points,
triangles and quantities into the one pickle, cutting `3 + N_quantities` files per
partition to 1 (~147k → 18.4k files at 18,400 ranks — the dominant cost on
metadata-bound Lustre/GPFS). `sequential_distribute_load` detects a filename string vs an
inline array, so legacy multi-file dumps still load and `single_file=False` restores the
old layout. (2) Pickle is Python-only and version-fragile; that is acceptable for a
regenerable internal cache, but is precisely why the *mesh* artifact — which users may
inspect, archive and share — was kept on NetCDF.

---

## Structure inlets: level the stage, don't flatten the depth (issue #229)

`Structure_operator` wrote its transfer back as a **uniform depth** across each inlet
(`Inlet.set_depths()` → `stage = elevation + d̄`). Volume was exact, but on a sloping bed
that tilts a level water surface onto the bed, disturbing a lake at rest by half the bed
elevation range across the inlet, every timestep, at zero discharge (3.4e-02 m on a 1/50
bed) — in every compute mode.

**Two models were tried:**

1. **Surface shift** (`e899973f`, reverted in `f4e2ab11`): move every stage by
   `new_avg_depth − old_avg_depth`, clamp at the bed. Well balanced, but it removed the
   smoothing the culvert feedback loop relies on: discharging into an initially **dry**
   inlet, the shallow-water dynamics inside the inlet build a ridge, the enquiry reading
   swings, and the discharge oscillates and stalls
   (`test_pipe_culvert_cpu_gpu_velocity_match` failed on every CI platform; with the
   momentum write it even collapsed the global timestep).

2. **Stage leveling** (this decision): apply the volume change with the "well-mixed
   reservoir" model the operators assume — water finds its level. Adding volume raises
   the *lowest* stages to a common level; removing volume lowers the *highest* stages to
   a common level, clamping each cell at its bed. Exactly zero transfer is a bit-exact
   no-op (the lake-at-rest case). On a **flat bed** filling a level or dry inlet
   reproduces the old uniform-depth write exactly, which is why the flat-bed culvert
   tests are untouched by construction — the property attempt 1 lacked. This is the same
   model `Inlet.set_stages_evenly()` already uses for inlet-operator inflow.

**Momentum**: the physics produces ONE momentum per inlet. The old uniform write was
consistent with uniform depth; over leveled (non-uniform) depths a uniform momentum
gives a nearly-dry cell an enormous velocity and collapses the timestep. So
`Inlet.set_average_momenta()` writes it **depth-weighted** (`m·dᵢ/d̄`): the area-weighted
average momentum is exactly `m`, and the *velocity* field is uniform (`m/d̄`). On uniform
depth it reduces to the old write.

**Level-finding is bisection** (100 iterations → last bit of a double) on the monotone
volume(L) function — no sort, so the same loop runs in numpy
(`level_stages_to_average()`, `anuga/structures/inlet.py`) and inside one GPU team
(`culvert_level_inlet_surface()`, `gpu_culvert_operator.c`). The two must change
together. In parallel each rank levels its own share by the same per-area volume change
(global average passed in); the seam across a partitioned inlet is the same order as the
old write's.

**Non-idempotence**: the old `stage = bed + depth` write could be applied twice
harmlessly; leveling cannot. The mode-2 MPI path therefore does not scatter in PHASE 3b —
the batched PHASE 4 scatter already covers a non-master rank's own inlet.

### Forcing-function classes removed; file-driven fields kept via the operators (2026-09-18)

**Context:** P2.10 — the deprecated classes were the only way to drive wind /
pressure from an sww or tms file (`use_coordinates=False`, the `_fast` variants).
`Barometric_pressure_operator` explicitly refused `use_coordinates=False`.

**Decision:** Remove the classes (4.1 breaking change, announced in 4.0.0) but
carry the file pathway over rather than drop it: both operators take
`use_coordinates=False` with a `file_function` (centroids for wind, nodes for
pressure) and evaluate the file's precomputed time series for all points in one
vectorised step (`evaluate_file_function_all_points`), which is what the `_fast`
classes were hand-rolling. Rainfall/inflow tests assert on stage deltas after a
single operator call (`domain.timestep = 1.0; op()`), not on `explicit_update`,
so they no longer need a legacy-mode pin.

### Slope for the depth-slope / energy-slope closures comes from centroids, not edges (2026-09-18)

**Context:** `core_tau_b_over_rho` took `S` for `[T-7]`/`[T-7e]` from the bed /
stage EDGE values by the divergence theorem. The DE extrapolation rewrites those
every step as limited stage minus limited height, so the erosion stress was
whatever the hydrodynamic limiter left: zero along reflective walls (the ghost
mirrors the interior), three-quarters at a kink. Found by the sediment_erosion
validation case, whose wall cells never eroded.

**Decision:** least-squares gradient of the CENTROID values over the cell and
its neighbours, one-sided at boundaries (1D projection when the neighbours are
collinear). A plane gives the exact slope in every cell; a slope discontinuity
is smeared over one cell either side, which is inherent to any centroid stencil
and acceptable. Computed in a pass of its own into `sediment_slope_work` (n),
because the source loop writes the bed it would otherwise be reading through
its neighbours: done inline, mode 1 and mode 2 disagreed at 1e-3 from the first
step (`test_the_depth_slope_shear_closure_agrees` caught it). The bedload kernel
already had a read-only first pass, so it takes the slope inline. Not done: reading the user's original elevation vertex values
(would not track an evolving bed). Known and documented: `'depth_slope'` with
bed evolution on is self-amplifying (erosion -> rougher bed -> steeper local
slopes -> more erosion); measured 0.8 m of scour in 60 s on a 1e-3 slope in
still water. That is the closure, not the gradient estimator -- the old edge
slope diverged too -- and it is why the validation case holds the bed fixed.

### GPU mode-2 time-varying boundaries stay on the Python-orchestrated loop (2026-09-19)

**Context:** the single-call C RK loop pushes Python-evaluated boundaries to the
device once per step, so RK2/RK3 substeps reused the step-start value (4 mm
stage error on a rising tide, issue #170). PR #171 routed such domains to the
Python-orchestrated GPU loop, which refreshes per substep and bit-matches mode 1
at a measured cost of ≤4%.

**Decision:** keep that routing as the design ("option B") and close #170.
"Option A" (evaluate the boundary at each substep time in Python and pass a
per-substep buffer into `evolve_one_rk*_step_gpu`) would buy back a few percent
at the price of a device-mapped buffer and changes to every step function, and
cannot be validated on real GPU hardware in CI until #333 has a runner.

### Bounds for the depth-slope closure: a per-cell cap and a frozen slope (2026-09-20)

**Context:** `'depth_slope'` reads S from the evolving bed and feeds on the
roughness uneven erosion creates (KNOWN_ISSUES). anugaSed contains it with an
undocumented domain-global clamp `S <- min(S, mean(S)/2)` that also divides by a
mean cell size.

**Decision:** two explicit, documented bounds on `set_shear_closure` rather than
reproducing that clamp: `max_slope` caps S per cell (both slope closures), and
`freeze_slope=True` (depth-slope only) takes S from the bed at setup and keeps it
for the run. The frozen slope is computed on the host (`bed_slope_magnitude`, the
kernel's least-squares estimator in numpy) into `sediment_slope_work`, which is
now mapped `to:` the device so every interface rebuild re-uploads it; the kernel
skips its slope pre-pass when frozen. Computed eagerly at `set_shear_closure`
and again when the first grain size allocates the array (the later wins), so no
device sync is needed mid-run. The sediment_erosion validation case now evolves
its bed under a frozen slope against a per-cell ODE with depth feedback.

