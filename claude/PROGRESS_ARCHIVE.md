# ANUGA Progress Archive

Historical record of completed work. Active tracking: `claude/PROGRESS.md`.

---

## Code Improvement Actions (completed items)

Source: `docs/code_improvement_actions.md`
Generated: 2026-03-23

### Priority 1 — Quick wins ✅ Complete

#### 1.1 Fix mutable default arguments (~43 functions)

- [x] `anuga/caching/caching.py:145` *(2026-03-24)*
- [x] `anuga/file/sww.py:535` *(2026-03-24)*
- [x] `anuga/parallel/parallel_boyd_box_operator.py:22` *(2026-03-24)*
- [x] `anuga/abstract_2d_finite_volumes/ermapper_grids.py:8,88,203` *(2026-03-24)*
- [x] Full repo audit — also fixed parallel_structure_operator, parallel_boyd_pipe_operator, parallel_weir_orifice_trapezoid_operator, parallel_internal_boundary_operator, parallel_operator_factory, riverwall, util.py *(2026-03-24)*

#### 1.2 Replace bare `except:` with specific exception types

- [x] `anuga/utilities/system_tools.py` — already OK *(2026-03-24)*
- [x] `anuga/shallow_water/boundaries.py` — already OK *(2026-03-24)*
- [x] `anuga/caching/caching.py` — already OK *(2026-03-24)*
- [x] `anuga/abstract_2d_finite_volumes/tests/test_quantity.py` — already OK *(2026-03-24)*
- [x] `anuga/abstract_2d_finite_volumes/tests/test_generic_domain.py` — already OK *(2026-03-24)*
- [x] `anuga/file_conversion/dem2pts.py` — already OK *(2026-03-24)*

#### 1.3 Convert file operations to use `with` statements

- [x] `anuga/file/csv_file.py:47,196,206,216,224` *(2026-03-24)*
- [x] `anuga/file/ungenerate.py:16` *(2026-03-24)*
- [ ] `anuga/file/urs.py:29` — intentionally skipped: file handle stored as `self.mux_file` for iterator lifecycle
- [x] `anuga/utilities/system_tools.py:29` *(2026-03-24)*
- [x] Audit `anuga/file/` for remaining bare `open()` calls *(2026-04-03)*

#### 1.4 Fix invalid escape sequences in docstrings

- [x] `anuga/utilities/norms.py:15` *(2026-03-24)*
- [x] `python -W error::DeprecationWarning -c "import anuga"` — clean *(2026-03-24)*

#### 1.5 Delete large commented-out dead code

- [x] `anuga/file_conversion/dem2pts.py:164–281` — 118-line pre-vectorisation loop deleted *(2026-03-24)*
- [x] `anuga/abstract_2d_finite_volumes/neighbour_mesh.py:615–668` — 53-line disabled block deleted *(2026-03-24)*
- [x] Grep for large legacy comment blocks in `shallow_water/` and `operators/` *(2026-04-03)*

### Priority 2 — Correctness and stability ✅ Complete

- [x] 2.1 Fix silent error suppression in `set_quantity.py` — documented expected ValueError *(2026-03-24)*
- [x] 2.2 Log xarray import failures in `rate_operators.py` — `log.debug(...)` *(2026-03-24)*
- [x] 2.3 Address FIXME items — `boyd_box_operator.py`, `fit.py`, `polygon.py`, `rate_operators.py` *(2026-03-24)*

### Priority 3 — Test coverage (completed)

#### 3.1 Add tests for untested operator classes ✅ Complete

- [x] `Bed_shear_erosion_operator`, `Circular_erosion_operator`, `Flat_slice_erosion_operator`, `Flat_fill_slice_erosion_operator` *(2026-03-24)*
- [x] `Collect_max_quantities_operator`, `Collect_max_stage_operator` — `test_collect_operators.py` *(2026-03-24)*
- [x] `Elliptic_operator` — `test_elliptic_operator.py` *(2026-03-24)*
- [x] `Circular_rate_operator`, `Circular_set_quantity_operator`, `Circular_set_stage_operator` *(2026-03-24)*

#### 3.2 Add tests for untested structure classes ✅ Complete

- [x] `Structure_operator` base class — `test_structure_operator.py` *(2026-03-24)*
- [x] `Internal_boundary_operator` — `test_internal_boundary_operator.py` *(2026-03-24)*
- [x] `RiverWall` — `Test_riverwall_notebook` class (5 tests) *(2026-04-13, commit a62e9c96)*

### Priority 4 — API and code quality (completed)

- [x] 4.2 Standardise naming in `pmesh/mesh.py` — 39 methods renamed; camelCase kept as deprecated wrappers *(2026-03-24)*
- [x] 4.3 Deprecate camelCase `get_CFL`/`set_CFL` in `generic_domain.py` *(2026-03-24)*
- [x] 4.4 Add `__all__` to `anuga/__init__.py` and sub-package `__init__.py` files *(2026-03-24)*

### Priority 5 — Performance (completed)

- [x] 5.1 Vectorise loops — `fit.py:598`, `csv_file.py:136`, `util.py:786` *(2026-03-24)*

### Priority 6 — Documentation improvements ✅ Complete

- [x] 6.1 `boyd_box_operator.py`, `boyd_pipe_operator.py`, `weir_orifice_trapezoid_operator.py` — full NumPy-style docstrings *(2026-03-24)*
- [x] 6.2 `rate_operators.py`, `erosion_operators.py` — Returns sections added *(2026-03-24)*

---

## Documentation Improvement Actions ✅ All 20 complete

Source: `docs/doc_improvement_actions.md` — Generated: 2026-03-23

| # | Item | Done |
|---|------|------|
| 1 | Fill out `visualisation/use_domain_plotter.rst` | 2026-03-23 |
| 2 | Fix `reference/index.rst` navigation | 2026-03-23 |
| 3 | Fix `anuga_user_manual/version.txt` stale SVN variables | 2026-03-23 |
| 4 | Add `setup_anuga_script/checkpointing.rst` | 2026-03-23 |
| 5 | Add `reference/file_formats.rst` | 2026-03-23 |
| 6 | Add `troubleshooting.rst` | 2026-03-23 |
| 7 | Expand `setup_anuga_script/boundaries.rst` | 2026-03-23 |
| 8 | Add comparison table to `setup_anuga_script/operators.rst` | 2026-03-23 |
| 9 | Add descriptions to `examples/index.rst` notebooks | 2026-03-23 |
| 10 | Add MPI section to `install_anuga_developers.rst` | 2026-03-23 |
| 11 | Clarify OpenMP support in `use_parallel_openmp.rst` | 2026-03-23 |
| 12 | Soften QGIS version in `use_qgis.rst` | 2026-03-23 |
| 13 | Add parallel decision guide to `parallel/index.rst` | 2026-03-23 |
| 14 | Add annotated TOML example to `toml_scenario/index.rst` | 2026-03-23 |
| 15 | Add GPU/`multiprocessor_mode=2` note in parallel docs | 2026-03-23 |
| 16 | Standardise quantity names in `initial_conditions.rst` | 2026-03-23 |
| 17 | Reconcile Python version statements across install docs | 2026-03-23 |
| 18 | Port mathematical background into Sphinx | 2026-03-23 |
| 19 | Add cross-references from RST pages to user manual | 2026-03-23 |
| 20 | Add `reference/validation.rst` | 2026-03-23 |

---

## Additional Enhancements ✅ All 57 complete

| Item | Files | Done |
|------|-------|------|
| Suppress triangle library verbose output in pytest | `anuga/pmesh/mesh.py` | 2026-03-26 |
| Suppress General_mesh logging in test | `anuga/abstract_2d_finite_volumes/tests/test_pmesh_to_mesh.py` | 2026-03-26 |
| Replace `print_timestepping_statistics()` calls in tests with `pass` | `anuga/shallow_water/tests/test_sw_domain_openmp.py` | 2026-03-26 |
| Add `memory_stats()` and `print_memory_stats()` | `anuga/utilities/system_tools.py` | 2026-03-26 |
| Add memory usage to `timestepping_statistics()` output | `anuga/abstract_2d_finite_volumes/generic_domain.py` | 2026-03-26 |
| Export `memory_stats`, `print_memory_stats` from `anuga` | `anuga/__init__.py` | 2026-03-26 |
| Export `distribute_basic_mesh`, `distribute_basic_mesh_collaborative` from `anuga` | `anuga/__init__.py` | 2026-03-26 |
| Add `basic_mesh_from_mesh_file()` factory function | `anuga/abstract_2d_finite_volumes/basic_mesh.py` | 2026-03-26 |
| Export `basic_mesh_from_mesh_file` from `anuga` | `anuga/__init__.py` | 2026-03-26 |
| Fast/slow test infrastructure (`--run-fast` flag, `@pytest.mark.slow`) | `conftest.py`, `pyproject.toml` | 2026-03-26 |
| Mark 10 slow tests across 5 test files | Various test files | 2026-03-26 |
| Document `--run-fast` in developer install docs | `docs/source/installation/install_anuga_developers.rst` | 2026-03-26 |
| Update `CLAUDE.md` with `--run-fast` and slow marker info | `CLAUDE.md` | 2026-03-26 |
| Declare missing runtime deps in `pyproject.toml`; add `[parallel]`, `[data]`, `[dev]` extras; fix classifiers | `pyproject.toml` | 2026-03-26 |
| Add EPSG/CRS support to `Geo_reference` — `epsg` property, `is_located()`, non-UTM support via pyproj, `write/read_NetCDF`, fix pre-existing zone/hemisphere bug in `read_NetCDF` | `anuga/coordinate_transforms/geo_reference.py` | 2026-03-26 |
| 23 new tests for EPSG/CRS behaviour | `anuga/coordinate_transforms/tests/test_geo_reference.py` | 2026-03-26 |
| New CRS documentation page; `Geo_reference` API reference; cross-references | `docs/source/setup_anuga_script/coordinate_reference.rst`, `docs/source/reference/anuga.Geo_reference.rst` | 2026-03-26 |
| Create `claude/` session-continuity directory | `claude/` | 2026-03-26 |
| Incorporate Hydrata REFACTOR_PLAN.md into claude/ docs | `claude/PROGRESS.md`, `DECISIONS.md`, `KNOWN_ISSUES.md` | 2026-03-26 |
| Fix `sww_merge` not propagating `hemisphere`, `epsg`, and `timezone` | `anuga/utilities/sww_merge.py` | 2026-03-28 |
| Add `sww2vtu` converter — SWW → VTU + PVD for ParaView | `anuga/file_conversion/sww2vtu.py` | 2026-03-28 |
| GPU verbose flag — suppresses C printf output during pytest | `gpu_domain.h`, `gpu_domain_core.c`, `gpu_boundaries.c`, `sw_domain_gpu_ext.pyx` | 2026-04-01 |
| Fix pyproj DeprecationWarning for 1-element arrays (NumPy ≥ 2.0) | `redfearn.py`, `tif2point_values.py` | 2026-04-01 |
| Fix ReadTheDocs shallow-clone version showing `0.0.0+unknown` | `.readthedocs.yaml` | 2026-04-02 |
| Vectorise `get_flow_through_cross_section` | `anuga/shallow_water/shallow_water_domain.py` | 2026-04-03 |
| Add ruff linting config and fix all genuine violations | `pyproject.toml`, various `.py` files | 2026-04-03 |
| L1-L4 logging refactor: `TeeStream`, lazy log file, `set_logfile()`, `log.verbose()`, `log.file_only()` | `anuga/utilities/log.py`, `anuga/scenario/prepare_data.py`, scripts | 2026-04-05 |
| Add logging documentation page | `docs/source/setup_anuga_script/logging.rst` | 2026-04-05 |
| Archive CuPy/CUDA files out of `anuga/shallow_water/` into `archive/cupy_cuda/` | `archive/cupy_cuda/` | 2026-04-05 |
| Fix `test_sww2csv_multiple_files` stale-file pollution | `anuga/abstract_2d_finite_volumes/tests/test_gauge.py` | 2026-04-05 |
| CI: add `pytest-regressions` to all 13 conda environment YMLs | `environments/environment_*.yml` | 2026-04-05 |
| CI: drop Python 3.8/3.9; fix `list \| np.ndarray` PEP-604 annotation | `.github/workflows/conda-setup.yml`, `pyproject.toml` | 2026-04-05 |
| Fix NPY002 test recalibration | `anuga/geospatial_data/tests/test_geospatial_data.py` | 2026-04-05 |
| Propagate v3.3.0, v3.3.1, v3.3.2 tags/releases to GeoscienceAustralia remote | `ga` remote | 2026-04-05 |
| L5: 715 `log.critical()` → `log.info()` across 70+ production files | 70+ `anuga/**/*.py` | 2026-04-06 |
| Drop Python 3.9 | `pyproject.toml`, `.github/workflows/conda-setup.yml` | 2026-04-06 |
| **anuga_animate_sww_gui** — parallel frames, zoom, elev quantity, terrain colormap, Sphinx docs | `scripts/anuga_animate_sww_gui.py`, `anuga/utilities/animate.py`, `_animate_worker.py` | 2026-04-21 |
| **anuga_sww_gui** — Baked overlay generation (elev contours + mesh baked into PNG frames) | `scripts/anuga_sww_gui.py`, `anuga/utilities/animate.py`, `_animate_worker.py` | 2026-04-24 |
| **anuga_sww_gui** — Multi-point timeseries picking, tab10 palette, legend, CSV export, Clear button | `scripts/anuga_sww_gui.py` | 2026-04-24 |
| **anuga_sww_gui** — Save Frame / Export Frame time-selection dialog | `scripts/anuga_sww_gui.py` | 2026-04-24 |
| **anuga_sww_gui** — 3-tab ttk.Notebook UI reorganisation | `scripts/anuga_sww_gui.py` | 2026-04-24 |
| **anuga_sww_gui** — Basemap checkbox for mesh viewer and save dialog | `scripts/anuga_sww_gui.py` | 2026-04-24 |
| **anuga_sww_gui** — Updated in-app help and Sphinx RST for all new features; fresh screenshots | `scripts/anuga_sww_gui.py`, `docs/source/visualisation/use_sww_gui.rst` | 2026-04-24 |
| **P2.3 `create_riverwalls` refactor** — `_validate_riverwall_inputs`, `_match_edges_to_segments`, `_build_hydraulic_properties`; `create_riverwalls` reduced to ~50-line orchestrator | `anuga/structures/riverwall.py` | 2026-04-25 |
| **P2.2 `Generic_Domain.__init__` refactor** — `_init_mesh`, `_init_quantities`, `_init_parallel`, `_init_timestepping`; `__init__` reduced to ~25 lines | `anuga/abstract_2d_finite_volumes/generic_domain.py` | 2026-04-25 |
| **`test_shallow_water_domain.py` cleanup** — removed duplicate/unused imports, 66 debug prints, dead skeleton; net −101 lines | `anuga/shallow_water/tests/test_shallow_water_domain.py` | 2026-04-25 |
| **Split `test_shallow_water_domain.py` into 5 files** — `test_flux.py` (15), `test_boundaries_sw.py` (9), `test_extrapolation_sw.py` (14), `test_physics_sw.py` (21); registered in meson.build | `anuga/shallow_water/tests/` | 2026-04-25 |
| **Fix 383 pytest warnings** — `np.array(netcdf_var)` → `netcdf_var[:]` in animate.py; zero-timestep guard in rate_operators.py; message-based filterwarnings for 5 deprecated forcing classes | `anuga/utilities/animate.py`, `anuga/operators/rate_operators.py`, `pyproject.toml` | 2026-04-25 |
| **anuga_sww_gui** — Basemap wet/dry smooth rendering: `LinearTriInterpolator` + `imshow` path in `_animated_frame`; zoom-aware grid (600 px across visible span, min 200 px) | `anuga/utilities/animate.py` | 2026-05-01 |
| **anuga_sww_gui** — Live x/y/triangle hover readout: status bar restructured; `_trifinder` cached; `_on_hover` sets right-side `_coord_var` StringVar | `scripts/anuga_sww_gui.py` | 2026-05-01 |
| **anuga_sww_gui** — Zoom clipping for mesh/elevation overlays: vertex-based exclusion (any vertex outside xlim/ylim); contour levels computed from visible elevation range | `scripts/anuga_sww_gui.py` | 2026-05-01 |
| **anuga_sww_gui** — `↻ Regenerate Frames` button prompt after zoom set/reset; reverts to `Generate Frames` at startup and when generation starts | `scripts/anuga_sww_gui.py` | 2026-05-01 |
| **anuga_sww_gui** — 14 new CLI parameters: `--vmin`, `--vmax`, `--cmap`, `--cmap-reverse`, `--mindepth`, `--flat-view`, `--outdir`, `--dpi`, `--stride`, `--alpha`, `--epsg`, `--basemap`/`--no-basemap`, `--basemap-provider` | `scripts/anuga_sww_gui.py` | 2026-05-01 |
| **anuga_sww_gui** — TOML config file support: `--config FILE.toml` CLI arg; Save Config / Load Config buttons; sectioned `[render]`/`[generate]`/`[file]` layout; CLI args override config | `scripts/anuga_sww_gui.py` | 2026-05-01 |
| **anuga_sww_gui** — Cross-section discharge panel: pick two points, compute Q(t) via `get_flow_through_cross_section`; cyan overlay markers + line on animation; vertical cursor synced with playback | `scripts/anuga_sww_gui.py` | 2026-05-01 |
| **anuga_sww_gui** — Cross-section panel repositioned below animation canvas (matching timeseries panel layout) | `scripts/anuga_sww_gui.py` | 2026-05-01 |

---

## Hydrata Refactor Plan ✅ Phases 0–4 complete

Source: [Hydrata/anuga_core REFACTOR_PLAN.md](https://github.com/Hydrata/anuga_core/blob/anuga-4.0-refactor-plan/REFACTOR_PLAN.md)

### Phase 0 — Test Infrastructure ✅

- [x] **0.1** Fix test isolation — `tempfile.mktemp` → `mkstemp`, `set_datadir('.')` → `mkdtemp()` *(2026-04-03)*
- [x] **0.2** Add test markers — `@pytest.mark.slow`, `--run-fast` flag *(2026-03-26)*
- [x] **0.3** Golden-master snapshots — 6 `pytest-regressions` tests *(2026-04-04)*
- [x] **0.4** Coverage baseline — `.coveragerc` with `branch=true, fail_under=55` *(2026-04-03)*
- [x] **0.5** CI test matrix — PRs: `--run-fast`; pushes to main/develop: full suite *(2026-04-03)*

### Phase 1 — Dependency Consolidation ✅

- [x] **1.1** Declare runtime deps in `pyproject.toml`; add `[parallel]`, `[data]`, `[dev]` extras *(2026-03-26)*
- [x] **1.2** Remove dead deps — GDAL fully removed; NPY002 fixes *(2026-04-04)*
- [x] **1.3** Delete `setup.py` — already absent *(2026-03-26)*
- [x] **1.4** Fix classifiers *(2026-03-26)*

### Phase 2 — Linting & Code Quality ✅

- [x] **2.1** Add ruff configuration *(2026-04-03)*
- [x] **2.2** Pre-commit hooks — `.pre-commit-config.yaml` with ruff *(2026-04-03)*
- [x] **2.3** CI enforcement — `.github/workflows/lint.yml` *(2026-04-03)*

### Phase 3 — Code Deduplication ✅

- [x] **3.1** Unify quantity kernels — single `quantity_openmp_ext.pyx` *(commit 5c191dc7)*
- [x] **3.2** Consolidate parallel operator wrappers — 3 helpers extracted; −125 lines net *(2026-04-12)*
- [x] **3.3** Merge duplicate culvert classes *(merged via PR #118)*
- [x] **3.4** Clean up `system_tools.py` — 335 lines removed *(2026-04-13, commit f083ad29)*

### Phase 4 — Expanded Test Coverage ✅

- [x] **4.1** Modernise test patterns — deferred to opportunistic pass
- [x] **4.2** Integrate validation tests — 33 `validate_*.py` scripts *(2026-04-10)*
- [x] **4.3** Coverage targets — extended `.coveragerc` omit rules; `fail_under=52` *(2026-04-10)*
- [x] **4.4** Push coverage to 63% — systematic new-test pass across 10 files *(2026-04-13)*
- [x] **4.5** Scenario module tests — 3 new test files, 33 tests *(2026-04-14)*

---

## Riverwall Throughflow ✅ Complete

Full plan: `claude/archive/RIVERWALL_THROUGHFLOW_PLAN.md`

- [x] **RW1** Add `Cd_through` to `hydraulic_variable_names` and `default_riverwallPar` *(2026-04-04)*
- [x] **RW2** Add `gpu_adjust_edgeflux_with_throughflow()` to `gpu_device_helpers.h` *(2026-04-04)*
- [x] **RW3** Call new function in `core_kernels.c` after existing weir call *(2026-04-04)*
- [x] **RW4** No separate CPU path needed — `core_kernels.c` shared via include *(2026-04-04)*
- [x] **RW5** Tests: 6 new tests *(2026-04-04)*
- [x] **RW6** Update docstring and user docs *(2026-04-04)*

---

## Quantity Memory Reduction ✅ Complete

Full plan: `claude/archive/QUANTITY_MEMORY_PLAN.md`
Target achieved: ~54% memory reduction (800 MB → ~368 MB for 10-quantity 1M-triangle domain).

- [x] **QM1** Introduce `qty_type` concept *(2026-04-09)*
- [x] **QM2** Lazy `vertex_values` property on all quantity types *(2026-04-09)*
- [x] **QM3** Strip update arrays from `elevation` *(2026-04-09)*
- [x] **QM4** Strip all arrays except `centroid_values` from `friction` *(2026-04-09)*
- [x] **QM5** Reduce `height`, `xvelocity`, `yvelocity` to centroid + edge only *(2026-04-09)*
- [x] **QM6** Make `x_gradient`, `y_gradient`, `phi` lazy for ALL types *(2026-04-10)*
- [x] **QM7** Shared gradient workspace on domain *(2026-04-13, commit 22559a5b)*

---

## Domain Work Array Memory Reduction ✅ Complete

~740 MB saved at N=2.25M triangles across three improvements.

- [x] **DM1** Defer all C-extension work arrays from `__init__` to first evolve step — 9 dead arrays removed, only 3 live arrays remain *(2026-04-15)*
- [x] **DM2** `edge_flux_type`/`edge_river_wall_counter` lazy for non-riverwall simulations *(2026-04-15)*
- [x] **DM3** `domain_memory_stats`, `print_domain_memory_stats`, `domain_struct_stats`, `print_domain_struct_stats` added to `system_tools.py` *(2026-04-15)*

### Measured end-to-end result (issue #33 benchmark, re-run 2026-07-06)

`mpiexec -np 2 /usr/bin/time -f "…%M…" python -u examples/parallel/run_parallel_rectangular.py --sqrtN N`,
`OMP_NUM_THREADS=1`, process-0 Max RSS. "Before" = the original figures in issue #33.

| sqrtN | triangles | RSS before | RSS now | reduction |
|------:|----------:|-----------:|--------:|:---------:|
| 250 | 250,000 | 710 MB | 511 MB | −28% |
| 500 | 1,000,000 | 2.5 GB | 1.37 GB | −45% |
| 750 | 2,250,000 | 5.1 GB | 2.74 GB | **−46%** |

The saving grows with N (fixed interpreter/import overhead dominates at small N; domain
arrays dominate at large N), so RSS is roughly halved at 2.25M triangles.

Hardware-independent per-rank domain arrays (`print_domain_memory_stats`, rank 0 after
`distribute`): quantities 55.7 / 222.0 / 499.0 MB and total numpy 102.7 / 409.6 / 920.5 MB
for local N = 125,706 / 501,465 / 1,127,216. At sqrtN=750 the breakdown is geometry 224,
connectivity 146, quantities 499, work arrays 26, **river wall 0.00** (DM2 lazy), other ~26,
total 920.5 MB. Remaining lever: rank 0's peak building the *full* domain before `distribute`.
Posted to issue #33.

---

## Benchmark Suite ✅ Complete

- [x] **B1** Single-process benchmark — `benchmarks/run_benchmarks.py` + `compare_benchmarks.py` *(2026-04-07)*
- [x] **B2** MPI distribution benchmark — `benchmarks/distribute_benchmarks.py` + `run_benchmark_grid.py` *(2026-04-07)*

---

## Bug Fixes ✅ Complete

- [x] **BF1** `Basic_mesh.reorder()` stale neighbours — ghost triangle count fix *(2026-04-07)*
- [x] **BF2** GPU test tolerances — relaxed to `atol=0.02` for real GPU hardware *(2026-04-11)*
- [x] **BF3** Mannings operator RuntimeWarning — `safe_h = maximum(height, 1e-15)` *(2026-04-11)*
- [x] **BF4** Rate_operator empty-check for numpy array — `hasattr(..., '__len__') and len(...) == 0` *(2026-04-11)*
- [x] **BF5** GPU_AWARE_MPI segfault — host staging buffers added in `gpu_halo.c` *(2026-04-11)*
- [x] **BF6** Rate_operator parallel false CPU-only — empty-indices operators marked `_gpu_initialized=True` *(2026-04-11)*
- [x] **BF7** Double `get_triangle_containing_point` call in parallel inlet enquiry *(2026-04-12)*
- [x] **BF8** Threshold-triggered spatial index — `MeshQuadtree` after 5 calls *(2026-04-12)*

---

## GPU / OpenMP Offloading — Phases 1–3 ✅ Complete

Full plan: `claude/archive/GPU_DEVELOPMENT_PLAN.md`

### Phase 1 — Correctness and test coverage ✅

- [x] **G1.1** File_boundary GPU support *(2026-04-09)*
- [x] **G1.2** Device memory check *(2026-04-09)*
- [x] **G1.3** Slot limit assertions → dynamic heap growth (superseded by G3.3) *(2026-04-07)*
- [x] **G1.4** End-to-end regression test; multi-rank halo exchange test; culvert test *(2026-04-07/09)*
- [x] **G1.5** SSP-RK3 GPU support *(2026-04-09)*

### Phase 2 — Performance validation ✅

- [x] **G2.1** GPU benchmark suite — `benchmarks/run_gpu_benchmarks.py` *(2026-04-10)*
- [x] **G2.2** GPU-aware MPI validation — runtime detection via `MPIX_Query_*` *(2026-04-10)*
- [x] **G2.3** NVTX/OMPT profiling hooks — `gpu_nvtx.h`, 10 kernel markers *(2026-04-10)*
- [x] **G2.4** Weak scaling scripts — `benchmarks/run_weak_scaling.py`, `scripts/hpc/weak_scaling.slurm` *(2026-04-10)*

### Phase 3 — Feature parity ✅

- [x] **G3.1** Gate/weir operators on GPU *(2026-04-10)*
- [x] **G3.2** Riverwall GPU support *(2026-04-10)*
- [x] **G3.3** Dynamic operator slot limits *(2026-04-10)*
- [x] **G3.4** GPU documentation page *(2026-04-10)*

---

## Kinematic Viscosity Parallelisation ✅ Complete (session 27, 2026-04-27)

- [x] **KV1** Remove Apple OpenMP guards from 4 C files (`sparse.c`, `kinematic_viscosity_operator.c`, `cg.c`, `fitsmooth.c`) — plain `#include "omp.h"` now that conda-forge llvm-openmp supports macOS *(2026-04-27)*
- [x] **KV2** Serial path: `parabolic_solve` routed through C CG (`cg_solve_c_precon`) with Jacobi preconditioner; `_build_parabolic_csr()` builds n×n parabolic matrix via vectorised numpy *(2026-04-27)*
- [x] **KV3** MPI parallel path (Option B distributed CG): `_exchange_ghost_vector` (non-blocking Irecv/Isend, tag 198), `_distributed_dot` (Allreduce SUM), `_parabolic_matvec_distributed` (ghost exchange before SpMV, n_full-length result), `_parabolic_solve_distributed` (standard CG loop on owned triangles only). `parallel_safe()` returns True. *(2026-04-27)*
- [x] **KV4** Tests: `run_parallel_kv_operator.py` + `test_parallel_kv_operator.py` (serial-vs-3proc xvelocity comparison, max diff 8.6×10⁻⁶); `run_parallel_kv_unit_tests.py` + `test_parallel_kv_unit_tests.py` (4 in-process MPI assertions: ghost exchange global-index round-trip, distributed dot Allreduce, matvec identity at dt=0, CG self-consistency). Bug fix: `test_select_alpha_degenerate_falls_back_to_default` was platform-dependent on Windows py3.10/3.11/3.13 due to numpy gradient differences — now uses `return_curve=True` to branch on actual kappa. Commits `61418742`, `5498f98d`. All CI passed. *(2026-04-27)*

---

## Mode-2 GPU-build test triage + isolated runner tooling ✅ Complete (session 40, 2026-06-17)

- [x] **M40.1** Pin 11 mode-1-only `shallow_water` tests to legacy so the unified default is green on a GPU build (Session 39's all-green was the CPU build, where device memory == host memory). 9 white-box tests call `compute_forcing_terms()`/`compute_fluxes()` and assert on host `semi_implicit_update`/`explicit_update` (mode-2 GPU computes on-device, never syncs back): in `test_forcing.py`, `test_friction.py`, `test_physics_sw.py` (Manning cases), `test_data_manager.py::test_sww_extrema`. 2 numerical tests compare against legacy-recorded references and diverge ~1e-6 from mode-2's reduction/eval order: `test_regression_snapshots.py::test_dam_break_DE1_stage_snapshot` (helpers pinned so the whole file is deterministic), `test_sww_interrogate.py::test_get_maximum_inundation_de0`. No-op for the legacy default. Result: 408 pass / 2 skip under `-cm unified` on the GPU build. Commit `0c50947d`. *(2026-06-17)*
- [x] **M40.2** Install the per-test process-isolation harness as a console command: moved `anuga/shallow_water/tests/run_isolated_tests.py` → `scripts/anuga_run_isolated_tests.py`, registered in `meson.build` (`configure_file` → bindir, matching the other `anuga_*` scripts). Made it install-safe: importlib-resolved default target, cwd-seeded rootdir, `_abs_nodeid` passes absolute/`--pyargs` ids through. Commit `37eccc6d`. *(2026-06-17)*
- [x] **M40.3** Added `-cm`/`--compute-mode {legacy,unified}` to set `ANUGA_DEFAULT_COMPUTE_MODE` for every child process (omit to inherit; banner prints the resolved mode). Verified via `test_default_is_legacy` (skips under unified, passes under legacy). Commit `34401cde`. *(2026-06-17)*
- [x] **M40.4** Docs: green-run note in `KNOWN_ISSUES.md` (commit `f57d0532`); compute-mode testing story across `CLAUDE.md`, new `CONVENTIONS.md` → "Compute mode in tests", `DECISIONS.md` rationale entry, `SESSION_GUIDE.md` Session 40 summary + quick-reference rows (commit `905cb1f2`). *(2026-06-17)*

## Documentation Overhaul ✅ Complete (session 47, 2026-07-07)

Merged to `anuga-community/develop` as PRs **#157–#164** (admin-merged by number).

- [x] **D47.1** Restructured the docs into a standard-user **Contents** vs advanced **Appendices** split, and added a landing-page **quick start** (install + minimal runnable model). Moved Parallelisation into Contents and the TOML/ANUGA-Viewer/QGIS pages into the standard sections; led the appendices with the developer + new GPU-install pages (PR #157).
- [x] **D47.2** New **Conventions & units** primer (`setup_anuga_script/conventions.rst`): coordinate system + `Geo_reference`, quantity units and `depth = stage - elevation`, wet/dry thresholds, `anuga.config` physical defaults, time conventions (PR #157).
- [x] **D47.3** New evolve **"Stability and blow-ups"** section: reading `delta t`/`steps`, warning signs, common causes, and a live `Domain.diagnose_timestep` snippet linked to troubleshooting (PR #157).
- [x] **D47.4** Converted every narrative "Reference" block into compact `autosummary` tables that link into the **API Reference** (expanded API consolidated there), and made each class page's method/attribute summary **link to per-method signatures** via `autodoc_default_options={'members':True}` + `sphinx.ext.napoleon` (PR #157).
- [x] **D47.5** Reframed the **ANUGA Viewer** as the recommended fast viewer for large `.sww` files (dropped "legacy"), linked to `anuga-viewer.readthedocs.io` (PR #157).
- [x] **D47.6** New GPU-install appendix (`appendices/install_gpu.rst`): NVIDIA HPC SDK / `nvc` setup, `install_anuga_nvc.sh` vs manual build, `gpu_arch` table, build-dir-switch warning (PR #157).
- [x] **D47.7** Content-review fixes: typos, heading-level consistency, a missing `pip install` path + version smoke-test, quantity units, `../../examples` link cleanup (PR #157).
- [x] **D47.8** New **Citing**, **Contributing**, and **Glossary** meta-pages, surfacing `CITATION.cff`, the user-manual DOI, Apache-2.0, `CONTRIBUTING.rst`, and 15 domain terms (stage/elevation/depth, yieldstep/finaltime/duration/outputstep, CFL, riverwall, DE algorithm, …) (PR #157).
- [x] **D47.9** Issue #32 "Make riverwalls transmissive": verified `Cd_through` throughflow is already active in **both** compute modes (legacy `_openmp_compute_fluxes_central` delegates to the shared `core_compute_fluxes_central`); confirmed bit-identical empirically; documented on the issue and in the docs.
- [x] **D47.10** Fixed malformed-RST source docstrings surfaced by autodoc `members`: `Quantity`/`Domain` ×6 (PR #158) and `internal_boundary_operator` ×2 (PR #160).
- [x] **D47.11** `install_anuga_nvc.sh` now builds into an already-activated conda env if present, else falls back to `$HOME/miniforge3` + `anuga_env_$PY` (PR #158).
- [x] **D47.12** Warning-free local build: fixed `html_static_path` (`_static` dir) and the `Geo_reference.epsg` duplicate — napoleon rendered a class-docstring *Attribute* **and** the real property, resolved with `napoleon_use_ivar = True` (PR #159).
- [x] **D47.13** Read the Docs `develop` build made **clean**: it had surfaced 62 warnings the local build hid — 56 `ipython3` Pygments-lexer (RTD lacks IPython → added `ipython` to `docs/requirements.txt`, PR #161) + the operator docstrings. Verified via RTD build 33476518 (only the harmless MPI-less `Could not import mpi4py` remains).
- [x] **D47.14** Added a Contributing **"Building the documentation"** note: reproduce the RTD build from a clean `docs/requirements.txt` env; a stray IPython or `nbsphinx_execute=never` can hide warnings RTD shows (PR #162).
- [x] **D47.15** Recorded the **branch policy** in `ROADMAP.md` (and memory): do not merge `develop` → `main` until the team cuts v4.0.0 (PR #163). Session 47 summary added to `SESSION_GUIDE.md` (PR #164).

## GPU mode-2 Time_boundary substep fix ✅ Complete (option B; session 47, 2026-07-08)

- [x] **GB.1** Root-caused a mode-1 vs mode-2 divergence with time-varying
  boundaries. A rising-tide (`Time_boundary`) flood diverged **~4.3e-3 m** (GPU
  vs CPU-legacy); reflective/steady boundaries agree to ~1e-12–1e-15. Isolated
  it by RK substep count: **DE0 / DE_ader2** (single substep) match to machine
  precision; **DE1 (rk2) / DE2 (rk3)** diverge — and only with **Python-evaluated
  boundaries** (Time/File/Field, wave, Flather, transmissive-set-stage). Cause:
  the single-call **C RK loop** (`_evolve_one_rk*_step_c`) sets those boundaries
  on the device **once per timestep**, reusing that value for every RK substep,
  whereas mode-1 calls `update_boundary()` **before each substep** — an O(dt)
  boundary-forcing error. Verified by call-time instrumentation (mode-1 evaluates
  the boundary 2× per rk2 step, mode-2 1×) and by `use_c_rk_loop=False`
  (Python-orchestrated loop) → 1e-15.
- [x] **GB.2** Fix (**option B**, PR #171): `_has_python_evaluated_gpu_boundaries()`
  gates the mode-2 RK2/RK3 dispatch — domains with any Python-evaluated boundary
  route to the Python-orchestrated GPU loop (refreshes the boundary per substep →
  bit-matches mode-1, ~4e-3 → ~1e-15). Reflective/steady keep the fast C loop;
  single-substep DE0/DE_ader2 untouched. Benchmarked GPU cost **≤ ~4%** (within
  noise, 14k–640k triangles).
- [x] **GB.3** Regression test `Test_GPU_TimeBoundarySubstep` in
  `test_DE_gpu_omp.py` (DE1/DE2 mode-1==mode-2 to atol 1e-6 + routing checks);
  full GPU file green (69/69) via the isolated runner. No rebuild needed
  (pure-Python change).
- [ ] **Option A** — proper per-substep evaluation inside the C RK loop: tracked
  as **issue #170** (see Remaining / Deferred in `PROGRESS.md`).

## Fractional-step operator evaluation timing ✅ Complete (session 47, 2026-07-08)

Fractional-step operators are applied by the evolve loop *before* it advances
`relative_time` from t to t+dt, so they should evaluate forcing at the pre-step
time **t** (the mode-2 code documents that t+dt is "one step too far"). First
confirmed operators/structures are applied **every inner timestep** in both modes
(13/13 and 46/46 evals ≫ yieldsteps — not a yieldstep-only issue). Then found and
fixed an operator-*time* bug in the RK schemes:

- [x] **OT.1 — DE1 (rk2), PR #174.** The mode-1 rk2 body advanced `relative_time`
  mid-step (for the substep-2 boundary) and never restored it, so its operators
  evaluated forcing at **t+dt**, diverging from mode-2 (which uses t) by ~4e-4 for
  a time-varying `Rate_operator`/`Inlet`. Restore the pre-step time at the end of
  the mode-1 rk2 body → DE1 matches DE0/DE2/DE_ader2 and mode-2 (0.0). Added
  `Test_GPU_OperatorTimeAlignment` (cross-mode, all algorithms).
- [x] **OT.2 — CPU regression test, PR #175.** The cross-mode guard is GPU-only
  (skips on standard builds), so added `test_operator_timing.py` — a mode-1-only
  test that runs anywhere: a time-varying operator must be evaluated at t
  (last inner step's eval < finaltime). Reverting a fix fails it. Registered in
  `tests/meson.build`.
- [x] **OT.3 — DE2 (rk3), PR #177 (closes #176).** Subtler than DE1: **all three**
  rk3 paths left time advanced (mode-1 body, mode-2 C loop, mode-2 GPU loop), so
  DE2 was post-step in *both* modes — self-consistent (mode-1 == mode-2), so the
  cross-mode check missed it. Restore the pre-step time in the mode-1 body and
  mode-2 GPU loop; drop the advance in the mode-2 C loop (matching the rk2 C loop).
  DE2 now evaluates operators at t in both modes. Added
  `test_rk3_operator_evaluated_at_pre_step_time`.

Result: all four flow algorithms evaluate fractional-step operators at the
pre-step time t, mode-1 == mode-2. No prior test used a time-varying operator, so
the bug was invisible (coverage gap now closed). CPU suite 2610 pass, GPU 73/73.
