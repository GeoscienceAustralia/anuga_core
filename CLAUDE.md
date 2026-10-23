# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Session Continuity

The `claude/` directory contains documents to orient a new session quickly:

| File | Contents |
|------|---------|
| `claude/SESSION_GUIDE.md` | Branches, common commands, recent session summaries |
| `claude/ROADMAP.md` | Release plan — v3.3.2 shipped, v4.0.0 (SC26) in progress |
| `claude/PROGRESS.md` | Overview totals and remaining/deferred items |
| `claude/PROGRESS_ARCHIVE.md` | Full history of all completed work with dates |
| `claude/FUTURE_WORK.md` | Prioritised recommendations for future sessions (P1–P3) |
| `claude/DECISIONS.md` | Key design choices with rationale |
| `claude/CONVENTIONS.md` | Coding style, naming, testing, import conventions |
| `claude/KNOWN_ISSUES.md` | Surprises, gotchas, non-obvious behaviour |

**Start here when picking up existing work:** read `claude/SESSION_GUIDE.md` then `claude/FUTURE_WORK.md`.

## What is ANUGA

ANUGA is a Python package (with C/Cython extensions) for simulating shallow water equations — primarily for tsunami and flood modeling. The core algorithm is a finite volume method on unstructured triangular meshes.

## Build System

ANUGA uses **meson-python** (not setuptools). Build requirements: `meson`, `ninja`, `Cython`, `pybind11`, `numpy >= 2.0.0`.

```bash
# Standard install (requires compilers in environment)
pip install --no-build-isolation -v -e .

# Recommended: use a conda environment first
conda env create --name anuga_env --file environments/environment_3.10.yml
conda activate anuga_env
pip install --no-build-isolation -v -e .
```

The `--no-build-isolation` flag is recommended (not strictly required in all environments) to ensure meson-python uses the existing numpy/Cython from the conda environment.

OpenMP is enabled conditionally: Linux/macOS use gcc/clang with OpenMP; Windows requires mingw compilers.

## Running Tests

```bash
# Run all tests (~1 600 tests, ~3 minutes)
pytest --pyargs anuga

# Fast run — skips parallel (MPI) and other marked slow tests (~1 500 tests, ~40 s)
pytest --pyargs anuga --run-fast

# Run only slow tests
pytest --pyargs anuga -m slow

# Run a single test file
pytest anuga/shallow_water/tests/test_shallow_water_domain.py

# Run a specific test
pytest anuga/shallow_water/tests/test_shallow_water_domain.py::TestCase::test_name

# CI runs tests from sandpit/ directory with OMP_NUM_THREADS=1
cd sandpit && OMP_NUM_THREADS=1 pytest -rs --pyargs anuga
```

### Testing a GPU-offload (nvc) build

On the standard **CPU build** (`-Dgpu_offload=false`, gcc) everything above runs
as-is, in one process — including the mode-2 GPU tests (`test_DE_gpu_omp.py`),
whose `omp target` regions execute on the host.

On a **GPU-offload build** (`-Dgpu_offload=true`, nvc) the NVHPC OpenMP-target
runtime aborts the process once many mode-2 GPU domains have been created in it,
so the full suite cannot run in one process. Do this instead:

```bash
# 1. Bulk suite — legacy default (domains stay on the CPU), GPU file excluded:
OMP_NUM_THREADS=1 ANUGA_DEFAULT_COMPUTE_MODE=legacy \
  pytest anuga/shallow_water/tests/ \
  --ignore=anuga/shallow_water/tests/test_DE_gpu_omp.py

# 2a. GPU tests — one fresh process per CLASS (fast). Do NOT use --forked
#     (CUDA is fork-unsafe). test_DE_gpu_omp.py auto-skips in a normal in-process
#     run on a GPU build; this runner sets ANUGA_GPU_TESTS_ISOLATED=1 to opt in:
bash anuga/shallow_water/tests/run_gpu_tests_isolated.sh

# 2b. Maximum isolation — one fresh process per TEST FUNCTION, with a per-test
#     timeout so a genuine hang is reported (not stuck). Works on any target.
#     Installed as the `anuga_run_isolated_tests` command (scripts/, via meson);
#     in a source checkout run scripts/anuga_run_isolated_tests.py directly:
anuga_run_isolated_tests                                  # the GPU file
anuga_run_isolated_tests --pyargs anuga.shallow_water -k riverwall --timeout 120
anuga_run_isolated_tests --pyargs anuga.shallow_water

# 2c. Whole anuga.shallow_water set under the unified default — GREEN on a GPU
#     build *only via the isolated runner* (each test in a fresh process, so the
#     NVHPC abort never accumulates). -cm/--compute-mode sets the default per-domain
#     compute path (ANUGA_DEFAULT_COMPUTE_MODE) for every child:
anuga_run_isolated_tests --pyargs anuga.shallow_water -cm unified   # 408 pass, 2 skip
```

`-cm legacy|unified` (omit to inherit the environment) is the ergonomic way to
flip compute mode; it is equivalent to exporting `ANUGA_DEFAULT_COMPUTE_MODE`.

Do **not** run `ANUGA_DEFAULT_COMPUTE_MODE=unified` over the full suite **in one
process** on a GPU build — every default domain then offloads and the suite
aborts early. Use the isolated runner (2c) on a GPU build, or validate the
unified default in one process on the **CPU build** (the documented all-green
run). See `claude/KNOWN_ISSUES.md` for the root cause.

**Compute mode and tests:** a handful of white-box tests probe mode-1-only host
state — they call `compute_forcing_terms()` / `compute_fluxes()` and assert on
the host `semi_implicit_update` / `explicit_update` arrays, which mode-2
('unified') computes on-device and never syncs back. Such tests pin themselves
to legacy with `domain.set_compute_mode('legacy')` right after constructing the
domain (see the existing pins in `test_forcing.py`, `test_friction.py`,
`test_physics_sw.py`). Follow that pattern for any new test that inspects those
internal Python update arrays. See `claude/CONVENTIONS.md` → "Compute mode in
tests".

Tests are marked slow with `@pytest.mark.slow`. All tests under `anuga/parallel/tests/`
are automatically treated as slow (they spawn MPI subprocesses). To mark a new slow test:

```python
import pytest

@pytest.mark.slow
def test_my_expensive_computation(self):
    ...
```

Test files follow the pattern `anuga/*/tests/test_*.py`. There are ~126 test files.

Validation tests (against analytical solutions and experimental data) live in `validation_tests/` and are run separately via `run_auto_validation_tests.py`.

## Code Quality

```bash
ruff check anuga/path/to/module.py   # linting (configured in pyproject.toml)
ruff check --fix anuga/path/to/module.py   # auto-fix safe violations
```

Pre-commit hooks run `ruff check --fix` and `ruff-format` automatically on staged files.

## Architecture

### Core Data Flow

```
User script
  → anuga.Domain (shallow_water/domain.py)
      → Mesh (abstract_2d_finite_volumes/mesh.py)
      → Quantities (abstract_2d_finite_volumes/quantity.py)
      → Operators/Structures (operators/, structures/)
  → .sww file output (NetCDF format)
```

### Key Packages

**`abstract_2d_finite_volumes/`** — The mathematical core. Contains `Domain` base class, `Mesh`, `Quantity`, and the finite volume update machinery. The shallow water domain inherits from this.

**`shallow_water/`** — The main solver. `Domain` here (inherits from abstract domain) implements the shallow water equations. The time-stepping and flux computations are in Cython/C extensions (`sw_domain_ext.c`, `sw_domain_openmp_ext.c`).

**`structures/`** — Culverts, weirs, inlets modeled as operators that transfer flow between mesh regions. `Structure_operator` is the base; `Boyd_box_operator`, `Boyd_pipe_operator`, `Weir_orifice_trapezoid_operator` are main implementations.

**`operators/`** — Domain operators applied each timestep: `Rate_operator` (rainfall/extraction), `Kinematic_viscosity_operator`, `Bed_shear_erosion_operator`, etc.

**`fit_interpolate/`** — Fitting point cloud data onto unstructured meshes, and interpolating mesh quantities to arbitrary points. Computationally intensive; has C extensions.

**`geometry/`** — Polygon operations (inside/outside tests, clipping, alpha shapes). Core to mesh generation and region definitions.

**`parallel/`** — MPI-based domain decomposition via `pmesh2domain`. Requires `mpi4py` and `pymetis`. Parallel domain inherits from the serial `Domain`.

**`file/` and `file_conversion/`** — Reading/writing SWW (NetCDF), DEM, ESRI grid, URS (tsunami source) formats.

**`coordinate_transforms/`** — UTM ↔ lat/lon conversions. `Geo_reference` class tracks the origin offset used to keep coordinates numerically stable.

### Public API

`anuga/__init__.py` exports the full public API (~1000 lines). Key entry points:
- `anuga.Domain` — create simulation domain
- `anuga.create_domain_from_regions()` — mesh + domain from polygon regions
- `anuga.file_function()` — interpolate SWW data as a time function
- Boundary classes: `Reflective_boundary`, `File_boundary`, `Dirichlet_boundary`, etc.
- `anuga.Boyd_box_operator`, `anuga.Inlet_operator`, etc. — structure operators

### Global Configuration

`anuga/config.py` contains physical and numerical constants used throughout:
- `g = 9.8` (gravity)
- `epsilon = 1.0e-6` (wet/dry threshold)
- `manning = 0.03` (default friction)
- `minimum_allowed_height = 1.0e-05`

### C/Cython Extensions

Each package with performance-critical code has a `*.pyx` (Cython) or `*_ext.c` source. The `meson.build` files in each subdirectory compile these. The generated `.c` files (e.g., `sw_domain_openmp_ext.c`) are in `.gitignore` but appear as untracked in `git status` — they are build artifacts.

### GPU/Parallel Development

The `develop_gpu`/`develop_cupy` and `sp26` branches explore GPU acceleration. CUDA examples are in `examples/cuda/`. The parallel module uses MPI; GPU work uses CuPy as a NumPy replacement.
