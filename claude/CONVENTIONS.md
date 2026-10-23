# Coding Conventions

Conventions observed in or established for the ANUGA codebase.

---

## Python Style

- **PEP 8** required — use `autopep8` to auto-fix before committing
- **pyflakes** must produce no warnings
- `numpy` imported as `num` throughout (not `np`) — this is a project-wide convention
- `import logging; log = logging.getLogger(__name__)` or via `import anuga.utilities.log as log`
- No f-strings in C extension interface code (plain `%` formatting preferred there)

---

## Docstrings

- **NumPy style** for all new docstrings (not Sphinx `:param:` style)
- Include `Parameters`, `Returns`, `Raises`, `Notes`, `Examples` sections as appropriate
- Types in docstrings use Python type names (`float`, `int`, `ndarray`, `Domain`, etc.)
- Units should be stated in the type field or description: `rate : float — Flow rate in m³/s`

---

## Testing

- Test files: `anuga/*/tests/test_*.py`
- Use `unittest.TestCase` (existing pattern) or plain `pytest` functions (new tests)
- Slow tests (>5s): decorate with `@pytest.mark.slow` or add `pytestmark = pytest.mark.slow`
  at module level if all tests in the file are slow
- Parallel tests live in `anuga/parallel/tests/` — auto-marked slow by `conftest.py`
- Run all tests: `pytest --pyargs anuga`
- Run fast tests only: `pytest --pyargs anuga --run-fast`
- Run only slow tests: `pytest --pyargs anuga -m slow`

---

## Compute mode in tests

Domains default to the `legacy` compute path (mode 1); `unified` (mode 2) routes
the step through the C `gpu_ext` kernels and, on a GPU build, offloads to the
device. The default is selectable per process via `ANUGA_DEFAULT_COMPUTE_MODE`
(`legacy`|`unified`), so tests should not assume one mode unless they pin it.

- **Pin mode-1-only white-box tests to legacy.** A test that calls
  `compute_forcing_terms()` / `compute_fluxes()` directly and asserts on the host
  `quantities[...].semi_implicit_update` / `explicit_update` arrays is exercising
  mode-1 internals: mode-2 computes those on-device and never syncs the host
  arrays back, so they read stale zeros. Pin such a test right after constructing
  the domain:

  ```python
  domain = Domain(points, vertices)
  domain.set_compute_mode('legacy')   # white-box: reads host *_update arrays
  ```

  Existing examples: `test_forcing.py`, `test_friction.py`, `test_physics_sw.py`
  (Manning-friction cases), `test_data_manager.py::test_sww_extrema`. Numerical
  snapshot/reference tests recorded under legacy (e.g.
  `test_regression_snapshots.py`) are pinned the same way — mode-2 differs at the
  ~1e-6 level from a different reduction/eval order.
- **End-to-end (evolve-based) tests need no pin** — they run correctly in either
  mode; mode-1-vs-mode-2 equivalence is covered separately in `test_DE_gpu_omp.py`.
- **Run a suite under a chosen mode** with the isolated runner's flag:
  `anuga_run_isolated_tests --pyargs anuga.shallow_water -cm unified` (see
  `CLAUDE.md` → "Testing a GPU-offload (nvc) build").

---

## Naming

- New public methods: `snake_case`
- camelCase methods in `pmesh/mesh.py` are deprecated — use their snake_case equivalents
- Quantity names: `stage`, `elevation`, `xmomentum`, `ymomentum` (not `xmom`/`ymom`)
- File format: prefer `.sww` extension for output files
- Domain variable typically named `domain`; mesh variable `mesh`

---

## Mutable Defaults

Never use mutable default arguments. Always use `None` and create inside the function body:

```python
# Wrong:
def f(x=[]):
    pass

# Correct:
def f(x=None):
    if x is None:
        x = []
```

---

## None Checks

Always use `is None` / `is not None` — never `== None` (can raise `ValueError` for ndarrays).

---

## File Operations

Use `with open(...)` context managers. The only exception is when a file handle must persist
across multiple method calls (e.g., `urs.py` iterator pattern).

---

## Imports

- Standard library imports first, then third-party, then anuga-internal
- Lazy imports (inside functions) are used throughout for optional dependencies
  (e.g., `import psutil`, `import xarray`) — keep this pattern
- `from anuga.utilities.system_tools import ...` for utility functions

---

## Build System

- **meson-python** — not setuptools
- Always use `pip install --no-build-isolation -e .` (the `--no-build-isolation` flag is required)
- C extensions: `.pyx` files compiled by Cython; generated `.c` files are in `.gitignore`
- OpenMP enabled conditionally per platform in `meson.build`

---

## Parallel Code

- MPI parallel code requires `mpi4py` and `pymetis`
- Parallel domain inherits from serial `Domain`
- `anuga.distribute()` for full domain; `anuga.distribute_basic_mesh()` for basic mesh
- Test parallel code with `mpiexec -np N python script.py`

---

## Error Handling

- Do not use bare `except:` — always specify exception type(s)
- Silent `except ImportError: pass` for optional dependencies → use `log.debug("... not available: %s", e)`
- Expected exceptions with no action needed: add a comment explaining why it is expected
