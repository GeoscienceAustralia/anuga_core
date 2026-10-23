# ANUGA Release Roadmap

Last updated: 2026-07-08

> **Branch policy (decided 2026-07-07):** `develop` will **not** be merged into
> `main` until the team explicitly decides to cut **v4.0.0**. Until then, all
> work and upstream syncs stay on `develop` (PR from the `stoiver` fork →
> `anuga-community/anuga_core:develop`, admin-merged by PR number). `main` is
> held at the current released line — do not open or merge develop→main PRs.

> **Tagging convention (decided 2026-07-08):** release tags are bare version
> numbers (e.g. `3.3.7`, **no** `v` prefix) and must be created as **annotated**
> tags — `git tag -a 3.3.8 -m "ANUGA 3.3.8"` — not lightweight. The `3.x` tags to
> date are lightweight, which is why a plain `git describe` skips them and falls
> back to the old annotated `1.3.1`; the version machinery (`_git_version.py`)
> still works because it uses `git describe --tags`, but annotating future tags
> makes plain `git describe` (and any tooling that relies on it) correct too.

> **Validation reports (decided 2026-10-01):** publishing a GitHub release
> also builds the validation reports (DE0, DE1, DE_ader2), attaches them to the
> release and, for final releases, opens a Zenodo draft to review and publish.
> Setup and per-release steps: `claude/RELEASE_VALIDATION_REPORTS.md`.

---

## Version 3.3.0 — Imminent

**Branch:** `develop` → `main` on [anuga-community/anuga_core](https://github.com/anuga-community/anuga_core)

### Release steps

1. Merge all active `develop_*` feature branches into `develop`
2. Run full test suite: `pytest --pyargs anuga`
3. Merge `develop` → `main`
4. Tag the release commit on `main` with an **annotated** tag (bare version, no
   `v` prefix): `git tag -a 3.3.0 -m "ANUGA 3.3.0"` (see Tagging convention above)
5. Publish to PyPI (triggered by tag via CI or manual `twine upload`)
6. Submit conda-forge feedstock PR to bump version and hash

### What's in 3.3.0 (develop branch improvements)

- TOML/Excel scenario configuration (`develop_excel`)
- EPSG coordinate reference system support (`develop_epsg`)
- Code quality improvements (mutable defaults, docstrings, API `__all__`, naming)
- New tests for operators and structures
- Fast/slow test infrastructure (`--run-fast`, `@pytest.mark.slow`)
- Memory reporting in timestepping statistics
- New API functions: `memory_stats`, `basic_mesh_from_mesh_file`, `distribute_basic_mesh`
- Documentation: 20 Sphinx improvements, MPI section, troubleshooting, validation

### Pre-release checklist

- [ ] All `develop_*` branches merged into `develop`
- [ ] Full test suite passes cleanly (`pytest --pyargs anuga`)
- [ ] `pyproject.toml` version bumped to `3.3.0`
- [ ] `CHANGELOG` / release notes updated
- [ ] Sphinx docs build without warnings (`cd docs && make html`)
- [ ] `python -W error::DeprecationWarning -c "import anuga"` — no warnings

---

## Version 4.0.0 — Future

**Foundation:** `sp26` branch → `develop` → `main`

### sp26 branch

The `sp26` branch is a **research project** forming the basis of a paper submitted to
**Supercomputing 2026 (SC26)**. It implements:

- GPU acceleration via OpenMP target-offloading (`multiprocessor_mode=2`)
- High-performance parallel solver improvements
- The next major capability jump for ANUGA

### v4.0.0 development plan

1. SC26 paper submitted / accepted
2. Merge `sp26` into `develop` as foundation for v4.0.0
3. Apply Hydrata refactor plan phases (see `PROGRESS.md`):
   - **Phase 3** — deduplicate quantity kernels, consolidate parallel operators (coordinate with sp26 merge)
   - **Phase 2** — ruff linting, pre-commit hooks
   - **Phase 4** — lift coverage to 65%, automate validation tests
4. **Only once the team decides to release 4.0.0:** stabilise, merge `develop` →
   `main`, and add an **annotated** tag `git tag -a 4.0.0 -m "ANUGA 4.0.0"`
   (bare version, no `v` prefix — see Tagging convention above). Until that
   decision `main` stays put and `develop` is not merged into it (see the branch
   policy at the top of this file).

### Coordination notes

- **Quantity kernel unification (H3.1)** is high-risk because sp26 likely modifies
  `quantity_ext_openmp.pyx`. Do this work *after* sp26 is merged into develop, not before.
- **Parallel operator consolidation (H3.2)** similarly — sp26 may change the parallel
  domain structure. Coordinate carefully.
- **Linting (H2)** and **dependency declaration (H1)** are safe to do any time — no
  conflict risk with sp26.

---

## PyPI / conda-forge packaging

### PyPI

ANUGA is published to PyPI as the `anuga` package. Publication is triggered by tagging
a release on `main`. The `pyproject.toml` `[project]` section controls the metadata.

**Known gap:** `pyproject.toml` currently only declares `numpy>=2.0.0` as a runtime
dependency. Before or alongside the 3.3.0 release, add the missing deps:
`scipy`, `netCDF4`, `matplotlib`, `meshpy`, `dill`, `pymetis`, `pyproj`, `affine`.
(See `KNOWN_ISSUES.md` and `PROGRESS.md` item H1.1.)

### conda-forge

The conda-forge feedstock is a separate repository. To release a new version:
1. Fork/update the `anuga` feedstock recipe
2. Bump `version`, update the source SHA256 hash
3. Verify that all runtime dependencies are listed in `meta.yaml` `requirements/run`
4. Submit PR to conda-forge/staged-recipes or the feedstock repo

**Note:** conda-forge `meta.yaml` and `pyproject.toml` dependency lists must be kept
in sync. The conda-forge recipe is the authoritative source for the conda package.
