# Ubuntu PPA for ANUGA — feasibility

Assessment for issue #25 ("Possibility of creating PPA for ANUGA"), 2026-08-23.
Facts below were checked against the Launchpad API and packages.ubuntu.com on
that date, not assumed.

**Summary: technically possible, but it means maintaining three source packages,
one of which has a licence problem — and the users it would serve are already
served by pip wheels, conda-forge and the Docker images.** Recommend not doing
it unless someone specifically wants `apt install python3-anuga`.

---

## What a PPA actually requires

Launchpad builds binaries from a Debian *source* package you upload with `dput`.
So it needs: a Launchpad account and GPG key, a `debian/` directory
(`control`, `rules`, `changelog`, `copyright`, `watch`), and a rebuild+upload for
**every Ubuntu series × every release**. Launchpad does not accept a wheel.

## Dependency audit

Ubuntu already carries almost everything:

| ANUGA needs | in Ubuntu? | version seen |
|---|---|---|
| `numpy>=2.0` | yes | 2.4.6 |
| `scipy>=1.11` | yes | 1.17.1 |
| `matplotlib>=3.7` | yes | 3.8.3 |
| `netCDF4>=1.6` | yes | 1.7.4 |
| `pyproj>=3.6` | yes | 3.7.2 |
| `xarray` | yes | 2026.04.0 |
| `dill>=0.3.7` | yes | 0.4.1 |
| build: `meson`, `python3-mesonpy`, `cython3`, `pybind11-dev` | yes | 1.7.0 / 0.19.0 / 3.1.6 / 3.0.1 |
| **`meshpy>=2022.1`** | **no** | — |
| **`pymetis>=2023.1`** | **no** | — |

So a PPA is not one package; it is **three**: `pymetis`, `meshpy`, `anuga`.

## The two gaps are not the same size

**`pymetis`** — merely unpackaged. METIS is Apache-2.0, so this is ordinary
packaging work. ANUGA imports it *lazily* (only in
`anuga/parallel/partitioning.py`, at call time), so it can be `Recommends:`
rather than `Depends:` — serial ANUGA works without it.

**`meshpy`** — a licence problem, not an effort problem. MeshPy wraps Jonathan
Shewchuk's **Triangle**, which is free for non-commercial use only. Ubuntu ships
`triangle` in **multiverse**, the component for licence-restricted software,
which is almost certainly why `python3-meshpy` does not exist in main or
universe. ANUGA imports meshpy at **import time**
(`anuga/mesh_engine/mesh_engine.py`, reached by `import anuga`), so it is a hard
dependency.

There is a fallback in that file — `except ImportError: import triangle` — but
the `triangle` PyPI package wraps the same Shewchuk code, so it carries the same
restriction. `python3-triangle` is not in Ubuntu either.

A PPA *may* still distribute this (Triangle permits free redistribution; it is
sale and commercial inclusion that need permission), but it is a question to
answer deliberately rather than discover after publishing.

## Series and Python versions

`requires-python = ">=3.10, <3.15"`, so viable targets are:

| series | Python | viable |
|---|---|---|
| jammy 22.04 LTS | 3.10 | yes |
| noble 24.04 LTS | 3.12 | yes |
| resolute 26.04 | 3.14 | yes |
| bionic/focal and older | < 3.10 | no |

Three series × a rebuild per release, plus the two dependency packages, is the
ongoing cost.

## What already serves these users

Issue #25 asked for "an installation possibility for Ubuntu Linux users". Since
it was filed, three arrived:

* **pip wheels** — manylinux2014 x86_64 for cp310–cp314, so
  `pip install anuga` works on jammy, noble and resolute with the system Python.
  Verified against 4.0.0rc1 in a clean venv.
* **conda-forge** — `conda install -c conda-forge anuga`.
* **Docker images** — CPU, GPU and GPU+MPI, published to GHCR.

A PPA adds one thing those do not: `apt install python3-anuga`, integrated with
system packages and upgraded by `apt upgrade`. That is a real convenience for
sysadmin-managed machines and teaching labs, and nothing else provides it.

## Recommendation

Don't, unless someone asks for the apt integration specifically. The cost is
three source packages across three series, indefinitely, and the meshpy licence
question has to be settled first. If it is wanted anyway, the order is:

1. Settle the Triangle/meshpy redistribution question — it gates everything.
2. Package `pymetis` (easy, Apache-2.0). Useful on its own.
3. Package `meshpy`, or make ANUGA's meshpy import lazy so the PPA can ship
   ANUGA with meshpy as a `Suggests:` and let users get it from pip.
4. Package `anuga` itself; `debian/rules` is thin over `pybuild` + meson-python.

Step 3's alternative — making the meshpy import lazy — is worth doing regardless
of the PPA: it would let `import anuga` succeed without a mesh generator
present, which matters for any environment where Triangle's licence is awkward.
