<!--
Draft announcement for the 4.0.0rc1 soak, to Ole, Rudy, Petar, david and Jorge.
Written 2026-08-23; RC published to TestPyPI the same day (21 artifacts,
cp310-cp314, Linux/macOS/Windows + sdist).

Kept here so the wording and the per-person asks survive the session. If the
RC slips or is rebuilt, check the version and the Towradgi figures still match
before sending. See claude/RELEASE_PLAN_4.0.0.md Phase 2.
-->

Subject: ANUGA 4.0.0rc1 available for testing on TestPyPI

Hi Ole, Rudy, Petar, David, Jorge,

A release candidate for ANUGA 4.0.0 is up on TestPyPI and I'd appreciate a look
before we tag. This is the first release off the develop line since 3.3.10 —
993 commits — so it's worth a few pairs of eyes on real models rather than just
the test suite.

Installing (the second index is needed — TestPyPI doesn't carry numpy, scipy,
netCDF4 etc., so without it the install fails on dependencies):

    pip install --index-url https://test.pypi.org/simple/ \
                --extra-index-url https://pypi.org/simple/ anuga==4.0.0rc1

Wheels are there for Python 3.10–3.14 on Linux, macOS (arm64 and x86_64) and
Windows, plus an sdist.


What's in it
------------

* A unified compute mode. The solver and its operators now share one C
  implementation that runs CPU-multicore, and offloads to a GPU on a build made
  with the NVIDIA HPC SDK. It's opt-in and not GPU-only:

      domain.set_compute_mode('unified')      # CPU-multicore by default
      anuga.set_gpu_offload(True)             # separate, needs an nvc build

  The default stays 'legacy', so existing scripts run unchanged.

* DE_ader2, a new timestepping method — about 1.75x faster than DE1 at
  equivalent accuracy.

* Rainfall: raster rainfall input and Australian Rainfall & Runoff support
  (Raster_rate_operator, ARR_rate_operator and friends).

* Container images (CPU, GPU, GPU+CUDA-aware-MPI) published to GHCR.

* anuga.culvert_flows is REMOVED — it was superseded by the structure operators
  years ago. The forcing classes (Rainfall, Inflow, Wind_stress,
  Barometric_pressure) are deprecated and go in 4.1.


The one thing that changes results
----------------------------------

Structure operators (culverts, weirs) used to write their transfer back as a
uniform DEPTH across the inlet. On a sloping bed that tilts a level water
surface onto the bed, so a lake at rest was disturbed by roughly half the bed
elevation range across the inlet — every timestep, even with no flow through
the structure. They now level the water surface instead.

  * Flat bed under the inlet: results are bit-identical. Nothing changes.
  * Sloping bed: results move, in proportion to the bed range across the inlet.

On the Towradgi catchment (Petar's case study — 22 culverts on real terrain,
median inlet bed spread 0.77 m) the peak stage difference was 0.32 m.

If you have a calibrated model with structures on sloping ground, this is the
thing to re-run and compare. Three of our four flat-bed validation cases were
bit-identical, which is exactly why this went unnoticed for so long — the
validation suite couldn't see it.

The upgrade guide has a short snippet that measures how exposed a given model
is (it prints the bed range across each inlet; near zero means unaffected).


What I'd find most useful
-------------------------

* Rudy and David — run a model you know well and tell me if anything moved that
  shouldn't have. Structures on sloping ground are the place to look first
  (see above). A view on whether the removals and deprecations are the right
  call for a 4.0 would also be welcome.
* Petar — two things, and the first matters more. Towradgi is the model I used
  to measure the structure change above, and it came out as the largest effect
  we've seen: 0.32 m peak stage difference, 2340 cells shifted by more than a
  millimetre. That's because its inlets sit on real terrain — median bed spread
  across an inlet is 0.77 m, worst is 4.95 m across 12 triangles at Collins St.
  I believe the new behaviour is the correct one, but you know that catchment
  far better than I do, so a view on whether the changed result looks physically
  sensible would be the single most valuable piece of feedback here.
  Second, the GPU path on your RTX card if you have time — it's had a lot of
  work and only limited hardware coverage.
* Jorge — anything on the parallel/MPI side.
* Ole — anything in general.

Please don't use this for production work; it's a release candidate and 4.0.0rc1
sorts before 4.0.0, so a plain "pip install anuga" won't pick it up.


Known issue
-----------

Windows CI is currently pinning the conda-forge mingw toolchain to an older
build. conda-forge's newer packages produce binaries that abort at startup
("stack smashing detected") — it hit several projects, not just us, and is being
sorted upstream (conda-forge/m2w64-sysroot-feedstock#23). It affects how we
build on Windows, not the wheels you'd install.


Full release notes and the upgrade guide are in the repo:
RELEASE_NOTES.md and UPGRADING_TO_4.0.md on the develop branch.

Feedback in the next few days would be ideal — happy to hold the tag if
something turns up.

Thanks,
Steve
