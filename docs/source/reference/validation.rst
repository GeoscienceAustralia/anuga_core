.. _anuga-validation:

.. currentmodule:: anuga

Validation Test Suite
=====================

ANUGA ships a validation test suite in the ``validation_tests/`` directory of
the repository.  The tests verify the correctness of the numerical solver
against analytical solutions, published experimental data, and benchmark
problems.  They are separate from the unit tests (run with ``pytest``) and are
intended to be run periodically, or whenever a significant change is made to
the solver.

.. contents:: Contents
   :local:
   :depth: 2


Running the validation tests
-----------------------------

Change into the ``validation_tests/`` directory and run the suite runner:

.. code-block:: bash

   cd validation_tests
   python run_auto_validation_tests.py

The runner searches all subdirectories recursively for scripts whose names
start with ``validate_`` and end with ``.py``, then executes them in
alphabetical order.  Each script uses Python's ``unittest`` framework to check
that the numerical solution meets a specified tolerance.

**Selecting a flow algorithm**

The ``-alg`` argument selects the flow algorithm:

.. code-block:: bash

   python run_auto_validation_tests.py -alg DE_ader2

The default is ``DE1``. Any algorithm accepted by
``Domain.set_flow_algorithm`` can be given, for example ``DE0``, ``DE1``,
``DE2`` or ``DE_ader2``.

**Long cases**

The HEC-RAS bridge and weir behaviour cases take much longer than the rest and
are skipped unless ``-l`` (``--long``) is given:

.. code-block:: bash

   python run_auto_validation_tests.py -alg DE1 -l

**Verbose output**

.. code-block:: bash

   python run_auto_validation_tests.py -v

**Running a single test**

Each ``validate_*.py`` script can also be run in isolation from its own
directory:

.. code-block:: bash

   cd validation_tests/analytical_exact/dam_break_dry
   python validate_dam_break_dry.py

**Producing a PDF report**

Each test directory contains a ``produce_results.py`` script that runs the
simulation, plots the results and typesets that case's ``report.pdf``. To
produce the combined report for all cases (excluding the large case
studies):

.. code-block:: bash

   cd validation_tests/reports
   python validations_produce_results.py -alg DE1 -l

This writes ``validations_report_alg_DE1.pdf``. It needs ``pdflatex`` and
``bibtex`` with the ``datatool``, ``placeins`` and ``titlesec`` packages (on
Debian or Ubuntu, ``texlive-latex-extra``). It exits non-zero, after still
typesetting what it could, if any case fails to produce its results.

Add ``-nr`` (``--no-run``) to rebuild the plots and the report from the
``.sww`` files already in the case directories, without running the
simulations again.

To build the combined report locally, run each algorithm and collect its
figures before running the next (each run overwrites the case directories),
then combine:

.. code-block:: bash

   cd validation_tests/reports
   for alg in DE0 DE1 DE_ader2; do
       python validations_produce_results.py -alg $alg -l
       python validations_combine_algs.py collect -alg $alg
   done
   python validations_combine_algs.py combine

This writes ``validations_report_combined.pdf``. New cases go in
``validations_report_body.tex``, which both reports share.

Validation reports for each release
------------------------------------

The PDF reports for every ANUGA release are built from the tagged commit, for
the ``DE0``, ``DE1`` and ``DE_ader2`` flow algorithms, by the
``validation-report.yml`` workflow. The main one is the **combined report**
(``anuga-<version>-validation-combined.pdf``): every case once, a summary
table of each case's outcome and run time per algorithm, and each figure
that differs between the algorithms shown as a row of panels, one per
algorithm. The full report for each algorithm on its own is published
alongside it. They are attached to the release on the
`GitHub releases page <https://github.com/anuga-community/anuga_core/releases>`_,
and for final releases (not release candidates) they are also archived on
Zenodo as a version of the *ANUGA validation reports* record, with its own
DOI. The Zenodo record is linked from the GitHub release. Cite the version
that matches the ANUGA release you used (see :doc:`../citing`).

The workflow is a release gate: if any case fails to produce its results, no
reports are attached.


Nightly run in CI
-----------------

The suite also runs every night in GitHub Actions (``validation.yml``), on
Linux with one Python, against both ``main`` and ``develop``, including the
long HEC-RAS behaviour cases. The unit tests say the code does what the code
says; this run says the solver still reproduces the benchmarks it is
validated against, which a limiter or wet/dry change can break while every
unit test stays green.

Each run's job summary shows the ``VALIDATION SUMMARY`` table, and the full
log plus the plots each case wrote are kept as an artifact for two weeks. A
failing scheduled run opens, or comments on, an issue titled *Validation
tests failing*. The workflow can also be started by hand from the Actions
tab for any branch, algorithm (``-alg``) and with or without the long cases,
which is the way to validate a feature branch before merging it.

Test suite structure
---------------------

The tests are organised into six categories:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Directory
     - Contents
   * - ``analytical_exact/``
     - Comparisons against closed-form analytical solutions to the shallow
       water equations.  These tests are fast and provide the most rigorous
       correctness check.
   * - ``experimental_data/``
     - Comparisons against physical laboratory experiments.  Results are
       validated against measured water-surface time series and runup
       extents.
   * - ``case_studies/``
     - Larger, real-world simulations (e.g. patong beach tsunami).  These
       can take many hours and require 8–16 GB of RAM.  They are excluded
       from the standard automated run.
   * - ``behaviour_only/``
     - Tests that verify qualitative behaviour against results from a
       previous ANUGA run (no independent analytical reference).
   * - ``other_references/``
     - Comparisons against results published in the literature other than
       the experimental datasets above.
   * - ``sediment/``
     - Sediment transport: closed-form solutions for settling, entrainment,
       a settling basin, equilibrium load and a migrating bed hump, and four
       laboratory flumes with measured data.


Benchmark descriptions
-----------------------

Dry avalanche — analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``analytical_exact/avalanche_dry/``

**Physical scenario:** A mass of water on a dry, frictionless inclined plane
released from rest.  The exact solution for the runup and rundown of the
water front is known analytically.

**What is checked:** The simulated water-surface profile and front position at
several output times are compared against the analytical solution.  The
normalised RMS error must be below a set tolerance.

**Key techniques:** Rectangular cross mesh, ``Dirichlet_boundary`` inflow,
wet/dry front tracking.

**Run time:** < 1 minute.

Wet avalanche — analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``analytical_exact/avalanche_wet/``

**Physical scenario:** Same inclined-plane geometry as the dry avalanche but
with a pre-existing thin layer of water on the slope.  The analytical solution
accounts for the interaction between the incoming and resident water.

**What is checked:** Water-surface profile versus analytical solution at
multiple times.

**Run time:** < 1 minute.

Dry dam break — analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``analytical_exact/dam_break_dry/``

**Physical scenario:** Instantaneous removal of a dam separating still water
from a dry bed.  The Riemann solution provides exact depth and velocity
profiles for the resulting shock and rarefaction wave.

**What is checked:** Depth and velocity profiles at a fixed output time are
compared against the Ritter (1892) analytical solution.  Peak error must be
below tolerance.

**Key techniques:** Rectangular cross mesh, zero-friction flat bed, initial
depth step as ``Dirichlet_boundary``.

**Run time:** < 1 minute.

Wet dam break — analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``analytical_exact/dam_break_wet/``

**Physical scenario:** Same as the dry dam break but with a shallow layer of
water on the downstream side.  The exact solution includes a reflected shock
wave in addition to the rarefaction fan.

**What is checked:** Depth and velocity profiles versus the exact Stoker (1957)
wet dam-break solution.

**Run time:** < 1 minute.

Okushiri Island tsunami runup — experimental data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``experimental_data/okushiri/``

**Physical scenario:** The 1993 Okushiri Island tsunami (Benchmark 2 from the
Third International Workshop on Long-Wave Runup Models, Catalina Island 2004).
A physical scale model was used to measure water-surface time series at gauges
5, 7, and 9 and maximum runup around the island.

**What is checked:** Simulated water-surface elevation time series at the three
gauge locations are compared against the experimental measurements.  The
normalised error at each gauge must be below a set tolerance.

**Key techniques:** Unstructured mesh from ``Benchmark_2.msh``, bathymetry
from ``Benchmark_2_Bathymetry.asc``, incoming wave from
``Benchmark_2_input.tms`` via
``Transmissive_n_momentum_zero_t_momentum_set_stage_boundary``, gauge
extraction via ``file_function``.

**Data source:** `<http://www.cee.cornell.edu/longwave/>`_ (Third IWLWRM
benchmark problem 2).

**Run time:** 10–30 minutes depending on mesh resolution and hardware.


Sediment settling in still water — analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``sediment/sediment_settling/``

**Physical scenario:** A tank of still water 1 m deep carries a uniform
volumetric concentration of 0.01 of 100 micron sediment. There is no flow, so
nothing is entrained; the sediment settles onto the bed at the deposition
rate :math:`D = d^{*} v_s c` with a constant near-bed profile factor, and the
bed rises by the deposited volume over the packing fraction. The depth falls
as the bed rises, which feeds back on the depth-averaged concentration; the
reference solution integrates that coupled system to machine precision, and
reduces to a plain exponential when the feedback is dropped.

**What is checked:** The cell-mean concentration and bed rise against the
reference at every output time (relative :math:`L^1` error below 1%); that
every cell follows the same curve; that the free surface does not move; and
that the sediment mass in the water column plus the bed stays at its initial
value. It exercises the deposition term, the sediment mass balance and the
bed update of the sediment transport operator.

**Key techniques:** ``add_sediment_fraction`` with ``initial_concentration``,
``set_deposition(near_bed='constant')``, reading the tracer and the
time-varying bed back from the SWW file.

**Run time:** a few seconds.


Sediment entrainment into still water — analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``sediment/sediment_erosion/``

**Physical scenario:** A tank of still water 1 m deep over a plane bed of
slope :math:`10^{-3}`, with 0.5 mm sand and no initial load. There is no
flow, so the case selects the depth-slope shear closure with its slope
frozen at the setup bed, under which the bed shear stress of a cell depends
only on its depth. The Smith-McLean entrainment law :math:`E = v_s E^{*}`
fires and the concentration relaxes towards :math:`c_{eq} = E/(d^{*} v_s)`
while the bed lowers by the entrained volume over the packing fraction; as
the depth grows so does the stress, and the reference integrates that
coupled per-cell system to machine precision. A second fraction of the same
grain with a critical Shields stress of 10 is the threshold control and
must stay at zero. The slope is frozen because the closure otherwise feeds
on the roughness that uneven erosion creates, leaving no closed-form
reference.

**What is checked:** Every cell's concentration and bed lowering against its
own curve at every output time (relative :math:`L^1` errors below 1%,
worst cell included, so the slope estimate must reproduce the plane along
the walls and in the corners); that the control fraction is never
entrained; that the bed feedback is visible against the fixed-bed
exponential; that sediment mass is conserved between water column and bed;
and that the free surface and the momenta do not move. It exercises the
entrainment law and its threshold, the depth-slope closure, the
entrainment-deposition balance in the sediment mass equation, and the bed
update under erosion.

**Key techniques:** ``set_shear_closure('depth_slope', freeze_slope=True)``,
``set_bed_material('noncohesive')``, a second ``add_sediment_fraction``
with a high ``tau_c_star`` as the control, reading the tracer and the
time-varying bed back from the SWW file.

**Run time:** a few seconds.

Settling basin with flow — analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``sediment/sediment_settling_basin/``

**Physical scenario:** A straight channel 300 m long carrying a uniform
steady flow (flat frictionless bed, 1 m deep, 0.5 m/s, held exact by
Dirichlet boundaries at both ends). Sediment-laden water enters at the
inflow with 100 micron silt at a fixed concentration; entrainment is
switched off and the load settles at the well-mixed deposition rate. At
steady state the concentration decays along the channel as
:math:`c_0 e^{-x/L_s}` with the settling length :math:`L_s = q/(d^{*}v_s)`
(62.5 m here), and the bed rises where the sediment lands at
:math:`d^{*} v_s c(x)/(1-\lambda)`. The rising bed lowers the depth and
raises the velocity together, so :math:`q` and :math:`L_s` are unchanged.

**What is checked:** The steady concentration profile cell by cell
(relative :math:`L^1` error below 1%, measured 5e-4); the bed-rise rate
over the last 100 s (below 1%, measured 5e-4); that the flow stays the
uniform steady state; that the profile is steady and uniform across the
channel; and the sediment budget, net boundary inflow equals water column
plus packed bed volume, to 1e-3 (measured 1e-14). It exercises the tracer
advection scheme with an inflow boundary concentration, the deposition
term, the bed update under deposition, and the tracer conservation
accounting.

**Key techniques:** ``set_tracer_boundary`` for the inflow concentration,
``add_sediment_fraction(tau_c_star=0.0)`` to switch entrainment off,
``check_tracer_conservation`` for the budget, ``Dirichlet_boundary`` at both
ends of a frictionless channel for an exact uniform flow.

**Run time:** a few seconds.

Equilibrium suspended load in normal flow — analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``sediment/sediment_equilibrium_flow/``

**Physical scenario:** A plane channel 300 m long of slope 2e-4 with
Manning roughness 0.03 carries uniform flow at its normal depth of 1 m,
held by Dirichlet boundaries at both ends. Under the quadratic-drag shear
closure with the Manning friction factor the bed stress in normal flow is
exactly :math:`\rho g h S`, so the Shields stress of the 100 micron sand,
the entrainment rate and the equilibrium concentration
:math:`c_{eq} = E/(d^{*} v_s)` are closed-form. Clear water enters at the
inflow and the load relaxes to the equilibrium along the channel as
:math:`c_{eq}(1 - e^{-x/L_s})`, :math:`L_s = q/(d^{*} v_s)` (59 m here).
The bed is held fixed so the flow and the reference stay exact.

**What is checked:** The steady concentration profile cell by cell
(relative :math:`L^1` error below 1%, measured 2e-4); the outflow
concentration against the reference (1e-3); that the flow stays the normal
flow (depth and discharge within 1e-3); that the profile is steady and
uniform across the channel; and that the fixed bed does not move. It
exercises the quadratic-drag closure with the Manning friction factor,
the non-cohesive entrainment law in flowing water, its balance with
deposition under advection, and ``bed_evolution=False``.

**Key techniques:** ``set_shear_closure('quadratic_drag')`` (the default),
``initialize_sediment_operator(bed_evolution=False)``,
``set_tracer_boundary`` for the clear-water inflow, ``Dirichlet_boundary``
at both ends of a Manning channel for exact normal flow.

**Run time:** a few seconds.

Migrating bed hump under bedload — analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``sediment/sediment_bed_hump/``

**Physical scenario:** The classic Exner test of Hudson and Sweby (2003).
A frictionless channel 1000 m long carries uniform subcritical flow, 10 m
deep at 1 m/s, held by Dirichlet boundaries at both ends, over a flat bed
with a low sin-squared hump. Bedload follows Grass's law
:math:`q_b = A_g u^3` and the bed moves by its divergence; the inflow and
outflow are open to bedload. With the free surface and discharge frozen,
Exner reduces to a scalar conservation law whose solution is carried along
characteristics: the hump migrates downstream, faster the higher the bed,
and its front steepens towards a shock. The run stops at half the shock
time.

**What is checked:** The bed against the characteristic solution cell by
cell at the end of the run (relative :math:`L^1` error below 10% of the
hump volume, measured 5.5% at 5 m cells and 9.9% at 10 m: the bedload
flux is first order); the crest position (within two cells); that the bed
volume is conserved, the flat reaches upstream and downstream do not move,
and the flow stays close to uniform. It exercises the Grass law, the
bedload divergence and bed update, and the open bedload boundaries. This
case is what exposed the instability of the earlier bare centred bedload
flux.

**Key techniques:** ``set_bedload('grass', K=A_g, open_boundaries=(...))``,
``initialize_sediment_operator(bed_evolution=True)``, a registered fraction
whose grain size does not enter (Grass is total load).

**Run time:** about a minute.

van Rijn's sediment pick-up flume — experimental data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``sediment/van_rijn_pickup_flume/``

**Physical scenario:** A Delft Hydraulics flume (van Rijn 1986). Uniform flow
0.25 m deep at 0.67 m/s, initially free of sediment, runs from a rigid bed
onto a bed of 230 micron sand and picks up sand until the suspended load
reaches equilibrium. The suspension is strongly stratified (Rouse number
1.15). The data are the depth-integrated suspended loads measured at four
stations. The bed is held fixed, as in the experiment.

**What is checked:** That the flow stays at normal depth and discharge and
reaches a steady state; that the suspended load increases along the flume;
and, scored the way van Rijn scores a transport formula, that the
geometric mean of the predicted-to-measured load ratio over the four
stations lies between 0.5 and 2 with no station outside 1/3 to 3. It
exercises the de Leeuw et al. (2020) entrainment law and the Rouse near-bed
ratio.

**Key techniques:** ``set_bed_material('noncohesive', entrainment='de_leeuw')``,
``set_deposition(law='d_star', near_bed='rouse')``, a clear-water inflow
through ``set_tracer_boundary``.

**Run time:** about a minute and a half.


van Rijn's migrating trench — experimental data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``sediment/van_rijn_trench/``

**Physical scenario:** Delft Hydraulics flume experiments (van Rijn 1986).
A flow 0.39 m deep at 0.51 m/s, carrying 160 micron sand in equilibrium,
crosses a trench cut across the flume. The flow slows over the trench, the
suspended sand settles on the downstream slope and the floor, and the
trench fills and migrates downstream. The data are the bed profiles
measured after 15 hours for side slopes of 1:3, 1:7 and 1:10. As in van
Rijn's own model, the entrainment constant is calibrated so that the
approach flow is in equilibrium at the measured supply.

**What is checked:** The bed level at the measured points (RMS error below
3 cm, against a fill 8--10 cm deep), and the position of the deepest point
of the fill (within 1.5 m of the measured one). It exercises the settling
lag, deposition on a non-uniform flow, bedload and the Exner coupling.

**Key techniques:** ``set_bedload('wong_parker_eq24', open_boundaries=...)``,
an evolving bed with a morphological factor.

**Run time:** the 1:10 trench takes a few minutes. All three trenches run
only with ``-l`` (``--long``).


Wang and Ribberink's settling flume — experimental data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``sediment/wang_ribberink_settling/``

**Physical scenario:** A Delft flume (Wang and Ribberink 1986) whose only
sediment process is deposition. A uniform flow 0.2155 m deep at 0.558 m/s
carries 100 micron sand in equilibrium onto 16 m of perforated plate.
Every grain that reaches the bed falls through the plate and is not picked
up again. The depth-averaged concentration decays along the test section.
The data are the concentration profiles measured at six or seven stations
in each of two runs.

**What is checked:** The concentration along the test section from 1 m on,
as the RMS of the log of the predicted-to-measured ratio (below 0.15). The
case isolates the deposition closure and its adaptation from everything
else: no entrainment, no bedload, no bed change. It scores the two-layer
closure; the instantaneous closure gives 0.73 and 0.58 on the two runs.

**Key techniques:** ``set_deposition(law='d_star', near_bed='rouse',
adaptation='two_layer')``, ``initialize_sediment_operator(bed_evolution=False)``.

**Run time:** under a minute. Run 2 is included only with ``-l``
(``--long``).


Dam break over a sand bed around a pier — experimental data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Directory:** ``sediment/zaragoza_pier_dambreak/``

**Physical scenario:** The Zaragoza Kinect experiments (Segovia-Burillo et
al. 2026; data open on Zenodo, DOI 10.5281/zenodo.17777387). An 8 cm
dam break runs down a 24 cm wide flume over a 5 cm layer of sand, past a
pier on the centreline: a 3 cm cylinder (case P1) or a rounded rectangle
(P2). The dam break is repeated three times over the evolving bed. A Kinect
sensor over the pier records the bed after each event and the water
surface during it. The flow is supercritical and the transport is bedload.

**What is checked:** The bed change after the first dam break, as an RMS
error over the Kinect window (below 6 mm, against a scatter in the
measured field of about 4 mm). The scour depth at the pier is reported, not
tested. It exercises the Exner coupling under a transient, a rigid base
under a finite sand layer, the angle of repose around a scour hole, and a
wetting front over an erodible bed.

**Key techniques:** a pier as a mesh hole, ``set_erodible_base``,
``set_angle_of_repose``, ``set_bedload('wong_parker_eq24', ...)``.

**Run time:** case P2 is included only with ``-l`` (``--long``).


Adding a new validation test
------------------------------

1. Create a subdirectory under the appropriate category
   (``analytical_exact/``, ``experimental_data/``, etc.).
2. Write a numerical simulation script (e.g. ``numerical_mytest.py``) that
   saves an SWW file.
3. Write a ``validate_mytest.py`` script using Python's ``unittest`` that
   loads the SWW output, extracts the relevant quantities, and asserts that
   errors are within tolerance.
4. Optionally write ``produce_results.py`` (calling
   ``anuga.validation_utilities.produce_report``) and a ``report.tex``
   template for the PDF report generator.

The runner will pick up ``validate_mytest.py`` automatically on the next run.
