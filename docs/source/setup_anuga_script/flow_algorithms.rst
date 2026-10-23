
.. currentmodule:: anuga

Choosing a Flow Algorithm
=========================

ANUGA's flow algorithm controls how the shallow water equations are advanced in
time.  All available algorithms use the same 2nd-order discontinuous-elevation
spatial reconstruction (the "DE" family), but differ in their time-stepping
scheme and the aggressiveness of their slope limiters.

.. code-block:: python

    domain.set_flow_algorithm('DE1')   # set before evolve


Quick comparison
----------------

.. list-table::
   :header-rows: 1
   :widths: 14 12 12 12 14 36

   * - Algorithm
     - Time order
     - Flux calls / step
     - CFL
     - Slope limiter β
     - When to use
   * - ``DE0``
     - 1st
     - 1
     - 0.9
     - 0.5 (conservative)
     - Fast exploratory runs; very large or very dry domains.  **Not for wave
       propagation** — see the warning below
   * - ``DE1``
     - 2nd
     - 2
     - 1.0
     - 1.0 (full)
     - **General-purpose default — robust and accurate**
   * - ``DE_ader2``
     - 2nd
     - 1
     - 1.0
     - 1.0 (full)
     - DE1 accuracy at DE0's cost (ADER predictor); **future default — still being validated**
   * - ``DE2``
     - 3rd
     - 3
     - 1.0
     - 1.0 (full)
     - Rarely needed; reserved for highly transient flows

Recommendation
--------------

* **Start with DE1.** It is 2nd-order in both space and time, uses the full
  slope limiter, and is the most widely validated algorithm in ANUGA.  It is
  roughly half the speed of DE0 (two flux evaluations per step), but produces
  materially better results in flows with steep gradients, bores, or rapidly
  varying inundation fronts.

* **Use DE0 when speed matters more than temporal accuracy.** It is
  approximately twice as fast as DE1 per simulated second, but uses a more
  conservative slope limiter (β = 0.5) and first-order time stepping.  Results
  are still 2nd-order in space, so spatial gradients are well resolved; the
  1st-order time stepping introduces a mild diffusivity that is often acceptable
  for long-duration flood studies on large, mildly varying domains.

  .. warning::

     **DE0 damps propagating waves.**  On the ``deep_wave`` validation case,
     where a sinusoid is forced at one edge of the domain, DE0 retains 66 % of
     the forced amplitude by the far edge, against 97 % for DE1 and 101 % for
     ``DE_ader2``.

     Do not use DE0 for tsunami propagation, seiching, tidal resonance or any
     other problem where a wave has to travel a long way without losing
     amplitude.  Use DE1, or ``DE_ader2`` at DE0's cost per step.

     β = 0.5 also costs DE0 accuracy on steady flow down a slope: on
     ``rundown_mild_slope_coarse`` it carries 88 % of the prescribed discharge
     where DE1 carries 100 %.  β < 1 cannot reconstruct a linear stage exactly
     — the limiter scales the gradient by β on genuinely linear data — so the
     bed the solver reconstructs (``z = w - h``) comes out flatter than the
     real one and gravity drives the flow less hard.

  .. note::

     **Do not "fix" this by raising DE0's β.**  The reduced β is not a
     conservative preference, it is what keeps forward Euler stable.  Run with
     β = 1.0, DE0:

     * goes unstable on ``deep_wave`` — the stage reaches 25–29 against a
       forced amplitude of 1.0, and the momentum 900 against 32;
     * loses monotonicity at the dry dam break — total variation 112 against
       the exact 10.0, with a 0.43 m overshoot where β = 0.5 has none, and
       twelve times the RMS error;
     * gets worse on the dry avalanche, RMS 0.29 against 0.051.

     If you need both stability and wave fidelity, that is what DE1 and
     ``DE_ader2`` are for: they pair β = 1.0 with a second-order time
     integrator that can carry it.

* **DE_ader2 is the future default.** It achieves the same 2nd-order temporal
  accuracy as DE1 but requires only one flux call per timestep (the same as
  DE0) by using an ADER Cauchy–Kovalewski predictor.  The first timestep
  bootstraps with a plain Euler step to establish the initial CFL dt; from the
  second step onwards it runs at full 2nd-order accuracy.  It is still
  being validated across the full range of ANUGA's test cases, so DE1 remains
  the safer choice for critical production work at present.

* **Avoid DE2 in routine use.** The SSP-RK3 scheme adds a third flux call per
  step (50% more expensive than DE1) for only marginal accuracy gains in the
  smooth-solution regime typical of shallow water problems.  Reserve it for
  specialised studies where 3rd-order temporal accuracy is demonstrably needed.


Algorithm details
-----------------

DE0 — Euler, CFL 0.9
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    domain.set_flow_algorithm('DE0')

* **Time stepping:** forward Euler — one extrapolation and one flux call per
  step.
* **CFL factor:** 0.9 (10% safety margin; slightly more conservative than the
  other algorithms).
* **Slope limiter:** β = 0.5, applied to both wet and transitional cells.  The
  reduced β damps spurious oscillations near wet/dry fronts and makes DE0 more
  robust on very heterogeneous terrain.  It also damps genuine waves: see the
  warning above before using DE0 for wave propagation.
* **Minimum allowed height:** 1 × 10\ :sup:`−12` m (effectively zero); water
  can persist in very thin layers.
* **Cost:** ~1× (baseline).


DE1 — Runge–Kutta 2, CFL 1.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    domain.set_flow_algorithm('DE1')

* **Time stepping:** 2nd-order Runge–Kutta (Heun's predictor–corrector) —
  two flux calls per step.
* **CFL factor:** 1.0.
* **Slope limiter:** β = 1.0 (full reconstruction, no artificial limiting in
  smooth regions).
* **Minimum allowed height:** 1 × 10\ :sup:`−5` m.
* **Cost:** ~2× relative to DE0.


DE_ader2 — ADER-2 Cauchy–Kovalewski, CFL 1.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    domain.set_flow_algorithm('DE_ader2')

* **Time stepping:** 2nd-order ADER scheme using a Cauchy–Kovalewski (C-K)
  predictor.  After an initial Euler bootstrap step, each subsequent step
  requires only a single flux call: edge values are advanced to the midpoint
  time level Q\ :sup:`n+1/2` in-place using local SWE derivatives, and the
  flux is computed from those predicted edges.
* **CFL factor:** 1.0.
* **Slope limiter:** β = 1.0 (same as DE1).
* **Minimum allowed height:** 1 × 10\ :sup:`−5` m (same as DE1).
* **Cost:** ~1× (same as DE0 after the bootstrap step).
* **Note:** The predictor reuses the timestep from the previous step
  (``_ader2_prev_dt``), which is the global-minimum CFL dt across all MPI
  ranks.  This ensures consistency in parallel runs.


DE2 — SSP-RK3, CFL 1.0
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    domain.set_flow_algorithm('DE2')

* **Time stepping:** strong-stability-preserving 3rd-order Runge–Kutta — three
  flux calls per step.
* **CFL factor:** 1.0.
* **Slope limiter:** β = 1.0.
* **Minimum allowed height:** 1 × 10\ :sup:`−5` m.
* **Cost:** ~3× relative to DE0.


Algorithm and GPU compatibility
--------------------------------

All four algorithms run on both CPU and GPU (OpenMP target offloading).  See
:doc:`../parallel/use_gpu_offloading` for a full compatibility table and
instructions on enabling GPU acceleration.


Low-Froude correction
---------------------

The standard Kurganov–Noelle–Petrova (KNP) central-upwind flux scheme applies
a numerical diffusion term proportional to the local wave speed.  In
subcritical or nearly-still-water flows — estuaries, tidal flats, ponded
floodplains — the velocity is much smaller than the wave speed
(:math:`\mathrm{Fr} = |\mathbf{u}| / \sqrt{gh} \ll 1`) and this diffusion can damp slow-moving features
more than the physical solution warrants.

The ``low_froude`` setting reduces the momentum diffusion term by a factor
``local_fr`` computed from the local Froude number at each edge:

.. code-block:: python

    domain.set_low_froude(0)   # default — no correction
    domain.set_low_froude(1)   # aggressive clamping for strongly subcritical flows
    domain.set_low_froude(2)   # smooth scaling, recommended for tidal / estuary

.. list-table::
   :header-rows: 1
   :widths: 10 20 70

   * - Value
     - Name
     - Effect on momentum diffusion scaling factor ``local_fr``
   * - ``0``
     - ``LOW_FROUDE_OFF``
     - ``local_fr = 1.0`` — no correction; standard KNP diffusion
   * - ``1``
     - ``LOW_FROUDE_1``
     - ``local_fr = sqrt(clamp(Fr², 0.001, 1.0))`` — clamps between ~0.032
       and 1.0; aggressive reduction for strongly subcritical flows
   * - ``2``
     - ``LOW_FROUDE_2``
     - ``local_fr = sqrt(clamp(Fr², 0.01, 1.0))`` with smooth floor —
       stays at 0.1 for Fr < 0.01 then ramps smoothly to 1.0;
       gentler and recommended for mixed-regime problems

Only the **momentum** diffusion terms are scaled; the **mass** flux is never
modified, so conservation is preserved exactly.

**When to use it**

* **0 (default):** suitable for most coastal, flood, and dam-break problems
  where the flow regularly passes through or near critical (Fr ≈ 1).

* **1:** strongly subcritical flows where Fr stays well below 0.1 throughout
  (e.g., estuarine hydraulics, slow tidal inundation).  The aggressive clamping
  provides maximum diffusion reduction but can introduce mild oscillations near
  steep gradients.

* **2:** the recommended non-zero choice for problems with mixed flow regimes
  (some fast, some very slow reaches).  The smoother ramp avoids the step-change
  behaviour of mode 1 near Fr ≈ 0.001.

The setting is independent of the flow algorithm and can be combined with any
of DE0, DE1, DE_ader2, or DE2:

.. code-block:: python

    domain.set_flow_algorithm('DE1')
    domain.set_low_froude(2)


.. seealso::

   :ref:`api_domain` — full API of ``set_flow_algorithm``,
   ``get_flow_algorithm`` and ``set_low_froude``.
