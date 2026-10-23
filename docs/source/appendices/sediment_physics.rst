.. _sediment_physics:

.. currentmodule:: anuga

Sediment physics: choosing the laws
===================================

.. note::

   **You can skip this page to begin with.** ``add_sediment_fraction`` picks a
   working set of laws for a sand bed, and :ref:`sediment` shows how to run a
   model with them. This appendix is for when you need to say *which* physics,
   rather than accept the defaults.

This page is a short review of the formulations ANUGA implements: for each
term in the governing equations, what the published alternatives are, what they
assume, and how to tell which one your problem wants. Sources are cited by
label -- [FG21]_, [RDy26]_ and so on -- and collected in `References`_ at the
end.

.. note::

   **Two kinds of reference appear on this page, and only one of them is
   something you can look up.**

   *Citations* -- [FG21]_, [aSM16]_, [Rou37]_ and the rest -- are published
   papers, listed in full under `References`_. Those are the sources for the
   physics, and they are where to go to check a formulation.

   *Spec numbers* -- "spec 4.1.1", "§9.5" -- point into
   :doc:`physics_spec`, the specification the implementation was written
   from, which is published alongside this page. Its section numbering is the
   one used in the source comments too, so the code, this page and the
   specification all refer to a term the same way.

   *Bracketed labels* -- :spec:`E-1`, :spec:`T-7` and so on -- name a **term
   in the physics**, not a reference. Each one links to the equation that
   defines it in the specification, and the same labels appear in the source
   comments and in the output of ``sediment_summary()``, so a term can be
   traced from the code to the equation and on to the paper. Every label is
   also an index entry, grouped under **physics label** in the
   :ref:`genindex`, which is the reliable way to look one up: the site search
   splits on the bracket and hyphen, so it will not find :speclit:`T-7` as
   typed. The verification evidence is in
   ``anuga/shallow_water/tests/test_sediment_*.py``.

Where the sources disagree, they disagree about physics rather than notation,
and the page says so. Erosion is the clearest case: the cohesive and
non-cohesive routes are not competing fits to the same data but descriptions of
different bed material, so choosing between them is a statement about the bed.

.. _sediment_notation:

Notation
--------

The symbols are those of the specification, and its notation table
(:doc:`physics_spec`, section 1) is the full list. Three things about the
notation are worth knowing before reading on, because they are easy to trip
over and the specification's table does not spell them out.

**Grain-size percentiles.** :math:`D_{50}` is the median grain diameter: half
the bed material by mass is finer. :math:`D_{84}` is the diameter that 84% is
finer than, and so on for any :math:`D_{xx}`. The coarser percentiles are used
where the roughness is set by the largest grains present rather than the
typical one -- which is why the ``'wilson'`` closure below asks for
:math:`D_{50}` on a sand bed but :math:`D_{84}` on gravel or boulders.

.. warning::

   Three symbols are overloaded, by long convention in this literature, and all
   three appear on this page.

   :math:`D` is **grain diameter** in the shear and bedload relations, and the
   **deposition flux** in the mass balance. :math:`m` is the **conserved
   variable** :math:`h\,c` in the transport equation, and the **exponent** in
   the bedload power law :math:`q_b^{*} = K\tau_x^{\,m}`. :math:`S` is the
   **excess-stress ratio** in the Smith & McLean entrainment :spec:`E-1`, and
   the **water-surface slope** in the depth-slope closures :spec:`T-7` and
   :spec:`T-7e`.

   Which is meant is unambiguous from the equation, but they are worth
   flagging.

**Sign convention.** :math:`E` is positive from the bed into the water column,
:math:`D` is positive from the water column onto the bed, and the bed rises
under net deposition.


.. _sediment_governing_equations:

The equations being solved
--------------------------

Everything on this page is a choice of closure for one of the terms below. It is
worth reading the equations first: most of the parameters name a term here.

A sediment fraction is a tracer, so it starts from :ref:`the tracer transport
equation <tracer_transport_equation>` -- the conserved variable is mass per unit
area, :math:`m_s = h\,c_s`, for fraction :math:`s = 1 \dots N_s`. The state
vector the solver carries is

.. math::

   \mathbf{U} = \begin{bmatrix} h & uh & vh & m_1 & \dots & m_{N_s}\end{bmatrix}^{T}

**Suspended transport** [RDy26]_, which is [DL09]_ written per fraction. The
tracer equation with a source: what the bed gives up and what settles out of
the water column.

.. math::

   \frac{\partial m_s}{\partial t}
   + \frac{\partial (u\,m_s)}{\partial x}
   + \frac{\partial (v\,m_s)}{\partial y}
   = E_s - D_s + S_{m_s}
   \qquad \text{[G-3]}

:math:`E_s` is entrainment from the bed and :math:`D_s` deposition onto it, both
per fraction. :math:`S_{m_s}` is an optional external supply -- hillslope yield,
a tributary load, rainfall washoff -- and is zero unless you set one.

**Bed evolution.** What leaves the water column arrives at the bed, and the bed
moves by the volume it gains, allowing for pore space:

.. math::

   \frac{\partial z}{\partial t}
   = \frac{1}{1 - \lambda}\sum_{s=1}^{N_s} \bigl(D_s - E_s\bigr)
   \qquad \text{[G-4]}

with :math:`\lambda` the bed porosity, since a deposited volume
:math:`(1-\lambda)\,dz` of grains fills a bed volume :math:`dz`. This is the
Exner equation [Exn25]_, in the form used by [P14]_ and [FG21]_.

A morphological acceleration factor :math:`M` (``morphological_factor``,
Delft3D's MORFAC) multiplies the bed change of every step, here and in the
bedload contribution [G-5], with the erodible-base limiter [L-5] scaled to
match; the water column is left alone. The suspension adapts in seconds and
the bed in hours, so :math:`M` steps of bed change per hydrodynamic step
reaches a morphological time :math:`M` times the simulated one at the same
cost, while the bed changes little over one hydrodynamic adjustment time.
Default 1; the bed and water-column budgets then differ by exactly
:math:`M`.

The sum matters once there is more than one fraction: there is **one** bed,
and every fraction exchanges with it. The kernel accumulates a single
:math:`dz` per cell over :math:`s` and applies it once, so fractions can
offset each other -- sand entraining while gravel deposits leaves the bed
still, though neither process has stopped. By default there is no
per-fraction bed: nothing records which fraction the last millimetre came
from, which is why the erodible base :spec:`L-5` shares a shortage between
fractions proportionally rather than by any order of priority. Bed
composition, below, records it.

**Bed composition.** With :meth:`~anuga.Domain.set_bed_composition` the bed of
each cell is an active layer of thickness :math:`H_A` over a substrate, after
[Hir71]_ and [Par04]_ (see [Blo08]_ for the alternatives). The state is the
solid volume per unit bed area of each fraction in each layer, :math:`A_s` and
:math:`B_s` [m], so that

.. math::

   \frac{1}{1-\lambda}\sum_s \left(A_s + B_s\right) = z - z_\text{base}.

Entrainment of fraction :math:`s` is scaled by its share of the active layer
at the start of the step,

.. math::

   F_s = \frac{A_s}{\sum_j A_j}, \qquad E_s = F_s\, E^{*}_s ,

where :math:`E^{*}_s` is the entrainment law's rate. After the limits
:spec:`L-1` and :spec:`L-2`, the net exchange :math:`S_s = E_s - D_s` is
capped at what the cell holds of the fraction,

.. math::

   S_s \le \frac{A_s + B_s}{\Delta t\, M},

which replaces the shared proportional scale of :spec:`L-5`. The exchange passes
through the active layer, :math:`A_s \leftarrow A_s - M\,\Delta t\,S_s` (a
fraction driven below zero takes the rest from its substrate), and the layer
is then restored to its solid content :math:`L_A = H_A(1-\lambda)`. On
aggradation the overflow moves down at the active composition: each
:math:`A_s` gives the fraction :math:`(\sum A - L_A)/\sum A` of itself to
:math:`B_s`. On degradation the deficit is refilled at the substrate
composition: each :math:`B_s` gives the fraction
:math:`\min\!\left(1, (L_A - \sum A)/\sum B\right)` of itself to :math:`A_s`.
Every transfer moves an amount of one fraction between layers, so each
fraction's total (suspended and both layers) is conserved to round-off.

Under a flow that lifts only a fine fraction, with :math:`F^{sub}` its share
of the substrate, the active layer's fine content :math:`a` obeys

.. math::

   \frac{da}{dt} = -\left(1 - F^{sub}\right) \frac{E^{*}}{L_A}\, a ,

so the surface coarsens and fine entrainment decays exponentially: the bed
armours. The Caltech Wax Lake sediment model [Wan23]_ uses the same scheme with :math:`H_A` = 1 cm; it stores fractions
rather than amounts and has no porosity. Bedload and angle-of-repose
relaxation still move the bed as one material.

**Bedload.** When bedload is switched on it moves the bed too, by the divergence
of the bedload transport vector :math:`\mathbf{q}_b`:

.. math::

   \frac{\partial z}{\partial t}
   = -\frac{1}{1 - \lambda}\,\nabla \cdot
     \sum_{s=1}^{N_s} \mathbf{q}_{b,s}
   \qquad \text{[G-5]}

Both act on the same :math:`z`, and when both are active their contributions sum.

.. note::

   **What is not coupled.** The bed feeds back into the flow -- a moving :math:`z`
   changes the bed slope in the shallow water source term -- but the sediment does
   **not** feed back into momentum: the fluid density is held at :math:`\rho`
   regardless of :math:`c_s`, and the momentum equations are unchanged by the
   presence of sediment.

   Neither [FG21]_ nor [RDy26]_ does that coupling either. [LM15]_ is the
   contrasting approach: a fully coupled system in which the sediment is part
   of the hyperbolic problem rather than a tracer advected by an
   already-computed flux. It is the assumption most likely to fail for very energetic flows,
   where the mixture starts to behave like a debris flow.

   Whether the bed moves at all is itself a choice --
   ``set_sediment_parameters(bed_evolution=False)`` holds :math:`z` fixed and
   solves only the transport equation above.

The rest of this page is the closures: what :math:`E_s`, :math:`D_s`,
:math:`\mathbf{q}_b` and the bed shear stress they all depend on are taken to be.

--------------

Bed shear stress
----------------

Every erosion and deposition rate below depends on the bed shear stress, so this
is the term to get right first. [FG21]_, [RDy26]_ and [P14]_ all specify the
same quadratic drag law, in different notation:

.. math::

   \tau_b = \rho\, f_c\, |\mathbf{v}|^2 \qquad \text{[T-1]}

where :math:`f_c` is the Darcy-Weisbach friction factor divided by eight,
:math:`f_c = f/8`. With a Manning closure :math:`f_c = g n^2 h^{-1/3}`.

From it follow the shear velocity, the dimensionless (Shields) stress, and the
excess stress that drives transport:

.. math::

   u_* = \sqrt{\tau_b/\rho} = |\mathbf{v}|\sqrt{f_c}
   \qquad \text{[T-2]}

.. math::

   \tau^{*} = \frac{\tau_b}{(\rho_s - \rho)\, g\, D} = \frac{u_*^2}{R\,g\,D}
   \qquad \text{[T-3]}

.. math::

   \tau_x = \tau^{*} - \tau_c^{*} \qquad \text{[T-4]}

with :math:`R = \rho_s/\rho - 1` the submerged specific gravity and
:math:`\tau_c^{*}` the critical Shields stress for that grain size.

Velocity is recovered from momentum with ANUGA's depth-limiting form [MR12]_,
so that a vanishing depth does not produce a divergent velocity. [RDy26]_
adopts the same form:

.. math::

   u = \frac{(uh)\,h}{h^2 + h_\epsilon^2}
   \qquad \text{[T-5]}

with :math:`h_\epsilon` a small regularisation depth: it leaves :math:`u`
unchanged where the water is deep and drives it smoothly to zero as the cell
dries, instead of dividing by a depth approaching zero.

The alternative depth-slope closure, :math:`\tau_b = \rho\,g\,h\,S`, is the one
[aSM16]_ used, and is kept for reproducing anugaSed's results.


Two independent choices feed :math:`\tau_b`: how the stress is formed, and what
friction factor goes into it.

.. _71-set_shear_closure----how:

``set_shear_closure`` -- how
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   domain.set_shear_closure('quadratic_drag')   # default
   domain.set_shear_closure('depth_slope')
   domain.set_shear_closure('energy_slope')

.. list-table::
   :header-rows: 1
   :widths: 26 50 24

   * - value
     - expression
     - label
   * - ``'quadratic_drag'``
     - :math:`\tau_b = \rho\, f_c\, |\mathbf{v}|^2`
     - :spec:`T-1`
   * - ``'depth_slope'``
     - :math:`\tau_b = \rho\, g\, h\, S`, :math:`S` from the bed
     - :spec:`T-7`
   * - ``'energy_slope'``
     - :math:`\tau_b = \rho\, g\, h\, S`, :math:`S` from the free surface
     - :spec:`T-7e`

``'quadratic_drag'`` is the default and the right choice for unsteady or
rapidly varying flow -- dam breaks, floods, anything with significant
inertia.

``'depth_slope'`` assumes locally uniform flow, where friction balances gravity.
It is what anugaSed uses, so choose it when reproducing their results
(divergence **D1** in the spec). It degrades where that balance does not hold.

:spec:`T-7e` is :spec:`T-7` with the equilibrium assumption dropped. Under the
shallow-water assumption the free surface *is* the energy grade line, so where
the bed slope is a poor proxy for it -- backwater, a pool-riffle sequence, a
flat bed drawing down, a dam break -- the free-surface slope is what actually
drives the flow. It is also what the older ``Bed_shear_erosion_operator`` used,
which makes it the closure to pick when reproducing a model built on that
operator; see :ref:`coming_from_erosion_operators`.

The three are interchangeable by construction: the kernel returns
:math:`\tau_b/\rho`,
so everything downstream is unchanged by the choice.

.. _72-set_sediment_friction----what:

``set_sediment_friction`` -- what
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   domain.set_sediment_friction('constant')    # default
   domain.set_sediment_friction('wilson', bed='gravel', grain_size=0.02)
   domain.set_sediment_friction('larsen_lamb', k_s=0.05, r_d=2.0, r_br=2.0)

``'wilson'`` [Wil04]_ and ``'larsen_lamb'`` [LL16]_ are not callable with the
mode alone -- they
require a length scale and refuse without one, rather than inventing a
default:

- ``'wilson'`` needs ``grain_size > 0`` (:math:`D_{50}` for sand,
  :math:`D_{84}` for gravel or boulder -- see :ref:`sediment_notation`);
- ``'larsen_lamb'`` needs either ``k_s`` or ``sigma_br``. There is no universal
  ``sigma_br``: it is site-measured, and [LL16]_ report about 5 m at Moses
  Coulee.

.. list-table::
   :header-rows: 1
   :widths: 24 30 46

   * - mode
     - spec
     - when
   * - ``'constant'``
     - :spec:`T-6`
     - default: ``f_c`` from the domain's Manning ``n``. Ordinary flood and channel work.
   * - ``'wilson'``
     - :spec:`T-8` to :spec:`T-10`
     - depth-dependent, from grain size. Shallow flow over coarse beds, where relative submergence matters.
   * - ``'larsen_lamb'``
     - :spec:`T-13` to :spec:`T-15`
     - partitions total stress into grain and form drag. Bedforms or roughness elements, where only the grain part drives sediment.

``bed`` is ``'sand'``, ``'gravel'`` or ``'boulder'`` -- one curve each,
:spec:`T-8` to :spec:`T-10`. Which grain-size percentile ``grain_size`` should
carry depends on it: :math:`D_{50}` for ``'sand'``, :math:`D_{84}` for
``'gravel'`` and ``'boulder'``.

.. warning::

   ``bed`` selects the *curve*; ``grain_size`` sets the relative submergence
   :math:`h/D` that curve is evaluated at. They are independent inputs, and
   nothing ties them together, so a mismatched pair runs without error and
   quietly gives the wrong friction.

   The gravel and boulder relations are logarithmic in :math:`h/D`, so too
   small a ``grain_size`` inflates the submergence and collapses
   :math:`f_c`. At :math:`h = 1` m, ``bed='boulder'`` gives
   :math:`f_c = 0.0016` at ``grain_size=2e-4`` against :math:`0.031` at a
   plausible ``0.5`` -- a factor of 19 in :math:`f_c`, and therefore in
   :math:`\tau_b`. On a test channel that under-predicted scour four-fold.

   ANUGA warns when ``grain_size`` looks implausible for the chosen ``bed``
   (roughly Wentworth, widened: sand 6e-5 to 2e-3 m, gravel 2e-3 to 0.25 m,
   boulder 0.05 to 10 m). It warns rather than refuses -- an unusual bed is a
   legitimate choice -- but check the pairing before ignoring it.

   Note also that ``grain_size`` is the roughness length scale of the **bed
   surface**, not the diameter passed to :meth:`add_sediment_fraction`.

``grain_size`` (m) is the roughness length
scale; ``k_s`` (m) is the roughness height; ``r_d`` and ``r_br`` (default 2.0) are
Larsen-Lamb's drag partitioning ratios.

This affects **only** the sediment source term. The hydrodynamic friction
operator is untouched, so momentum still sees the domain's Manning ``n``
whatever you choose here.

--------------


Choosing ``tau_c_star``
~~~~~~~~~~~~~~~~~~~~~~~

0.04 is a reasonable default for sand. It is the threshold at which grains
begin to move, and erosion is zero below it, so it sets *when* the bed becomes
active, not how fast. If the bed does not erode when you expect it to, check
this against the Shields curve for your grain size before adjusting anything
else.



--------------

Erosion: naming the bed material
--------------------------------

The three routes are different physics, not different tunings of one law. Which
you want is a statement about the bed material.

**Non-cohesive** -- sand, gravel, boulders. A saturating near-bed reference
concentration after [SM77]_ and [Par98]_, as used by [FG21]_, recovered as a
flux by the settling velocity. The threshold is a Shields stress [Shi36]_:

.. math::

   E^{*} = \frac{0.65\, \gamma_0\, S}{1 + \gamma_0\, S},
   \qquad S = \frac{\tau^{*}}{\tau_c^{*}} - 1
   \qquad \text{[E-1]}

.. math::

   E = v_s\, E^{*} \qquad \text{[E-2]}

with :math:`\gamma_0 = 0.0024` empirical. The 0.65 is the maximum packing
fraction, so :math:`E^{*}` saturates rather than growing without bound.

A second non-cohesive relation, from [dL20]_, is driven by the
*skin-friction* shear velocity over the settling velocity and by the Froude
number, and carries no Shields threshold; it is the one the Delta-X Wax Lake
Delta sediment model uses for sand and for mud transported as flocculated bed
material load [Wan23]_:

.. math::

   \begin{aligned}
   E^{*} &= \frac{A\, X^{\beta}}{1 + 3 A\, X^{\beta}}, \\
   X &= \left(\frac{u_{*,\mathrm{sk}}}{v_s}\right)^{\alpha}
               \mathrm{Fr} - 0.015,
   \qquad \mathrm{Fr} = \frac{U}{\sqrt{g h}}
   \end{aligned}
   \qquad \text{[E-6]}

with :math:`E = v_s E^{*}` as before and :math:`E^{*}` the near-bed
concentration at :math:`0.1\,h`, capped at 1/3. The skin-friction shear
velocity follows from Manning-Strickler on a roughness :math:`k_s`:
:math:`U / u_{*,\mathrm{sk}} = 8.1\,(H_{sk}/k_s)^{1/6}`,
:math:`u_{*,\mathrm{sk}} = \sqrt{g H_{sk} S}`, with the friction slope
:math:`S = \tau_b / (\rho g h)` from the active shear closure and
:math:`H_{sk} \le h`. Two constant sets are provided:

.. list-table::
   :header-rows: 1
   :widths: 22 14 14 14 36

   * - ``de_leeuw_fit``
     - :math:`A`
     - :math:`\alpha`
     - :math:`\beta`
     - source
   * - ``de_leeuw_2020``
     - 4.74e-4
     - 1.5
     - 1.18
     - [dL20]_ Eq 26a, sand and gravel (default)
   * - ``nghiem_2022``
     - 7.04e-4
     - 0.945
     - 1.81
     - as the Delta-X Wax Lake model uses it [Wan23]_, after Nghiem et al.
       (2022); mud as flocs and sand

.. code-block:: python

   domain.set_bed_material('noncohesive', entrainment='de_leeuw',
                           de_leeuw_fit='nghiem_2022',
                           skin_roughness=None)   # default k_s = 3 d per fraction

``skin_roughness`` fixes :math:`k_s` in metres for every fraction ([dL20]_
use :math:`3 D_{84}` of the bed); ``A``, ``alpha``, ``beta`` and ``threshold``
override the chosen constants. Pair it with the Rouse near-bed profile and a
reference height of :math:`0.1\,h`, which is where :math:`E^{*}` is defined.

**Cohesive** -- silt, clay, cohesive bank material; the regime of the Rio
Puerco field data [P13]_ that anugaSed was built for. An excess *dimensional*
stress law, calibrated by jet test [HS01]_ and specified for ANUGA by
[aSM16]_:

.. math::

   E = K_e\,(\tau_b - \tau_c) \qquad \text{[E-3]}

.. math::

   K_e = \frac{0.2 \times 10^{-6}}{\sqrt{\tau_c}}
   \quad [\mathrm{m^3\,N^{-1}\,s^{-1}}],
   \qquad \tau_c = \tau_c^{*}\,(\rho_s - \rho)\, g\, D_{50}
   \qquad \text{[E-5]}

**Partheniades** [Par65]_ -- the form [RDy26]_ uses, also on dimensional
stress, normalised by the threshold:

.. math::

   E = K_e \, \frac{\tau_b - \tau_c}{\tau_c}
   \quad \text{for } \tau_b > \tau_c, \text{ else } 0
   \qquad \text{[E-4]}

Note that :spec:`E-1` is written in Shields stress while :spec:`E-3` and :spec:`E-4`
are in dimensional stress; that is the usual source of confusion between them.


.. code-block:: python

   domain.set_bed_material('noncohesive')   # default
   domain.set_bed_material('cohesive', tau_crit=0.088, K_e=6.742e-7)
   domain.set_bed_material('partheniades', tau_crit=0.088, K_e=...)

The argument is the **material**, not the formula, because these describe
different sediment rather than competing
descriptions of the same sediment. Picking the wrong one is a physics error.

.. list-table::
   :header-rows: 1
   :widths: 26 37 37

   * - material
     - law
     - when
   * - ``noncohesive``
     - :spec:`E-1`/:spec:`E-2` Shields, Smith & McLean (default), or
       :spec:`E-6` de Leeuw (``entrainment='de_leeuw'``)
     - sand, gravel, boulders; mud as flocculated bed material load with
       :spec:`E-6`
   * - ``cohesive``
     - :spec:`E-3` Hanson & Simon
     - clay, silt, consolidated mud
   * - ``partheniades``
     - :spec:`E-4` Partheniades
     - cohesive, where you have a site-calibrated ``K_e``

``tau_crit`` (Pa, default 0.088) and ``K_e`` (m3/N/s) apply to the cohesive
routes only; the non-cohesive route takes its threshold per class from
``tau_c_star`` instead. The default ``K_e = 6.742e-7`` is anugaSed's.

The choice is not a small correction. On the same channel over 30 s, the
non-cohesive route scours 3-6 cm while the cohesive route accretes about a
millimetre -- the sign of the bed change reverses. See
``examples/sediment/README.md``.

--------------

.. _sediment_deposition:

Deposition
----------

Deposition is the settling flux out of the water column: the near-bed
concentration times the settling velocity, after [P14]_ and [FG21]_.

.. math::

   D = c_b\, v_s = d^{*}(Z)\, c\, v_s \qquad \text{[D-1]}

What :math:`c_b` is depends on the adaptation closure below. Taking it as
:math:`d^{*} c`, the equilibrium stratification of the moment, was the
default up to ANUGA 4.0 and is still available as
``set_deposition(adaptation='none')``; since 4.1 the default carries
:math:`c_b` as a state (:spec:`D-5`), because three flume datasets show the
instantaneous form depositing two to four times too fast.

The alternative, from [RDy26]_, is a threshold form with a critical
*deposition* stress, which switches deposition off in flow strong enough to
keep grains suspended:

.. math::

   D = v_s\, c \left(1 - \frac{\tau_b}{\tau_d}\right)
   \quad \text{for } \tau_b < \tau_d, \text{ else } 0
   \qquad \text{[D-2]}

Setting :math:`\tau_d = 0` disables deposition entirely, which is how the
passive-transport benchmarks are run.

Both forms take the near-bed concentration to be the equilibrium one for the
local flow at every instant. The vertical profile actually adjusts over a time
of order :math:`h/(\alpha w_s)` [GV85]_: grains entrained at the bed must
diffuse up the column before they are carried, and grains high in the column
must settle through it before deposition is felt. Their depth-integrated model
relaxes the load toward the same equilibrium at that rate,

.. math::

   \begin{aligned}
   E - D &= \alpha\, v_s\, (c_{eq} - c), \\
   \frac{1}{\alpha} &= \frac{a}{h} + \left(1 - \frac{a}{h}\right)
   \exp\!\left[-1.5\left(\frac{a}{h}\right)^{-1/6} \frac{w_s}{u_*}\right]
   \end{aligned}
   \qquad \text{[D-3]}

with :math:`\alpha` in the closed form of [ADS88]_ and
:math:`c_{eq} = E^{*}/d^{*}`. Both :math:`E` and :math:`D` are scaled by
:math:`\alpha/d^{*}`, so every equilibrium is unchanged and only the transient
slows; :math:`\alpha \to 1` in the well-mixed limit and :math:`h/a` when fully
stratified. It is off by default. The case for it is van Rijn's pick-up flume,
which reaches equilibrium in about 15 depths without it against more than 40
measured, and his migrating trench, which fills about 25 % too fast.

The rate scaling has no memory and is symmetric, and in the trench it
over-corrects. The alternative [D-4] keeps the stratification of the
suspension as a state, the ratio :math:`r_b = c_b/c` carried with the flow
and relaxed toward :math:`d^{*}` over the settling time
:math:`T_D = (z_c - a)/w_s` from the centroid :math:`z_c` of the equilibrium
profile, only when :math:`d^{*}` has risen above it:

.. math::

   \frac{\partial (h r_b)}{\partial t} + \nabla\cdot(h r_b \mathbf{u})
   = h\,\frac{d^{*} - r_b}{T_D}, \qquad D = v_s\, r_b\, c
   \qquad \text{[D-4]}

so a parcel entering slower water deposits first at the stratification it
brought with it. It never acts in uniform or quickening flow, and it does not
represent the slow filling of the upper column over a bed loading clear
water, which needs a layered suspension.

The layered suspension is [D-5]: a near-bed layer :math:`h_1 = a` and the
rest :math:`h_2 = h - a`, the fraction's tracer still carrying the total
:math:`m = hc` and a second tracer the upper layer's :math:`m_2 = h_2 c_2`,

.. math::

   \begin{aligned}
   \frac{\partial m_2}{\partial t} + \nabla\cdot(m_2 \mathbf{u})
   &= K\,(c_1 - c_2) - v_s\, c_2, \qquad D = v_s\, c_1, \\
   K &= v_s \frac{\rho}{1-\rho}, \qquad \rho = \frac{1/d^{*} - a/h}{1 - a/h}
   \end{aligned}
   \qquad \text{[D-5]}

with the exchange :math:`K` chosen so that the two-layer equilibrium
reproduces the Rouse ratio :math:`d^{*}` exactly. Settling moves sediment
down and the exchange up, so a parcel entering slower water settles its
upper-layer load out over about :math:`h_2/w_s`, and a bed loading clear
water fills the near-bed layer first and the depth-averaged load only as
sediment is exchanged up. This is the form that lengthens the adaptation
in both flume cases.

The settling velocity itself is [FC04]_, smooth across the Stokes-to-turbulent
transition and branch-free. [Die82]_ is the more accurate polynomial fit for
natural irregular grains, at the cost of a branchy evaluation:

.. math::

   v_s = \frac{R\,g\,d^2}{C_1 \nu + \sqrt{0.75\, C_2\, R\, g\, d^3}}
   \qquad \text{[S-1]}

with :math:`C_1 = 18`, :math:`C_2 = 0.4` for smooth spheres, and
:math:`1.0`, :math:`1.1` for natural irregular grains.


.. code-block:: python

   domain.set_deposition(law='d_star', tau_d=0.0, near_bed='constant',
                         reference_height_floor=0.01)

``law``
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 37 37

   * - value
     - expression
     - when
   * - ``'d_star'``
     - :math:`D = d^{*}\, c\, v_s`, :spec:`D-1`
     - default; always deposits
   * - ``'threshold'``
     - :spec:`D-2`, deposition only where :math:`\tau_b < \tau_d`
     - when you need deposition suppressed under strong flow

``tau_d`` (Pa) is the threshold for ``'threshold'`` and is ignored otherwise.

.. _near_bed_d_star:

.. _62-near_bed----the-d-ratio:

``near_bed`` -- the ``d*`` ratio
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`d^{*}` is the ratio of near-bed to depth-averaged concentration, and it
is what makes :spec:`D-1` a near-bed law rather than a depth-averaged one. It is a
function of the Rouse number [Rou37]_

.. math::

   Z = \frac{v_s}{\kappa\, u_*}, \qquad \kappa = 0.41
   \qquad \text{[S-2]}

At low :math:`Z` the grain is well mixed through the column and
:math:`d^{*} \to 1` (washload); at high :math:`Z` it concentrates near the bed
and :math:`d^{*} \gg 1`.

``'constant'`` takes :math:`d^{*} = 1` -- a uniform suspension. ``'rouse'``
evaluates the defining ratio given by [aSM16]_ after [DL09]_, obtained by
requiring that the sediment discharge be the depth integral of concentration
times velocity over a Rouse-Vanoni concentration profile and a logarithmic
velocity profile:

.. math::

   d^{*} = \frac{\displaystyle\int_a^h \ln(z/z_0)\, dz}
                {\displaystyle\int_a^h
                   \left(\frac{h-z}{h-a}\cdot\frac{a}{z}\right)^{Z}
                   \ln(z/z_0)\, dz}
   \qquad \text{[S-4]}

with :math:`a` the reference height, in the sense of [vR84]_. ANUGA evaluates a fitted form of this
rather than the integral itself, for the reasons given below.


Deposition is driven by the concentration *at the bed*, but the transported
quantity is depth-averaged. ``d* = c_b/c`` bridges them.

.. list-table::
   :header-rows: 1
   :widths: 23 77

   * - value
     - meaning
   * - ``'constant'``
     - :math:`d^{*}` is whatever each fraction was given
       (:math:`d^{*}=1` is well-mixed). Default.
   * - ``'rouse'``
     - :math:`d^{*}` from the Rouse profile :spec:`S-4`, recomputed per cell
       per step

``'constant'`` with ``d* = 1`` is the well-mixed assumption: fine sediment,
vigorous mixing, shallow flow. It is also what the analytic decay solutions
assume, so use it when comparing against them.

``'rouse'`` is the physical choice when the profile is stratified -- coarser
grains, or deeper and slower flow, where near-bed concentration genuinely
exceeds the mean. It costs an evaluation of the fitted ``d*(Z, a/h)``
polynomial per cell per fraction per step (§9.5 of the spec; 28 terms, maximum
error 0.82% over ``Z`` in [0.01, 2.5], ``a/h`` in [1e-3, 0.15], clamped at the
edges rather than extrapolated).

``reference_height_floor`` (default 0.01) is the floor on ``a/h``. The Rouse
profile is singular as the reference height approaches the bed, so ``a/h`` is
never allowed below this. Lowering it admits more stratification and more
near-bed concentration; it is a numerical guard, not a physical parameter, and
0.01 sits comfortably inside the fit range.

Near-bed concentration is bounded by ``c_pack`` :spec:`L-4` regardless. That bound
exists because equilibrium Rouse ``d*`` at vanishing shear will otherwise
deposit the entire water column in under a second.

.. _adaptation_lag:

``adaptation`` -- the lag of the near-bed concentration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both deposition laws take the near-bed concentration to be the *equilibrium*
one for the local flow at every instant, so the exchange
:math:`E - D = d^{*} v_s (c_{eq} - c)` responds at once to a change in the
flow. The vertical profile actually adjusts over a time of order
:math:`h/(\alpha w_s)`: grains entrained at the bed must diffuse up the column
before they are carried, and grains high in the column must settle through it
before deposition is felt. ``adaptation`` switches on the depth-integrated lag
of [GV85]_, :spec:`D-3` above, with :math:`\alpha` in the closed form of
[ADS88]_. Both erosion and deposition are scaled by :math:`\alpha/d^{*}`, so
every equilibrium concentration is exactly what it was and only the transient
slows: the load adapts over :math:`h/(\alpha w_s)` instead of
:math:`h/(d^{*} w_s)`.

.. code-block:: python

   domain.set_deposition(law='d_star', near_bed='rouse',
                         reference_height_floor=0.1, adaptation='armanini')

.. list-table::
   :header-rows: 1
   :widths: 23 77

   * - value
     - meaning
   * - ``'none'``
     - no lag; the instantaneous exchange. Default.
   * - ``'armanini'``
     - :math:`\alpha(w_s/u_*, a/h)` per cell per fraction from the closed
       form; 1 in the well-mixed limit, :math:`h/a` when fully stratified
   * - ``'constant'``
     - :math:`\alpha` = ``adaptation_alpha`` everywhere
   * - ``'carried'``
     - :spec:`D-4`: the near-bed ratio :math:`r_b = c_b/c` is carried per
       fraction (tracer ``<name>_nearbed_ratio``, registered at the first
       ``evolve``) and relaxed toward :math:`d^{*}` over the settling time
       :math:`(z_c - a)/w_s` when the flow has slowed; deposition uses
       :math:`r_b c`. One-sided, with memory; changes nothing in uniform or
       quickening flow
   * - ``'two_layer'``
     - :spec:`D-5`: a near-bed layer :math:`a` thick and the rest of the
       column, each fraction carrying its upper-layer mass (tracer
       ``<name>_upper``, registered at the first ``evolve``); settling
       down, an exchange up set to reproduce :math:`d^{*}` at equilibrium;
       deposition :math:`v_s c_1`. Lags in both directions; the equilibrium
       is unchanged. ``layer_fraction`` sets the near-bed layer thickness
       (default the reference height) and ``exchange_factor`` scales the
       rate of the partition's relaxation; neither moves the equilibrium.
       ``velocity_profile=True`` advects each layer at its log-law mean
       velocity (:math:`u_1/\bar u = 1 + \ln f_1/L`,
       :math:`u_2/\bar u = 1 - f_1 \ln f_1/((1-f_1)L)`,
       :math:`L = \kappa/\sqrt{f_c}`), so the sediment-rich near-bed
       layer lags the flow and the transport is :math:`\int u c\,dz`
       rather than :math:`\bar u h c`; the equilibrium concentration is
       still unchanged

``rouse_beta='van_rijn'`` divides the Rouse number of the ``'rouse'`` fit
by van Rijn's (1984b) :math:`\beta = 1 + 2 (w_s/u_*)^2` (at most 2), his
allowance for sediment being mixed more strongly than momentum, and
``rouse_scale`` multiplies it; both flatten the profile and lower
:math:`d^{*}`. Van Rijn's trench profiles (1986b, Fig. 17) sit at an
effective Rouse number near 0.5 where the plain value is 0.79; on the bed
profiles themselves the correction did not help, see the case's notes.

Pair ``'armanini'`` with ``near_bed='rouse'``: the lag is the difference
between the equilibrium stratification :math:`d^{*}` and the effective
exchange :math:`\alpha`, and with the well-mixed ``d* = 1`` the option
*speeds up* a stratified suspension rather than slowing it. ``a`` is the
fraction's reference height with the same floor as ``'rouse'``.

**The default.** ``'two_layer'`` with ``layer_fraction=0.2`` and
``velocity_profile=True`` is what :meth:`~anuga.Domain.set_deposition` selects if you say
nothing; ``adaptation='none'`` restores the instantaneous exchange that
was the default up to ANUGA 4.0. The near-bed tracers it registers at the
first ``evolve`` are the only change to a script's tracer list, and they
sit after any tracer of your own.

**Why it exists.** Against van Rijn's flume measurements
(``validation_tests/sediment/van_rijn_*``) the instantaneous
exchange reaches its equilibrium load within about 15 depths of a clear-water
inflow where the flume took more than 40, and fills a dredged trench about
25 % too fast: in a decelerating flow the near-bed concentration is not yet
the equilibrium one because the grains high in the column have not settled
through it. The lag is off by default so that existing results are unchanged.

--------------

Bedload
-------

Bedload moves grains along the bed rather than through the water column.
Following [FG21]_, the transport rate is a power law in the excess Shields
stress, made dimensional by the grain size:

.. math::

   q_b^{*} = K\, \tau_x^{\,m} \qquad \text{[K-1]}

.. math::

   q_b = q_b^{*}\, \sqrt{\left(\tfrac{\rho_s}{\rho} - 1\right) g}\; D^{3/2}
   \qquad \text{[K-2]}

and it moves the bed by its divergence:

.. math::

   \frac{\partial z}{\partial t}
     = -\frac{1}{1-\lambda}\, \nabla \cdot \mathbf{q}_b
   \qquad \text{[K-3]}

:math:`q_b` is a magnitude, so it is directed along the flow to give the
vector :math:`\mathbf{q}_b`, following [Par98]_.

.. list-table::
   :header-rows: 1
   :widths: 34 14 12 14 26

   * - Parameter set
     - :math:`K`
     - :math:`m`
     - :math:`\tau_c^{*}`
     - Notes
   * - [WP06]_ Eq 24, correcting [MPM48]_
     - 3.97
     - 1.5
     - 0.0495
     - the default; bedload only
   * - .. _spec-k-5:

       [EH67]_, as :speclit:`K-5`
     - :math:`0.05/f_c`
     - 2.5
     - 0
     - **total load**, no threshold
   * - .. _spec-k-6:

       Grass, as :speclit:`K-6`
     - :math:`A_g` (given)
     - 3
     - 0
     - **total load**, in velocity, not stress

Grass's law is the one the classic Exner test cases are written in, and the
only one with a closed-form reference: it is dimensional as it stands,

.. math::

   \mathbf{q}_b = A_g\, |\mathbf{u}|^{m-1}\, \mathbf{u} \qquad \text{[K-6]}

so :speclit:`K-2` does not apply, no grain size enters, and :math:`A_g`
(s\ :sup:`2`/m, between 0 and 1) is a calibration with no default.

.. warning::

   Engelund & Hansen is a **total load** relation -- it already includes
   suspended transport. Selecting it therefore *replaces* the suspended source
   rather than supplementing it, and running both would double-count. Note also
   that its :math:`K` is friction-dependent, not a constant, and that it has no
   threshold: subtracting a :math:`\tau_c^{*}` from it would be a different
   model.


.. code-block:: python

   domain.set_bedload('wong_parker_eq24')   # K=3.97, m=1.5, tau_c*=0.0495
   domain.set_bedload('wong_parker_eq23')   # K=4.93, m=1.6, tau_c*=0.0470
   domain.set_bedload('engelund_hansen')
   domain.set_bedload('grass', K=0.001)     # A_g = 0.001, m = 3
   domain.set_bedload('off')                # default
   # bedload passes through these boundaries (zero gradient); all others closed
   domain.set_bedload('grass', K=0.001, open_boundaries=('inflow', 'outflow'))

Bedload :spec:`K-1`-:spec:`K-4` transports sediment along the bed rather than in
suspension, and drives its own bed evolution term :spec:`G-5`. It is **off by
default**: it is a separate transport mode, not a refinement of suspension,
and enabling it changes what the model represents.

Enable it when the grain size is coarse enough that a significant fraction of
the load moves without going into suspension -- sand and gravel beds under
moderate flow. Leave it off for fine, fully suspended sediment.

The two Wong & Parker variants are their corrected Meyer-Peter-Muller fits;
Eq 24 is the default. ``K``, ``m`` and ``tau_c_star`` override the formula's
constants if you have a calibration.

Bedload only redistributes: it moves sediment between cells and conserves the
total exactly. The flux across each edge is centred with Rusanov dissipation
scaled by a bound on the bed-wave speed, :math:`|\partial(\mathbf{q}_b
\cdot \mathbf{n})/\partial z|`; both sides of an edge form the same flux
with opposite sign, which is what makes it antisymmetric and therefore
conservative, and the dissipation is what keeps a migrating bed form from
growing oscillations (see the ``sediment_bed_hump`` validation case). The
scheme is first order in the bed; see ``test_sediment_bedload.py``.

Boundary edges carry no bedload unless their tag is named in
``open_boundaries``, across which the flux is the cell's own
:math:`\mathbf{q}_b \cdot \mathbf{n}` (zero gradient): an outflow carries
bedload away at the rate it arrives, an inflow supplies it at the rate the
first cell carries it off. Name the inflow and outflow of a reach, or the
inflow cell exports and never imports and digs a hole that travels
downstream at the bed-wave speed. Walls stay closed, which is what keeps a
closed domain exactly conservative.

A boundary can instead carry a **prescribed** bedload inflow,
``supply={tag: q_b}`` in m\ :sup:`2`/s volumetric per unit width (the mass
rate divided by the grain density); the tag is opened if it is not already.
Use it where the supply is known, as in a flume fed at a set rate. The
zero-gradient import equals the inflow cell's own export, so a cell that
aggrades under a fixed inflow stage sees its transport and hence its import
rise, a feedback that runs away within hours in the van Rijn trench case; a
prescribed supply does not depend on the cell at all. Give it only to inflow
tags: an outflow edge with a supply set still carries the supply in.

.. code-block:: python

   domain.set_bedload('wong_parker_eq24', open_boundaries=['outflow'],
                      supply={'inflow': 0.01 / 2650.0})   # 0.01 kg/s/m of quartz


.. _sediment_references:

References
----------

Cited by label throughout the page. The short labels are also the ones used in
the source comments and in :doc:`physics_spec`, so a term can be traced from
the code to the paper it comes from. This is the one list for both pages: the
specification cites the same labels and refers here rather than keeping its
own.

.. [ADS88] Armanini, A. and Di Silvio, G. (1988). A one-dimensional model for
   the transport of a sediment mixture in non-equilibrium conditions.
   *Journal of Hydraulic Research*, 26(3), 275-292.

.. [Blo08] Blom, A. (2008). Different approaches to handling vertical and
   streamwise sorting in modeling river morphodynamics. *Water Resources
   Research*, 44, W03415. doi:10.1029/2006WR005474

.. [DL09] Davy, P. and Lague, D. (2009). Fluvial erosion/transport equation of
   landscape evolution models revisited. *Journal of Geophysical Research:
   Earth Surface*, 114, F03007. doi:10.1029/2008JF001146

.. [dL20] de Leeuw, J., Lamb, M. P., Parker, G., Moodie, A. J., Haught, D.,
   Venditti, J. G. and Nittrouer, J. A. (2020). Entrainment and suspension of
   sand and gravel. *Earth Surface Dynamics*, 8, 485-504.
   doi:10.5194/esurf-8-485-2020

.. [Die82] Dietrich, W. E. (1982). Settling velocity of natural particles.
   *Water Resources Research*, 18(6), 1615-1626.

.. [EH67] Engelund, F. and Hansen, E. (1967). *A Monograph on Sediment
   Transport in Alluvial Streams.* Teknisk Forlag, Copenhagen.
   Scanned copy held by the TU Delft repository
   (``uuid:81101b08-04b5-4082-9121-861949c336c9``); redistribution is
   restricted, so access it there rather than circulating a copy.

.. [Exn25] Exner, F. M. (1925). Über die Wechselwirkung zwischen Wasser und
   Geschiebe in Flüssen. *Sitzungsberichte der Akademie der Wissenschaften*,
   Vienna.

.. [FC04] Ferguson, R. I. and Church, M. (2004). A simple universal equation
   for grain settling velocity. *Journal of Sedimentary Research*, 74(6),
   933-937.

.. [GV85] Galappatti, G. and Vreugdenhil, C. B. (1985). A depth-integrated
   model for suspended sediment transport. *Journal of Hydraulic Research*,
   23(4), 359-377.

.. [FG21] Fassett, C. I. and Goudge, T. A. (2021). Modeling the hydrodynamics,
   sediment transport, and valley incision of outlet-forming floods from
   Martian crater lakes. *Journal of Geophysical Research: Planets*, 126,
   e2021JE006979. doi:10.1029/2021JE006979

.. [HS01] Hanson, G. J. and Simon, A. (2001). Erodibility of cohesive
   streambeds in the loess area of the midwestern USA. *Hydrological
   Processes*, 15(1), 23-38.

.. [Hir71] Hirano, M. (1971). River bed degradation with armouring.
   *Proceedings of the Japan Society of Civil Engineers*, 195, 55-65.

.. [LL16] Larsen, I. J. and Lamb, M. P. (2016). Progressive incision of the
   Channeled Scablands by outburst floods. *Nature*, 538(7624), 229-232.
   doi:10.1038/nature19817. Uses ANUGA.
   Uses ANUGA. The finding relevant here: including erosion thresholds gives
   flood discharges five to ten times smaller than full-to-the-brim estimates.

.. [LM15] Liu, X., Mohammadian, A., Kurganov, A. and Infante Sedano, J. A.
   (2015). Well-balanced central-upwind scheme for a fully coupled shallow
   water system modeling flows over erodible bed. *Journal of Computational
   Physics*, 300, 202-218. doi:10.1016/j.jcp.2015.07.043
   A fully coupled hyperbolic system -- shallow water, non-equilibrium
   suspended sediment and Exner bed evolution -- on triangular grids, the
   contrasting approach to the decoupled tracer path taken here.

.. [MPM48] Meyer-Peter, E. and Müller, R. (1948). Formulas for bed-load
   transport. *Proceedings of the 2nd Meeting of the IAHR*, Stockholm, 39-64.

.. [MR12] Mungkasi, S. and Roberts, S. G. (2012). Approximations of the
   Carrier-Greenspan periodic solution to the shallow water wave equations for
   flows on a sloping beach. *International Journal for Numerical Methods in
   Fluids*, 69(4), 763-780. Source of the velocity regularisation ANUGA uses.

.. [P13] Perignon, M. C., Tucker, G. E., Griffin, E. R. and Friedman, J. M.
   (2013). Effects of riparian vegetation on topographic change during a large
   flood event, Rio Puerco, New Mexico. *Journal of Geophysical Research:
   Earth Surface*, 118(3), 1193-1209. doi:10.1002/jgrf.20073

.. [P14] Perignon, M. C. (2014). *A Rolling Stone Gathers No Moss.* PhD thesis,
   University of Colorado Boulder.

.. [Wan23] Wang, D., Salter, G. and Lamb, M. P. (2023). Delta-X: Matlab
   Model for Wax Lake Delta Land Accretion. ORNL DAAC, Oak Ridge, Tennessee,
   USA. doi:10.3334/ORNLDAAC/2309

.. [Par04] Parker, G. (2004). *1D Sediment Transport Morphodynamics with
   Applications to Rivers and Turbidity Currents.* E-book, University of
   Illinois.

.. [Par65] Partheniades, E. (1965). Erosion and deposition of cohesive soils.
   *Journal of the Hydraulics Division, ASCE*, 91(1), 105-139.

.. [Par98] Parker, G. (1998). *1D Sediment Transport Morphodynamics with
   Applications to Rivers and Turbidity Currents.* e-book.

.. [RDy26] Feng, D., Tan, Z., Xu, D., Johnson, J. and Bisht, G. (2026).
   *RDycore-sediment v1.0.* EGUsphere preprint 2026-4859.
   doi:10.5194/egusphere-2026-4859

.. [Rou37] Rouse, H. (1937). Modern conceptions of the mechanics of fluid
   turbulence. *Transactions of the ASCE*, 102, 463-505.

.. [SM77] Smith, J. D. and McLean, S. R. (1977). Spatially averaged flow over a
   wavy surface. *Journal of Geophysical Research*, 82(12), 1735-1746.

.. [Shi36] Shields, A. (1936). *Anwendung der Ähnlichkeitsmechanik und der
   Turbulenzforschung auf die Geschiebebewegung.* Mitteilungen der
   Preussischen Versuchsanstalt für Wasserbau und Schiffbau, Berlin.

.. [Wil04] Wilson, L., Ghatan, G. J., Head, J. W. and Mitchell, K. L. (2004).
   Mars outflow channels: a reappraisal of the estimation of water flow
   velocities from water depths, regional slopes, and channel floor
   properties. *Journal of Geophysical Research: Planets*, 109, E09003.
   doi:10.1029/2004JE002281
   Equations 3-4 relate Darcy-Weisbach to Manning and 13-15 give
   :math:`f` from grain size and depth. Not to be confused with [Wil66]_.

.. [Wil66] Wilson, K. C. (1966). Bed-load transport at high shear stress.
   *Journal of the Hydraulics Division*, 92(6), 49-59.
   doi:10.1061/JYCEAJ.0001562. Source of the erosion coefficient cited by
   [P14]_ equation 3.8. Not the same Wilson as [Wil04]_.

.. [WP06] Wong, M. and Parker, G. (2006). Reanalysis and correction of bed-load
   relation of Meyer-Peter and Müller using their own database. *Journal of
   Hydraulic Engineering*, 132(11), 1159-1168.
   Equation 24. The corrected rates are at most half the original
   Meyer-Peter and Müller values, which is why [FG21]_ use a lower
   :math:`K`.

.. [aSM16] Perignon, M. C. (2016). *Using the Sediment Transport and Vegetation
   Operators in ANUGA.* anugaSed manual,
   https://github.com/mperignon/anugaSed/blob/master/docs/anugaSed_manual.pdf.
   The authoritative specification for the original anugaSed code, which is
   at https://github.com/mperignon/anugaSed.

.. [vR84] van Rijn, L. C. (1984). Sediment transport, part II: suspended load
   transport. *Journal of Hydraulic Engineering*, 110(11), 1613-1641.

.. note::

   Two things in the source that look like citations but are not. ``A22`` and
   ``A23`` are equation numbers in [RDy26]_'s appendix, and ``FP64`` in the GPU
   sources is double-precision floating point.
