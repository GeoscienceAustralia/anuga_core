.. currentmodule:: anuga

.. _sediment:

Sediment transport
==================

This documents the user-facing interface: every parameter, its units, its
default, and how to choose between the alternative methods. It assumes you
know ANUGA, and it does not derive the physics.

Labels like :spec:`E-1` name a term in the physics; each links to the
equation that defines it in :doc:`/appendices/physics_spec`, the
specification the implementation was written from, which is published with
these docs. Section numbers such as "spec 2.6" point into the same page.

A sediment fraction **is** a tracer with settling parameters attached, so
:ref:`tracers` covers the transport, boundary and conservation machinery that
this page builds on.

Verification evidence for these terms is in
``anuga/shallow_water/tests/test_sediment_*.py``; runnable
examples are in ``examples/sediment/``.

--------------

The shortest useful program
---------------------------

.. code-block:: python

   import anuga

   domain = anuga.rectangular_cross_domain(40, 10, len1=100.0, len2=25.0)
   domain.set_flow_algorithm('DE0')
   domain.set_quantity('elevation', lambda x, y: -0.01 * x)
   domain.set_quantity('stage', lambda x, y: -0.01 * x + 0.5)   # 0.5 m deep
   domain.set_quantity('friction', 0.03)
   domain.set_boundary({'left':   anuga.Dirichlet_boundary([0.8, 1.0, 0.0]),
                        'right':  anuga.Dirichlet_boundary([-0.5, 0.0, 0.0]),
                        'top':    anuga.Reflective_boundary(domain),
                        'bottom': anuga.Reflective_boundary(domain)})

   domain.add_sediment_fraction('sand', diameter=2.0e-4)   # <- the only new line

   for t in domain.evolve(yieldstep=1.0, finaltime=30.0):
       pass

   print('peak concentration : %.3f' % domain.get_tracer('sand').max())

After 30 s this scours 5-22 cm from the bed and carries a peak concentration of
about 0.145.

.. note::

   The stage follows the bed (``-0.01 * x + 0.5``) and the boundaries drive
   flow through the channel. A *flat* stage over a sloping bed is hydrostatic
   -- the water does not move, the bed shear stress never reaches the Shields
   threshold, and nothing erodes. Sediment needs flow; if a model reports zero
   concentration everywhere, check the hydrodynamics before the sediment
   settings.

:meth:`~anuga.Domain.add_sediment_fraction` is the entry point. One call gives you a transported
concentration, erosion, deposition, the settling velocity, the bed exchange,
and the limiters, with defaults chosen for a sand bed. It creates the sediment
operator and registers the fractional step for you, so there is nothing else to
wire up.

Everything below is about changing those defaults.

Print what you configured
-------------------------

.. code-block:: python

   print(domain.sediment_summary())

This is the single most useful call in the interface. It reports the active
configuration -- every law selected, every scalar in force, and each grain
size's derived settling velocity -- as text. For the program above, every line
but the first is a default:

::

   sediment configuration
     fractions          : 1  --  sand (d=0.0002 m)
     erosion            : Shields / Smith-McLean, non-cohesive (sand, gravel), gamma0=0.0024   [E-1]
     deposition         : D = c_b v_s with c_b carried   [D-1] under [D-5]
     near-bed d*        : constant, per fraction
     shear closure      : quadratic drag, tau_b = rho f_c |v|^2   [T-1]
     friction closure   : constant n, from the domain friction quantity
     bedload            : off
     bed evolution      : True  (Phase 4, evolving)
     porosity lambda    : 0.3
     morphological M    : 1
     c_max      [L-2]   : 0.3
     c_pack     [L-4]   : 0.65
     rho_w              : 1000 kg/m3
     adaptation [D-5]   : two-layer suspension, near-bed layer 0.2 h, exchange set to the Rouse equilibrium, rate x 1 (tracers registered at evolve)
     velocity profile   : layers advected at their log-law mean speeds
     erodible base [L-5]: none (unlimited depth)
     angle of repose    : off
     erodible region    : whole domain
     ([E-1] and the like are cross-references to the term in the physics;
      see the Sediment physics appendix -- the description before each
      label is the whole story.)
     per fraction:
       sand       d=0.0002 m  v_s=2.6695e-02 m/s  R=1.65  tau_c*=0.04

Settling velocity in particular is *derived*, not set: if ``v_s`` is not what
you expected, the diameter or the fluid properties are not what you thought.
Print this at the top of every run.

Two of those lines are the default suspension closure rather than anything
the program asked for: the two-layer suspension with the log-law velocity
split (:spec:`D-5`), which is also why the deposition line says the near-bed
concentration is *carried* instead of the plain :math:`d^{*} c` of
:spec:`D-1`. :meth:`set_deposition(adaptation='none')
<anuga.Domain.set_deposition>` returns to the instantaneous exchange that was
the default before ANUGA 4.1.

Lines appear only when they apply -- configure the Rouse near-bed profile, for
instance, and an ``a/h floor`` line joins the list. Everything below this point
is about changing one of these defaults, and the summary is how you check that
the change took.

The interface at a glance
-------------------------

Choices are made by naming the **physics**, never by setting a flag. Each
call below links to its entry in the :class:`~anuga.Domain` reference, where
the docstring lists every argument, the values it accepts and what each one
means:

.. list-table::
   :header-rows: 1
   :widths: 58 32 10

   * - call
     - chooses
     - spec
   * - :meth:`initialize_sediment_operator(...) <anuga.Domain.initialize_sediment_operator>`
     - sediment transport on, and the domain-wide parameters
     - 2.2
   * - :meth:`add_sediment_fraction(name, diameter, ...) <anuga.Domain.add_sediment_fraction>`
     - a sediment fraction to carry, and its own properties
     - 2.2
   * - :meth:`set_bed_material(material, ...) <anuga.Domain.set_bed_material>`
     - the erosion law
     - 4.1.1
   * - :meth:`set_deposition(law, near_bed, ...) <anuga.Domain.set_deposition>`
     - the deposition law and near-bed ratio
     - 4.4
   * - :meth:`set_shear_closure(closure) <anuga.Domain.set_shear_closure>`
     - how ``tau_b`` is formed
     - 3.2
   * - :meth:`set_sediment_friction(mode, ...) <anuga.Domain.set_sediment_friction>`
     - the friction factor feeding ``tau_b``
     - 3.3
   * - :meth:`set_bedload(formula, ...) <anuga.Domain.set_bedload>`
     - bedload transport, or off
     - 5
   * - :meth:`set_sediment_parameters(...) <anuga.Domain.set_sediment_parameters>`
     - the scalar physical properties
     - 2.4, 6
   * - :meth:`set_erodible_base(...) <anuga.Domain.set_erodible_base>`
     - the depth below which nothing erodes
     - 4.5
   * - :meth:`set_erodible_region(...) <anuga.Domain.set_erodible_region>`
     - where erosion may act at all
     - 4.5
   * - :meth:`set_bed_composition(...) <anuga.Domain.set_bed_composition>`
     - per-class bed composition (active layer over a substrate)
     - 4.5
   * - :meth:`set_angle_of_repose(...) <anuga.Domain.set_angle_of_repose>`
     - relaxation of over-steep bed slopes
     - 7
   * - :meth:`set_tracer_source(name, values) <anuga.Domain.set_tracer_source>`
     - an external source
     - 2.6
   * - :meth:`set_tracer_boundary(name, tag, value) <anuga.Domain.set_tracer_boundary>`
     - inflow concentration, per boundary tag
     - 2.5

.. seealso::

   :ref:`sediment_physics`
      What each of those calls is choosing between -- the erosion, deposition,
      shear and bedload laws, and how to tell which one your problem wants.
      The defaults are a working sand-bed configuration, so you can leave them
      alone until you need to say otherwise.

   :doc:`/appendices/physics_spec`
      Where the bracketed labels are defined. :spec:`E-1`, :spec:`T-1` and the
      rest name a term in the physics; they appear in this page, in the source,
      and in the output of ``sediment_summary()``, and each links to its
      equation there.

Order does not matter, with one exception noted under :ref:`grain_sizes`
below: call them before ``evolve()``, in whatever order reads best.

**Anything not in that table is internal.** The domain carries roughly fifteen
``sediment_*`` arrays (``sediment_qbx``, ``sediment_settling_velocity``,
``sediment_erosion_mode``, ...) that exist to be handed to the C kernel. Setting
them directly can leave the GPU mapping stale, and no validation runs. Use the
setters; they invalidate the device mapping for you.

--------------

.. _sediment_decision_chart:

Which options does my problem need?
-----------------------------------

Most of the choices follow from a few physical questions about the bed and
the flow. Work down the chart; anything you skip keeps its default.

.. figure:: sediment_decision_chart.*
   :alt: Decision chart: seven questions about the bed and the flow, each
         pointing to the setter and option it selects.
   :width: 100%

   Choosing the sediment options. The numbered list below says the same
   with the reasons.

1. **What is the bed made of?** This picks the erosion law, and it is a
   statement about the material, not a tuning choice (spec 4.1.1). Sand,
   gravel, and mud that is transported as flocs are
   :meth:`set_bed_material('noncohesive') <anuga.Domain.set_bed_material>`.
   The default Smith--McLean law suits a sand or gravel bed with a
   critical Shields stress you trust; ``entrainment='de_leeuw'`` is the
   choice for sand-bedded rivers and deltas (the van Rijn flume cases
   validate it), with ``de_leeuw_fit='nghiem_2022'`` for mud and sand as
   in the Delta-X Wax Lake model. Clay or silt that binds is
   ``set_bed_material('cohesive', tau_crit=...)``, with ``tau_crit`` in
   pascals, ideally from a jet test.

2. **How does it move?** The Rouse number :math:`Z = w_s/(\kappa u_*)`,
   with :math:`\kappa = 0.41`, compares settling with turbulent lifting.
   Above about 2.5 the grains stay near the bed: switch bedload on with
   :meth:`set_bedload('wong_parker_eq24') <anuga.Domain.set_bedload>` (and
   give its ``open_boundaries``). From about 0.8 to 2.5 both modes matter:
   bedload as well as the suspension. Below about 0.8 it is suspension
   and the default (no bedload) is right. Estimate :math:`u_*` from the
   flow you expect, for example :math:`u_* = U n \sqrt{g} / h^{1/6}` from
   Manning, and :math:`w_s` from
   :meth:`~anuga.Domain.settling_velocity`. The table gives :math:`Z` for
   quartz (Ferguson--Church, smooth-sphere constants); a dash means the
   grains do not move at all (Shields stress below 0.04).

   .. list-table::
      :header-rows: 1
      :widths: 22 18 15 15 15 15

      * - grain
        - :math:`w_s` (m/s)
        - :math:`u_*` = 0.02
        - 0.05
        - 0.1
        - 0.2 m/s
      * - silt, 20 µm
        - 0.0004
        - 0.04
        - 0.02
        - 0.01
        - 0.005
      * - very fine sand, 100 µm
        - 0.008
        - 0.98
        - 0.39
        - 0.20
        - 0.10
      * - fine sand, 250 µm
        - 0.038
        - 4.6
        - 1.8
        - 0.92
        - 0.46
      * - coarse sand, 1 mm
        - 0.18
        - --
        - 9.0
        - 4.5
        - 2.3
      * - gravel, 2 mm
        - 0.30
        - --
        - 15
        - 7.3
        - 3.7

3. **Is the suspension stratified?** Deposition goes with the
   concentration near the bed, and that exceeds the depth average as soon
   as :math:`Z` is more than about 0.1: the Rouse ratio :math:`d^*` is
   1.3 to 1.6 there, 2 to 4 at :math:`Z = 0.3`, and 5 to 20 at 0.8,
   depending on the reference height. Only fine silt in fast flow is well
   mixed. Use
   :meth:`set_deposition(near_bed='rouse') <anuga.Domain.set_deposition>`,
   which computes the ratio per cell; otherwise keep the default
   :math:`d^* = 1`. Keep the default ``adaptation='two_layer'`` either way:
   it lets the near-bed concentration lag a change of flow, which is what
   the Wang and Ribberink flume measured.

4. **What sets the bed shear?** For rivers, floodplains and dam breaks keep
   the defaults: ``'quadratic_drag'`` with the domain's own Manning
   :math:`n`. Where there is vegetation, the stems take most of the drag;
   :meth:`~anuga.Domain.set_vegetation_drag` with ``sediment_shear='bed'``
   gives the sediment only the bed's share (see :doc:`vegetation`). The
   other closures, ``'depth_slope'`` and the ``'wilson'`` and
   ``'larsen_lamb'`` friction modes, exist to reproduce anugaSed and
   megaflood studies, not for ordinary floods.

5. **How much of the bed can erode?** By default the bed is bottomless.
   A finite layer over rock, a lining or an apron is
   :meth:`~anuga.Domain.set_erodible_base`; scour confined to one reach
   or a structure that must not scour is
   :meth:`~anuga.Domain.set_erodible_region`; banks or scour holes that
   should slump rather than stand vertical are
   :meth:`set_angle_of_repose(35.0) <anuga.Domain.set_angle_of_repose>`
   (a numerical device; read its notes first).

6. **Several grain sizes, and the bed sorts them?** Mud that should only
   come from where the bed holds mud, a sand channel through a muddy marsh,
   a bed that armours as it loses its fines: use
   :ref:`bed composition <sediment_bed_composition>`, after the erodible
   base. With one grain size, or a bed that is the same mixture everywhere,
   leave it off.

7. **Is the bed much slower than the flow?** Morphology over days to years
   with a flow that adjusts in minutes:
   ``set_sediment_parameters(morphological_factor=M)``, with :math:`M`
   small enough that the flow still sees the bed change slowly. For a
   check against an analytic solution with a fixed bed,
   ``bed_evolution=False``.

Then print :meth:`~anuga.Domain.sediment_summary` and check it says what
you meant.

The chart leads to a *sensible* configuration, not a *calibrated* one. The
coefficients still need data: the van Rijn trench case, for example,
calibrates the de Leeuw constant on the approach flow before it tests
anything (see :ref:`anuga-validation`). The recipes under
:ref:`sediment_choosing_a_configuration` are worked answers to the chart for
common cases.

--------------

Setting sediment up
-------------------

Sediment transport is configured through two calls, and the split between them
is the point: a parameter belongs to exactly one of the two, so there is never
a question of which call wins.

.. code-block:: python

   domain.initialize_sediment_operator(porosity=0.28, rho_w=1000.0)  # the run
   domain.add_sediment_fraction('sand', diameter=2.0e-4)                    # a grain
   domain.add_sediment_fraction('silt', diameter=2.0e-5, tau_c_star=0.11)   # another

:meth:`~anuga.Domain.initialize_sediment_operator` takes what describes the **run**;
:meth:`~anuga.Domain.add_sediment_fraction` takes what describes **one fraction**. Neither accepts the
other's parameters -- passing ``diameter=`` to the first, or ``rho_w=`` to the
second, is a ``TypeError`` rather than a silently ignored argument.

The two may be called in either order, and :meth:`~anuga.Domain.initialize_sediment_operator` is
optional: :meth:`~anuga.Domain.add_sediment_fraction` creates the operator with default domain-wide
parameters if none exists, which is why the program at the top of this page is
a single line.

.. _sediment_operator_init:

``initialize_sediment_operator``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   domain.initialize_sediment_operator(
       porosity=None, c_max=None, c_pack=None, bed_evolution=None, rho_w=None,
       description=None, label=None, logging=False, verbose=False)

Switches sediment transport on and returns the operator. The physical
parameters are the same ones :ref:`set_sediment_parameters <scalar_parameters>`
takes, and are documented there; passing them here is a convenience, and
:meth:`~anuga.Domain.set_sediment_parameters` can still change them afterwards.

**One operator per domain.** Calling it again returns the same operator,
applying any parameters given the second time. That is not just tidiness: the
kernel makes a single pass over every registered fraction, so a second
operator in the fractional-step list would apply the bed exchange twice per
timestep.

.. code-block:: python

   op = domain.initialize_sediment_operator(porosity=0.28)
   op is domain.initialize_sediment_operator(rho_w=1025.0)   # True

The return value is the operator, which you need only for controlling operator
order -- see :ref:`operator_order` below.

.. _grain_sizes:

``add_sediment_fraction``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   domain.add_sediment_fraction(
       name, diameter,
       rho_s=2650.0, tau_c_star=0.04, d_star=1.0, beta=None,
       initial_concentration=0.0, reference_height=None, **settling_kwargs)

.. list-table::
   :header-rows: 1
   :widths: 24 10 12 54

   * - parameter
     - units
     - default
     - meaning
   * - ``name``
     - --
     - required
     - label; also the tracer name
   * - ``diameter``
     - m
     - required
     - grain diameter :math:`d`; sets :math:`v_s` via :spec:`S-1`
   * - ``rho_s``
     - kg/m3
     - 2650
     - sediment density (quartz); enters as
       :math:`R = \rho_s/\rho_w - 1`
   * - ``tau_c_star``
     - --
     - 0.04
     - critical Shields stress :math:`\tau_c^{*}`, :spec:`E-1`
   * - ``d_star``
     - --
     - 1.0
     - near-bed ratio :math:`d^{*} = c_b/c`; 1.0 is well-mixed
   * - ``initial_concentration``
     - --
     - 0.0
     - volumetric :math:`c`, uniform
   * - ``beta``
     - --
     - domain's
     - edge reconstruction limiter, shared by all tracers
   * - ``reference_height``
     - m
     - ``None``
     - Rouse reference height :math:`a` in :spec:`S-4`; see
       :ref:`near_bed_d_star`

Multiple fractions are independent: each has its own concentration, settling
velocity and critical stress, and each exchanges with the same bed. Call it
once per fraction.

There is no ``rho_w`` here. Water density is a property of the fluid, and there
is one fluid, so it lives on :meth:`~anuga.Domain.initialize_sediment_operator` and
:meth:`~anuga.Domain.set_sediment_parameters`. Changing it afterwards recomputes :math:`R` and
:math:`v_s` for every fraction already registered.


What belongs to a fraction, and what to the run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only the grain's own properties are per fraction; every physics choice, and
the bed itself, is shared.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Per fraction
     - Set by
   * - ``diameter``
     - ``add_sediment_fraction``
   * - settling velocity :math:`v_s`
     - computed from ``diameter``, ``rho_s`` and any ``settling_kwargs``
   * - ``tau_c_star``
     - ``add_sediment_fraction``
   * - ``d_star``, ``reference_height``
     - ``add_sediment_fraction``
   * - submerged specific gravity :math:`R`
     - computed from ``rho_s`` and the shared ``rho_w``

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Shared by every fraction
     - Set by
   * - ``porosity``, ``c_max``, ``c_pack``, ``rho_w``, ``bed_evolution``
     - ``initialize_sediment_operator`` or ``set_sediment_parameters``
   * - erosion law, ``tau_crit``, ``K_e``, ``K_partheniades``
     - :ref:`set_bed_material <sediment_physics>`
   * - deposition law, ``tau_d``, near-bed :math:`d^{*}` mode
     - ``set_deposition``
   * - shear closure; friction mode and its parameters
     - ``set_shear_closure``, ``set_sediment_friction``
   * - bedload mode, ``K``, ``m``, bedload's own ``tau_c_star``, the
       boundaries open to it and any prescribed ``supply`` across them
     - ``set_bedload``
   * - ``beta``
     - one ``beta_tracer`` for every tracer on the domain
   * - **the bed**
     - one ``elevation``, which all fractions erode and deposit onto

Because the shared settings live on the **domain**, not on the operator, order
does not matter: :meth:`set_bed_material('cohesive') <anuga.Domain.set_bed_material>` called after both fractions
are registered applies to both.

Fractions occupy tracer slots in call order, so fraction ``s`` is tracer
``s``. **The one ordering rule**: do not interleave :meth:`~anuga.Domain.add_tracer` and
:meth:`~anuga.Domain.add_sediment_fraction` on the same domain if you rely on that correspondence.


.. _operator_order:

Controlling operator order
~~~~~~~~~~~~~~~~~~~~~~~~~~

Fractional-step operators run in the order they are created. If another
operator must run before the sediment one -- an external source that the bed
exchange then consumes, say -- create the sediment operator at the point you
want it in the sequence, and add the fractions afterwards:

.. code-block:: python

   My_source_operator(domain, ...)          # runs first
   domain.initialize_sediment_operator()    # then this

   for nm, d50 in grain_sizes:
       domain.add_sediment_fraction(nm, diameter=d50)

:meth:`~anuga.Domain.add_sediment_fraction` will not displace an operator that already exists, so the
order established here survives however many fractions follow. Ordinary
models do not need this: calling :meth:`~anuga.Domain.add_sediment_fraction` straight away, as
everywhere else on this page, puts the operator in a sensible place by itself.


Choosing a name
~~~~~~~~~~~~~~~

``name`` is a label you choose, not a value from a list. Nothing is inferred
from it: ``diameter`` and the settling parameters do all the work, so a grain
size called ``'boulder'`` with ``diameter=2e-5`` is silt, and behaves like
silt.

.. code-block:: python

   domain.add_sediment_fraction('fine_sand', diameter=1.5e-4)
   domain.add_sediment_fraction('mud',       diameter=2.0e-5)

   domain.get_sediment_names()        # ['fine_sand', 'mud']
   domain.get_tracer('fine_sand')     # its concentration, per cell

It is also the **tracer** name, which is what constrains it. It must be a
non-empty string, unique on the domain, and it may not be

* the name of a quantity -- ``stage``, ``elevation``, ``friction``,
  ``xmomentum``, ``ymomentum``, ``height``, ``x``, ``y``, ``xvelocity``,
  ``yvelocity``. Both a quantity and a tracer are written to the sww as
  ``<name>_c``, so a fraction called ``stage`` would overwrite the stage in
  the output.
* anything beginning ``max_``, which is reserved for the running maxima
  ``Collect_max_quantities_operator`` writes.

Both are refused rather than allowed to corrupt the output. The name also
becomes the sww variable ``<name>_c``, so pick something you will recognise
when you open the file six months later.

.. warning::

   Do not confuse this with ``set_sediment_friction(bed='sand')``. **That**
   ``'sand'`` is one of a fixed set -- ``'sand'``, ``'gravel'``,
   ``'boulder'`` -- selecting a roughness closure, and has nothing to do with
   what you called your fraction. ``name`` is the only sediment argument
   that is free text; every other choice below names a physics option from a
   fixed vocabulary.


If you leave a fraction out
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``name`` and ``diameter`` describe one fraction, so they travel together:
both are required positional parameters, and omitting either is a
``TypeError`` at the call rather than a partial registration.

.. code-block:: python

   domain.add_sediment_fraction('sand')
   # TypeError: add_sediment_fraction() missing 1 required positional argument: 'diameter'

Setting sediment up with **no** fraction at all is legal, and is what
:meth:`~anuga.Domain.initialize_sediment_operator` on its own does -- see :ref:`operator_order`.
It is a configured run waiting for its sediment, not an error. But it is only
useful as a step on the way to :meth:`~anuga.Domain.add_sediment_fraction`: an operator with no grain
sizes transports nothing, and evolving that way warns rather than completing
silently with the bed untouched.

.. note::

   **This includes bedload.** Moving the bed without carrying anything in
   suspension is a reasonable thing to want, and it does not remove the need
   for a fraction: the diameter and :math:`R` that set the Shields stress
   -- and so the transport vector :math:`\mathbf{q}_b` -- live on a grain
   size. With none registered, the bedload kernel returns immediately and the
   bed does not move.

   For bed evolution with no suspended sediment, register the fraction and
   choose the total-load formula, which turns the suspended exchange off:

   .. code-block:: python

      domain.add_sediment_fraction('sand', diameter=2.0e-4)
      domain.set_bedload('engelund_hansen')    # [K-5], total load

   On the channel at the top of this page that scours about 14 cm of bed while
   the concentration stays at zero throughout.


.. _42-choosing-tau_c_star:

Initial and boundary concentrations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   domain.set_tracer('sand', 0.001)          # uniform, or an array of centroids
   domain.set_tracer_boundary('sand', 'inflow', 0.02)   # entering across 'inflow'

Both are volumetric concentration ``c`` (dimensionless), not ``h*c``. The
conserved quantity is ``m = h*c``; the interface works in ``c`` throughout.

--------------

.. _scalar_parameters:

Scalar parameters
-----------------

.. code-block:: python

   domain.set_sediment_parameters(porosity=0.30, c_max=0.30, c_pack=0.65,
                                  bed_evolution=True, rho_w=1000.0,
                                  morphological_factor=1.0)

All optional; only what you pass is changed. All are validated.

.. list-table::
   :header-rows: 1
   :widths: 37 10 15 37

   * - parameter
     - units
     - default
     - meaning
   * - ``porosity``
     - --
     - 0.30
     - bed porosity ``lambda``, :spec:`G-4`. Sediment volume leaving suspension is ``(1-lambda) dz``; the rest is pore space filled from the water column. LM15 use 0.28.
   * - ``c_max``
     - --
     - 0.30
     - :spec:`L-2`, ceiling on depth-averaged concentration (FG21; anugaSed use 0.20).
   * - ``c_pack``
     - --
     - 0.65
     - :spec:`L-4`, maximum packing bounding *near-bed* ``c_b = d* c``. Only bites when ``d* != 1``.
   * - ``bed_evolution``
     - --
     - ``True``
     - whether the bed moves
   * - ``rho_w``
     - kg/m3
     - 1000
     - fluid density used to form dimensional ``tau_b``
   * - ``morphological_factor``
     - --
     - 1
     - morphological acceleration ``M`` (Delft3D's MORFAC): every step's bed change from :spec:`G-4` and :spec:`G-5` is multiplied by ``M``, the erodible base :spec:`L-5` respected, the water column untouched. ``M`` hydrodynamic steps of bed change per step: a 15 h flume becomes 1.5 h with ``M = 10``. Valid while the bed changes little over one hydrodynamic adjustment time; a tidal case wants ``M`` such that the tidal cycles still average out. Bed and water-column sediment budgets then differ by exactly ``M``.

.. _91-bed_evolution-is-the-coupling-stage:

``bed_evolution`` is the coupling stage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``False`` gives a **fixed bed**: sediment is entrained and deposited, and
concentration evolves, but elevation never changes. Choose it when

- comparing against analytic solutions, which assume constant depth;
- comparing against RDycore v1.0, which is configured this way;
- isolating a transport question from a morphology question.

``True`` (the default) evolves the bed through :spec:`G-4` and :spec:`G-5`. Choose it
for any real morphological problem. Both bed terms are applied in a single
fractional step.

--------------

The non-erodible base
---------------------

.. code-block:: python

   domain.set_erodible_base(depth=0.5)          # 0.5 m of erodible material
   domain.set_erodible_base(elevation=z_rock)   # or an absolute surface, (n,)
   domain.set_erodible_base()                   # remove it again

By default the bed is **bottomless**: erosion lowers it for as long as the flow
has the strength to. That is right for a deep alluvial channel and wrong
wherever the erodible layer is finite -- a reach floored by an outcrop, a lined
culvert, a dam apron, a soil layer of known depth over rock. :spec:`L-5` gives it a
floor.

The base is a **per-centroid field**, because bedrock is a surface. ``depth=``
measures down from the elevation set so far and is recorded as an elevation at
the moment of the call, so later changes to the elevation quantity do not drag
it around. ``elevation=`` gives the surface directly, in the domain's datum.
Scalars broadcast; arrays must be ``(n,)``. Give exactly one.

A base above the bed is rejected rather than silently accepted -- it would mean
negative erodible thickness, which is a mistake, not a configuration.

.. code-block:: python

   domain.erodible_thickness()   # (n,) metres remaining; 0 means bedrock

:meth:`~anuga.Domain.sediment_summary` reports the range and how many cells have reached bedrock.

What it guarantees
~~~~~~~~~~~~~~~~~~

The limit is applied to the **source**, not by clamping elevation. Erosion is
scaled back to what the remaining thickness can supply, so the sediment that is
not eroded never enters the water column and the budget still closes to machine
precision. Clamping ``z`` afterwards would leave suspended sediment that came
from nowhere.

Where several classes compete for the last of the material they are scaled by
one shared proportional factor, not served in registration order: by default
the bed is one shared material, so no class has a better claim, and the answer
must not depend on the order you registered the fractions. (With
:ref:`bed composition <sediment_bed_composition>` each class is instead capped
at what the cell holds of it.) Deposition is never scaled -- it is what
replenishes the bed.

The two transport routes give **different strengths of guarantee**, and it is
worth knowing which you are relying on:

.. list-table::
   :header-rows: 1
   :widths: 33 33 33

   * - route
     - floor is
     - why
   * - suspended exchange :spec:`G-4`
     - **exact**
     - the limit is on the exchange term itself
   * - bedload :spec:`G-5`
     - within one step's flux
     - bedload is a divergence; see below

Bedload only redistributes, and stays exactly conservative with a base present,
because the limit is applied to the transport vector and to whole edges -- both
of which the two cells sharing an edge evaluate identically. The price is that
the floor is not exact: closing an edge for a cell that cannot pay also cancels
its neighbour's inflow, so the deficit migrates. Measured overshoot is
5.1e-6 m on a 1.0e-2 m layer. If you need bedload's floor exact, that is a
known limitation with a known fix (iterating the exhaustion flag to a fixed
point), not a mystery.

Restricting erosion to a region
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The base says how *deep* erosion may go; a region says *where* it may happen
at all.

.. code-block:: python

   domain.set_erodible_region(polygon=breach)                 # ONLY here erodes
   domain.set_erodible_region(polygon=apron, erodible=False)  # everywhere BUT here
   domain.set_erodible_region(center=[x, y], radius=25.0)     # a circle
   domain.set_erodible_region(indices=ids)                    # triangles directly
   domain.set_erodible_region(my_region)                      # a Region object
   domain.set_erodible_region()                               # remove it

The keyword arguments are the ones the region-based operators
(``Erosion_operator`` and friends) already take, resolved by the same ``Region``
class, so a polygon that selects a set of cells there selects the same set
here. A region that selects no cells is rejected rather than silently doing
nothing -- that is almost always a polygon in the wrong coordinates.

**Passing a Region is the general form.** ``Region`` understands more than the
keywords above -- ``line=``, ``poly=``, ``expand_polygon=`` -- so build one and hand
it over when you need those:

.. code-block:: python

   from anuga.abstract_2d_finite_volumes.region import Region

   domain.set_erodible_region(Region(domain, line=thalweg))
   domain.set_erodible_region(Region(domain, polygon=reach, expand_polygon=True))

Both do reach through: a line selects the cells it crosses (83 of 960 on a test
mesh), and ``expand_polygon`` genuinely changes the selection (504 cells against
480), because it intersects on vertices rather than centroids.

It must be a ``Region`` built on **this** domain -- one built elsewhere carries
triangle indices that mean nothing here, and is refused rather than silently
mis-selecting. A bare list of points passed positionally is refused too, with a
message pointing at ``polygon=``: taking it as a region would select every cell
and look like it worked.

**Locked means unscourable, not inert.** A locked cell is held at the
elevation it has when you call this, by giving it zero erodible thickness --
the restriction is :spec:`L-5` with the layer set to nothing, not a separate
mechanism. Sediment may still settle onto it, which is what a concrete apron
or a rock bar does in the field, and that new material is erodible again
because it now sits above the base. Under genuinely erosive flow such a cell
sits at exactly net zero: the limiter scales erosion back until it just
cancels deposition, so nothing piles up on a scoured apron.

The two compose, in either order, and neither discards the other:

.. code-block:: python

   domain.set_erodible_base(depth=0.4)        # 0.4 m of erodible material
   domain.set_erodible_region(polygon=reach)  # but only inside this reach

Where they disagree the stricter wins. :meth:`~anuga.Domain.sediment_summary` reports both, and
the thickness range it prints covers only the erodible cells -- locked ones
carry zero thickness and would otherwise drag the minimum to zero whatever the
layer is.

Cost
~~~~

None when unset. With no base configured the kernels take the path they took
before the feature existed, and produce bitwise identical results -- which is
asserted, not assumed, in ``test_sediment_erodible_base.py`` check E1.

--------------

.. _sediment_bed_composition:

Bed composition
---------------

.. code-block:: python

   domain.add_sediment_fraction('mud', diameter=2.0e-5)
   domain.add_sediment_fraction('sand', diameter=1.2e-4)
   domain.set_erodible_base(depth=1.0)
   domain.set_bed_composition({'mud': fmud, 'sand': 1.0 - fmud},
                              active_layer=0.01)

By default the bed is one shared material: every fraction can be entrained
from it wherever the shear allows, so a fine fraction keeps coming out of a
coarse bed for as long as the bed lasts.
:meth:`~anuga.Domain.set_bed_composition` tracks how much of each fraction the
bed holds instead, in a thin **active layer** over a **substrate** (the Hirano
active-layer model). A fraction is then only entrained where the bed holds it,
and a bed that is losing its fines **armours**: its surface coarsens and the
fine entrainment dies away.

When to use it
~~~~~~~~~~~~~~

Use it whenever more than one fraction can be entrained and the bed is not the
same mixture everywhere -- sandy channels beside muddy marsh, a gravel bed with
a fine matrix -- or a run is long enough for armouring to matter. Without it,
a mud fraction is picked up from deep, high-shear channel beds that are really
sand. In the Delta-X Wax Lake case that ran away within a day: suspended mud
reached hundreds of g/L against 0.12 g/L in the river. With composition the
eroded area after one day fell from 37 km² to 2.4 km², and the median suspended
mud stayed at the inflow concentration.

With a single fraction it changes nothing: the share is always 1 and the cap is
the whole thickness, so the answer is the shared-bed one.

Setting it up
~~~~~~~~~~~~~

Call it after every :meth:`~anuga.Domain.add_sediment_fraction`, after
:meth:`~anuga.Domain.set_sediment_parameters` (it uses the porosity) and after
:meth:`~anuga.Domain.set_erodible_base`, whose thickness it divides between the
two layers. A fraction added afterwards is refused, and so is a call with no
erodible base.

``fractions``
   Fraction name to volume fraction of the active layer: a scalar, a
   per-centroid array, or a function ``f(x, y)`` of the centroid coordinates
   (as in ``set_quantity``). Fractions not named are 0. They must sum to 1 in
   every cell, to within 10\ :sup:`-6`. ``None`` switches composition off.
``active_layer``
   The active-layer thickness, in metres of bed (pores included). Default
   0.01 m, the Caltech Wax Lake model's value. For sand and gravel a few times
   the coarse grain size is the usual choice; a thinner layer armours faster.
``substrate``
   The substrate's composition, in the same form. Defaults to ``fractions``.

A cell with no erodible thickness holds nothing, whatever it is given.

**A starting composition from a map.** Arrays are the usual route, for example
the nearest value of a gridded mud-fraction map:

.. code-block:: python

   from scipy.spatial import cKDTree

   xy_map, fmud_map = ...                     # map points (absolute) and values
   cc = domain.get_centroid_coordinates(absolute=True)
   fmud = fmud_map[cKDTree(xy_map).query(cc)[1]]
   domain.set_bed_composition({'mud': fmud, 'sand': 1.0 - fmud},
                              active_layer=0.01)

What happens each step
~~~~~~~~~~~~~~~~~~~~~~

The layers hold **amounts**, the solid volume of each fraction per unit bed
area, not fractions, so each fraction is conserved exactly. In each cell:

1. The entrainment law gives a potential rate; it is multiplied by the
   fraction's share of the active layer. An empty active layer supplies
   nothing.
2. A fraction may not take more out of the cell in one step than the cell
   holds of it, in both layers together. This replaces the shared proportional
   scale of the erodible base.
3. The exchange passes through the active layer, which is then restored to its
   thickness. When the bed builds up, the overflow goes down into the
   substrate with the **active** layer's mix. When it erodes, the active layer
   is topped up from the substrate with the **substrate's** mix.

The equations are in the :ref:`Sediment physics <sediment_physics>` appendix.

Reading the composition
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   F = domain.get_bed_composition('active')     # {'mud': (n,), 'sand': (n,)}
   Fs = domain.get_bed_composition('substrate')

Cells that hold nothing give 0 for every fraction. The host copy is current at
each yieldstep; on a GPU it is copied back with the bed. The composition is not
written to the SWW file, so save what you need at yieldsteps, for example with
``numpy.savez``. :meth:`~anuga.Domain.sediment_summary` reports the active
layer and the mean fractions.

Interactions and limits
~~~~~~~~~~~~~~~~~~~~~~~

- **Compute modes.** The update is in the source kernel shared by both compute
  modes, so legacy and unified (CPU or GPU) give the same answer. It costs
  under 3% of a sediment step.
- **Parallel.** The update is cell-local, and each rank's ghost cells take
  their owners' bed after every step, so a partitioned run conserves each
  fraction and matches the serial one.
- **Morphological factor.** The layers move by the factor times the physical
  exchange, as the bed does.
- **Bedload and angle-of-repose relaxation** still move the bed as one
  material, so the layers and the elevation drift apart where they act.
  Combining them with composition warns.
- **One substrate layer.** There is no stratigraphy below the active layer:
  what goes down is mixed into one substrate.
- **Bed roughness** is still taken per fraction (3 × the grain size), not from
  the active layer's mixture.

--------------

Angle-of-repose relaxation
--------------------------

.. code-block:: python

   domain.set_angle_of_repose(35.0)     # degrees; FG21 use 35
   domain.set_angle_of_repose(None)     # off again (the default)

Where the centroid-to-centroid bed slope exceeds the critical angle, bed
material is moved downslope until it does not. Without it, scour will cut a
vertical face that in the field would collapse.

**It is off by default, and that is a considered default.** FG21, whose
formulation this is, are explicit that it is *a numerical heuristic, not
physics* -- real slope failures are advective. It exists to stop the rest of
the model breaking on over-steep slopes, and it has a side effect worth
knowing before you switch it on: it limits the steepness of canyon walls and
knickpoints, and so suppresses knickpoint retreat that may be real.

**Mass is conserved.** Material removed from an over-steep cell is deposited on
its neighbours, never discarded (measured drift 2.3e-13 m3 on 4.0e2 m3 of bed).
This is the sharpest difference from ANUGA's ``sanddune_erosion_operator``, which
lowers an over-steep cell and lets the material vanish.

It respects :spec:`L-5`: a cell cannot slump away material it is not allowed to
lose, so a locked cell or one already at its base stays put and its neighbours
relax around it.

The sweep count, which will surprise you
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is an explicit diffusion solve, and convergence from a badly over-steep
bed is slow. An over-steep cone (36.8 degrees) needed **793 sweeps** to reach a
30 degree limit from cold.

That is not what the per-timestep cap is sized for. In a running model the bed
is already near-relaxed and each step needs a handful of sweeps. The cap
(``max_sweeps``, default 50) exists for the pathological case, and **hitting it
is not fatal** -- progress carries over, so the bed keeps relaxing on
subsequent steps. It is reported so you know relaxation is lagging rather than
finished:

::

   Sediment_transport_operator: angle-of-repose relaxation hit its 50-sweep cap at
   t = 0.3 s; the bed may still exceed 30.0 degrees.

If you *start* from a bed steeper than the critical angle, expect that on the
first few steps. Either let it settle, or raise ``max_sweeps`` for that run.
``operator.repose_sweeps`` and ``operator.repose_cap_hits`` are available if you
want to watch it, and ``timestepping_statistics()`` prints the sweep count.

``relax`` defaults to 1.0. That is the fastest **stable** setting, not an
aggressive one -- the kernel already divides by the edge count for stability,
and 1.0 converged the cone in 793 sweeps against 2400+ at 0.3. Lower it only if
you see something pathological.

In parallel
~~~~~~~~~~~

The kernel sweeps internally so it stays on the device, and elevation is
exchanged once per timestep afterwards. Spec 7 asks for an exchange per sweep.
The consequence is that relaxation crosses a subdomain boundary one sweep per
*timestep* rather than one per sweep, so a slump spanning a boundary relaxes
more slowly there than it would in serial. Serial and single-subdomain results
are unaffected. See the specification, section 7.1.

--------------

External sources
----------------

.. code-block:: python

   domain.set_tracer_source('sand', values)   # array over centroids, or a scalar

Adds a source to the tracer equation directly (spec 2.6), in units of ``m``
(that is, ``h*c``) per second. This is what the manufactured-solution tests use
to impose an analytic forcing, and it is the hook for anything the bed
exchange does not describe -- a lateral inflow, a point discharge, a
prescribed release.

--------------

Running on the GPU
------------------

.. code-block:: python

   domain.set_compute_mode('unified')

``'legacy'`` is the CPU/OpenMP path; ``'unified'`` is the shared path that can
run on the device. **Both paths share the same kernel file, so the physics is
the same code** -- ``core_kernels.c`` -- and ``test_sediment_gpu.py`` holds
them to agreement.

Nothing about the sediment configuration changes between modes: set it up the
same way and switch the mode.

``'unified'`` selects the *unified* code path; whether that path actually
offloads to a device is a property of the **build**, not of this call. A build
without offload compiles the same kernels under ``CPU_ONLY_MODE`` and runs them
on the host. :meth:`set_compute_mode('unified') <anuga.Domain.set_compute_mode>` therefore does not fail on a
machine with no GPU, but it also does not report one:
``domain.get_compute_mode()`` reads ``'unified'`` either way. To find out what
you are actually running on, ask the build:

.. code-block:: python

   import anuga
   anuga.gpu_offload_enabled()    # True if this build offloads

To confirm kernels are reaching the device on a run, set
``NVCOMPILER_ACC_NOTIFY=1`` in the environment. Polling ``nvidia-smi`` is
unreliable for this -- the sampling interval misses short kernel bursts.

Call :meth:`~anuga.Domain.set_compute_mode` **after** the sediment setup. Each setter
invalidates the device mapping, so configuring sediment after selecting
``'unified'`` simply forces the mapping to be rebuilt.

--------------

.. _coming_from_erosion_operators:

Coming from the erosion operators
---------------------------------

ANUGA has carried a family of erosion operators for a long time --
``Polygonal_erosion_operator``, ``Circular_erosion_operator``,
``Bed_shear_erosion_operator``, ``Flat_slice_erosion_operator``,
``Flat_fill_slice_erosion_operator`` and ``Sanddune_erosion_operator``. They
still work and nothing is scheduled for removal, but they are the expensive
way to evolve a bed: each runs in Python on the host every timestep, and none
is GPU-safe, so under compute mode ``'unified'`` each one forces a
GPU-to-host sync on every RK step.

Measured on 115,200 triangles, the overhead an erosion operator adds above a
plain run is **12.9x that of the entire sediment transport module on the GPU**,
and 3.3x on the CPU -- while modelling less and conserving nothing.

Structural equivalents
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - erosion operator argument
     - sediment equivalent
   * - ``base=``
     - :meth:`set_erodible_base` (``elevation=``), :spec:`L-5`
   * - ``polygon=``, ``center=``/``radius=``
     - :meth:`set_erodible_region`
   * - Sanddune's repose relaxation
     - :meth:`set_angle_of_repose`
   * - eroded material simply disappears
     - ``set_deposition(law='threshold', tau_d=0.0)`` suppresses redeposition

``Bed_shear_erosion_operator``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This one maps almost exactly. It forms ``1000 * 9.81 * d * EN_slope`` -- that
is :math:`\rho g h S` with :math:`S` the **energy** slope -- and then erodes
``de = tau_b / shear_factor * dt``. That is the cohesive Hanson & Simon law
:spec:`E-3`, :math:`E = K_e(\tau_b - \tau_c)`, with no threshold and
:math:`K_e = 1/\texttt{shear\_factor}`:

.. code-block:: python

   domain.add_sediment_fraction('sand', diameter=2.0e-4)
   domain.set_shear_closure('energy_slope')            # tau_b = rho g h S   [T-7e]
   domain.set_bed_material('cohesive', tau_crit=1e-9,
                           K_e=0.5 / shear_factor)     # E = K_e tau_b
   domain.set_deposition(law='threshold', tau_d=0.0)   # no redeposition
   domain.set_erodible_base(elevation=base)
   domain.set_erodible_region(polygon=polygon)

Use ``'energy_slope'`` rather than ``'depth_slope'``: the old operator used the
free-surface slope, and choosing the bed slope instead is what costs the
agreement. For both, :math:`S` is the least-squares gradient of the centroid
values over the cell and its neighbours (one-sided at boundaries), not the
edge values the limiter rebuilds each step, so a wall cell sees the same
slope as its neighbours. ``'depth_slope'`` on an evolving bed feeds back on
itself -- erosion roughens the bed and steeper local slopes erode faster --
which is the closure's own property; anugaSed contains it with the global
clamp described in the specification. ANUGA offers two stated bounds
instead: :meth:`set_shear_closure('depth_slope', max_slope=0.05) <anuga.Domain.set_shear_closure>` caps the slope
per cell, and :meth:`set_shear_closure('depth_slope', freeze_slope=True) <anuga.Domain.set_shear_closure>` takes
the slope from the bed at setup and keeps it, so the closure sees the reach
slope it was written for while the bed evolves under it
(:meth:`~anuga.Domain.bed_slope_magnitude` returns that slope). Prefer the
other two closures when the bed moves and neither bound is wanted. On a sloping channel over 30 s, correlation of the bed-change field
against the original operator:

.. list-table::
   :header-rows: 1
   :widths: 60 40

   * - closure
     - correlation with the old operator
   * - ``'depth_slope'`` -- bed slope, :spec:`T-7`
     - 0.815
   * - ``'energy_slope'`` -- free surface, :spec:`T-7e`
     - **0.964**
   * - ``'energy_slope'`` with ``K_e`` calibrated
     - **0.988**

``K_e`` needs calibrating: the two gradient reconstructions differ, so
``1/shear_factor`` is the right form but not the right constant. Halving it
matched the case above. Calibrate against a run you trust.

The others
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - operator
     - notes
   * - ``Polygonal``, ``Circular``
     - Erode at ``de = |momentum| * dt``. There is no sediment law of that
       form, and it is dimensionally inconsistent -- ``|momentum|`` is
       m\ :sup:`2`/s, so ``de`` is not a length. Treat it as a tuning knob,
       not a rate, and move to the Shields route :spec:`E-1` (the default)
       rather than trying to reproduce it.
   * - ``Flat_slice``, ``Flat_fill_slice``
     - Not erosion: they set elevation to a target value. Use
       ``set_quantity('elevation', ...)`` or ``Set_elevation_operator``.
       Sediment has no equivalent because these are not scour models.
   * - ``Sanddune``
     - Erosion plus repose. The repose half is
       :meth:`set_angle_of_repose`; the erosion half is another excess-shear
       law, as for ``Bed_shear`` above.

.. warning::

   Nothing prevents enabling an erosion operator **and** sediment transport on
   the same domain. Both write ``elevation`` and the bed changes simply add.
   One conserves mass and the other does not, so the sum is unlikely to mean
   anything -- pick one.

--------------

.. _sediment_choosing_a_configuration:

Choosing a configuration
------------------------

If you do not know where to start:

- **Sand bed, flood or dam break, morphology wanted.** Defaults, plus one
  sediment fraction:
  ``domain.add_sediment_fraction('sand', diameter=2e-4)``. Add
  :meth:`set_bedload('wong_parker_eq24') <anuga.Domain.set_bedload>` if the grains are coarse enough to move
  along the bed.
- **Fine cohesive sediment, muddy estuary.** :meth:`set_bed_material('cohesive') <anuga.Domain.set_bed_material>`
  with a ``tau_crit`` you trust, ``d*`` left at 1.0.
- **Deep, slow, stratified flow.** :meth:`set_deposition(near_bed='rouse') <anuga.Domain.set_deposition>`, and give
  each class a ``reference_height``.
- **Shallow flow over gravel.** :meth:`set_sediment_friction('wilson', bed='gravel', grain_size=...) <anuga.Domain.set_sediment_friction>`.
- **Reproducing anugaSed.** :meth:`set_bed_material('cohesive') <anuga.Domain.set_bed_material>` and
  :meth:`set_shear_closure('depth_slope') <anuga.Domain.set_shear_closure>`; see ``examples/sediment/``.
- **Comparing against an analytic solution.**
  :meth:`set_sediment_parameters(bed_evolution=False) <anuga.Domain.set_sediment_parameters>` and leave ``d*`` at 1.0.
- **A finite erodible layer over rock.** :meth:`set_erodible_base(depth=...) <anuga.Domain.set_erodible_base>`, and
  check :meth:`~anuga.Domain.erodible_thickness` afterwards to see where it bit.
- **Scour confined to one structure or reach.** :meth:`set_erodible_region(polygon=...) <anuga.Domain.set_erodible_region>`,
  or ``erodible=False`` to lock an apron while the rest of the domain erodes.
- **A dune or a steep bank that should collapse rather than stand vertical.**
  :meth:`set_angle_of_repose(35.0) <anuga.Domain.set_angle_of_repose>`, and read section 11.1 first.

Then print :meth:`~anuga.Domain.sediment_summary` and check it says what you meant.

--------------

Validation
----------

The tests in ``anuga/shallow_water/tests/test_sediment_*.py`` are
verification: they show the equations are solved correctly. Validation is the
separate claim that they are the right equations for a real flow. The cases for
that are in ``validation_tests/sediment/``. They run nightly with the rest of the
validation suite, and :ref:`anuga-validation` describes each one and how to run
it.

Five cases compare against closed-form solutions:

- ``sediment_settling``: deposition, the sediment mass balance and the bed
  update, for sediment settling out of still water.
- ``sediment_erosion``: the Smith--McLean entrainment law, its threshold, the
  depth-slope closure and the bed lowering, for a bed eroding under still
  water.
- ``sediment_settling_basin``: the steady concentration profile, the bed-rise
  rate and the sediment budget of a settling basin with flow.
- ``sediment_equilibrium_flow``: clear water picking up sediment in normal flow,
  approaching the equilibrium concentration.
- ``sediment_bed_hump``: the classic Exner test, a bed hump migrating under
  Grass bedload.

Four compare against laboratory measurements:

- ``van_rijn_pickup_flume``: clear water picking up sand from a bed, against
  the measured suspended load.
- ``van_rijn_trench``: suspended sand settling into a dredged trench, against
  the measured bed profiles.
- ``wang_ribberink_settling``: a suspension settling out through a bed that
  takes everything that reaches it, which tests deposition alone.
- ``zaragoza_pier_dambreak``: a dam break over a sand bed around a pier,
  against the bed measured by a Kinect sensor.

No field-scale case has been built yet. The notebook
:doc:`/examples/notebook_sediment_transport` walks through the settling case.
