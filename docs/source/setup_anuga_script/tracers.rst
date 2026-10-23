.. currentmodule:: anuga

.. _tracers:

Passive Tracers
===============

A **tracer** is a depth-averaged concentration carried along by the water: a
salinity, a dye release, a pollutant, a temperature proxy. It is *passive* —
it is transported by the flow but does not affect it, so adding one never
changes the hydrodynamics.

Tracers are the foundation suspended sediment is built on: a sediment grain
size **is** a tracer with settling parameters attached, so everything on this
page applies to sediment too. See :ref:`sediment`.


The shortest useful program
---------------------------

.. code-block:: python

   import numpy as np

   import anuga

   domain = anuga.rectangular_cross_domain(20, 20, len1=100.0, len2=100.0)
   domain.set_quantity('elevation', 0.0)
   domain.set_quantity('stage', 1.0)
   domain.set_quantity('xmomentum', 2.0)          # 2 m/s of flow to carry it

   Br = anuga.Reflective_boundary(domain)
   Bt = anuga.Transmissive_boundary(domain)
   domain.set_boundary({'left': anuga.Dirichlet_boundary([1.0, 2.0, 0.0]),
                        'top': Br, 'bottom': Br, 'right': Bt})

   # A tracer, released as a patch in the upstream half
   domain.add_tracer('dye')
   x = domain.centroid_coordinates[:, 0]
   domain.set_tracer('dye', np.where(x < 50.0, 1.0, 0.0))

   # The water flowing in at 'left' is clean. Say so: an unset inflow
   # boundary would also bring c = 0, but silently -- see below.
   domain.set_tracer_boundary('dye', 'left', 0.0)

   domain.set_name('dye_release')
   for t in domain.evolve(yieldstep=1.0, finaltime=60.0):
       print(domain.timestepping_statistics())

   change, flux, discrepancy = domain.check_tracer_conservation('dye')
   print('mass change %g, boundary flux %g, discrepancy %g'
         % (change, flux, discrepancy))

The plume advects downstream and leaves through the open right-hand boundary,
so the mass falls and the boundary flux accounts for it exactly::

   mass change -4999.9, boundary flux -4999.9, discrepancy 8.19e-12

The run writes ``dye_release.sww`` containing a ``dye_c`` variable, one value
per cell per timestep.


What is actually conserved
--------------------------

The quantity the solver integrates is not the concentration but

.. math::

   m = h\,c

the mass per unit area, with :math:`h` the water depth. The concentration
:math:`c` is *derived* from :math:`m` each substep, exactly as height is
derived from stage.

This matters in two places.

* A **dry cell carries no concentration.** When :math:`h` falls below
  ``minimum_allowed_height`` the derived :math:`c` is zero, not the last wet
  value.

* :meth:`set_tracer` reads the **current** depth to form :math:`m = h c`, so
  set ``stage`` and ``elevation`` *before* the tracer, not after. Setting them
  afterwards leaves :math:`m` inconsistent with the depth it was formed from.

Tracers are deliberately **not** :class:`Quantity` objects. They live in one
contiguous ``(n_tracers, N)`` block so the flux kernel can stride them and the
GPU can map them in a single transfer. The practical consequence is that they
do not appear in ``domain.quantities`` and cannot be named in
``quantities_to_be_stored``; the accessors below are the interface.


Registering and setting
-----------------------

.. code-block:: python

   domain.add_tracer('salinity', initial_value=0.035)
   domain.add_tracer('dye')                       # defaults to zero

   domain.set_tracer('dye', 0.5)                  # uniform
   domain.set_tracer('dye', array_of_len_N)       # one value per cell

   domain.get_tracer('dye')                       # the concentration, per cell
   domain.get_tracer_names()                      # in registration order

Order matters. A tracer's slot is fixed at registration, and a sediment grain
size
occupies the tracer slot of the same index, so do not interleave
:meth:`add_tracer` and :meth:`~anuga.Domain.add_sediment_fraction` if you rely on that
correspondence.

A tracer may not be named after a quantity. Both are written to the sww as
``<name>_c``, so a tracer called ``stage`` would overwrite the stage in the
output; :meth:`~anuga.Domain.add_tracer` refuses those names, and names beginning ``max_``, which
are reserved for running maxima.

Registering a tracer **reallocates** every tracer array, so a reference held to
one beforehand goes stale. Add all the tracers first, then seed their values.


Reconstruction order
--------------------

``beta`` controls the edge reconstruction, analogous to ``beta_w`` for stage:
``0`` selects first order, ``> 0`` a limited second order.

**One value is shared by every tracer** — the C domain carries a single
``beta_tracer`` scalar — so passing a second, different value is an error
rather than a silent last-writer-wins:

.. code-block:: python

   domain.add_tracer('a', beta=1.0)
   domain.add_tracer('b', beta=0.0)     # ValueError: beta is shared


Boundary concentrations
-----------------------

What a tracer brings in across a boundary is set per boundary tag:

.. code-block:: python

   domain.set_tracer_boundary('salinity', 'ocean', 0.035)
   domain.set_tracer_boundary('salinity', 'river', 0.0)

   # time varying: re-evaluated each timestep
   domain.set_tracer_boundary('salinity', 'tide',
                              lambda t: 0.035 + 0.001 * math.sin(t / 4000.0))

   # or one value per edge of that tag, ordered as domain.tag_boundary_cells[tag]
   domain.set_tracer_boundary('salinity', 'river', per_edge_array)

**The value is only used where water flows IN.** The flux kernel picks the
upwind concentration edge by edge from the sign of the water flux:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Flow at the edge
     - Concentration used
   * - inflow (:math:`n \cdot U` into the domain)
     - the boundary value set here
   * - outflow (:math:`n \cdot U` out of the domain)
     - the interior edge value

So there is nothing to set for an outflow, and nothing that *can* be set —
prescribing a concentration on an outflow would over-determine the advection,
and the kernel ignores it. A boundary that only ever lets water out needs no
call at all, and one that alternates gets Dirichlet-on-inflow /
transmissive-on-outflow automatically: the characteristic condition, with no
switch to throw.

.. warning::

   **An unset boundary brings in** :math:`c = 0`. The array is zero-filled at
   :meth:`add_tracer`, so with no call an inflow carries clean water.

   That is a modelling assumption, not a neutral default. For salinity or
   suspended sediment it is usually wrong, and it is invisible in the output —
   the run completes and the numbers look reasonable. Set it explicitly
   wherever water enters.

On a distributed sub-domain a tag this rank owns no part of is ignored, so the
same call works on every rank. In serial an unknown tag is an error, where it
can only be a typo.


.. _tracers_at_inlets:

Tracers at inlets, structures and rate operators
------------------------------------------------

An :class:`~anuga.Inlet_operator` adds or removes water inside the domain, and
the water it moves carries tracer:

.. code-block:: python

   # river inflow at a known sediment concentration (constant or f(t))
   anuga.Inlet_operator(domain, river_line, Q=river_discharge,
                        tracer_concentrations={'mud': 4.8e-5,
                                               'sand': lambda t: sand_rating(t)})

   # an outlet needs nothing: the water it pumps out takes its tracer with it
   anuga.Inlet_operator(domain, outlet_line, Q=-500.0)

* **Inflow.** The added water carries the concentration given for each tracer,
  a constant or a function of time (averaged over each step, as ``Q`` is).
  Tracers not named come in at :math:`c = 0`, the same modelling assumption as
  an unset boundary.
* **Extraction.** The operator treats its region as one level pool, so the
  water leaves at the pool's mean concentration and the cells left behind keep
  that concentration. Draining the region dry removes all of its tracer. In
  parallel the pool spans every rank holding part of the region.

Both are exact: the domain's tracer changes by precisely what the water carried
across the inlet, in either compute mode and in parallel. Each operator keeps
the running totals in ``op.tracers.total_in`` and ``op.tracers.total_out``,
per tracer name. (Before 4.1 an outlet removed only water, so its tracer
stayed behind and concentrated without bound as the cells drained.)

The inlet exchange is not a boundary flux, so a domain with inlets does not
satisfy the identity in the next section; add the inlets' net transfers to the
boundary flux integral to close the budget.

**Structures** -- culverts (:class:`~anuga.Boyd_box_operator`,
:class:`~anuga.Boyd_pipe_operator`), weirs, bridges and internal boundaries --
carry tracer through with the water they move, with nothing to set. Each cell
of the inflow region that loses water keeps its concentration, and the tracer
it gives up arrives in the outflow region with the water, spread over its
cells in proportion to the water each gains. Tracer is conserved exactly, a
structure that moves no water leaves the tracer untouched, and in parallel it
is carried across sub-domains when the two ends are on different ranks. The
running total is in ``op._structure_tracers.total_moved`` on the host paths.
(Before 4.1 a culvert moved water only: the tracer stayed upstream and clean
water arrived downstream.)

Structures move water inside the domain, so unlike inlets they do not change
the domain's tracer mass and the identity below still holds.

**Rate operators** add water over a region (rain, a distributed inflow) or
take it away (a negative rate). What the tracers do is yours to choose:

.. code-block:: python

   # rain carrying a pollutant; clean rain is the default
   anuga.Rate_operator.rainfall(domain, rate=20.0,
                                tracer_concentrations={'nitrate': 2.0e-6})

   # a drain or pump: the tracer leaves with the water
   anuga.Rate_operator(domain, rate=-1.0e-4, polygon=drain,
                       tracer_extraction='carry')

   # evaporation: water leaves, salt stays (the default, 'retain')
   anuga.Rate_operator(domain, rate=-2.0e-8, tracer_extraction='retain')

   # one choice per tracer: dissolved nitrate infiltrates, sediment is filtered
   anuga.Rate_operator(domain, rate=-5.0e-6,
                       tracer_extraction={'nitrate': 'carry', 'mud': 'retain'})

With ``'carry'`` a cell keeps its concentration as it loses water, so draining
it dry removes all of its tracer. With ``'retain'`` the tracer stays however
much water leaves; a cell drained dry keeps its tracer, which reappears when
the cell rewets. Both are exact per cell and work in both compute modes. A
rate operator with no concentration given and every tracer retaining does
exactly what it did before 4.1, at no extra cost.


Checking conservation
---------------------

The tracer counterpart of the water balance. For a domain with no source
terms,

.. math::

   \text{mass}(t) - \text{mass}(0) = \text{boundary flux integral}(t)

to within the timestepping error.

.. code-block:: python

   domain.get_tracer_mass('dye')                   # integral of m = h*c
   domain.get_tracer_boundary_flux_integral('dye') # net across the boundary
   change, flux, discrepancy = domain.check_tracer_conservation('dye')

``discrepancy`` is the conservation error. These return absolute quantities
rather than a relative error: for a tracer that is zero almost everywhere a
relative measure is meaningless, and you know the scale that matters.

Sign follows the water balance: **positive is into the domain**. Nothing needs
setting up — :meth:`add_tracer` registers the operator that accumulates the
boundary flux, and the baseline is taken from the initial condition.

All three are **collective** in parallel, like
:meth:`~anuga.Domain.get_water_volume`: every rank must call them, or the ones
that do will block.


Running maxima
--------------

``Collect_max_quantities_operator`` tracks the running maximum of every
registered tracer alongside stage, depth and speed, stored in the sww as
``max_<name>_c``:

.. code-block:: python

   from anuga.operators.collect_max_quantities_operator import (
       Collect_max_quantities_operator)

   domain.add_tracer('dye')
   Collect_max_quantities_operator(domain, store_to_sww=True)

Register the tracers **before** creating the operator: the names are fixed at
construction, and an sww variable cannot appear after the file is created.

.. note::

   The maxima do not sample the initial condition. Collection happens on the
   fractional-step call, which first runs after the first update, so a cell
   whose maximum is its *starting* value reports the first post-step value
   instead. This is long-standing behaviour for stage, depth and speed; the
   tracers follow it so that every ``max_*`` variable in an sww means the same
   thing.


Output and plotting
-------------------

Each tracer is written as one dynamic centroid variable ``<name>_c``. They are
centroid quantities, so nothing is interpolated to vertices that the solver
never computed.

The file names its own tracers in a ``tracer_names`` attribute — ``salinity_c``
is otherwise indistinguishable from ``stage_c`` — and the readers use it:

.. code-block:: python

   from anuga.utilities import plot_utils as util

   util.get_tracer_names('dye_release.sww')     # ['dye']

   p = util.get_centroids('dye_release.sww')
   p.tracers['dye']                             # (ntimes, ncells)

``tracers`` is a dict rather than an attribute per tracer, so a tracer can
never shadow ``stage`` or ``vel``. ``get_output`` exposes the same dict.

Tracers also appear in the sww GUI's quantity menu as ``tracer_<name>``, and
their maxima as ``max_tracer_<name>``, once a file containing them is loaded.

Set ``domain.store_tracers = False`` to keep them out of the sww.


Parallel and reordered domains
------------------------------

Tracers are carried through :func:`~anuga.distribute` and survive
:meth:`~anuga.Domain.reorder`. Register them and set their values **before** distributing:

.. code-block:: python

   if myid == 0:
       domain = build_domain()
       domain.add_tracer('salinity', initial_value=0.035)
   else:
       domain = None

   domain = anuga.distribute(domain)

   # boundary concentrations AFTER distribute: they are sized by each
   # sub-domain's own boundary
   domain.set_tracer_boundary('salinity', 'ocean', 0.035)

The per-timestep halo exchange is handled for you. Each rank counts only the
boundary edges of cells it owns, so the conservation figures above are
whole-domain values, not per-rank ones.


Compute modes
-------------

Tracers work in both the ``legacy`` and ``unified`` compute modes, and on the
GPU. The two modes agree to round-off, which is the property the GPU tests use
as their oracle.

In ``unified`` mode the tracer state lives on the device between yield steps.
That is invisible for ordinary use, but if you reach into the arrays directly
— ``domain.tracer_centroid_values`` and friends — the host copy will be stale.
Use the accessors, or pin the domain with ``domain.set_compute_mode('legacy')``.
