.. _glossary:

Glossary
========

.. glossary::

   stage
      The absolute water-surface elevation (m), measured from the same vertical
      datum as :term:`elevation`. Water :term:`depth` is ``stage - elevation``,
      so stage is *not* the water depth.

   elevation
      The bed elevation / bathymetry (m).

   depth
      Water depth — the derived quantity ``stage - elevation``. A cell is dry
      where the depth is (near) zero; see :term:`wetting and drying`.

   quantity
      A physical field stored on the mesh (``stage``, ``elevation``,
      ``friction``, ``xmomentum``, ``ymomentum``). Accessed via
      ``domain.quantities[...]``.

   operator
      An object applied each timestep that modifies domain quantities — for
      example rainfall, extraction, culverts or inlets. See
      :doc:`setup_anuga_script/operators`.

   riverwall
      An infinitely-thin internal wall (levee, weir or embankment) aligned with
      triangle edges via :term:`breaklines <breakline>`. It blocks or controls
      flow until overtopped. See :doc:`setup_anuga_script/riverwalls`.

   breakline
      A polyline that the mesh generator forces triangle edges to follow — used
      to align the mesh with :term:`riverwalls <riverwall>`, channels or roads.

   yieldstep
      The interval (s) at which :meth:`~anuga.Domain.evolve` returns control to
      the Python loop (and, by default, writes output). The internal
      :term:`CFL`-limited timestep is usually much smaller.

   finaltime
      The simulation time (s) at which :meth:`~anuga.Domain.evolve` stops. Give
      either ``finaltime`` or :term:`duration`, not both.

   duration
      The length of time (s) to evolve for, measured from the current simulation
      time — an alternative to :term:`finaltime`.

   outputstep
      An optional coarser interval (s) at which the state is written to the
      :term:`SWW file`, while :meth:`~anuga.Domain.evolve` still yields every
      :term:`yieldstep`. Use it to interact frequently while keeping the ``.sww``
      file small; it should be a whole multiple of the yieldstep.

   CFL
      The Courant–Friedrichs–Lewy stability condition, which bounds the explicit
      timestep by the mesh size and wave speed. ANUGA chooses the inner timestep
      automatically to satisfy it.

   wetting and drying
      ANUGA's handling of moving shorelines: a cell is treated as dry once its
      :term:`depth` falls below ``minimum_allowed_height``. See
      :doc:`setup_anuga_script/conventions`.

   SWW file
      The NetCDF output file (``.sww``) holding the time series of
      :term:`quantities <quantity>` on the mesh. Viewed with the
      :ref:`ANUGA Viewer <use_anuga_viewer>`, ``SWW_plotter``, or GIS.

   DE algorithm
      The "discontinuous-elevation" family of flow algorithms (``DE0``, ``DE1``,
      ``DE_ader2``, ``DE2``). See :doc:`setup_anuga_script/flow_algorithms`.
