Examples
========

.. currentmodule:: anuga

The examples below range from a minimal introductory script through to
real-world flood and tsunami validation cases.  Each notebook can be
downloaded and run locally; the script example is also available in the
``examples/simple_examples/`` directory of the repository.

.. only:: html

.. toctree::
   :maxdepth: 1

   script_simple_example
   notebook_simple_example
   notebook_create_domain_from_regions
   notebook_create_domain_with_riverwalls
   notebook_flooding_example
   notebook_tsunami_benchmark
   notebook_sediment_transport


.. list-table::
   :header-rows: 1
   :widths: 30 50 20

   * - Example
     - Description
     - Level
   * - :doc:`Simple script <script_simple_example>`
     - Wave runup on a linearly sloping bed using a rectangular mesh.
       Introduces the core ANUGA workflow: mesh creation, ``set_quantity``,
       boundary conditions, and the evolve loop.  No external data files
       required.
     - Beginner
   * - :doc:`Simple notebook <notebook_simple_example>`
     - The same runup scenario presented as an interactive Jupyter notebook
       with inline animation via ``Domain_plotter``.  A good starting point
       for exploring ANUGA interactively.
     - Beginner
   * - :doc:`Domain from regions <notebook_create_domain_from_regions>`
     - Builds an unstructured mesh from a bounding polygon with multiple
       tagged boundary segments using ``create_domain_from_regions``.
       Demonstrates mesh refinement with ``maximum_triangle_area`` and
       multi-tag boundary assignment.
     - Intermediate
   * - :doc:`Domain with riverwalls <notebook_create_domain_with_riverwalls>`
     - Adds infinitely-thin levee/riverwall structures to a mesh using
       ``breaklines`` (to align triangle edges with the wall) and
       ``create_riverwalls``.  Compares flow with and without the wall in
       place.
     - Intermediate
   * - :doc:`Flooding example <notebook_flooding_example>`
     - Real-world urban flood model of Merewether, Newcastle (NSW) driven
       by DEM topography, an ``Inlet_operator`` hydrograph, and building
       exclusion polygons.  Results are compared against benchmark
       observational data using ``SWW_plotter``.
     - Advanced
   * - :doc:`Tsunami benchmark <notebook_tsunami_benchmark>`
     - Okushiri Island tsunami runup (IWTS Benchmark 2) with incoming wave
       specified via
       ``Transmissive_n_momentum_zero_t_momentum_set_stage_boundary``
       and spatially interpolated bathymetry.  Simulated water-surface
       time series are compared against physical model gauge data.
     - Advanced
   * - :doc:`Sediment transport <notebook_sediment_transport>`
     - Suspended sediment with ``add_sediment_fraction``: scour of a sand
       bed in a sloping channel, with the concentration and bed change
       plotted on the mesh, then sediment settling out of still water
       compared against its reference solution (the same case as the
       ``sediment_settling`` validation test).
     - Intermediate

.. only:: html
