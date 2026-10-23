anuga.Quantity
==============

.. currentmodule:: anuga

.. autoclass:: Quantity

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::
   
      ~Quantity.__init__
      ~Quantity.backup_centroid_values
      ~Quantity.bound_vertices_below_by_constant
      ~Quantity.bound_vertices_below_by_quantity
      ~Quantity.compute_gradients
      ~Quantity.compute_local_gradients
      ~Quantity.extrapolate_first_order
      ~Quantity.extrapolate_second_order
      ~Quantity.extrapolate_second_order_and_limit_by_edge
      ~Quantity.extrapolate_second_order_and_limit_by_vertex
      ~Quantity.get_beta
      ~Quantity.get_extremum_index
      ~Quantity.get_gradients
      ~Quantity.get_integral
      ~Quantity.get_interpolated_values
      ~Quantity.get_maximum_index
      ~Quantity.get_maximum_location
      ~Quantity.get_maximum_value
      ~Quantity.get_minimum_index
      ~Quantity.get_minimum_location
      ~Quantity.get_minimum_value
      ~Quantity.get_name
      ~Quantity.get_values
      ~Quantity.get_vertex_values
      ~Quantity.interpolate
      ~Quantity.interpolate_from_edges_to_vertices
      ~Quantity.interpolate_from_vertices_to_edges
      ~Quantity.interpolate_old
      ~Quantity.limit
      ~Quantity.limit_edges_by_all_neighbours
      ~Quantity.limit_edges_by_neighbour
      ~Quantity.limit_vertices_by_all_neighbours
      ~Quantity.maximum
      ~Quantity.minimum
      ~Quantity.plot_quantity
      ~Quantity.save_centroid_data_to_csv
      ~Quantity.save_data_to_dem
      ~Quantity.save_to_array
      ~Quantity.saxpy_centroid_values
      ~Quantity.set_beta
      ~Quantity.set_boundary_values
      ~Quantity.set_boundary_values_from_edges
      ~Quantity.set_name
      ~Quantity.set_values
      ~Quantity.set_values_from_array
      ~Quantity.set_values_from_constant
      ~Quantity.set_values_from_file
      ~Quantity.set_values_from_function
      ~Quantity.set_values_from_geospatial_data
      ~Quantity.set_values_from_lat_long_grid_file
      ~Quantity.set_values_from_points
      ~Quantity.set_values_from_quantity
      ~Quantity.set_values_from_tif_file
      ~Quantity.set_values_from_utm_grid_file
      ~Quantity.set_values_from_utm_raster
      ~Quantity.set_vertex_values
      ~Quantity.smooth_vertex_values
      ~Quantity.update
   
   

   
   
   .. rubric:: Attributes

   .. autosummary::

      ~Quantity.counter
      ~Quantity.vertex_values




.. _quantity-memory-layout:

Memory layout (``qty_type``)
----------------------------

Each ``Quantity`` allocates a subset of its internal arrays based on its
*type*.  The type is resolved in the following priority order:

1. The explicit ``qty_type`` keyword argument passed to ``Quantity.__init__``.
2. ``domain._quantity_type_map[name]`` — the domain's per-name override table.
3. ``'evolved'`` — the default, which allocates every array for full backward
   compatibility.

The four built-in types and their memory footprint (assuming ``float64``,
8 bytes per element, *N* triangles) are:

.. list-table::
   :header-rows: 1
   :widths: 25 60 15

   * - ``qty_type``
     - Arrays allocated
     - Bytes / triangle
   * - ``'evolved'``
     - centroid, edge, explicit update, semi-implicit update, centroid backup.
       x/y-gradient and phi are **lazy** (allocated on first access).
     - 56 N eager
   * - ``'edge_diagnostic'``
     - centroid, edge.  x/y-gradient and phi are lazy.
     - 32 N
   * - ``'centroid_only'``
     - centroid only
     - 8 N
   * - ``'coordinate'``
     - centroid, edge, vertex (eager)
     - 56 N

``vertex_values`` is **lazily allocated** for all types except
``'coordinate'``: the backing array is created on first access and is
transparent to existing code.  Gradient arrays (``x_gradient``,
``y_gradient``) and ``phi`` are also lazy for all types — they are
allocated on first access, which occurs only when explicitly called (e.g.
by the erosion operators), not during a normal DE timestep.

**Default types used by the shallow-water Domain**:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Quantity
     - Type
     - Reason
   * - stage, xmomentum, ymomentum
     - ``'evolved'``
     - Time-stepped conserved quantities — need update arrays
   * - elevation
     - ``'edge_diagnostic'``
     - Static field; gradients lazy (allocated only if erosion operators
       call ``compute_local_gradients``)
   * - height, xvelocity, yvelocity
     - ``'edge_diagnostic'``
     - Derived fields; edge values needed for output, no update machinery
   * - friction
     - ``'centroid_only'``
     - Scalar parameter; never edge-interpolated
   * - x, y
     - ``'coordinate'``
     - Mesh coordinates; vertex values required immediately at construction

For a 10-quantity shallow-water domain the optimised layout uses **~368
bytes/triangle** versus **800 bytes/triangle** if every quantity were
``'evolved'`` — a saving of ~54 %.

**Example — adding a custom tracer with minimal memory:**

.. code-block:: python

   import anuga
   from anuga.abstract_2d_finite_volumes.quantity import Quantity

   domain = anuga.rectangular_cross_domain(100, 100)

   # Option 1: pass qty_type directly
   tracer = Quantity(domain, name='tracer', register=True,
                     qty_type='centroid_only')
   tracer.set_values(0.0)

   # Option 2: register via the domain's type map, then create with Quantity()
   # Note: set_quantity() and add_quantity() operate on *existing* quantities
   # only — they do not create a new one.  Use Quantity(..., register=True).
   domain._quantity_type_map['salinity'] = 'edge_diagnostic'
   salinity = Quantity(domain, name='salinity', register=True)
   salinity.set_values(0.0)
