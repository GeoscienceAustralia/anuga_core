
.. currentmodule:: anuga

Creating a Domain
=================

The first step in running an ANUGA model is to create a domain. This is done by
creating a mesh and then creating a domain from that mesh.

`rectangular_cross_domain`
--------------------------

The domain (mesh) can be created in a number of ways. The simplest way is to 
create a simple rectangular
domain using the :func:`rectangular_cross_domain` function.

For instance the following code creates a 1m  by 1m rectangular domain, with 
a 10 by 10 mesh, with the bottom left corner at (0,0).

.. code-block:: python

    domain = anuga.rectangular_cross_domain(10, 10)

See the :func:`rectangular_cross_domain` API for the full argument list, and the
:doc:`../examples/script_simple_example` or
:doc:`../examples/notebook_simple_example` for worked examples.

`create_domain_from_regions`
----------------------------

The usual method for creating a domain for practical problems is to create a
domain by defining a boundary polygon and then a set of regions within the domain to 
define area of different refinement levels and holes in the domain.
This is done using the :func:`create_domain_from_regions` 
function. 

The regions are defined by a list of polygons. 
Each polygon is defined by a list of points. The  most important polygon is the 
boundary polygon. This is the outer polygon that defines the boundary of the domain.
The segments of the boundary need to be tagged with boundary tags which will allow 
different boundary conditions to be applied to different segments.



Other polygons are the interior polygons that define the regions within the domain. 
These other polygons can be used to define regions with different refinement
levels (see `Refinement regions`_) and holes in the domain.

The following example creates a domain with a rectangular boundary 20m by 10m with 
boundary tags on the 4 sides of the rectangle, and the mesh having a maximum triangle 
area of 0.2 m^2. 

.. code-block:: python

   import anuga

   bounding_polygon = [[0.0, 0.0],
                    [20.0, 0.0],
                    [20.0, 10.0],
                    [0.0, 10.0]]

   boundary_tags={'bottom': [0],
                'right': [1],
                'top': [2],
                'left': [3]}

   domain = anuga.create_domain_from_regions(bounding_polygon,
                               boundary_tags, 
                               maximum_triangle_area = 0.2,
                               )



See the :func:`create_domain_from_regions` API for the full argument list, and the
:doc:`../examples/notebook_create_domain_from_regions` for a worked example.


Refinement regions
~~~~~~~~~~~~~~~~~~

A single ``maximum_triangle_area`` gives much the same resolution everywhere,
but a model usually needs fine triangles only in places: a river channel, a
town, the ground around a structure. ``interior_regions`` is a list of
``(polygon, maximum area)`` pairs. Each sets a smaller maximum triangle area
inside its polygon, and ``maximum_triangle_area`` still applies everywhere else:

.. code-block:: python

   town = [[8.0, 3.0], [12.0, 3.0], [12.0, 7.0], [8.0, 7.0]]

   domain = anuga.create_domain_from_regions(bounding_polygon,
                               boundary_tags,
                               maximum_triangle_area = 0.2,
                               interior_regions = [(town, 0.01)])

Triangles are now at most 0.01 m² in the town and 0.2 m² elsewhere, giving about
4,300 triangles against about 1,500 for the mesh above.

An equilateral triangle of edge length *h* has area 0.43 *h*\ ², so 0.01 m² means
edges of about 0.15 m. Halving the edge length takes four times as many
triangles. The smallest triangles also set the timestep of the whole run, so
refine only where the flow needs it.

The polygons have to follow a few rules:

* **Inside the boundary.** A polygon that is not fully inside the bounding
  polygon raises an error. With ``fail_if_polygons_outside = False`` it is
  ignored instead.
* **Nested or separate.** A polygon can sit inside another, in any order, and
  the inner area applies inside it, so ``[(town, 0.05), (street, 0.01)]`` works
  when ``street`` is inside ``town``. The polygons may share vertices and
  edges, and an inner region or hole may also touch the bounding polygon.
* **No crossing outlines.** Crossing outlines raise ``PolygonError``, because
  the mesh would be wrong: triangles larger than any area you asked for, and a
  finer region that does not get its resolution. The same goes for
  ``interior_holes``, and for a hole that crosses a region.
* **Not too close.** Outlines close to each other, or to the boundary, force
  small triangles between them.


Resolution from a size field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The outline of a refinement region becomes part of the mesh, as a ring of
segments that the triangulation honours exactly. That is what you
want when the outline is a physical feature. Often it is not. A region such
as "below 12 m elevation" or "within 500 m of the town" only says where the
mesh should be fine, and its outline can do harm: a sharp corner, or two
outlines passing close to each other, forces very small triangles, and since
ANUGA steps at the smallest ``radii / wave speed`` in the mesh, a single such
triangle sets the timestep of the whole run.

A *size field* gives the resolution without adding any segments. It is the
maximum triangle area as a function of position, and a triangle is refined
whenever its area exceeds the field at its centroid.
:func:`size_field_from_regions` builds one from the same
``(polygon, maximum area)`` list that ``interior_regions`` takes, so the town
above becomes:

.. code-block:: python

   field = anuga.size_field_from_regions(bounding_polygon,
                               [(town, 0.01)],
                               maximum_triangle_area = 0.2)

   domain = anuga.create_domain_from_regions(bounding_polygon,
                               boundary_tags,
                               maximum_triangle_area = 0.2,
                               size_field = field,
                               relax_passes = 30)

Unlike refinement regions, the polygons given to
:func:`size_field_from_regions` may touch, overlap (the finest wins), have sharp
corners or extend beyond the bounding polygon.
The field is graded: the target edge length holds inside each polygon and
grows away from it by ``grade`` metres per metre (0.5 by default), so the
mesh coarsens smoothly instead of stepping at the outline.

``size_field`` can also be a :class:`Size_field` built from your own raster
of target areas, or any function ``f(x, y)`` returning the maximum area at
absolute coordinates ``x, y``. A function is called a few times for every
triangle in the mesh, so for a large mesh rasterise it into a
:class:`Size_field`. Because the field is only looked at where triangle
centroids fall, it should vary smoothly.

A size field works together with the other arguments. Use ``breaklines`` and
``interior_holes`` for features that must be in the mesh, and
``interior_regions`` where the outline is physical; a triangle is refined
if either its region's maximum area or the field calls for it. The field
needs the ``meshpy`` triangulation backend.

``relax_passes`` is independent of the size field and can be used with any
mesh. Triangulation bounds the largest triangle and the smallest angle, but
not the smallest triangle, and leaves a few cells much smaller than the rest.
Relaxation moves the vertices of those cells, for up to the given number of
passes, to enlarge them. The vertices you gave (polygon corners, breakline
and hole vertices) do not move, so every outline keeps its exact shape; a
point the triangulation inserted part-way along a straight segment may slide
along it, which helps most beside short building sides. No triangles are
added or removed. On meshes of 0.5 to 4 million triangles built from a size
field it has roughly doubled the smallest ``radii``, and with it the
timestep.


Meshes from a DEM
~~~~~~~~~~~~~~~~~

Where water goes is written in the elevation model, so the DEM itself can
say where the mesh should be fine. :func:`create_domain_from_dem` routes the
DEM by D8 (each cell drains to its steepest downhill neighbour, after the
depressions are filled), which gives every cell the area draining through
it: the streams of a catchment, or the streets and overland flow paths of a
town. It then meshes fine along those paths and coarse on the slopes between
them, with a graded size field (so no segments are added), relaxes the worst
cells, and sets each triangle's elevation to the mean of the DEM over it:

.. code-block:: python

   domain = anuga.create_domain_from_dem('topography.asc',
                               maximum_triangle_area = 8.0,
                               bounding_polygon = bounding_polygon,
                               boundary_tags = boundary_tags,
                               streams = [(5000.0, 2.0), (1000.0, 4.0)],
                               interior_holes = houses)

``streams`` is a list of ``(drainage area, maximum triangle area)``: here
cells that 5000 m² or more drains through get triangles of up to 2 m², and
those draining 1000 m² up to 4 m². The target holds for ``band`` (1.5) edge
lengths either side of the flow path and then grows by ``grade`` (0.5) metres
per metre, as in :func:`size_field_from_regions`. Water is routed around
``interior_holes`` and leaves across the bounding polygon.

``floodplain`` refines the valley floors, the ground a stream spreads over
when it rises: ``(height, drainage area, maximum triangle area)`` refines the
cells no more than ``height`` above the stream they drain to (HAND, height
above nearest drainage), where that stream drains at least ``drainage
area``. ``regions`` adds polygons for what the DEM cannot show, such as a
town centre or a structure, again without adding their segments.

The DEM can be an ``.asc``, ``.grd``, ``.dem`` or ``.tif`` file, or a raster
tuple ``(x, y, Z)`` with ``Z[i, j]`` the elevation at ``(x[i], y[j])``. Use
:func:`size_field_from_dem` to build the field alone (to combine it, or to
plot it, with its ``drainage_area``), :func:`create_pmesh_from_dem` for the
mesh without a domain, and :func:`drainage_area_from_dem` to look at the
routing.

On the Merewether case (:doc:`../examples/notebook_create_domain_from_dem`)
this meshes the town with a quarter of the triangles of a uniform 1 m² mesh
and twice the timestep, and the flood levels at the observation points agree
to 2 cm.

**Pits.** On a coarse mesh a triangle's mean elevation can fill the gully it
drained by, leaving the triangle upstream lower than all its neighbours: a
pit that holds water for good. In a large catchment these can store a large
part of the rain. :func:`remove_mesh_pits` removes them on the mesh, from the
outlets you name: ``method='fill'`` raises each pit to its spill level, and
``method='breach'`` instead lowers the cells along its lowest way out, filling
only where the cut would be deeper than ``max_depth``. A small cell is more
likely a real depression than an artefact, so pass ``where`` to limit the
change to the coarse cells. ``create_domain_from_dem`` does this when given
``pits='fill'`` or ``'breach'`` with ``pit_outlets``; by default it leaves
the elevation as sampled.


.. seealso::

   :doc:`flow_algorithms`
      How to choose between DE0, DE1, DE_ader2, and DE2 — trade-offs between
      cost, accuracy, and robustness.

   :doc:`coordinate_reference`
      How to attach a coordinate reference system (UTM zone, national grid,
      or arbitrary local CRS) to a domain via :class:`Geo_reference`.

   `ANUGA User Manual — Chapter 7: The Domain <https://github.com/anuga-community/anuga_user_manual>`_
      Covers domain construction in depth, including mesh generation from
      polygon regions, geo-referencing, flow algorithm choices, and domain
      attributes.

.. seealso::

   :ref:`api_domain` — full API of the domain-creation functions
   (``rectangular_cross_domain``, ``create_domain_from_regions``,
   ``create_domain_from_dem``) and the
   ``Domain`` methods.