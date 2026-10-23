
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
levels and holes in the domain.

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


Resolution from a size field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each polygon in ``interior_regions`` becomes part of the mesh: its outline is
a ring of segments that the triangulation honours exactly. That is what you
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
``(polygon, maximum area)`` list that ``interior_regions`` takes:

.. code-block:: python

   town = [[8.0, 3.0], [12.0, 3.0], [12.0, 7.0], [8.0, 7.0]]

   field = anuga.size_field_from_regions(bounding_polygon,
                               [(town, 0.01)],
                               maximum_triangle_area = 0.2)

   domain = anuga.create_domain_from_regions(bounding_polygon,
                               boundary_tags,
                               maximum_triangle_area = 0.2,
                               size_field = field,
                               relax_passes = 30)

The polygons given to :func:`size_field_from_regions` may touch, overlap
(the finest wins), have sharp corners or extend beyond the bounding polygon.
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
passes, to enlarge them. Vertices on the boundary, on breaklines and on
region outlines do not move, and no triangles are added or removed. On
meshes of 0.5 to 4 million triangles built from a size field it has roughly
doubled the smallest ``radii``, and with it the timestep.


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
   (``rectangular_cross_domain``, ``create_domain_from_regions``) and the
   ``Domain`` methods.