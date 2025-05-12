
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

Check here for full documentation of the :func:`rectangular_cross_domain` function 
and check the :doc:`../../examples/script_simple_example` 
or :doc:`../../examples/notebook_simple_example` for examples of how to use it.

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
These other polygons can be used to define regions with different refinenment 
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



Check here for full documentation of the :func:`create_domain_from_regions` function 
and check the :doc:`../../examples/notebook_create_domain_from_regions` for an example of how to use it.


Reference
---------

.. autosummary::
   :toctree:  
   
   rectangular_cross_domain
   create_domain_from_regions