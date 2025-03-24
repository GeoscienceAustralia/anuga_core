
.. currentmodule:: anuga

Creating a Domain
=================

The first step in running an ANUGA model is to create a domain. This is done by
creating a mesh and then creating a domain from that mesh.

`rectangular_cross_domain` creates a simple rectangular domain.
----------------------------------------------------------------

The domain (mesh) can be created in a number of ways. The simplest way is to 
create a simple rectangular
domain using the :func:`rectangular_cross_domain` function.

.. code-block:: python

    from anuga import rectangular_cross_domain

    domain = rectangular_cross_domain(10, 10)


`create_domain_from_regions` creates a simple rectangular domain.
-----------------------------------------------------------------

The usual method for creating a domain for practical problems is to create a
domain from a set of regions. This is done using the :func:`create_domain_from_regions` 
function. 

The regions are defined by a list of polygons. 
Each polygon is defined by a list of points. The  most important polygon is the 
boundary polygon. This is the outer polygon that defines the boundary of the domain.
The segments of the boundary need to be tagged with boundary tags which will allow 
different boundary conditions to be applied to different segments.



Other polygons are the interior polygons that define the regions within the domain. 
These other polygons can be used to define regions with different refinenment 
levels and holes in the domain. 

The following example creates a domain with a rectangular boundary and a square hole in the middle. 

.. code-block:: python

    from anuga import create_domain_from_regions

    boundary_polygon = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    boundary_tags = {'left': [0], 'right': [1], 'top': [2], 'bottom': [3]}
    hole_polygon = [(4.0, 4.0), (6.0, 4.0), (6.0, 6.0), (4.0, 6.0)]
    hole_pt = (5.0, 5.0)

    domain = create_domain_from_regions(boundary_polygon = boundary_polygon, 
                                       boundary_tags = boundadry_tags,
                                        region_polygons  =[hole_polygon]
                                        hole_pts = [hole_pt])




.. autosummary::
   :toctree:  
   
   rectangular_cross_domain
   create_domain_from_regions