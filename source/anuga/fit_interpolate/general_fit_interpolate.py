"""Least squares interpolation.

   Implements a least-squares interpolation.

   Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
   Geoscience Australia, 2004.

DESIGN ISSUES
* what variables should be global?
- if there are no global vars functions can be moved around alot easier

* The public interface
__init__
interpolate
interpolate_block

"""

import time
import os
from warnings import warn

from Numeric import zeros, array, Float, Int, dot, transpose, concatenate, \
     ArrayType, allclose, take, NewAxis, arange

from anuga.caching.caching import cache
from anuga.pyvolution.neighbour_mesh import Mesh
from anuga.utilities.sparse import Sparse, Sparse_CSR
from anuga.utilities.cg_solve import conjugate_gradient, VectorShapeError
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.pyvolution.quad import build_quadtree
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.polygon import in_and_outside_polygon
from anuga.geospatial_data.geospatial_data import Geospatial_data, \
     ensure_absolute
from search_functions import search_tree_of_vertices



class FitInterpolate:
    
    def __init__(self,
                 vertex_coordinates,
                 triangles,
                 mesh_origin=None,
                 verbose=False,
                 max_vertices_per_cell=30):


        """ Build interpolation matrix mapping from
        function values at vertices to function values at data points

        Inputs:

          vertex_coordinates: List of coordinate pairs [xi, eta] of
	      points constituting a mesh (or an m x 2 Numeric array or
              a geospatial object)
              Points may appear multiple times
              (e.g. if vertices have discontinuities)

          triangles: List of 3-tuples (or a Numeric array) of
              integers representing indices of all vertices in the mesh.

          mesh_origin: A geo_reference object or 3-tuples consisting of
              UTM zone, easting and northing.
              If specified vertex coordinates are assumed to be
              relative to their respective origins.

          max_vertices_per_cell: Number of vertices in a quad tree cell
          at which the cell is split into 4.

          Note: Don't supply a vertex coords as a geospatial object and
              a mesh origin, since geospatial has its own mesh origin.
        """
        
        #Convert input to Numeric arrays
        triangles = ensure_numeric(triangles, Int)
        vertex_coordinates = ensure_absolute(vertex_coordinates,
                                         geo_reference = mesh_origin)

        #Don't pass geo_reference to mesh.  It doesn't work.
        self.mesh = Mesh(vertex_coordinates, triangles)
        self.mesh.check_integrity()
        self.root = build_quadtree(self.mesh,
                              max_points_per_cell = max_vertices_per_cell)
        #print "self.root",self.root.show() 
        
        
