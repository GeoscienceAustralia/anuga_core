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

from anuga.caching.caching import cache
from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh

from anuga.utilities.sparse import Sparse, Sparse_CSR
from anuga.utilities.cg_solve import conjugate_gradient, VectorShapeError
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.polygon import in_and_outside_polygon
from anuga.utilities.quad import build_quadtree

from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geospatial_data.geospatial_data import Geospatial_data, \
     ensure_absolute
from anuga.fit_interpolate.search_functions import set_last_triangle

import numpy as num


# tests fail if 2 is used
MAX_VERTICES_PER_CELL = 13 # A value of 8 or lower can cause problems,
                           # if a vert has 9 triangles.
                           
build_quadtree_time = 0

def get_build_quadtree_time():
    return build_quadtree_time

class FitInterpolate:
    
    def __init__(self,
                 vertex_coordinates=None,
                 triangles=None,
                 mesh=None,
                 mesh_origin=None,
                 verbose=False,
                 max_vertices_per_cell=None):


        """ Build interpolation matrix mapping from
        function values at vertices to function values at data points

        Pass in a mesh instance or vertex_coordinates and triangles
        and optionally mesh_origin
        
        Inputs:

          vertex_coordinates: List of coordinate pairs [xi, eta] of
	      points constituting a mesh (or an m x 2 numeric array or
              a geospatial object)
              Points may appear multiple times
              (e.g. if vertices have discontinuities)

          triangles: List of 3-tuples (or a numeric array) of
              integers representing indices of all vertices in the mesh.

        mesh: A mesh instance describing the mesh.

          mesh_origin: A geo_reference object or 3-tuples consisting of
              UTM zone, easting and northing.
              If specified vertex coordinates are assumed to be
              relative to their respective origins.

          max_vertices_per_cell: Number of vertices in a quad tree cell
          at which the cell is split into 4.

          Note: Don't supply a vertex coords as a geospatial object and
              a mesh origin, since geospatial has its own mesh origin.
        """
        global build_quadtree_time
        if max_vertices_per_cell == None:
            max_vertices_per_cell = MAX_VERTICES_PER_CELL
        if mesh is None:
            if vertex_coordinates is not None and  triangles is not None:
                # Fixme (DSG) Throw errors if triangles or vertex_coordinates
                # are None
            
                #Convert input to numeric arrays
                triangles = ensure_numeric(triangles, num.int)
                vertex_coordinates = ensure_absolute(vertex_coordinates,
                                                 geo_reference = mesh_origin)

                if verbose: print 'FitInterpolate: Building mesh'        
                self.mesh = Mesh(vertex_coordinates, triangles)
                #self.mesh.check_integrity() # Time consuming
            else:
                self.mesh = None
        else:
            self.mesh = mesh

        if self.mesh is not None:
            if verbose: print 'FitInterpolate: Building quad tree'
            #This stores indices of vertices
            t0 = time.time()
            #print "self.mesh.get_extent(absolute=True)", \
            #self.mesh.get_extent(absolute=True)
            self.root = build_quadtree(self.mesh,
                                       max_points_per_cell = max_vertices_per_cell)
            #print "self.root",self.root.show()
        
            build_quadtree_time =  time.time()-t0
            set_last_triangle()
        
    def __repr__(self):
        return 'Interpolation object based on: ' + repr(self.mesh)
