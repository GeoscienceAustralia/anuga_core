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

from builtins import object
import time
import os
from warnings import warn

from anuga.caching.caching import cache
from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh

from anuga.utilities.sparse import Sparse, Sparse_CSR
from anuga.utilities.cg_solve import conjugate_gradient, VectorShapeError
from anuga.utilities.numerical_tools import ensure_numeric

from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geospatial_data.geospatial_data import Geospatial_data, \
     ensure_absolute
from anuga.pmesh.mesh_quadtree import MeshQuadtree
import anuga.utilities.log as log

import numpy as num

                           
build_quadtree_time = 0

def get_build_quadtree_time():
    return build_quadtree_time

class FitInterpolate(object):
    
    def __init__(self,
                 vertex_coordinates=None,
                 triangles=None,
                 mesh=None,
                 mesh_origin=None,
                 verbose=False):


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

          Note: Don't supply a vertex coords as a geospatial object and
              a mesh origin, since geospatial has its own mesh origin.
        """

        # NOTE PADARN: The Fit_Interpolate class now uses a the c based
        # quad tree to store triangles, rather than the python based tree.
        # The tree is still stored at self.root. However, the subtrees of
        # the new quad tree can not be directly accessed by python as 
        # was previously possible.
        # Most of the previous functionality has been preserved.

        global build_quadtree_time
        if mesh is None:
            if vertex_coordinates is not None and triangles is not None:
                # Fixme (DSG) Throw errors if triangles or vertex_coordinates
                # are None

                # Convert input to numeric arrays
                triangles = ensure_numeric(triangles, int)
                vertex_coordinates = ensure_absolute(vertex_coordinates,
                                                 geo_reference=mesh_origin)

                if verbose:
                    log.critical('FitInterpolate: Building mesh')
					

                self.mesh = Mesh(vertex_coordinates, triangles)

                #self.mesh.check_integrity() # Time consuming
            else:
                self.mesh = None
        else:
            self.mesh = mesh

        if self.mesh is not None:
            if verbose:
                log.critical('FitInterpolate: Building quad tree')
            #This stores indices of vertices
            t0 = time.time()

            self.root = MeshQuadtree(self.mesh, verbose=verbose)
            build_quadtree_time = time.time() - t0
            # Padarn Note 06/12/12: Do I need this?
            #self.root.set_last_triangle()

    def build_quad_tree(self,verbose=False):
        self.root = MeshQuadtree(self.mesh, verbose=verbose)

    def __repr__(self):
        return 'Interpolation object based on: ' + repr(self.mesh)
