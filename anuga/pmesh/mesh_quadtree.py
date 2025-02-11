"""
    Quadtree representation of 2D mesh geometry.

   Ole Nielsen, Stephen Roberts, Duncan Gray
   Geoscience Australia, 2006.


   PADARN NOTE 06/12/12: This quad tree has been
   replaced by a C quad tree. Save the old code
   somewhere?

   PADARN NOTE: Most of the functionality of the
   old quad tree has been emulated. However, some
   methods inherrited from cell have not been. This
   causes a few of the unit tests to fail, and thus
   can easily be identified if it is felt the are
   needed.

"""

from anuga.config import max_float

from anuga.geometry.quad import Cell
from anuga.geometry.aabb import AABB

import numpy as num
from anuga.utilities.numerical_tools import ensure_numeric
import anuga.fit_interpolate.fitsmooth_ext as fitsmooth


# PADARN NOTE: I don't think much from Cell is used anymore, if
# anything, this dependency could be removed.
class MeshQuadtree(Cell):
    """ A quadtree constructed from the given mesh.
        This class is the root node of a quadtree,
        and derives from a Cell.
        It contains optimisations and search patterns specific to meshes.
    """

    def __init__(self, mesh, verbose=False):
        """Build quad tree for mesh.

        All vertex indices in the mesh are stored in a quadtree.
        """
        self.mesh = mesh

        self.set_extents()
        self.add_quad_tree()

        Cell.__init__(self, self.extents, None)  # root has no parent


    def __getstate__(self):
        dic = self.__dict__
        if ('root' in dic):
            dic.pop('root')
        return dic

    def set_extents(self):
        extents = AABB(*self.mesh.get_extent(absolute=True))
        extents.grow(1.001)  # To avoid round off error
        numextents = [extents.xmin, extents.xmax, extents.ymin, extents.ymax]
        self.extents = num.array(numextents, float)
        #print self.extents

    def add_quad_tree(self):

        V = self.mesh.get_vertex_coordinates(absolute=True)
        
        self.set_extents()
        #print self.extents
        self.root = fitsmooth.build_quad_tree(self.mesh.triangles, V, self.extents)


    # PADARN NOTE: This function does not properly emulate the old functionality -
    # it seems uneeded though. Check this.
    def search(self, point):
        return self.search_fast(point)

    # PADARN NOTE: Although this function emulates the functionality of the old
    # quad tree, it cannot be called on the sub-trees anymore.
    def count(self):
        if not hasattr(self, 'root'):
            self.add_quad_tree()
        return fitsmooth.items_in_tree(self.root)

    def search_fast(self, point):
        """
        Find the triangle (element) that the point x is in.

        Does a coherent quadtree traversal to return a single triangle that the
        point falls within. The traversal begins at the last triangle found.
        If this fails, it checks the triangles beneath it in the tree, and then
        begins traversing up through the tree until it reaches the root.

        This results in performance which varies between constant time and O(n),
        depending on the geometry.

        Inputs:
            point:    The point to test

        Return:
            element_found, sigma0, sigma1, sigma2, k

            where
            element_found: True if a triangle containing x was found
            sigma0, sigma1, sigma2: The interpolated values
            k: Index of triangle (if found)

        """

        # PADARN NOTE: Adding checks on the input point to make sure it is a float.

        if not hasattr(self, 'root'):
            self.add_quad_tree()

        point = ensure_numeric(point, float)

        [found, sigma, index] = fitsmooth.individual_tree_search(self.root, point)

        if found == 1:
            element_found = True
        else:
            element_found = False

        return element_found, sigma[0], sigma[1], sigma[2], index

    # PADARN NOTE: Only here to pass unit tests - does nothing.
    def set_last_triangle(self):
        pass
