"""
General functions used in fit and interpolate.

   Ole Nielsen, Stephen Roberts, Duncan Gray
   Geoscience Australia, 2006.

"""

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.config import max_float

from anuga.geometry.quad import Cell
from anuga.geometry.aabb import AABB

from anuga.utilities import compile as compile_c
if compile_c.can_use_C_extension('polygon_ext.c'):
    # Underlying C implementations can be accessed
    from polygon_ext import _is_inside_triangle        
else:
    MESSAGE = 'C implementations could not be accessed by %s.\n ' % __file__
    MESSAGE += 'Make sure compile_all.py has been run as described in '
    MESSAGE += 'the ANUGA installation guide.'
    raise Exception(MESSAGE)

import numpy as num
 

LAST_TRIANGLE = [[[-1, num.array([[max_float, max_float],
                               [max_float, max_float],
                               [max_float, max_float]]),
                      num.array([[max_float, max_float],
                               [max_float, max_float],
                               [max_float, max_float]])], -10]]



class MeshQuadtree(Cell):
    """ A quadtree constructed from the given mesh.
        This class is the root node of a quadtree,
        and derives from a Cell.
        It contains optimisations and search patterns specific to meshes.
    """
    def __init__(self, mesh):
        """Build quad tree for mesh.

        All vertex indices in the mesh are stored in a quadtree.
        """
        
        extents = AABB(*mesh.get_extent(absolute=True))   
        extents.grow(1.001) # To avoid round off error
        Cell.__init__(self, extents, None)  # root has no parent

        self.last_triangle = None        
        N = len(mesh)
        self.mesh = mesh
        self.set_last_triangle()  

        # Get x,y coordinates for all vertices for all triangles
        V = mesh.get_vertex_coordinates(absolute=True)
        
        normals = mesh.get_normals()
        
        # Check each triangle
        for i in range(N):
            i3 = 3*i
            x0, y0 = V[i3, :]
            x1, y1 = V[i3+1, :]
            x2, y2 = V[i3+2, :]

            node_data = [i, V[i3:i3+3, :], normals[i, :]]

            # insert a tuple with an AABB, and the triangle index as data
            self.insert_item((AABB(min([x0, x1, x2]), max([x0, x1, x2]), \
                             min([y0, y1, y2]), max([y0, y1, y2])), \
                             node_data))

    def search_fast(self, point):
        """
        Find the triangle (element) that the point x is in.

        Inputs:
            point:    The point to test
        
        Return:
            element_found, sigma0, sigma1, sigma2, k

            where
            element_found: True if a triangle containing x was found
            sigma0, sigma1, sigma2: The interpolated values
            k: Index of triangle (if found)

        """
        
        point = ensure_numeric(point, num.float)
                 
        # check the last triangle found first
        element_found, sigma0, sigma1, sigma2, k = \
                   self._search_triangles_of_vertices(self.last_triangle, point)
        if element_found:
            return True, sigma0, sigma1, sigma2, k

        branch = self.last_triangle[0][1]

        # test neighbouring tris
        tri_data = branch.test_leaves(point)          
        element_found, sigma0, sigma1, sigma2, k = \
                    self._search_triangles_of_vertices(tri_data, point)
        if element_found:
            return True, sigma0, sigma1, sigma2, k       

        # search rest of tree
        element_found = False
        next_search = [branch]
        while branch:               
            for sibling in next_search:
                tri_data = sibling.search(point)         
                element_found, sigma0, sigma1, sigma2, k = \
                            self._search_triangles_of_vertices(tri_data, point)
                if element_found:
                    return True, sigma0, sigma1, sigma2, k
            
            next_search = branch.get_siblings()                            
            branch = branch.parent
            if branch:
                tri_data = branch.test_leaves(point)     
                element_found, sigma0, sigma1, sigma2, k = \
                            self._search_triangles_of_vertices(tri_data, point)
                if element_found:
                    return True, sigma0, sigma1, sigma2, k      

        return element_found, sigma0, sigma1, sigma2, k


    def _search_triangles_of_vertices(self, triangles, point):
        """Search for triangle containing x among triangle list

        This is called by search_tree_of_vertices once the appropriate node
        has been found from the quad tree.
        
        Input check disabled to speed things up. 
        
        point is the point to test 
        triangles is the triangle list
        return the found triangle and its interpolation sigma.
        """  

        for node_data in triangles:             
            if bool(_is_inside_triangle(point, node_data[0][1], \
                        int(True), 1.0e-12, 1.0e-12)):
                normals = node_data[0][2]      
                n0 = normals[0:2]
                n1 = normals[2:4]
                n2 = normals[4:6]          
                xi0, xi1, xi2 = node_data[0][1]

                sigma0 = num.dot((point-xi1), n0)/num.dot((xi0-xi1), n0)
                sigma1 = num.dot((point-xi2), n1)/num.dot((xi1-xi2), n1)
                sigma2 = num.dot((point-xi0), n2)/num.dot((xi2-xi0), n2)

                # Don't look for any other triangles in the triangle list
                self.last_triangle = [node_data]
                return True, sigma0, sigma1, sigma2, node_data[0][0] # tri index
        return False, -1, -1, -1, -10

        
        
    def set_last_triangle(self):
        """ Reset last triangle.
            The algorithm is optimised to find nearby triangles to the
            previously found one. This is called to reset the search to
            the root of the tree.
        """
        self.last_triangle = LAST_TRIANGLE
        self.last_triangle[0][1] = self # point at root by default          

    

                
