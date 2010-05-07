"""
General functions used in fit and interpolate.

   Ole Nielsen, Stephen Roberts, Duncan Gray
   Geoscience Australia, 2006.

"""
import time

from anuga.utilities.numerical_tools import get_machine_precision
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.config import max_float
from anuga.geometry.quad import Cell
from anuga.geometry.aabb import AABB

from anuga.utilities import compile
if compile.can_use_C_extension('polygon_ext.c'):
    # Underlying C implementations can be accessed
    from polygon_ext import _is_inside_triangle        
else:
    msg = 'C implementations could not be accessed by %s.\n ' %__file__
    msg += 'Make sure compile_all.py has been run as described in '
    msg += 'the ANUGA installation guide.'
    raise Exception, msg

import numpy as num


initial_search_value = 'uncomment search_functions code first'#0
search_one_cell_time = initial_search_value
search_more_cells_time = initial_search_value

# FIXME(Ole): Could we come up with a less confusing structure?
# FIXME(James): remove this global var
LAST_TRIANGLE = [[-10,
                   (num.array([[max_float, max_float],
                               [max_float, max_float],
                               [max_float, max_float]]),
                    (num.array([1.,1.]),      
                     num.array([0.,0.]),      
                     num.array([-1.1,-1.1])))]]

last_triangle = LAST_TRIANGLE					 
					 
def search_tree_of_vertices(root, mesh, x):
    """
    Find the triangle (element) that the point x is in.

    Inputs:
        root: A quad tree of the vertices
        mesh: The mesh which the quad tree indexes into
        x:    The point being placed
    
    Return:
        element_found, sigma0, sigma1, sigma2, k

        where
        element_found: True if a triangle containing x was found
        sigma0, sigma1, sigma2: The interpolated values
        k: Index of triangle (if found)

    """
    global search_one_cell_time
    global search_more_cells_time

    # Search the last triangle first
    element_found, sigma0, sigma1, sigma2, k = \
        _search_triangles_of_vertices(last_triangle, x)
                   
    if element_found is True:
        return element_found, sigma0, sigma1, sigma2, k

    
    # Get triangles in the cell that the point is in.
    tri_indices = root.search(x[0], x[1])
    triangles = _trilist_from_indices(mesh, tri_indices)
    
    element_found, sigma0, sigma1, sigma2, k = \
                   _search_triangles_of_vertices(triangles, x)

    is_more_elements = True
    
    # while not element_found and is_more_elements:
        # triangles, branch = root.expand_search()
        # if branch == []:
            # # Searching all the verts from the root cell that haven't
            # # been searched.  This is the last try
            # element_found, sigma0, sigma1, sigma2, k = \
                           # _search_triangles_of_vertices(triangles, x)
            # is_more_elements = False
        # else:
            # element_found, sigma0, sigma1, sigma2, k = \
                       # _search_triangles_of_vertices(triangles, x)
                       
        
    return element_found, sigma0, sigma1, sigma2, k


def _search_triangles_of_vertices(triangles, x):
    """Search for triangle containing x amongs candidate_vertices in triangles

    This is called by search_tree_of_vertices once the appropriate node
    has been found from the quad tree.
    

    This function is responsible for most of the compute time in
    fit and interpolate.
    """
    global last_triangle

    x = ensure_numeric(x, num.float)     
    
    # These statments are needed if triangles is empty
    sigma2 = -10.0
    sigma0 = -10.0
    sigma1 = -10.0
    k = -10
    # For all vertices in same cell as point x
    element_found = False    
    for k, tri_verts_norms in triangles:
        tri = tri_verts_norms[0]
        tri = ensure_numeric(tri)        
        # k is the triangle index
        # tri is a list of verts (x, y), representing a triangle
        # Find triangle that contains x (if any) and interpolate
        
        # Input check disabled to speed things up.    
        if bool(_is_inside_triangle(x, tri, int(True), 1.0e-12, 1.0e-12)):
            
            n0, n1, n2 = tri_verts_norms[1]        
            sigma0, sigma1, sigma2 =\
                compute_interpolation_values(tri, n0, n1, n2, x)
                
            element_found = True    

            # Don't look for any other triangles in the triangle list
            last_triangle = [[k, tri_verts_norms]]
            break            
            
    return element_found, sigma0, sigma1, sigma2, k


def _trilist_from_indices(mesh, indices):
    """return a list of lists. For the inner lists,
    The first element is the triangle index,
    the second element is a list.for this list
       the first element is a list of three (x, y) vertices,
       the following elements are the three triangle normals.

    """

    ret_list = []
    for i in indices:
        vertices = mesh.get_vertex_coordinates(triangle_id=i, absolute=True)
        n0 = mesh.get_normal(i, 0)
        n1 = mesh.get_normal(i, 1)
        n2 = mesh.get_normal(i, 2) 
        ret_list.append([i, [vertices, (n0, n1, n2)]])
    return ret_list
                
            
def compute_interpolation_values(triangle, n0, n1, n2, x):
    """Compute linear interpolation of point x and triangle.
    
    n0, n1, n2 are normal to the tree edges.
    """

    # Get the three vertex_points of candidate triangle k
    xi0, xi1, xi2 = triangle

    sigma0 = num.dot((x-xi1), n0)/num.dot((xi0-xi1), n0)
    sigma1 = num.dot((x-xi2), n1)/num.dot((xi1-xi2), n1)
    sigma2 = num.dot((x-xi0), n2)/num.dot((xi2-xi0), n2)

    return sigma0, sigma1, sigma2

def set_last_triangle():
    global last_triangle
    last_triangle = LAST_TRIANGLE
    
def search_times():

    global search_one_cell_time
    global search_more_cells_time

    return search_one_cell_time, search_more_cells_time

def reset_search_times():

    global search_one_cell_time
    global search_more_cells_time
    search_one_cell_time = initial_search_value
    search_more_cells_time = initial_search_value


class MeshQuadtree(Cell):
    """ A quadtree constructed from the given mesh.
    """
    def __init__(self, mesh):
        """Build quad tree for mesh.

        All vertex indices in the mesh are stored in a quadtree.
        """
        
        extents = AABB(*mesh.get_extent(absolute=True))   
        extents.grow(1.001) # To avoid round off error
        Cell.__init__(self, extents)
        
        N = len(mesh)

        # Get x,y coordinates for all vertices for all triangles
        V = mesh.get_vertex_coordinates(absolute=True)
        
        # Check each triangle
        for i in range(N):
            x0, y0 = V[3*i, :]
            x1, y1 = V[3*i+1, :]
            x2, y2 = V[3*i+2, :]

            # insert a tuple with an AABB, and the triangle index as data
            self._insert((AABB(min([x0, x1, x2]), max([x0, x1, x2]), \
                             min([y0, y1, y2]), max([y0, y1, y2])), \
                             i))

