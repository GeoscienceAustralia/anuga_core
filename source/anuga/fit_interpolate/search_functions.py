"""
General functions used in fit and interpolate.

   Ole Nielsen, Stephen Roberts, Duncan Gray
   Geoscience Australia, 2006.

"""
import time

from anuga.utilities.polygon import is_inside_triangle
from anuga.utilities.numerical_tools import get_machine_precision
from anuga.config import max_float

import Numeric as num


initial_search_value = 'uncomment search_functions code first'#0
search_one_cell_time = initial_search_value
search_more_cells_time = initial_search_value

LAST_TRIANGLE = [[-10,
                   (num.array([[max_float,max_float],
                               [max_float,max_float],
                               [max_float,max_float]]),
                    (num.array([1,1], num.Int),      #array default#
                     num.array([0,0], num.Int),      #array default#
                     num.array([-1.1,-1.1])))]]

def search_tree_of_vertices(root, mesh, x):
    """
    Find the triangle (element) that the point x is in.

    Inputs:
        root: A quad tree of the vertices
        mesh: The underlying mesh
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

    # This will be returned if element_found = False
    element_found = False    
    sigma2 = -10.0
    sigma0 = -10.0
    sigma1 = -10.0
    k = -10.0

    # Search the last triangle first
    try:
        element_found, sigma0, sigma1, sigma2, k = \
            _search_triangles_of_vertices(mesh, last_triangle, x)
            
    except:
        print 'This should never happen:', last_triangle
        element_found = False
                   
    if element_found is True:
        return element_found, sigma0, sigma1, sigma2, k

    
    # Get triangles in the cell that the point is in.
    # Triangle is a list, first element triangle_id,
    # second element the triangle
    triangles = root.search(x[0], x[1])
    is_more_elements = True
    element_found, sigma0, sigma1, sigma2, k = \
                   _search_triangles_of_vertices(mesh,
                                                 triangles, x)
    #search_one_cell_time += time.time()-t0
    #print "search_one_cell_time",search_one_cell_time
    #t0 = time.time()
    while not element_found and is_more_elements:
        triangles, branch = root.expand_search()
        if branch == []:
            # Searching all the verts from the root cell that haven't
            # been searched.  This is the last try
            element_found, sigma0, sigma1, sigma2, k = \
                           _search_triangles_of_vertices(mesh, triangles, x)
            is_more_elements = False
        else:
            element_found, sigma0, sigma1, sigma2, k = \
                       _search_triangles_of_vertices(mesh, triangles, x)
        #search_more_cells_time += time.time()-t0
    #print "search_more_cells_time", search_more_cells_time
        
    return element_found, sigma0, sigma1, sigma2, k


def _search_triangles_of_vertices(mesh, triangles, x):
    """Search for triangle containing x amongs candidate_vertices in mesh

    This is called by search_tree_of_vertices once the appropriate node
    has been found from the quad tree.
    

    This function is responsible for most of the compute time in
    fit and interpolate.
    """
    global last_triangle
    
    # These statments are needed if triangles is empty
    element_found = False
    sigma2 = -10.0
    sigma0 = -10.0
    sigma1 = -10.0
    k = -10

    # For all vertices in same cell as point x
    for k, tri_verts_norms in triangles:
        tri = tri_verts_norms[0]

        # k is the triangle index
        # tri is a list of verts (x, y), representing a tringle
        # Find triangle that contains x (if any) and interpolate
        
        # Input check disabled to speed things up.
        if is_inside_triangle(x, tri, 
                              closed=True,
                              check_inputs=False):
            
            n0, n1, n2 = tri_verts_norms[1]        
            sigma0, sigma1, sigma2 =\
                compute_interpolation_values(tri, n0, n1, n2, x)
                
            element_found = True    

            # Don't look for any other triangles in the triangle list
            last_triangle = [[k, tri_verts_norms]]
            break            

            
    return element_found, sigma0, sigma1, sigma2, k

            
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
