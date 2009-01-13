"""
General functions used in fit and interpolate.

   Ole Nielsen, Stephen Roberts, Duncan Gray
   Geoscience Australia, 2006.

"""
import time

from anuga.utilities.numerical_tools import get_machine_precision
from anuga.config import max_float

import Numeric as num


initial_search_value = 'uncomment search_functions code first'#0
search_one_cell_time = initial_search_value
search_more_cells_time = initial_search_value

#FIXME test what happens if a 
LAST_TRIANGLE = [[-10,[(num.array([max_float,max_float]),
                        num.array([max_float,max_float]),
                        num.array([max_float,max_float])),
                       (num.array([1,1]),num.array([0,0]),num.array([-1.1,-1.1]))]]]

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

    #Find triangle containing x:
    element_found = False

    # This will be returned if element_found = False
    sigma2 = -10.0
    sigma0 = -10.0
    sigma1 = -10.0
    k = -10.0

    # Search the last triangle first
    try:
        element_found, sigma0, sigma1, sigma2, k = \
            _search_triangles_of_vertices(mesh, last_triangle, x)
    except:
        element_found = False
                   
    #print "last_triangle", last_triangle
    if element_found is True:
        #print "last_triangle", last_triangle
        return element_found, sigma0, sigma1, sigma2, k

    # This was only slightly faster than just checking the
    # last triangle and it significantly slowed down
    # non-gridded fitting
    # If the last element was a dud, search its neighbours
    #print "last_triangle[0][0]", last_triangle[0][0]
    #neighbours = mesh.get_triangle_neighbours(last_triangle[0][0])
    #print "neighbours", neighbours
    #neighbours = []
  #   for k in neighbours:
#         if k >= 0:
#             tri = mesh.get_vertex_coordinates(k,
#                                                    absolute=True)
#             n0 = mesh.get_normal(k, 0)
#             n1 = mesh.get_normal(k, 1)
#             n2 = mesh.get_normal(k, 2) 
#             triangle =[[k,(tri, (n0, n1, n2))]]
#             element_found, sigma0, sigma1, sigma2, k = \
#                            _search_triangles_of_vertices(mesh,
#                                                          triangle, x)
#             if element_found is True:
#                 return element_found, sigma0, sigma1, sigma2, k
            
    #t0 = time.time()
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
    
    # these statments are needed if triangles is empty
    #Find triangle containing x:
    element_found = False

    # This will be returned if element_found = False
    sigma2 = -10.0
    sigma0 = -10.0
    sigma1 = -10.0
    k = -10
    
    #For all vertices in same cell as point x
    for k, tri_verts_norms in triangles:
        tri = tri_verts_norms[0]
        n0, n1, n2 = tri_verts_norms[1]
        # k is the triangle index
        # tri is a list of verts (x, y), representing a tringle
        # Find triangle that contains x (if any) and interpolate
        element_found, sigma0, sigma1, sigma2 =\
                       find_triangle_compute_interpolation(tri, n0, n1, n2, x)
        if element_found is True:
            # Don't look for any other triangles in the triangle list
            last_triangle = [[k,tri_verts_norms]]
            break
    return element_found, sigma0, sigma1, sigma2, k


            
def find_triangle_compute_interpolation(triangle, n0, n1, n2, x):
    """Compute linear interpolation of point x and triangle k in mesh.
    It is assumed that x belongs to triangle k.max_float
    """

    # Get the three vertex_points of candidate triangle k
    xi0, xi1, xi2 = triangle

    # this is where we can call some fast c code.
      
    # Integrity check - machine precision is too hard
    # so we use hardwired single precision 
    epsilon = 1.0e-6
    
    if  x[0] > max(xi0[0], xi1[0], xi2[0]) + epsilon:
        # print "max(xi0[0], xi1[0], xi2[0])", max(xi0[0], xi1[0], xi2[0])
        return False,0,0,0
    if  x[0] < min(xi0[0], xi1[0], xi2[0]) - epsilon:
        return False,0,0,0
    if  x[1] > max(xi0[1], xi1[1], xi2[1]) + epsilon:
        return False,0,0,0
    if  x[1] < min(xi0[1], xi1[1], xi2[1]) - epsilon:
        return False,0,0,0
    
    # machine precision on some machines (e.g. nautilus)
    epsilon = get_machine_precision() * 2
    
    # Compute interpolation - return as soon as possible
    #  print "(xi0-xi1)", (xi0-xi1)
    # print "n0", n0
    # print "dot((xi0-xi1), n0)", dot((xi0-xi1), n0)
    
    sigma0 = num.dot((x-xi1), n0)/num.dot((xi0-xi1), n0)
    if sigma0 < -epsilon:
        return False,0,0,0
    sigma1 = num.dot((x-xi2), n1)/num.dot((xi1-xi2), n1)
    if sigma1 < -epsilon:
        return False,0,0,0
    sigma2 = num.dot((x-xi0), n2)/num.dot((xi2-xi0), n2)
    if sigma2 < -epsilon:
        return False,0,0,0
    
    # epsilon = 1.0e-6
    # we want to speed this up, so don't do assertions
    #delta = abs(sigma0+sigma1+sigma2-1.0) # Should be close to zero
    #msg = 'abs(sigma0+sigma1+sigma2-1) = %.15e, eps = %.15e'\
    #      %(delta, epsilon)
    #assert delta < epsilon, msg


    # Check that this triangle contains the data point
    # Sigmas are allowed to get negative within
    # machine precision on some machines (e.g. nautilus)
    #if sigma0 >= -epsilon and sigma1 >= -epsilon and sigma2 >= -epsilon:
    #    element_found = True
    #else:
    #    element_found = False 
    return True, sigma0, sigma1, sigma2

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
