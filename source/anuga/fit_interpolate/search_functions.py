"""
General functions used in fit and interpolate.

   Ole Nielsen, Stephen Roberts, Duncan Gray
   Geoscience Australia, 2006.

"""
from Numeric import dot

from anuga.utilities.numerical_tools import get_machine_precision

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
    #Find triangle containing x:
    element_found = False

    # This will be returned if element_found = False
    sigma2 = -10.0
    sigma0 = -10.0
    sigma1 = -10.0
    k = -10.0
            
    #Get triangles in the cell that the point is in.
    # Triangle is a list, first element triangle_id,
    # second element the triangle
    triangles = root.search(x[0], x[1])
    is_more_elements = True

    element_found, sigma0, sigma1, sigma2, k = \
                   _search_triangles_of_vertices(mesh,
                                                 triangles, x)
    while not element_found and is_more_elements:
        triangles, branch = root.expand_search()
        if branch == []:
            # Searching all the verts from the root cell that haven't
            # been searched.  This is the last try
            element_found, sigma0, sigma1, sigma2, k = \
                           _search_triangles_of_vertices(mesh,triangles, x)
            is_more_elements = False
        else:
            element_found, sigma0, sigma1, sigma2, k = \
                       _search_triangles_of_vertices(mesh,triangles, x)

    return element_found, sigma0, sigma1, sigma2, k


def _search_triangles_of_vertices(mesh, triangles, x):
    """Search for triangle containing x amongs candidate_vertices in mesh

    This is called by search_tree_of_vertices once the appropriate node
    has been found from the quad tree.
    

    This function is responsible for most of the compute time in
    fit and interpolate.
    """

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
            break
    return element_found, sigma0, sigma1, sigma2, k


            
def find_triangle_compute_interpolation(triangle, n0, n1, n2, x):
    """Compute linear interpolation of point x and triangle k in mesh.
    It is assumed that x belongs to triangle k.
    """

    # Get the three vertex_points of candidate triangle k
    xi0, xi1, xi2 = triangle

    # this is where we can call some fast c code.
      
    # Integrity check - machine precision is too hard
    # so we use hardwired single precision 
    epsilon = 1.0e-6
    
    if  x[0] > max(xi0[0], xi1[0], xi2[0]) + epsilon:
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
    
    sigma0 = dot((x-xi1), n0)/dot((xi0-xi1), n0)
    if sigma0 < -epsilon:
        return False,0,0,0
    sigma1 = dot((x-xi2), n1)/dot((xi1-xi2), n1)
    if sigma1 < -epsilon:
        return False,0,0,0
    sigma2 = dot((x-xi0), n2)/dot((xi2-xi0), n2)
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
