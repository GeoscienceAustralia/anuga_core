
from Numeric import concatenate, reshape, take, allclose
from Numeric import array, zeros, Int, Float, sqrt, sum

from anuga.coordinate_transforms.geo_reference import Geo_reference

class General_mesh:
    """Collection of 2D triangular elements

    A triangular element is defined in terms of three vertex ids,
    ordered counter clock-wise, each corresponding to a given node
    which is represented as a coordinate set (x,y).
    Vertices from different triangles can point to the same node.
    The nodes are implemented as an Nx2 Numeric array containing the
    x and y coordinates.


    To instantiate:
       Mesh(nodes, triangles)

    where

      nodes is either a list of 2-tuples or an Nx2 Numeric array of
      floats representing all x, y coordinates in the mesh.

      triangles is either a list of 3-tuples or an Mx3 Numeric array of
      integers representing indices of all vertices in the mesh.
      Each vertex is identified by its index i in [0, N-1].


    Example:

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        e = [2.0, 2.0]

        nodes = [a, b, c, e]
        triangles = [ [1,0,2], [1,2,3] ]   # bac, bce

        # Create mesh with two triangles: bac and bce    
        mesh = Mesh(nodes, triangles)



    Other:

      In addition mesh computes an Mx6 array called vertex_coordinates.
      This structure is derived from coordinates and contains for each
      triangle the three x,y coordinates at the vertices.

      See neighbourmesh.py for a specialisation of the general mesh class
      which includes information about neighbours and the mesh boundary.

      The mesh object is purely geometrical and contains no information
      about quantities defined on the mesh.

    """

    #FIXME: It would be a good idea to use geospatial data as an alternative
    #input
    def __init__(self, nodes, triangles,
                 geo_reference=None,                 
                 number_of_full_nodes=None,
                 number_of_full_triangles=None,                 
                 verbose=False):
        """Build triangular 2d mesh from nodes and triangle information

        Input:
        
          nodes: x,y coordinates represented as a sequence of 2-tuples or
                 a Nx2 Numeric array of floats.
                 
          triangles: sequence of 3-tuples or Mx3 Numeric array of
                     non-negative integers representing indices into
                     the nodes array.
       
          georeference (optional): If specified coordinates are
          assumed to be relative to this origin.


        number_of_full_nodes and number_of_full_triangles relate to
        parallelism when each mesh has an extra layer of ghost points and
        ghost triangles attached to the end of the two arrays.
        In this case it is usefull to specify the number of real (called full)
        nodes and triangles. If omitted they will default to all.
          
        """

        if verbose: print 'General_mesh: Building basic mesh structure' 

        self.triangles = array(triangles, Int)
        self.nodes = array(nodes, Float)


        # Register number of elements and nodes 
        self.number_of_triangles = N = self.triangles.shape[0]
        self.number_of_nodes = self.nodes.shape[0]        

        

        if number_of_full_nodes is None:
            self.number_of_full_nodes = self.number_of_nodes
        else:
            assert int(number_of_full_nodes)
            self.number_of_full_nodes = number_of_full_nodes            


        if number_of_full_triangles is None:
            self.number_of_full_triangles = self.number_of_triangles           
        else:
            assert int(number_of_full_triangles)            
            self.number_of_full_triangles = number_of_full_triangles
        
        
        #print self.number_of_full_nodes, self.number_of_nodes
        #print self.number_of_full_triangles, self.number_of_triangles
        
            

        # FIXME: this stores a geo_reference, but when coords are returned
        # This geo_ref is not taken into account!
        if geo_reference is None:
            self.geo_reference = Geo_reference() #Use defaults 
        else:
            self.geo_reference = geo_reference

        # Input checks
        msg = 'Triangles must an Mx3 Numeric array or a sequence of 3-tuples. '
        msg += 'The supplied array has the shape: %s'\
               %str(self.triangles.shape)
        assert len(self.triangles.shape) == 2, msg

        msg = 'Nodes must an Nx2 Numeric array or a sequence of 2-tuples'
        msg += 'The supplied array has the shape: %s'\
               %str(self.nodes.shape)
        assert len(self.nodes.shape) == 2, msg

        msg = 'Vertex indices reference non-existing coordinate sets'
        assert max(max(self.triangles)) <= self.nodes.shape[0], msg


        # FIXME: Maybe move to statistics?
        # Or use with get_extent
        xy_extent = [ min(self.nodes[:,0]), min(self.nodes[:,1]) ,
                      max(self.nodes[:,0]), max(self.nodes[:,1]) ]

        self.xy_extent = array(xy_extent, Float)


        # Allocate space for geometric quantities
        self.normals = zeros((N, 6), Float)
        self.areas = zeros(N, Float)
        self.edgelengths = zeros((N, 3), Float)

        # Get x,y coordinates for all triangles and store
        self.vertex_coordinates = V = self.compute_vertex_coordinates()


        # Initialise each triangle
        if verbose:
            print 'General_mesh: Computing areas, normals and edgelenghts'
            
        for i in range(N):
            if verbose and i % ((N+10)/10) == 0: print '(%d/%d)' %(i, N)
           
            x0, y0 = V[3*i, :]
            x1, y1 = V[3*i+1, :]
            x2, y2 = V[3*i+2, :]            

            # Area
            self.areas[i] = abs((x1*y0-x0*y1)+(x2*y1-x1*y2)+(x0*y2-x2*y0))/2

            msg = 'Triangle (%f,%f), (%f,%f), (%f, %f)' %(x0,y0,x1,y1,x2,y2)
            msg += ' is degenerate:  area == %f' %self.areas[i]
            assert self.areas[i] > 0.0, msg


            # Normals
            # The normal vectors
            #   - point outward from each edge
            #   - are orthogonal to the edge
            #   - have unit length
            #   - Are enumerated according to the opposite corner:
            #     (First normal is associated with the edge opposite
            #     the first vertex, etc)
            #   - Stored as six floats n0x,n0y,n1x,n1y,n2x,n2y per triangle

            n0 = array([x2 - x1, y2 - y1])
            l0 = sqrt(sum(n0**2))

            n1 = array([x0 - x2, y0 - y2])
            l1 = sqrt(sum(n1**2))

            n2 = array([x1 - x0, y1 - y0])
            l2 = sqrt(sum(n2**2))

            # Normalise
            n0 /= l0
            n1 /= l1
            n2 /= l2

            # Compute and store
            self.normals[i, :] = [n0[1], -n0[0],
                                  n1[1], -n1[0],
                                  n2[1], -n2[0]]

            # Edgelengths
            self.edgelengths[i, :] = [l0, l1, l2]

        
        # Build vertex list
        if verbose: print 'Building vertex list'         
        self.build_vertexlist()

            

    def __len__(self):
        return self.number_of_triangles
    

    def __repr__(self):
        return 'Mesh: %d vertices, %d triangles'\
               %(self.nodes.shape[0], len(self))

    def get_normals(self):
        """Return all normal vectors.
        Return normal vectors for all triangles as an Nx6 array
        (ordered as x0, y0, x1, y1, x2, y2 for each triangle)
        """
        return self.normals


    def get_normal(self, i, j):
        """Return normal vector j of the i'th triangle.
        Return value is the numeric array slice [x, y]
        """
        return self.normals[i, 2*j:2*j+2]


    def get_nodes(self, absolute=False):
        """Return all nodes in mesh.

        The nodes are ordered in an Nx2 array where N is the number of nodes.
        This is the same format they were provided in the constructor
        i.e. without any duplication.

        Boolean keyword argument absolute determines whether coordinates
        are to be made absolute by taking georeference into account
        Default is False as many parts of ANUGA expects relative coordinates.
        (To see which, switch to default absolute=True and run tests).        
        """

        N = self.number_of_full_nodes
        V = self.nodes[:N,:]
        if absolute is True:
            if not self.geo_reference.is_absolute():
                V = self.geo_reference.get_absolute(V)
                
        return V
        
        

    def get_vertex_coordinates(self, absolute=False):
        """Return vertex coordinates for all triangles. 
        
        Return all vertex coordinates for all triangles as a 3*M x 2 array
        where the jth vertex of the ith triangle is located in row 3*i+j and
        M the number of triangles in the mesh.

        Boolean keyword argument absolute determines whether coordinates
        are to be made absolute by taking georeference into account
        Default is False as many parts of ANUGA expects relative coordinates.
        """
        
        V = self.vertex_coordinates
        if absolute is True:
            if not self.geo_reference.is_absolute():
                V = self.geo_reference.get_absolute(V)
            
        return V



    def get_vertex_coordinate(self, i, j, absolute=False):
        """Return coordinates for vertex j of the i'th triangle.
        Return value is the numeric array slice [x, y]
        """

        V = self.get_vertex_coordinates(absolute=absolute)
        return V[3*i+j, :]


    def compute_vertex_coordinates(self):
        """Return all vertex coordinates for all triangles as a 3*M x 2 array
        where the jth vertex of the ith triangle is located in row 3*i+j.

        This function is used to precompute this important structure. Use
        get_vertex coordinates to retrieve the points.
        """

        M = self.number_of_triangles
        vertex_coordinates = zeros((3*M, 2), Float)

        for i in range(M):
            for j in range(3):
                k = self.triangles[i,j] # Index of vertex j in triangle i
                vertex_coordinates[3*i+j,:] = self.nodes[k]

        return vertex_coordinates



    def get_triangles(self, indices=None):
        """Get mesh triangles.

        Return Mx3 integer array where M is the number of triangles.
        Each row corresponds to one triangle and the three entries are
        indices into the mesh nodes which can be obtained using the method
        get_nodes()

        Optional argument, indices is the set of triangle ids of interest.
        """

        M = self.number_of_full_triangles

        if indices is None:
            indices = range(M)

        return take(self.triangles, indices)
    


    def get_disconnected_triangles(self):
        """Get mesh based on nodes obtained from get_vertex_coordinates.

        Return array Mx3 array of integers where each row corresponds to
        a triangle. A triangle is a triplet of indices into
        point coordinates obtained from get_vertex_coordinates and each
        index appears only once

        This provides a mesh where no triangles share nodes
        (hence the name disconnected triangles) and different
        nodes may have the same coordinates.

        This version of the mesh is useful for storing meshes with
        discontinuities at each node and is e.g. used for storing
        data in sww files.

        The triangles created will have the format

        [[0,1,2],
         [3,4,5],
         [6,7,8],
         ...
         [3*M-3 3*M-2 3*M-1]]          
        """

        M = len(self) # Number of triangles
        K = 3*M       # Total number of unique vertices
        T = reshape(array(range(K)).astype(Int), (M,3))

        return T     

    

    def get_unique_vertices(self,  indices=None):
        """FIXME(Ole): This function needs a docstring
        """
        triangles = self.get_triangles(indices=indices)
        unique_verts = {}
        for triangle in triangles:
           unique_verts[triangle[0]] = 0
           unique_verts[triangle[1]] = 0
           unique_verts[triangle[2]] = 0
        return unique_verts.keys()


    def build_vertexlist(self):
        """Build vertexlist indexed by vertex ids and for each entry (point id)
        build a list of (triangles, vertex_id) pairs that use the point
        as vertex.

        The vertex list will have length N, where N is the number of nodes
        in the mesh.

        Preconditions:
          self.nodes and self.triangles are defined

        Postcondition:
          self.vertexlist is built
        """

        vertexlist = [None]*self.number_of_nodes
        for i in range(self.number_of_triangles):

            a = self.triangles[i, 0]
            b = self.triangles[i, 1]
            c = self.triangles[i, 2]

            #Register the vertices v as lists of
            #(triangle_id, vertex_id) tuples associated with them
            #This is used for averaging multiple vertex values.
            for vertex_id, v in enumerate([a,b,c]):
                if vertexlist[v] is None:
                    vertexlist[v] = []

                vertexlist[v].append( (i, vertex_id) )

        self.vertexlist = vertexlist


    def get_extent(self, absolute=False):
        """Return min and max of all x and y coordinates

        Boolean keyword argument absolute determines whether coordinates
        are to be made absolute by taking georeference into account
        """
        


        C = self.get_vertex_coordinates(absolute=absolute)
        X = C[:,0:6:2].copy()
        Y = C[:,1:6:2].copy()

        xmin = min(X.flat)
        xmax = max(X.flat)
        ymin = min(Y.flat)
        ymax = max(Y.flat)

        return xmin, xmax, ymin, ymax


    def get_area(self):
    	"""Return total area of mesh
        """

        return sum(self.areas)

        
