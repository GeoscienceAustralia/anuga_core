import copy
import numpy as num

from anuga.coordinate_transforms.geo_reference import Geo_reference
import anuga.utilities.log as log

class General_mesh:
    """Collection of 2D triangular elements

    A triangular element is defined in terms of three vertex ids,
    ordered counter clock-wise, each corresponding to a given node
    which is represented as a coordinate set (x,y).
    Vertices from different triangles can point to the same node.
    The nodes are implemented as an Nx2 numeric array containing the
    x and y coordinates.


    To instantiate:
       Mesh(nodes, triangles)

    where

      nodes is either a list of 2-tuples or an Nx2 numeric array of
      floats representing all x, y coordinates in the mesh.

      triangles is either a list of 3-tuples or an Mx3 numeric array of
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

    # FIXME: It would be a good idea to use geospatial data as an alternative
    #        input
    def __init__(self,
                 nodes,
                 triangles,
                 geo_reference=None,
                 verbose=False):
        """Build triangular 2d mesh from nodes and triangle information

        Input:

          nodes: x,y coordinates represented as a sequence of 2-tuples or
                 a Nx2 numeric array of floats.

          triangles: sequence of 3-tuples or Mx3 numeric array of
                     non-negative integers representing indices into
                     the nodes array.

          georeference (optional): If specified coordinates are
          assumed to be relative to this origin.


        """

        if verbose: log.critical('General_mesh: Building basic mesh structure '
                                 'in ANUGA domain')

        self.triangles = num.array(triangles, num.int)
        self.nodes = num.array(nodes, num.float)

        # Register number of elements and nodes
        self.number_of_triangles = N = self.triangles.shape[0]
        self.number_of_nodes = self.nodes.shape[0]


        # FIXME: this stores a geo_reference, but when coords are returned
        # This geo_ref is not taken into account!
        if geo_reference is None:
            self.geo_reference = Geo_reference()    # Use defaults
        else:
            self.geo_reference = geo_reference

        # Input checks
        msg = ('Triangles must an Mx3 numeric array or a sequence of 3-tuples. '
               'The supplied array has the shape: %s'
               % str(self.triangles.shape))
        assert len(self.triangles.shape) == 2, msg

        msg = ('Nodes must an Nx2 numeric array or a sequence of 2-tuples'
               'The supplied array has the shape: %s' % str(self.nodes.shape))
        assert len(self.nodes.shape) == 2, msg

        msg = 'Vertex indices reference non-existing coordinate sets'
        assert num.max(self.triangles) < self.nodes.shape[0], msg

        # FIXME: Maybe move to statistics?
        # Or use with get_extent
        xy_extent = [min(self.nodes[:,0]), min(self.nodes[:,1]),
                     max(self.nodes[:,0]), max(self.nodes[:,1])]

        self.xy_extent = num.array(xy_extent, num.float)

        # Allocate space for geometric quantities
        self.normals = num.zeros((N, 6), num.float)
        self.areas = num.zeros(N, num.float)
        self.edgelengths = num.zeros((N, 3), num.float)

        # Get x,y coordinates for all triangle vertices and store
        self.vertex_coordinates = V = self.compute_vertex_coordinates()

        # Get x,y coordinates for all triangle edge midpoints and store
        self.edge_midpoint_coordinates  = self.compute_edge_midpoint_coordinates()

        # Initialise each triangle
        if verbose:
            log.critical('General_mesh: Computing areas, normals '
                         'and edgelengths')

        for i in range(N):
            if verbose and i % ((N+10)/10) == 0: log.critical('(%d/%d)' % (i, N))

            x0, y0 = V[3*i, :]
            x1, y1 = V[3*i+1, :]
            x2, y2 = V[3*i+2, :]


            i0 = self.triangles[i][0]
            i1 = self.triangles[i][1]
            i2 = self.triangles[i][2]

            assert x0 == self.nodes[i0][0]
            assert y0 == self.nodes[i0][1]

            assert x1 == self.nodes[i1][0]
            assert y1 == self.nodes[i1][1]

            assert x2 == self.nodes[i2][0]
            assert y2 == self.nodes[i2][1]
            
            # Area
            self.areas[i] = abs((x1*y0-x0*y1) + (x2*y1-x1*y2) + (x0*y2-x2*y0))/2

            msg = 'Triangle %g (%f,%f), (%f,%f), (%f, %f)' % (i,x0,y0,x1,y1,x2,y2)
            msg += ' is degenerate:  area == %f' % self.areas[i]
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
            n0 = num.array([x2-x1, y2-y1], num.float)
            l0 = num.sqrt(num.sum(n0**2))

            n1 = num.array([x0-x2, y0-y2], num.float)
            l1 = num.sqrt(num.sum(n1**2))

            n2 = num.array([x1-x0, y1-y0], num.float)
            l2 = num.sqrt(num.sum(n2**2))

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

        # Build structure listing which triangles belong to which node.
        if verbose: log.critical('Building inverted triangle structure')
        self.build_inverted_triangle_structure()

    def __len__(self):
        return self.number_of_triangles

    def __repr__(self):
        return ('Mesh: %d vertices, %d triangles'
                % (self.nodes.shape[0], len(self)))

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
        
    def get_edgelength(self, i, j):
        """Return length of j'th edge of the i'th triangle.
        Return value is the numeric array slice [x, y]
        """
        return self.edgelengths[i, j]
                

    def get_number_of_triangles(self):
        return self.number_of_triangles

    
    def get_number_of_nodes(self):
        return self.number_of_nodes


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

        N = self.number_of_nodes
        V = self.nodes[:N,:]
        if absolute is True:
            if not self.geo_reference.is_absolute():
                V = self.geo_reference.get_absolute(V)

        return V

    def get_node(self, i, absolute=False):
        """Return node coordinates for triangle i.

        Boolean keyword argument absolute determines whether coordinates
        are to be made absolute by taking georeference into account
        Default is False as many parts of ANUGA expects relative coordinates.
        (To see which, switch to default absolute=True and run tests).

        Note: This method returns a modified _copy_ of the nodes slice if 
              absolute is True.  If absolute is False, just return the slice.
              This is related to the ensure_numeric() returning a copy problem.
        """

        V = self.nodes[i,:]
        if absolute is True:
            if not self.geo_reference.is_absolute():
                # get a copy so as not to modify the internal self.nodes array
                V = copy.copy(V)
                V += num.array([self.geo_reference.get_xllcorner(),
                                self.geo_reference.get_yllcorner()], num.float)
        return V

    def get_vertex_coordinates(self, triangle_id=None, absolute=False):
        """Return vertex coordinates for all triangles.

        Return all vertex coordinates for all triangles as a 3*M x 2 array
        where the jth vertex of the ith triangle is located in row 3*i+j and
        M the number of triangles in the mesh.

        if triangle_id is specified (an integer) the 3 vertex coordinates
        for triangle_id are returned.

        Boolean keyword argument absolute determines whether coordinates
        are to be made absolute by taking georeference into account
        Default is False as many parts of ANUGA expects relative coordinates.
        """

        V = self.vertex_coordinates

        if triangle_id is None:
            if absolute is True:
                if not self.geo_reference.is_absolute():
                    V = self.geo_reference.get_absolute(V)
            return V
        else:
            i = triangle_id
            msg = 'triangle_id must be an integer'
            assert int(i) == i, msg
            assert 0 <= i < self.number_of_triangles

            i3 = 3*i
            if absolute is True and not self.geo_reference.is_absolute():
                offset=num.array([self.geo_reference.get_xllcorner(),
                                  self.geo_reference.get_yllcorner()], num.float)
                                  
                return V[i3:i3+3,:] + offset                                  
            else:
                return V[i3:i3+3,:]

    def get_vertex_coordinate(self, i, j, absolute=False):
        """Return coordinates for vertex j of the i'th triangle.
        Return value is the numeric array slice [x, y]
        """

        msg = 'vertex id j must be an integer in [0,1,2]'
        assert j in [0,1,2], msg

        V = self.get_vertex_coordinates(triangle_id=i, absolute=absolute)
        return V[j,:]



    def compute_vertex_coordinates(self):
        """Return all vertex coordinates for all triangles as a 3*M x 2 array
        where the jth vertex of the ith triangle is located in row 3*i+j.

        This function is used to precompute this important structure. Use
        get_vertex coordinates to retrieve the points.
        """

        M = self.number_of_triangles
        vertex_coordinates = num.zeros((3*M, 2), num.float)

        for i in range(M):
            for j in range(3):
                k = self.triangles[i,j] # Index of vertex j in triangle i
                vertex_coordinates[3*i+j,:] = self.nodes[k]

        return vertex_coordinates


    def get_edge_midpoint_coordinates(self, triangle_id=None, absolute=False):
        """Return edge midpoint coordinates for all triangles or from particular triangle.

        Return all edge midpoint coordinates for all triangles as a 3*M x 2 array
        where the jth midpoint of the ith triangle is located in row 3*i+j and
        M the number of triangles in the mesh.

        if triangle_id is specified (an integer) the 3 midpoint coordinates
        for triangle_id are returned.

        Boolean keyword argument absolute determines whether coordinates
        are to be made absolute by taking georeference into account
        Default is False as many parts of ANUGA expects relative coordinates.
        """

        E = self.edge_midpoint_coordinates
        
        if triangle_id is None:
            if absolute is True:
                if not self.geo_reference.is_absolute():
                    E = self.geo_reference.get_absolute(E)
            return E
        else:
            i = triangle_id
            msg = 'triangle_id must be an integer'
            assert int(i) == i, msg
            assert 0 <= i < self.number_of_triangles

            i3 = 3*i
            if absolute is True and not self.geo_reference.is_absolute():
                offset=num.array([self.geo_reference.get_xllcorner(),
                                  self.geo_reference.get_yllcorner()], num.float)

                return E[i3:i3+3,:] + offset                                  
            else:
                return E[i3:i3+3,:]


    def get_edge_midpoint_coordinate(self, i, j, absolute=False):
        """Return coordinates for edge midpoint j of the i'th triangle.
        Return value is the numeric array slice [x, y]
        """

        msg = 'edge midpoint id j must be an integer in [0,1,2]'
        assert j in [0,1,2], msg

        E = self.get_edge_midpoint_coordinates(triangle_id=i, absolute=absolute)
        return E[j,:] # Return (x, y) for edge mid point

    
    def compute_edge_midpoint_coordinates(self):
        """Return all edge midpoint coordinates for all triangles as a 3*M x 2 array
        where the jth edge midpoint of the ith triangle is located in row 3*i+j.

        This function is used to precompute this important structure. Use
        get_edge_midpoint_coordinates to retrieve the points.

        Assumes that vertex_coordinates have been computed
        """

        M = self.number_of_triangles
        E = num.zeros((3*M, 2), num.float)

        V = self.vertex_coordinates

        V0 = V[0:3*M:3, :]
        V1 = V[1:3*M:3, :]
        V2 = V[2:3*M:3, :]

        
        #print V.shape, V0.shape, V1.shape, V2.shape

        #print E.shape, E[0:3*M:3, :].shape, E[1:3*M:3, :].shape, E[2:3*M:3, :].shape
        E[0:3*M:3, :] = 0.5*(V1+V2)
        E[1:3*M:3, :] = 0.5*(V2+V0)
        E[2:3*M:3, :] = 0.5*(V0+V1)

        return E



    def get_triangles(self, indices=None):
        """Get mesh triangles.

        Return Mx3 integer array where M is the number of triangles.
        Each row corresponds to one triangle and the three entries are
        indices into the mesh nodes which can be obtained using the method
        get_nodes()

        Optional argument, indices is the set of triangle ids of interest.
        """


        if indices is None:
            return self.triangles

        return num.take(self.triangles, indices, axis=0)

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
        return num.reshape(num.arange(K, dtype=num.int), (M,3))

    def get_unique_vertices(self, indices=None):
        """FIXME(Ole): This function needs a docstring"""

        triangles = self.get_triangles(indices=indices)
        unique_verts = {}
        for triangle in triangles:
            unique_verts[triangle[0]] = 0
            unique_verts[triangle[1]] = 0
            unique_verts[triangle[2]] = 0
        return unique_verts.keys()

    def get_triangles_and_vertices_per_node(self, node=None):
        """Get triangles associated with given node.

        Return list of triangle_ids, vertex_ids for specified node.
        If node in None or absent, this information will be returned
        for all nodes in a list L where L[v] is the triangle
        list for node v.
        """

        triangle_list = []
        if node is not None:
            # Get index for this node
            first = num.sum(self.number_of_triangles_per_node[:node])

            # Get number of triangles for this node
            count = self.number_of_triangles_per_node[node]

            for i in range(count):
                index = self.vertex_value_indices[first+i]

                volume_id = index / 3
                vertex_id = index % 3

                triangle_list.append( (volume_id, vertex_id) )

            triangle_list = num.array(triangle_list, num.int)    #array default#
        else:
            # Get info for all nodes recursively.
            # If need be, we can speed this up by
            # working directly with the inverted triangle structure
            for i in range(self.number_of_nodes):
                L = self.get_triangles_and_vertices_per_node(node=i)
                triangle_list.append(L)

        return triangle_list

    def build_inverted_triangle_structure(self):
        """Build structure listing triangles belonging to each node

        Two arrays are created and store as mesh attributes

        number_of_triangles_per_node: An integer array of length N
        listing for each node how many triangles use it. N is the number of
        nodes in mesh.

        vertex_value_indices: An array of length M listing indices into
        triangles ordered by node number. The (triangle_id, vertex_id)
        pairs are obtained from each index as (index/3, index%3) or each
        index can be used directly into a flattened triangles array. This
        is for example the case in the quantity.c where this structure is
        used to average vertex values efficiently.

        Example:
        a = [0.0, 0.0] # node 0
        b = [0.0, 2.0] # node 1
        c = [2.0, 0.0] # node 2
        d = [0.0, 4.0] # node 3
        e = [2.0, 2.0] # node 4
        f = [4.0, 0.0] # node 5
        nodes = array([a, b, c, d, e, f])

        #                    bac,     bce,     ecf,     dbe
        triangles = array([[1,0,2], [1,2,4], [4,2,5], [3,1,4]])

        For this structure:
        number_of_triangles_per_node = [1 3 3 1 3 1]
        which means that node a has 1 triangle associated with it, node b
        has 3, node has 3 and so on.

        vertex_value_indices = [ 1  0  3 10  2  4  7  9  5  6 11  8]
        which reflects the fact that
        node 0 is used by triangle 0, vertex 1 (index = 1)
        node 1 is used by triangle 0, vertex 0 (index = 0)
                   and by triangle 1, vertex 0 (index = 3)
                   and by triangle 3, vertex 1 (index = 10)
        node 2 is used by triangle 0, vertex 2 (index = 2)
                   and by triangle 1, vertex 1 (index = 4)
                   and by triangle 2, vertex 1 (index = 7)
        node 3 is used by triangle 3, vertex 0 (index = 9)
        node 4 is used by triangle 1, vertex 2 (index = 5)
                   and by triangle 2, vertex 0 (index = 6)
                   and by triangle 3, vertex 2 (index = 11)
        node 5 is used by triangle 2, vertex 2 (index = 8)

        Preconditions:
          self.nodes and self.triangles are defined

        Postcondition:
          self.number_of_triangles_per_node is built
          self.vertex_value_indices is built
        """

        # Count number of triangles per node
        number_of_triangles_per_node = num.zeros(self.number_of_nodes,
                                                 num.int)       #array default#
        for volume_id, triangle in enumerate(self.get_triangles()):
            for vertex_id in triangle:
                number_of_triangles_per_node[vertex_id] += 1

        # Allocate space for inverted structure
        number_of_entries = num.sum(number_of_triangles_per_node)
        vertex_value_indices = num.zeros(number_of_entries, num.int) #array default#

        # Register (triangle, vertex) indices for each node
        vertexlist = [None] * self.number_of_nodes
        for volume_id in range(self.number_of_triangles):
            a = self.triangles[volume_id, 0]
            b = self.triangles[volume_id, 1]
            c = self.triangles[volume_id, 2]

            for vertex_id, node_id in enumerate([a, b, c]):
                if vertexlist[node_id] is None:
                    vertexlist[node_id] = []
                vertexlist[node_id].append((volume_id, vertex_id))

        # Build inverted triangle index array
        k = 0
        for vertices in vertexlist:
            if vertices is not None:
                for volume_id, vertex_id in vertices:
                    vertex_value_indices[k] = 3*volume_id + vertex_id
                    k += 1

        # Save structure
        self.number_of_triangles_per_node = number_of_triangles_per_node
        self.vertex_value_indices = vertex_value_indices

    def get_extent(self, absolute=False):
        """Return min and max of all x and y coordinates

        Boolean keyword argument absolute determines whether coordinates
        are to be made absolute by taking georeference into account
        """

        C = self.get_vertex_coordinates(absolute=absolute)
        X = C[:,0:6:2].copy()
        Y = C[:,1:6:2].copy()

        xmin = num.min(X)
        xmax = num.max(X)
        ymin = num.min(Y)
        ymax = num.max(Y)

        return xmin, xmax, ymin, ymax

    def get_areas(self):
        """Get areas of all individual triangles."""

        return self.areas

    def get_area(self):
        """Return total area of mesh"""

        return num.sum(self.areas)

    def set_georeference(self, g):
        self.geo_reference = g

    def get_georeference(self):
        return self.geo_reference

