#from builtins import str
#from builtins import range
#from builtins import object
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
        c = [2.0, 0.0]
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
                 use_inscribed_circle=False,
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

        self.verbose = verbose

        if verbose: log.critical('General_mesh: Building basic mesh structure')

        self.use_inscribed_circle = use_inscribed_circle

        self.triangles = num.array(triangles, int)

        if verbose:
            log.timingInfo("numTriangles, " + str(self.triangles.shape[0]))

        self.nodes = num.array(nodes, float)

        # Register number of elements and nodes
        self.number_of_triangles = N = int(self.triangles.shape[0])
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

        self.xy_extent = num.array(xy_extent, float)

        # Allocate space for geometric quantities
        self.normals = num.zeros((N, 6), float)
        self.areas = num.zeros(N, float)
        self.edgelengths = num.zeros((N, 3), float)

        # Get x,y coordinates for all triangle vertices and store
        self.centroid_coordinates = num.zeros((N, 2), float)

        #Allocate space for geometric quantities
        self.radii = num.zeros(N, float)

        # Get x,y coordinates for all triangle vertices and store
        self.vertex_coordinates = V = self.compute_vertex_coordinates()

        # Get x,y coordinates for all triangle edge midpoints and store
        self.edge_midpoint_coordinates  = self.compute_edge_midpoint_coordinates()

        # Initialise each triangle
        if verbose:
            log.critical('General_mesh: Computing areas, normals, '
                         'edgelengths, centroids and radii')


        # Calculate Areas
        V0 = V[0:3*N:3, :]
        V1 = V[1:3*N:3, :]
        V2 = V[2:3*N:3, :]


        # Area
        x0 = V0[:,0]
        y0 = V0[:,1]
        x1 = V1[:,0]
        y1 = V1[:,1]
        x2 = V2[:,0]
        y2 = V2[:,1]

        self.areas[:] = -((x1*y0-x0*y1) + (x2*y1-x1*y2) + (x0*y2-x2*y0))/2.0

        #areas = -((x0-x1)*(y2-y1) - (y0-y1)*(x2-x1))/2.0

        #assert num.allclose(self.areas, areas)

        ind = num.where(self.areas <= 0.0)
        msg = 'Degenerate Triangle(s) '+str(ind[0])
        assert num.all(self.areas > 0.0), msg


        #print V.shape, V0.shape, V1.shape, V2.shape

#        #print E.shape, E[0:3*M:3, :].shape, E[1:3*M:3, :].shape, E[2:3*M:3, :].shape
#        E[0:3*M:3, :] = 0.5*(V1+V2)
#        E[1:3*M:3, :] = 0.5*(V2+V0)
#        E[2:3*M:3, :] = 0.5*(V0+V1)

        i0 = self.triangles[:,0]
        i1 = self.triangles[:,1]
        i2 = self.triangles[:,2]

        assert num.allclose( x0, self.nodes[i0,0] )
        assert num.allclose( y0, self.nodes[i0,1] )

        assert num.allclose( x1, self.nodes[i1,0] )
        assert num.allclose( y1, self.nodes[i1,1] )

        assert num.allclose( x2, self.nodes[i2,0] )
        assert num.allclose( y2, self.nodes[i2,1] )


        xn0 = x2-x1
        yn0 = y2-y1
        l0 = num.sqrt(xn0**2 + yn0**2)

        xn0 /= l0
        yn0 /= l0

        xn1 = x0-x2
        yn1 = y0-y2
        l1 = num.sqrt(xn1**2 + yn1**2)

        xn1 /= l1
        yn1 /= l1

        xn2 = x1-x0
        yn2 = y1-y0
        l2 = num.sqrt(xn2**2 + yn2**2)

        xn2 /= l2
        yn2 /= l2

        # Compute and store

        self.normals[:,0] =  yn0
        self.normals[:,1] = -xn0

        self.normals[:,2] =  yn1
        self.normals[:,3] = -xn1

        self.normals[:,4] =  yn2
        self.normals[:,5] = -xn2

        self.edgelengths[:,0] = l0
        self.edgelengths[:,1] = l1
        self.edgelengths[:,2] = l2

        self.centroid_coordinates[:,0] = (x0 + x1 + x2) / 3
        self.centroid_coordinates[:,1] = (y0 + y1 + y2) / 3



        if self.use_inscribed_circle == False:
            #OLD code. Computed radii may exceed that of an
            #inscribed circle

            #Midpoints
            xm0 = (x1 + x2) / 2
            ym0 = (y1 + y2) / 2

            xm1 = (x2 + x0) / 2
            ym1 = (y2 + y0) / 2

            xm2 = (x0 + x1) / 2
            ym2 = (y0 + y1) / 2


            #The radius is the distance from the centroid of
            #a triangle to the midpoint of the side of the triangle
            #closest to the centroid

            d0 = num.sqrt((self.centroid_coordinates[:,0] - xm0)**2 + (self.centroid_coordinates[:,1] - ym0)**2)
            d1 = num.sqrt((self.centroid_coordinates[:,0] - xm1)**2 + (self.centroid_coordinates[:,1] - ym1)**2)
            d2 = num.sqrt((self.centroid_coordinates[:,0] - xm2)**2 + (self.centroid_coordinates[:,1] - ym2)**2)


            self.radii[:] = num.minimum(num.minimum(d0, d1), d2)

        else:
            #NEW code added by Peter Row. True radius
            #of inscribed circle is computed

            a = num.sqrt((x0-x1)**2+(y0-y1)**2)
            b = num.sqrt((x1-x2)**2+(y1-y2)**2)
            c = num.sqrt((x2-x0)**2+(y2-y0)**2)

            self.radii[:] = 2.0 * self.areas / (a+b+c)



#        for i in range(N):
#            if verbose and i % ((N+10)/10) == 0: log.critical('(%d/%d)' % (i, N))
#
#            x0, y0 = V[3*i, :]
#            x1, y1 = V[3*i+1, :]
#            x2, y2 = V[3*i+2, :]
#
#
#            i0 = self.triangles[i][0]
#            i1 = self.triangles[i][1]
#            i2 = self.triangles[i][2]
#
##            assert x0 == self.nodes[i0][0]
##            assert y0 == self.nodes[i0][1]
##
##            assert x1 == self.nodes[i1][0]
##            assert y1 == self.nodes[i1][1]
##
##            assert x2 == self.nodes[i2][0]
##            assert y2 == self.nodes[i2][1]
#
##            # Area
##            self.areas[i] = abs((x1*y0-x0*y1) + (x2*y1-x1*y2) + (x0*y2-x2*y0))/2
##
##            msg = 'Triangle %g (%f,%f), (%f,%f), (%f, %f)' % (i,x0,y0,x1,y1,x2,y2)
##            msg += ' is degenerate:  area == %f' % self.areas[i]
##            assert self.areas[i] > 0.0, msg
#
#            # Normals
#            # The normal vectors
#            #   - point outward from each edge
#            #   - are orthogonal to the edge
#            #   - have unit length
#            #   - Are enumerated according to the opposite corner:
#            #     (First normal is associated with the edge opposite
#            #     the first vertex, etc)
#            #   - Stored as six floats n0x,n0y,n1x,n1y,n2x,n2y per triangle
#            n0 = num.array([x2-x1, y2-y1], float)
#            l0 = num.sqrt(num.sum(n0**2))
#
#            n1 = num.array([x0-x2, y0-y2], float)
#            l1 = num.sqrt(num.sum(n1**2))
#
#            n2 = num.array([x1-x0, y1-y0], float)
#            l2 = num.sqrt(num.sum(n2**2))
#
#            # Normalise
#            n0 /= l0
#            n1 /= l1
#            n2 /= l2
#
##            # Compute and store
##            self.normals[i, :] = [n0[1], -n0[0],
##                                  n1[1], -n1[0],
##                                  n2[1], -n2[0]]
#
#            # Edgelengths
#            #self.edgelengths[i, :] = [l0, l1, l2]
#
#
#
#            #Compute centroid
##            centroid = num.array([(x0 + x1 + x2)/3, (y0 + y1 + y2)/3], float)
###            self.centroid_coordinates[i] = centroid
##
##
##            if self.use_inscribed_circle == False:
##                #OLD code. Computed radii may exceed that of an
##                #inscribed circle
##
##                #Midpoints
##                m0 = num.array([(x1 + x2)/2, (y1 + y2)/2], float)
##                m1 = num.array([(x0 + x2)/2, (y0 + y2)/2], float)
##                m2 = num.array([(x1 + x0)/2, (y1 + y0)/2], float)
##
##                #The radius is the distance from the centroid of
##                #a triangle to the midpoint of the side of the triangle
##                #closest to the centroid
##                d0 = num.sqrt(num.sum( (centroid-m0)**2 ))
##                d1 = num.sqrt(num.sum( (centroid-m1)**2 ))
##                d2 = num.sqrt(num.sum( (centroid-m2)**2 ))
##
##                #self.radii[i] = min(d0, d1, d2)
##
##            else:
##                #NEW code added by Peter Row. True radius
##                #of inscribed circle is computed
##
##                a = num.sqrt((x0-x1)**2+(y0-y1)**2)
##                b = num.sqrt((x1-x2)**2+(y1-y2)**2)
##                c = num.sqrt((x2-x0)**2+(y2-y0)**2)
##
##                self.radii[i]=2.0*self.areas[i]/(a+b+c)


        # Build structure listing which triangles belong to which node.
        if verbose: log.critical('General Mesh: Building inverted triangle structure')
        self.build_inverted_triangle_structure()

        if verbose: log.timingInfo("aoi, '%s'" % self.get_area())


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
                                self.geo_reference.get_yllcorner()], float)
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
                                  self.geo_reference.get_yllcorner()], float)

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
        vertex_coordinates = num.zeros((3*M, 2), float)

        k0 = self.triangles[:,0]
        k1 = self.triangles[:,1]
        k2 = self.triangles[:,2]

#        I = num.arange(M,dtype=int)
#
#        V0 = V[0:3*M:3, :]
#        V1 = V[1:3*M:3, :]
#        V2 = V[2:3*M:3, :]

        vertex_coordinates[0:3*M:3,:] = self.nodes[k0,:]
        vertex_coordinates[1:3*M:3,:] = self.nodes[k1,:]
        vertex_coordinates[2:3*M:3,:] = self.nodes[k2,:]

#        for i in range(M):
#            for j in range(3):
#                k = self.triangles[i,j] # Index of vertex j in triangle i
#                vertex_coordinates[3*i+j,:] = self.nodes[k]

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
                                  self.geo_reference.get_yllcorner()], float)

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
        E = num.zeros((3*M, 2), float)

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
        return num.reshape(num.arange(K, dtype=int), (M,3))

    def get_unique_vertices(self, indices=None):
        """Return indices to vertices as a sorted list.
           FIXME (Ole): It may not be needed anymore
        """

        triangles = self.get_triangles(indices=indices)
        unique_verts = {}
        for triangle in triangles:
            unique_verts[triangle[0]] = 0
            unique_verts[triangle[1]] = 0
            unique_verts[triangle[2]] = 0
        res = list(unique_verts.keys())
        res.sort() # Ensure uniqueness
        return res

        # Note Padarn 27/11/12:
        # This function was modified, but then it was deicded it was not
        # needed. It should be restored if it is used elsewhere in the code
        # (it was being used in quantity.py in the _set_vertex_values function).
        # Note however, the function in the head of the code is very slow and
        # could be easily sped up many fold.
        #
        # Have we profiled it? (Ole 31/5/2020)
        
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
            #first = num.sum(self.number_of_triangles_per_node[:node])

            first = self.node_index[node]
            # Get number of triangles for this node
            count = self.number_of_triangles_per_node[node]

            for i in range(count):
                index = self.vertex_value_indices[first+i]

                # FIXME(Ole): This must be floor division ('//')
                # However, tests pass either way. Need to update tests.
                volume_id = index // 3  
                #print('volume_id', volume_id, index)
                
                vertex_id = index % 3

                triangle_list.append( (volume_id, vertex_id) )

            triangle_list = num.array(triangle_list, int)    #array default#
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
        index can be used directly into a flat triangles array. This
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
#        number_of_triangles_per_node = num.zeros(self.number_of_nodes,
#                                                 int)       #array default#
#        for volume_id, triangle in enumerate(self.get_triangles()):
#            for vertex_id in triangle:
#                number_of_triangles_per_node[vertex_id] += 1

        # Need to pad number_of_triangles_per_node in case lone nodes at end of list
        #number_of_triangles_per_node = num.zeros(self.number_of_nodes, int)


        number_of_triangles_per_node = num.bincount(self.triangles.flat).astype(int)
        number_of_lone_nodes = self.number_of_nodes - len(number_of_triangles_per_node)


        orphan_nodes = num.argwhere(number_of_triangles_per_node==0)
        number_of_orphan_nodes = len(orphan_nodes)

        if number_of_orphan_nodes > 0 and self.verbose:
            msg = 'Node(s) %d not associated to a triangle.' % orphan_nodes[0]
            print(msg)

        if number_of_lone_nodes > 0:
            number_of_triangles_per_node =  \
               num.append(number_of_triangles_per_node,num.zeros(number_of_lone_nodes,int))

        #assert num.allclose(number_of_triangles_per_node_new, number_of_triangles_per_node)

        # Allocate space for inverted structure
        number_of_entries = num.sum(number_of_triangles_per_node)

        assert number_of_entries == 3*self.number_of_triangles

        #vertex_value_indices = num.zeros(number_of_entries, int) #array default#

        # Array of vertex_indices (3*vol_id+vertex_id) sorted into contiguous
        # order around each node. Use with number_of_triangles_per_node to
        # find vertices associated with a node.
        # ie There are  number_of_triangles_per_node[i] vertices
        vertex_value_indices = num.argsort(self.triangles.flat).astype(int)
        #vertex_value_indices = num.argsort(self.triangles.flatten())

#        node_index = num.zeros((self.number_of_nodes)+1, dtype = int)
#        node_index[0] = 0
#        for i in xrange(self.number_of_nodes):
#            node_index[i+1] = node_index[i] + number_of_triangles_per_node[i]

        node_index = num.zeros((self.number_of_nodes)+1, dtype = int)
        node_index[1:] = num.cumsum(number_of_triangles_per_node)

        #assert num.allclose(node_index,node_index_new)




        # Save structures
        self.node_index = node_index
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
