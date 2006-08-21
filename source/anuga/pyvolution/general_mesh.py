
from Numeric import concatenate, reshape, take, allclose
from Numeric import array, zeros, Int, Float, sqrt, sum

from coordinate_transforms.geo_reference import Geo_reference

class General_mesh:
    """Collection of triangular elements (purely geometric)

    A triangular element is defined in terms of three vertex ids,
    ordered counter clock-wise,
    each corresponding to a given coordinate set.
    Vertices from different elements can point to the same
    coordinate set.

    Coordinate sets are implemented as an N x 2 Numeric array containing
    x and y coordinates.


    To instantiate:
       Mesh(coordinates, triangles)

    where

      coordinates is either a list of 2-tuples or an Mx2 Numeric array of
      floats representing all x, y coordinates in the mesh.

      triangles is either a list of 3-tuples or an Nx3 Numeric array of
      integers representing indices of all vertices in the mesh.
      Each vertex is identified by its index i in [0, M-1].


    Example:
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        e = [2.0, 2.0]

        points = [a, b, c, e]
        triangles = [ [1,0,2], [1,2,3] ]   #bac, bce
        mesh = Mesh(points, triangles)

        #creates two triangles: bac and bce

    Other:

      In addition mesh computes an Nx6 array called vertex_coordinates.
      This structure is derived from coordinates and contains for each
      triangle the three x,y coordinates at the vertices.


        This is a cut down version of mesh from pyvolution mesh.py
    """

    #FIXME: It would be a good idea to use geospatial data as an alternative
    #input
    def __init__(self, coordinates, triangles,
                 geo_reference=None,
                 verbose=False):
        """
        Build triangles from x,y coordinates (sequence of 2-tuples or
        Mx2 Numeric array of floats) and triangles (sequence of 3-tuples
        or Nx3 Numeric array of non-negative integers).

        origin is a 3-tuple consisting of UTM zone, easting and northing.
        If specified coordinates are assumed to be relative to this origin.
        """

        if verbose: print 'General_mesh: Building basic mesh structure' 

        self.triangles = array(triangles,Int)
        self.coordinates = array(coordinates,Float)
        

        # FIXME: this stores a geo_reference, but when coords are returned
        # This geo_ref is not taken into account!
        if geo_reference is None:
            self.geo_reference = Geo_reference() #Use defaults 
        else:
            self.geo_reference = geo_reference

        #Input checks
        msg = 'Triangles must an Nx3 Numeric array or a sequence of 3-tuples'
        assert len(self.triangles.shape) == 2, msg

        msg = 'Coordinates must an Mx2 Numeric array or a sequence of 2-tuples'
        assert len(self.coordinates.shape) == 2, msg

        msg = 'Vertex indices reference non-existing coordinate sets'
        assert max(max(self.triangles)) <= self.coordinates.shape[0], msg


        #Register number of elements (N)
        self.number_of_elements = N = self.triangles.shape[0]

        # FIXME: Maybe move to statistics?
        # Or use with get_extent
        xy_extent = [ min(self.coordinates[:,0]), min(self.coordinates[:,1]) ,
                      max(self.coordinates[:,0]), max(self.coordinates[:,1]) ]

        self.xy_extent = array(xy_extent, Float)


        #Allocate space for geometric quantities
        self.normals = zeros((N, 6), Float)
        self.areas = zeros(N, Float)
        self.edgelengths = zeros((N, 3), Float)

        #Get x,y coordinates for all triangles and store
        self.vertex_coordinates = V = self.compute_vertex_coordinates()


        #Initialise each triangle
        if verbose:
            print 'General_mesh: Computing areas, normals and edgelenghts'
            
        for i in range(N):
            if verbose and i % ((N+10)/10) == 0: print '(%d/%d)' %(i, N)
           

            x0 = V[i, 0]; y0 = V[i, 1]
            x1 = V[i, 2]; y1 = V[i, 3]
            x2 = V[i, 4]; y2 = V[i, 5]

            #Area
            self.areas[i] = abs((x1*y0-x0*y1)+(x2*y1-x1*y2)+(x0*y2-x2*y0))/2

            msg = 'Triangle (%f,%f), (%f,%f), (%f, %f)' %(x0,y0,x1,y1,x2,y2)
            msg += ' is degenerate:  area == %f' %self.areas[i]
            assert self.areas[i] > 0.0, msg


            #Normals
            #The normal vectors
            # - point outward from each edge
            # - are orthogonal to the edge
            # - have unit length
            # - Are enumerated according to the opposite corner:
            #   (First normal is associated with the edge opposite
            #    the first vertex, etc)
            # - Stored as six floats n0x,n0y,n1x,n1y,n2x,n2y per triangle

            n0 = array([x2 - x1, y2 - y1])
            l0 = sqrt(sum(n0**2))

            n1 = array([x0 - x2, y0 - y2])
            l1 = sqrt(sum(n1**2))

            n2 = array([x1 - x0, y1 - y0])
            l2 = sqrt(sum(n2**2))

            #Normalise
            n0 /= l0
            n1 /= l1
            n2 /= l2

            #Compute and store
            self.normals[i, :] = [n0[1], -n0[0],
                                  n1[1], -n1[0],
                                  n2[1], -n2[0]]

            #Edgelengths
            self.edgelengths[i, :] = [l0, l1, l2]

        
        #Build vertex list
        if verbose: print 'Building vertex list'         
        self.build_vertexlist()

            

    def __len__(self):
        return self.number_of_elements

    def __repr__(self):
        return 'Mesh: %d triangles, %d elements'\
               %(self.coordinates.shape[0], len(self))

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



    def get_vertex_coordinates(self, obj=False, absolute=False):
        """Return all vertex coordinates.
        Return all vertex coordinates for all triangles as an Nx6 array
        (ordered as x0, y0, x1, y1, x2, y2 for each triangle)

        if obj is True, the x/y pairs are returned in a 3*N x 2 array.
        FIXME, we might make that the default.
        FIXME Maybe use keyword: continuous = False for this condition?
	FIXME - Maybe use something referring to unique vertices?

        Boolean keyword argument absolute determines whether coordinates
        are to be made absolute by taking georeference into account
        Default is False as many parts of ANUGA expects relative coordinates.
        (To see which, switch to default absolute=True and run tests).
        """

        V = self.vertex_coordinates
        if absolute is True:
            if not self.geo_reference.is_absolute():
            
                V0 = self.geo_reference.get_absolute(V[:,0:2])
                V1 = self.geo_reference.get_absolute(V[:,2:4])
                V2 = self.geo_reference.get_absolute(V[:,4:6])

                # This does double the memory need 
                V = concatenate( (V0, V1, V2), axis=1 )

                
        if obj is True:

            N = V.shape[0]
            return reshape(V, (3*N, 2))
        else:    
            return V


    def get_vertex_coordinate(self, i, j, absolute=False):
        """Return coordinates for vertex j of the i'th triangle.
        Return value is the numeric array slice [x, y]
        """

        V = self.get_vertex_coordinates(absolute=absolute)
        return V[i, 2*j:2*j+2]
    
        ##return self.vertex_coordinates[i, 2*j:2*j+2]


    def compute_vertex_coordinates(self):
        """Return vertex coordinates for all triangles as an Nx6 array
        (ordered as x0, y0, x1, y1, x2, y2 for each triangle)
        """

        #FIXME (Ole): Perhaps they should be ordered as in obj files??
        #See quantity.get_vertex_values
        #FIXME (Ole) - oh yes they should

        N = self.number_of_elements
        vertex_coordinates = zeros((N, 6), Float)

        for i in range(N):
            for j in range(3):
                k = self.triangles[i,j]  #Index of vertex 0
                v_k = self.coordinates[k]
                vertex_coordinates[i, 2*j+0] = v_k[0]
                vertex_coordinates[i, 2*j+1] = v_k[1]

        return vertex_coordinates

    def get_vertices(self, indices=None):
        """Get connectivity
        indices is the set of element ids of interest
        """

        if (indices ==  None):
            indices = range(len(self))  #len(self)=number of elements

        return  take(self.triangles, indices)

    #FIXME - merge these two (get_vertices and get_triangles)
    def get_triangles(self, obj=False):
        """Get connetivity
        Return triangles (triplets of indices into point coordinates)
        
        If obj is True return structure commensurate with replicated
        points, allowing for discontinuities
        (FIXME: Need good name for this concept)
        """

        if obj is True:
            m = len(self)  #Number of triangles
            M = 3*m        #Total number of unique vertices
            T = reshape(array(range(M)).astype(Int), (m,3))
        else:
            T = self.triangles

        return T     

    

    def get_unique_vertices(self,  indices=None):
        triangles = self.get_vertices(indices=indices)
        unique_verts = {}
        for triangle in triangles:
           unique_verts[triangle[0]] = 0
           unique_verts[triangle[1]] = 0
           unique_verts[triangle[2]] = 0
        return unique_verts.keys()

    def build_vertexlist(self):
        """Build vertexlist index by vertex ids and for each entry (point id)
        build a list of (triangles, vertex_id) pairs that use the point
        as vertex.

        Preconditions:
          self.coordinates and self.triangles are defined

        Postcondition:
          self.vertexlist is built
        """

        vertexlist = [None]*len(self.coordinates)
        for i in range(self.number_of_elements):

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

        
