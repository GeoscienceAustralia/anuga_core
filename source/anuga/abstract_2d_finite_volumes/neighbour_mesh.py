"""Classes implementing general 2D triangular mesh with neighbour structure.

   This structure is purely geometrical. Anything relating to quantities
   or timestepping is implemented in subclass domain.py.

   Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
   Geoscience Australia, 2004
"""

from general_mesh import General_mesh
from math import pi, sqrt
        

class Mesh(General_mesh):
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


    Mesh takes the optional third argument boundary which is a
    dictionary mapping from (element_id, edge_id) to boundary tag.
    The default value is None which will assign the default_boundary_tag
    as specified in config.py to all boundary edges.
    """

    #FIXME: Maybe rename coordinates to points (as in a poly file)
    #But keep 'vertex_coordinates'

    #FIXME: Put in check for angles less than a set minimum


    def __init__(self, coordinates, triangles,
                 boundary=None,
                 tagged_elements=None,
                 geo_reference=None,
                 number_of_full_nodes=None,
                 number_of_full_triangles=None,
                 use_inscribed_circle=False,
                 verbose=False):
        """
        Build triangles from x,y coordinates (sequence of 2-tuples or
        Mx2 Numeric array of floats) and triangles (sequence of 3-tuples
        or Nx3 Numeric array of non-negative integers).
        """



        from Numeric import array, zeros, Int, Float, maximum, sqrt, sum

        General_mesh.__init__(self, coordinates, triangles,
                              number_of_full_nodes=\
                              number_of_full_nodes,
                              number_of_full_triangles=\
                              number_of_full_triangles,
                              geo_reference=geo_reference,
                              verbose=verbose)

        if verbose: print 'Initialising mesh'         

        N = len(self) #Number_of_triangles

        self.use_inscribed_circle = use_inscribed_circle

        #Allocate space for geometric quantities
        self.centroid_coordinates = zeros((N, 2), Float)

        self.radii = zeros(N, Float)

        self.neighbours = zeros((N, 3), Int)
        self.neighbour_edges = zeros((N, 3), Int)
        self.number_of_boundaries = zeros(N, Int)
        self.surrogate_neighbours = zeros((N, 3), Int)

        #Get x,y coordinates for all triangles and store
        V = self.vertex_coordinates # Relative coordinates

        #Initialise each triangle
        if verbose: print 'Mesh: Computing centroids and radii'        
        for i in range(N):
            if verbose and i % ((N+10)/10) == 0: print '(%d/%d)' %(i, N)

            x0, y0 = V[3*i, :]
            x1, y1 = V[3*i+1, :]
            x2, y2 = V[3*i+2, :]                        

            #x0 = V[i, 0]; y0 = V[i, 1]
            #x1 = V[i, 2]; y1 = V[i, 3]
            #x2 = V[i, 4]; y2 = V[i, 5]

            #Compute centroid
            centroid = array([(x0 + x1 + x2)/3, (y0 + y1 + y2)/3])
            self.centroid_coordinates[i] = centroid


            if self.use_inscribed_circle == False:
                #OLD code. Computed radii may exceed that of an
                #inscribed circle

                #Midpoints
                m0 = array([(x1 + x2)/2, (y1 + y2)/2])
                m1 = array([(x0 + x2)/2, (y0 + y2)/2])
                m2 = array([(x1 + x0)/2, (y1 + y0)/2])

                #The radius is the distance from the centroid of
                #a triangle to the midpoint of the side of the triangle
                #closest to the centroid
                d0 = sqrt(sum( (centroid-m0)**2 ))
                d1 = sqrt(sum( (centroid-m1)**2 ))
                d2 = sqrt(sum( (centroid-m2)**2 ))

                self.radii[i] = min(d0, d1, d2)

            else:
                #NEW code added by Peter Row. True radius
                #of inscribed circle is computed

                a = sqrt((x0-x1)**2+(y0-y1)**2)
                b = sqrt((x1-x2)**2+(y1-y2)**2)
                c = sqrt((x2-x0)**2+(y2-y0)**2)

                self.radii[i]=2.0*self.areas[i]/(a+b+c)


            #Initialise Neighbours (-1 means that it is a boundary neighbour)
            self.neighbours[i, :] = [-1, -1, -1]

            #Initialise edge ids of neighbours
            #In case of boundaries this slot is not used
            self.neighbour_edges[i, :] = [-1, -1, -1]


        #Build neighbour structure
        if verbose: print 'Mesh: Building neigbour structure'                
        self.build_neighbour_structure()

        #Build surrogate neighbour structure
        if verbose: print 'Mesh: Building surrogate neigbour structure'
        self.build_surrogate_neighbour_structure()

        #Build boundary dictionary mapping (id, edge) to symbolic tags
        if verbose: print 'Mesh: Building boundary dictionary'
        self.build_boundary_dictionary(boundary)

        #Build tagged element  dictionary mapping (tag) to array of elements
        if verbose: print 'Mesh: Building tagged elements dictionary'        
        self.build_tagged_elements_dictionary(tagged_elements)

        #Update boundary indices FIXME: OBSOLETE
        #self.build_boundary_structure()

        #FIXME check integrity?
        if verbose: print 'Mesh: Done'                

    def __repr__(self):
        return General_mesh.__repr__(self) + ', %d boundary segments'\
               %(len(self.boundary))


    def set_to_inscribed_circle(self,safety_factor = 1):
        #FIXME phase out eventually
        N = self.number_of_triangles
        V = self.vertex_coordinates

        #initialising min and max ratio
        i=0
        old_rad = self.radii[i]
        x0 = V[i, 0]; y0 = V[i, 1]
        x1 = V[i, 2]; y1 = V[i, 3]
        x2 = V[i, 4]; y2 = V[i, 5]
        a = sqrt((x0-x1)**2+(y0-y1)**2)
        b = sqrt((x1-x2)**2+(y1-y2)**2)
        c = sqrt((x2-x0)**2+(y2-y0)**2)
        ratio = old_rad/self.radii[i]
        max_ratio = ratio
        min_ratio = ratio

        for i in range(N):
            old_rad = self.radii[i]
            x0 = V[i, 0]; y0 = V[i, 1]
            x1 = V[i, 2]; y1 = V[i, 3]
            x2 = V[i, 4]; y2 = V[i, 5]
            a = sqrt((x0-x1)**2+(y0-y1)**2)
            b = sqrt((x1-x2)**2+(y1-y2)**2)
            c = sqrt((x2-x0)**2+(y2-y0)**2)
            self.radii[i]=self.areas[i]/(2*(a+b+c))*safety_factor
            ratio = old_rad/self.radii[i]
            if ratio >= max_ratio: max_ratio = ratio
            if ratio <= min_ratio: min_ratio = ratio
        return max_ratio,min_ratio

    def build_neighbour_structure(self):
        """Update all registered triangles to point to their neighbours.

        Also, keep a tally of the number of boundaries for each triangle

        Postconditions:
          neighbours and neighbour_edges is populated
          number_of_boundaries integer array is defined.
        """

        #Step 1:
        #Build dictionary mapping from segments (2-tuple of points)
        #to left hand side edge (facing neighbouring triangle)

        N = len(self) #Number_of_triangles
        neighbourdict = {}
        for i in range(N):

            #Register all segments as keys mapping to current triangle
            #and segment id
            a = self.triangles[i, 0]
            b = self.triangles[i, 1]
            c = self.triangles[i, 2]
            if neighbourdict.has_key((a,b)):
                    msg = "Edge 2 of triangle %d is duplicating edge %d of triangle %d.\n" %(i,neighbourdict[a,b][1],neighbourdict[a,b][0])
                    raise Exception, msg
            if neighbourdict.has_key((b,c)):
                    msg = "Edge 0 of triangle %d is duplicating edge %d of triangle %d.\n" %(i,neighbourdict[b,c][1],neighbourdict[b,c][0])
                    raise msg
            if neighbourdict.has_key((c,a)):
                    msg = "Edge 1 of triangle %d is duplicating edge %d of triangle %d.\n" %(i,neighbourdict[c,a][1],neighbourdict[c,a][0])
                    raise msg

            neighbourdict[a,b] = (i, 2) #(id, edge)
            neighbourdict[b,c] = (i, 0) #(id, edge)
            neighbourdict[c,a] = (i, 1) #(id, edge)


        #Step 2:
        #Go through triangles again, but this time
        #reverse direction of segments and lookup neighbours.
        for i in range(N):
            a = self.triangles[i, 0]
            b = self.triangles[i, 1]
            c = self.triangles[i, 2]

            self.number_of_boundaries[i] = 3
            if neighbourdict.has_key((b,a)):
                self.neighbours[i, 2] = neighbourdict[b,a][0]
                self.neighbour_edges[i, 2] = neighbourdict[b,a][1]
                self.number_of_boundaries[i] -= 1

            if neighbourdict.has_key((c,b)):
                self.neighbours[i, 0] = neighbourdict[c,b][0]
                self.neighbour_edges[i, 0] = neighbourdict[c,b][1]
                self.number_of_boundaries[i] -= 1

            if neighbourdict.has_key((a,c)):
                self.neighbours[i, 1] = neighbourdict[a,c][0]
                self.neighbour_edges[i, 1] = neighbourdict[a,c][1]
                self.number_of_boundaries[i] -= 1


    def build_surrogate_neighbour_structure(self):
        """Build structure where each triangle edge points to its neighbours
        if they exist. Otherwise point to the triangle itself.

        The surrogate neighbour structure is useful for computing gradients
        based on centroid values of neighbours.

        Precondition: Neighbour structure is defined
        Postcondition:
          Surrogate neighbour structure is defined:
          surrogate_neighbours: i0, i1, i2 where all i_k >= 0 point to
          triangles.

        """

        N = len(self) #Number of triangles
        for i in range(N):
            #Find all neighbouring volumes that are not boundaries
            for k in range(3):
                if self.neighbours[i, k] < 0:
                    self.surrogate_neighbours[i, k] = i #Point this triangle
                else:
                    self.surrogate_neighbours[i, k] = self.neighbours[i, k]



    def build_boundary_dictionary(self, boundary = None):
        """Build or check the dictionary of boundary tags.
         self.boundary is a dictionary of tags,
         keyed by volume id and edge:
         { (id, edge): tag, ... }

         Postconditions:
            self.boundary is defined.
        """

        from anuga.config import default_boundary_tag

        if boundary is None:
            boundary = {}
            for vol_id in range(len(self)):
                for edge_id in range(0, 3):
                    if self.neighbours[vol_id, edge_id] < 0:
                        boundary[(vol_id, edge_id)] = default_boundary_tag
        else:
            #Check that all keys in given boundary exist
            for vol_id, edge_id in boundary.keys():
                msg = 'Segment (%d, %d) does not exist' %(vol_id, edge_id)
                a, b = self.neighbours.shape
                assert vol_id < a and edge_id < b, msg

                #FIXME: This assert violates internal boundaries (delete it)
                #msg = 'Segment (%d, %d) is not a boundary' %(vol_id, edge_id)
                #assert self.neighbours[vol_id, edge_id] < 0, msg

            #Check that all boundary segments are assigned a tag
            for vol_id in range(len(self)):
                for edge_id in range(0, 3):
                    if self.neighbours[vol_id, edge_id] < 0:
                        if not boundary.has_key( (vol_id, edge_id) ):
                            msg = 'WARNING: Given boundary does not contain '
                            msg += 'tags for edge (%d, %d). '\
                                   %(vol_id, edge_id)
                            msg += 'Assigning default tag (%s).'\
                                   %default_boundary_tag

                            #FIXME: Print only as per verbosity
                            #print msg

                            #FIXME: Make this situation an error in the future
                            #and make another function which will
                            #enable default boundary-tags where
                            #tags a not specified
                            boundary[ (vol_id, edge_id) ] =\
                                      default_boundary_tag



        self.boundary = boundary


    def build_tagged_elements_dictionary(self, tagged_elements = None):
        """Build the dictionary of element tags.
         self.tagged_elements is a dictionary of element arrays,
         keyed by tag:
         { (tag): [e1, e2, e3..] }

         Postconditions:
            self.element_tag is defined
        """
        from Numeric import array, Int

        if tagged_elements is None:
            tagged_elements = {}
        else:
            #Check that all keys in given boundary exist
            for tag in tagged_elements.keys():
                tagged_elements[tag] = array(tagged_elements[tag]).astype(Int)

                msg = 'Not all elements exist. '
                assert max(tagged_elements[tag]) < len(self), msg
        #print "tagged_elements", tagged_elements
        self.tagged_elements = tagged_elements

    def build_boundary_structure(self):
        """Traverse boundary and
        enumerate neighbour indices from -1 and
        counting down.

        Precondition:
            self.boundary is defined.
        Post condition:
            neighbour array has unique negative indices for boundary
            boundary_segments array imposes an ordering on segments
            (not otherwise available from the dictionary)

        Note: If a segment is listed in the boundary dictionary
        it *will* become a boundary - even if there is a neighbouring triangle.
        This would be the case for internal boundaries
        """

        #FIXME: Now Obsolete - maybe use some comments from here in
        #domain.set_boundary

        if self.boundary is None:
            msg = 'Boundary dictionary must be defined before '
            msg += 'building boundary structure'
            raise msg


        self.boundary_segments = self.boundary.keys()
        self.boundary_segments.sort()

        index = -1
        for id, edge in self.boundary_segments:

            #FIXME: One would detect internal boundaries as follows
            #if self.neighbours[id, edge] > -1:
            #    print 'Internal boundary'

            self.neighbours[id, edge] = index
            index -= 1


    def get_boundary_tags(self):
        """Return list of available boundary tags
        """

        tags = {}
        for v in self.boundary.values():
            tags[v] = 1

        return tags.keys()


    def get_boundary_polygon(self, verbose=False):
        """Return bounding polygon for mesh (counter clockwise)

        Using the mesh boundary, derive a bounding polygon for this mesh.
        If multiple vertex values are present (vertices stored uniquely), 
	the algorithm will select the path that contains the entire mesh.

        All points are in absolute UTM coordinates
        """
        
        from Numeric import allclose, sqrt, array, minimum, maximum
        from anuga.utilities.numerical_tools import angle, ensure_numeric     


        # Get mesh extent
        xmin, xmax, ymin, ymax = self.get_extent(absolute=True)
        pmin = ensure_numeric([xmin, ymin])
        pmax = ensure_numeric([xmax, ymax])        


        # Assemble dictionary of boundary segments and choose starting point
        segments = {}
        inverse_segments = {}
        p0 = None
        mindist = sqrt(sum((pmax-pmin)**2)) # Start value across entire mesh
        for i, edge_id in self.boundary.keys():
            # Find vertex ids for boundary segment
            if edge_id == 0: a = 1; b = 2
            if edge_id == 1: a = 2; b = 0
            if edge_id == 2: a = 0; b = 1

            A = self.get_vertex_coordinate(i, a, absolute=True) # Start
            B = self.get_vertex_coordinate(i, b, absolute=True) # End


            # Take the point closest to pmin as starting point
            # Note: Could be arbitrary, but nice to have
            # a unique way of selecting
            dist_A = sqrt(sum((A-pmin)**2))
            dist_B = sqrt(sum((B-pmin)**2))

            # Find lower leftmost point
            if dist_A < mindist:
                mindist = dist_A
                p0 = A
            if dist_B < mindist:
                mindist = dist_B
                p0 = B


            # Sanity check
            if p0 is None:
                raise Exception('Impossible')


            # Register potential paths from A to B
            if not segments.has_key(tuple(A)):
                segments[tuple(A)] = [] # Empty list for candidate points

            segments[tuple(A)].append(B)                


        # Start with smallest point and follow boundary (counter clock wise)
        polygon = [list(p0)]# Storage for final boundary polygon
        point_registry = {} # Keep track of storage to avoid multiple runs 
	                    # around boundary. This will only be the case if 
			    # there are more than one candidate.
                            # FIXME (Ole): Perhaps we can do away with polygon
                            # and use only point_registry to save space.

        point_registry[tuple(p0)] = 0                            
                            
        while len(point_registry) < len(self.boundary):

            candidate_list = segments[tuple(p0)]
            if len(candidate_list) > 1:
                # Multiple points detected (this will be the case for meshes 
		# with duplicate points as those used for discontinuous 
		# triangles with vertices stored uniquely).
                # Take the candidate that is furthest to the clockwise 
		# direction, as that will follow the boundary.
		#
		# This will also be the case for pathological triangles
		# that have no neighbours.

                if verbose:
                    print 'Point %s has multiple candidates: %s'\
                          %(str(p0), candidate_list)

                # Check that previous are not in candidate list
                #for p in candidate_list:
                #    assert not allclose(p0, p)

                # Choose vector against which all angles will be measured
                if len(polygon) > 1:    
                    v_prev = p0 - polygon[-2] # Vector that leads to p0 
		                              # from previous point
                else:
                    # FIXME (Ole): What do we do if the first point has 
		    # multiple candidates?
                    # Being the lower left corner, perhaps we can use the
                    # vector [1, 0], but I really don't know if this is 
		    # completely watertight.
                    v_prev = [1.0, 0.0]
                    

                # Choose candidate with minimum angle    
                minimum_angle = 2*pi
                for pc in candidate_list:

                    vc = pc-p0  # Candidate vector (from p0 to candidate pt)
                   
                    # Angle between each candidate and the previous vector
                    # in [-pi, pi]
                    ac = angle(vc, v_prev)
                    if ac > pi:
                        # Give preference to angles on the right hand side 
			# of v_prev 
                        # print 'pc = %s, changing angle from %f to %f'\
			# %(pc, ac*180/pi, (ac-2*pi)*180/pi)
                        ac = ac-2*pi

                    # Take the minimal angle corresponding to the 
		    # rightmost vector
                    if ac < minimum_angle:
                        minimum_angle = ac
                        p1 = pc             # Best candidate 
                        

                if verbose is True:
                    print '  Best candidate %s, angle %f'\
                          %(p1, minimum_angle*180/pi)
                    
            else:
                p1 = candidate_list[0]

		
            if point_registry.has_key(tuple(p1)):
                # We have reached a point already visited. 
		
		if allclose(p1, polygon[0]):
		    # If it is the initial point, the polygon is complete. 
		    
                    if verbose is True:
		        print '  Stop criterion fulfilled at point %s' %p1
		        print polygon		    
			
                    # We have completed the boundary polygon - yeehaa
		    break
		else:    
		    # The point already visited is not the initial point
		    # This would be a pathological triangle, but the 
		    # algorithm must be able to deal with this
		    pass
   
            else:
	        # We are still finding new points on the boundary
                point_registry[tuple(p1)] = len(point_registry)
            
            polygon.append(list(p1)) # De-Numeric each point :-)
            p0 = p1


        return polygon


    def check_integrity(self):
        """Check that triangles are internally consistent e.g.
        that area corresponds to edgelengths, that vertices
        are arranged in a counter-clockwise order, etc etc
        Neighbour structure will be checked by class Mesh
        """

        from anuga.config import epsilon
        from anuga.utilities.numerical_tools import anglediff

        from Numeric import sort, allclose

        N = len(self)
        #Get x,y coordinates for all vertices for all triangles
        V = self.get_vertex_coordinates()
        #Check each triangle
        for i in range(N):

            x0, y0 = V[3*i, :]
            x1, y1 = V[3*i+1, :]
            x2, y2 = V[3*i+2, :]
            
            #Check that area hasn't been compromised
            area = self.areas[i]
            ref = abs((x1*y0-x0*y1)+(x2*y1-x1*y2)+(x0*y2-x2*y0))/2
            msg = 'Wrong area for vertex coordinates: %f %f %f %f %f %f'\
                  %(x0,y0,x1,y1,x2,y2)
            assert abs((area - ref)/area) < epsilon, msg

            #Check that points are arranged in counter clock-wise order
            v0 = [x1-x0, y1-y0]
            v1 = [x2-x1, y2-y1]
            v2 = [x0-x2, y0-y2]
            a0 = anglediff(v1, v0)
            a1 = anglediff(v2, v1)
            a2 = anglediff(v0, v2)

            msg = '''Vertices (%s,%s), (%s,%s), (%s,%s) are not arranged
            in counter clockwise order''' %(x0, y0, x1, y1, x2, y2)
            assert a0 < pi and a1 < pi and a2 < pi, msg

            #Check that normals are orthogonal to edge vectors
            #Note that normal[k] lies opposite vertex k

            normal0 = self.normals[i, 0:2]
            normal1 = self.normals[i, 2:4]
            normal2 = self.normals[i, 4:6]

            for u, v in [ (v0, normal2), (v1, normal0), (v2, normal1) ]:

                #Normalise
                l_u = sqrt(u[0]*u[0] + u[1]*u[1])
                l_v = sqrt(v[0]*v[0] + v[1]*v[1])                

                x = (u[0]*v[0] + u[1]*v[1])/l_u/l_v #Inner product
                
                msg = 'Normal vector (%f,%f) is not perpendicular to' %tuple(v)
                msg += ' edge (%f,%f) in triangle %d.' %(tuple(u) + (i,))
                msg += ' Inner product is %e.' %x                
                assert x < epsilon, msg

        self.lone_vertices = []
        #Check that all vertices have been registered
        for node, count in enumerate(self.number_of_triangles_per_node):
        
            #msg = 'Node %d does not belong to an element.' %node
            #assert count > 0, msg
            if count == 0:
                self.lone_vertices.append(node)



        #Check neighbour structure
        for i in range(N):
            # For each triangle
            
            for k, neighbour_id in enumerate(self.neighbours[i,:]):

                #Assert that my neighbour's neighbour is me
                #Boundaries need not fulfill this
                if neighbour_id >= 0:
                    edge = self.neighbour_edges[i, k]
                    msg = 'Triangle %d has neighbour %d but it does not point back. \n' %(i,neighbour_id)
                    msg += 'Only points to (%s)' %(self.neighbours[neighbour_id,:])
                    assert self.neighbours[neighbour_id, edge] == i ,msg



        #Check that all boundaries have
        # unique, consecutive, negative indices

        #L = len(self.boundary)
        #for i in range(L):
        #    id, edge = self.boundary_segments[i]
        #    assert self.neighbours[id, edge] == -i-1


        #NOTE: This assert doesn't hold true if there are internal boundaries
        #FIXME: Look into this further.
        #FIXME (Ole): In pyvolution mark 3 this is OK again
        #NOTE: No longer works because neighbour structure is modified by
        #      domain set_boundary.
        #for id, edge in self.boundary:
        #    assert self.neighbours[id,edge] < 0
        #
        #NOTE (Ole): I reckon this was resolved late 2004?
        #
        #See domain.set_boundary



        #Check integrity of inverted triangle structure

        V = self.vertex_value_indices[:] #Take a copy
        V = sort(V)
        assert allclose(V, range(3*N))

        assert sum(self.number_of_triangles_per_node) ==\
               len(self.vertex_value_indices)

        # Check number of triangles per node
        count = [0]*self.number_of_nodes
        for triangle in self.triangles:
            for i in triangle:
                count[i] += 1

        assert allclose(count, self.number_of_triangles_per_node)


        # Check integrity of vertex_value_indices
        current_node = 0
        k = 0 # Track triangles touching on node
        for index in self.vertex_value_indices:

            if self.number_of_triangles_per_node[current_node] == 0:
                # Node is lone - i.e. not part of the mesh
                continue
            
            k += 1
            
            volume_id = index / 3
            vertex_id = index % 3
            
            msg = 'Triangle %d, vertex %d points to %d. Should have been %d'\
                  %(volume_id, vertex_id, self.triangles[volume_id, vertex_id], current_node)
            assert self.triangles[volume_id, vertex_id] == current_node, msg
                        
            if self.number_of_triangles_per_node[current_node] == k:
                # Move on to next node
                k = 0
                current_node += 1


    def get_lone_vertices(self):
        """Return a list of vertices that are not connected to any triangles.

        Precondition
        FIXME(DSG - DSG) Pull the code out of check integrity that builds this
                         structure.
        check_integrity has to have been called. 
        """
        return self.lone_vertices

    def get_centroid_coordinates(self, absolute=False):
        """Return all centroid coordinates.
        Return all centroid coordinates for all triangles as an Nx2 array
        (ordered as x0, y0 for each triangle)

        Boolean keyword argument absolute determines whether coordinates
        are to be made absolute by taking georeference into account
        Default is False as many parts of ANUGA expects relative coordinates.
        """

        V = self.centroid_coordinates
        if absolute is True:
            if not self.geo_reference.is_absolute():
                V = self.geo_reference.get_absolute(V)
            
        return V

        
    def get_radii(self):
        """Return all radii.
        Return radius of inscribed cirle for all triangles
        """
        return self.radii	



    def statistics(self):
        """Output statistics about mesh
        """

        from Numeric import arange
        from anuga.utilities.numerical_tools import histogram, create_bins

        vertex_coordinates = self.vertex_coordinates # Relative coordinates
        areas = self.areas
        x = vertex_coordinates[:,0]
        y = vertex_coordinates[:,1]


        #Setup 10 bins for area histogram
        bins = create_bins(areas, 10)
        #m = max(areas)
        #bins = arange(0., m, m/10)
        hist = histogram(areas, bins)

        str =  '------------------------------------------------\n'
        str += 'Mesh statistics:\n'
        str += '  Number of triangles = %d\n' %len(self)
        str += '  Extent [m]:\n'
        str += '    x in [%f, %f]\n' %(min(x), max(x))
        str += '    y in [%f, %f]\n' %(min(y), max(y))
        str += '  Areas [m^2]:\n'
        str += '    A in [%f, %f]\n' %(min(areas), max(areas))
        str += '    number of distinct areas: %d\n' %(len(areas))        
        str += '    Histogram:\n'

        hi = bins[0]
        for i, count in enumerate(hist):
            lo = hi
            if i+1 < len(bins):
                #Open upper interval                
                hi = bins[i+1]
                str += '      [%f, %f[: %d\n' %(lo, hi, count)                
            else:
                #Closed upper interval
                hi = max(areas)
                str += '      [%f, %f]: %d\n' %(lo, hi, count)

        N = len(areas)
        if N > 10:
            str += '    Percentiles (10%):\n'
            areas = areas.tolist()
            areas.sort()

            k = 0
            lower = min(areas)
            for i, a in enumerate(areas):        
                if i % (N/10) == 0 and i != 0: #For every 10% of the sorted areas                
                    str += '      %d triangles in [%f, %f]\n' %(i-k, lower, a)
                    lower = a
                    k = i
                    
            str += '      %d triangles in [%f, %f]\n'\
                   %(N-k, lower, max(areas))                    
                
                      
        str += '  Boundary:\n'
        str += '    Number of boundary segments == %d\n' %(len(self.boundary))
        str += '    Boundary tags == %s\n' %self.get_boundary_tags()  
        str += '------------------------------------------------\n'
        

        return str

