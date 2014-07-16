"""Classes implementing general 2D triangular mesh with neighbour structure.

   This structure is purely geometrical. Anything relating to quantities
   or timestepping is implemented in subclass domain.py.

   Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
   Geoscience Australia, 2004
"""

from general_mesh import General_mesh
from anuga.caching import cache
import anuga.utilities.log as log

from math import pi, sqrt

import numpy as num


class Mesh(General_mesh):
    """Collection of triangular elements (purely geometric)

    A triangular element is defined in terms of three vertex ids,
    ordered counter clock-wise,
    each corresponding to a given coordinate set.
    Vertices from different elements can point to the same
    coordinate set.

    Coordinate sets are implemented as an N x 2 numeric array containing
    x and y coordinates.


    To instantiate:
       Mesh(coordinates, triangles)

    where

      coordinates is either a list of 2-tuples or an Mx2 numeric array of
      floats representing all x, y coordinates in the mesh.

      triangles is either a list of 3-tuples or an Nx3 numeric array of
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
                 use_inscribed_circle=False,
                 verbose=False):
        """
        Build Mesh

            Input x,y coordinates (sequence of 2-tuples or Mx2 numeric array of floats)
            triangles (sequence of 3-tuples or Nx3 numeric array of non-negative integers).
        """




        General_mesh.__init__(self, coordinates, triangles,
                              geo_reference=geo_reference,
                              use_inscribed_circle=use_inscribed_circle,
                              verbose=verbose)

        if verbose: log.critical('Mesh: Initialising')

        N = len(self) #Number_of_triangles

        # Allocate arrays for neighbour data

        self.neighbours = -1*num.ones((N, 3), num.int)
        self.neighbour_edges = -1*num.ones((N, 3), num.int)
        self.number_of_boundaries = num.zeros(N, num.int)
        self.surrogate_neighbours = num.zeros((N, 3), num.int)

        #Get x,y coordinates for all triangles and store
        V = self.vertex_coordinates # Relative coordinates

#        #Initialise each triangle
#        if verbose: log.critical('Mesh: Computing centroids and radii')
#        for i in range(N):
#            if verbose and i % ((N+10)/10) == 0: log.critical('(%d/%d)' % (i, N))
#
#            x0, y0 = V[3*i, :]
#            x1, y1 = V[3*i+1, :]
#            x2, y2 = V[3*i+2, :]
#
#            #x0 = V[i, 0]; y0 = V[i, 1]
#            #x1 = V[i, 2]; y1 = V[i, 3]
#            #x2 = V[i, 4]; y2 = V[i, 5]
#
#            #Compute centroid
#            centroid = num.array([(x0 + x1 + x2)/3, (y0 + y1 + y2)/3], num.float)
#            self.centroid_coordinates[i] = centroid
#
#
#            if self.use_inscribed_circle == False:
#                #OLD code. Computed radii may exceed that of an
#                #inscribed circle
#
#                #Midpoints
#                m0 = num.array([(x1 + x2)/2, (y1 + y2)/2], num.float)
#                m1 = num.array([(x0 + x2)/2, (y0 + y2)/2], num.float)
#                m2 = num.array([(x1 + x0)/2, (y1 + y0)/2], num.float)
#
#                #The radius is the distance from the centroid of
#                #a triangle to the midpoint of the side of the triangle
#                #closest to the centroid
#                d0 = num.sqrt(num.sum( (centroid-m0)**2 ))
#                d1 = num.sqrt(num.sum( (centroid-m1)**2 ))
#                d2 = num.sqrt(num.sum( (centroid-m2)**2 ))
#
#                self.radii[i] = min(d0, d1, d2)
#
#            else:
#                #NEW code added by Peter Row. True radius
#                #of inscribed circle is computed
#
#                a = num.sqrt((x0-x1)**2+(y0-y1)**2)
#                b = num.sqrt((x1-x2)**2+(y1-y2)**2)
#                c = num.sqrt((x2-x0)**2+(y2-y0)**2)
#
#                self.radii[i]=2.0*self.areas[i]/(a+b+c)
#
#
#            #Initialise Neighbours (-1 means that it is a boundary neighbour)
#            self.neighbours[i, :] = [-1, -1, -1]
#
#            #Initialise edge ids of neighbours
#            #In case of boundaries this slot is not used
#            self.neighbour_edges[i, :] = [-1, -1, -1]


        #Build neighbour structure
        if verbose: log.critical('Mesh: Building neigbour structure')
        self.build_neighbour_structure()

        #Build surrogate neighbour structure
        if verbose: log.critical('Mesh: Building surrogate neigbour structure')
        self.build_surrogate_neighbour_structure()

        #Build boundary dictionary mapping (id, edge) to symbolic tags
        if verbose: log.critical('Mesh: Building boundary dictionary')
        self.build_boundary_dictionary(boundary)

        #Update boundary_enumeration
        self.build_boundary_neighbours()

        #Build tagged element  dictionary mapping (tag) to array of elements
        if verbose: log.critical('Mesh: Building tagged elements dictionary')
        self.build_tagged_elements_dictionary(tagged_elements)

        # Build a list of vertices that are not connected to any triangles
        self.lone_vertices = []
        #Check that all vertices have been registered
        for node, count in enumerate(self.number_of_triangles_per_node):
            #msg = 'Node %d does not belong to an element.' %node
            #assert count > 0, msg
            if count == 0:
                self.lone_vertices.append(node)

        #Update boundary indices FIXME: OBSOLETE
        #self.build_boundary_structure()



        #FIXME check integrity?
        if verbose: log.critical('Mesh: Done')
        if verbose: log.timingInfo("finishMesh, '%s'" % log.CurrentDateTime())
        if verbose: log.resource_usage_timing(log.logging.INFO, "finishMesh_")

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
        a = num.sqrt((x0-x1)**2+(y0-y1)**2)
        b = num.sqrt((x1-x2)**2+(y1-y2)**2)
        c = num.sqrt((x2-x0)**2+(y2-y0)**2)
        ratio = old_rad/self.radii[i]
        max_ratio = ratio
        min_ratio = ratio

        for i in xrange(N):
            old_rad = self.radii[i]
            x0 = V[i, 0]; y0 = V[i, 1]
            x1 = V[i, 2]; y1 = V[i, 3]
            x2 = V[i, 4]; y2 = V[i, 5]
            a = num.sqrt((x0-x1)**2+(y0-y1)**2)
            b = num.sqrt((x1-x2)**2+(y1-y2)**2)
            c = num.sqrt((x2-x0)**2+(y2-y0)**2)
            self.radii[i]=self.areas[i]/(2*(a+b+c))*safety_factor
            ratio = old_rad/self.radii[i]
            if ratio >= max_ratio: max_ratio = ratio
            if ratio <= min_ratio: min_ratio = ratio
        return max_ratio,min_ratio



    def build_neighbour_structure_python(self):
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
                    raise Exception(msg)
            if neighbourdict.has_key((b,c)):
                    msg = "Edge 0 of triangle %d is duplicating edge %d of triangle %d.\n" %(i,neighbourdict[b,c][1],neighbourdict[b,c][0])
                    raise Exception(msg)
            if neighbourdict.has_key((c,a)):
                    msg = "Edge 1 of triangle %d is duplicating edge %d of triangle %d.\n" %(i,neighbourdict[c,a][1],neighbourdict[c,a][0])
                    raise Exception(msg)

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

    def build_neighbour_structure(self):
        """Update all registered triangles to point to their neighbours.

        Also, keep a tally of the number of boundaries for each triangle

        Postconditions:
          neighbours and neighbour_edges is populated
          number_of_boundaries integer array is defined.
        """

        import neighbour_table_ext
        
        N = self.number_of_nodes

        
        neighbour_table_ext.build_neighbour_structure(N,
                                            self.triangles,
                                            self.neighbours,
                                            self.neighbour_edges,
                                            self.number_of_boundaries)
 

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
#        for i in xrange(N):
#            #Find all neighbouring volumes that are not boundaries
#            for k in xrange(3):
#                if self.neighbours[i, k] < 0:
#                    self.surrogate_neighbours[i, k] = i #Point this triangle
#                else:
#                    self.surrogate_neighbours[i, k] = self.neighbours[i, k]

        tmp_range = num.arange(N)
        for k in xrange(3):
            self.surrogate_neighbours[:,k] = \
              num.where(self.neighbours[:,k]<0, tmp_range, self.neighbours[:, k])


    def build_boundary_dictionary(self, boundary=None):
        """Build or check the dictionary of boundary tags.
        self.boundary is a dictionary of tags,
        keyed by volume id and edge:
        { (id, edge): tag, ... }

        Postconditions:
        self.boundary is defined.
        """

        from anuga.config import default_boundary_tag

        #arr_neighbours = num.array(self.neighbours)

        
        if boundary is None:
            boundary = {}

        from neighbour_mesh_ext import boundary_dictionary_construct
        boundary = boundary_dictionary_construct(len(self), default_boundary_tag, self.neighbours, boundary)
        

        self.boundary = boundary
        self.boundary_length = len(self.boundary)



    def build_boundary_dictionary_old(self, boundary = None):
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
            for vol_id in xrange(len(self)):
                for edge_id in xrange(0, 3):
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
            for vol_id in xrange(len(self)):
                for edge_id in xrange(0, 3):
                    if self.neighbours[vol_id, edge_id] < 0:
                        if not boundary.has_key( (vol_id, edge_id) ):
                            msg = 'WARNING: Given boundary does not contain '
                            msg += 'tags for edge (%d, %d). '\
                                   %(vol_id, edge_id)
                            msg += 'Assigning default tag (%s).'\
                                   %default_boundary_tag

                            #FIXME: Print only as per verbosity

                            #FIXME: Make this situation an error in the future
                            #and make another function which will
                            #enable default boundary-tags where
                            #tags a not specified
                            boundary[ (vol_id, edge_id) ] =\
                                      default_boundary_tag



        self.boundary = boundary
        self.boundary_length = len(self.boundary)


    def build_tagged_elements_dictionary(self, tagged_elements = None):
        """Build the dictionary of element tags.
         self.tagged_elements is a dictionary of element arrays,
         keyed by tag:
         { (tag): [e1, e2, e3..] }

         Postconditions:
            self.element_tag is defined
        """

        if tagged_elements is None:
            tagged_elements = {}
        else:
            #Check that all keys in given boundary exist
            for tag in tagged_elements.keys():
                tagged_elements[tag] = num.array(tagged_elements[tag], num.int)

                msg = 'Not all elements exist. '
                assert max(tagged_elements[tag]) < len(self), msg
        self.tagged_elements = tagged_elements

    def get_tagged_elements(self):
        return self.tagged_elements

#    def build_boundary_structure(self):
#        """Traverse boundary and
#        enumerate neighbour indices from -1 and
#        counting down.
#
#        Precondition:
#            self.boundary is defined.
#        Post condition:
#            neighbour array has unique negative indices for boundary
#            boundary_segments array imposes an ordering on segments
#            (not otherwise available from the dictionary)
#
#        Note: If a segment is listed in the boundary dictionary
#        it *will* become a boundary - even if there is a neighbouring triangle.
#        This would be the case for internal boundaries
#        """
#
#        #FIXME: Now Obsolete - maybe use some comments from here in
#        #domain.set_boundary
#
#        if self.boundary is None:
#            msg = 'Boundary dictionary must be defined before '
#            msg += 'building boundary structure'
#            raise Exception(msg)
#
#
#        self.boundary_segments = self.boundary.keys()
#        self.boundary_segments.sort()
#
#        index = -1
#        for id, edge in self.boundary_segments:
#
#            #FIXME: One would detect internal boundaries as follows
#            #if self.neighbours[id, edge] > -1:
#            #    log.critical('Internal boundary')
#
#            self.neighbours[id, edge] = index
#
#            self.boundary_enumeration[id,edge] = index
#
#            index -= 1
#


    def build_boundary_neighbours(self):
        """Traverse boundary and
        enumerate neighbour indices from -1 and
        counting down.

        Precondition:
            self.boundary is defined.
        Post condition:
            neighbours array has unique negative indices for boundary
            boundary_segments array imposes an ordering on segments
            (not otherwise available from the dictionary)

        """

        if self.boundary is None:
            msg = 'Boundary dictionary must be defined before '
            msg += 'building boundary structure'
            raise Exception(msg)

        self.boundary_enumeration = {}



        X = self.boundary.keys()
        X.sort()

        #print 'X', X
        index = -1
        for id, edge in X:
            self.neighbours[id, edge] = index

            self.boundary_enumeration[id,edge] = -index -1 

            index -= 1

        # Now we know number of boundaries
        M = len(self.boundary_enumeration)
        self.boundary_cells = num.zeros((M,),num.int)
        self.boundary_edges = num.zeros((M,),num.int)

        for id, edge in X:
            j = self.boundary_enumeration[id,edge]
            self.boundary_cells[j] = id
            self.boundary_edges[j] = edge

        # For each tag create list of boundary edges
        self.tag_boundary_cells = {}

        tags = self.get_boundary_tags()

        #print tags

        for tag in tags:
            self.tag_boundary_cells[tag] = []


        for j in xrange(self.boundary_length):
            id  = self.boundary_cells[j]
            edge = self.boundary_edges[j]
            tag = self.boundary[id, edge]
            #print tag, id, edge
            self.tag_boundary_cells[tag].append(j)


        #print self.tag_boundary_cells

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

        from anuga.utilities.numerical_tools import angle, ensure_numeric

        # Get mesh extent
        xmin, xmax, ymin, ymax = self.get_extent(absolute=True)
        pmin = ensure_numeric([xmin, ymin])
        pmax = ensure_numeric([xmax, ymax])

        # Assemble dictionary of boundary segments and choose starting point
        segments = {}
        inverse_segments = {}
        p0 = None

        # Start value across entire mesh
        mindist = num.sqrt(num.sum((pmax-pmin)**2))
        for i, edge_id in self.boundary.keys():
            # Find vertex ids for boundary segment
            if edge_id == 0: a = 1; b = 2
            if edge_id == 1: a = 2; b = 0
            if edge_id == 2: a = 0; b = 1

            A = self.get_vertex_coordinate(i, a, absolute=True)    # Start
            B = self.get_vertex_coordinate(i, b, absolute=True)    # End

            # Take the point closest to pmin as starting point
            # Note: Could be arbitrary, but nice to have
            # a unique way of selecting
            dist_A = num.sqrt(num.sum((A-pmin)**2))
            dist_B = num.sqrt(num.sum((B-pmin)**2))

            # Find lower leftmost point
            if dist_A < mindist:
                mindist = dist_A
                p0 = A
            if dist_B < mindist:
                mindist = dist_B
                p0 = B

            # Sanity check
            if p0 is None:
                msg = 'Impossible: p0 is None!?'
                raise Exception(msg)

            # Register potential paths from A to B
            if not segments.has_key(tuple(A)):
                segments[tuple(A)] = []    # Empty list for candidate points

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
                    log.critical('Point %s has multiple candidates: %s'
                                 % (str(p0), candidate_list))

                # Check that previous are not in candidate list
                #for p in candidate_list:
                #    assert not allclose(p0, p)

                # Choose vector against which all angles will be measured
                if len(polygon) > 1:
                    v_prev = p0 - polygon[-2]    # Vector that leads to p0
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
                    vc = pc-p0    # Candidate vector (from p0 to candidate pt)

                    # Angle between each candidate and the previous vector
                    # in [-pi, pi]
                    ac = angle(vc, v_prev)
                    if ac > pi:
                        # Give preference to angles on the right hand side
                        # of v_prev
                        ac = ac-2*pi

                    # Take the minimal angle corresponding to the
                    # rightmost vector
                    if ac < minimum_angle:
                        minimum_angle = ac
                        p1 = pc             # Best candidate

                if verbose is True:
                    log.critical('  Best candidate %s, angle %f'
                                 % (p1, minimum_angle*180/pi))
            else:
                p1 = candidate_list[0]

            if point_registry.has_key(tuple(p1)):
                # We have reached a point already visited.
                if num.allclose(p1, polygon[0]):
                    # If it is the initial point, the polygon is complete.
                    if verbose is True:
                        log.critical('  Stop criterion fulfilled at point %s'
                                     % str(p1))
                        log.critical(str(polygon))

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

            polygon.append(list(p1))    # De-numeric each point :-)
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

        N = len(self)

        # Get x,y coordinates for all vertices for all triangles
        V = self.get_vertex_coordinates()

#        # Check each triangle
#        for i in xrange(0):
#
#            x0, y0 = V[3*i, :]
#            x1, y1 = V[3*i+1, :]
#            x2, y2 = V[3*i+2, :]
#
#            # Check that area hasn't been compromised
#            area = self.areas[i]
#            ref = -((x1*y0-x0*y1)+(x2*y1-x1*y2)+(x0*y2-x2*y0))/2
#            msg = 'Triangle %i (%f,%f), (%f,%f), (%f, %f)' % (i, x0,y0,x1,y1,x2,y2)
#            msg += 'Wrong area: %f  %f'\
#                  %(area, ref)
#            assert abs((area - ref)/area) < epsilon, msg
#
#            msg = 'Triangle %i (%f,%f), (%f,%f), (%f, %f)' % (i, x0,y0,x1,y1,x2,y2)
#            msg += ' is degenerate:  area == %f' % self.areas[i]
#            assert area > 0.0, msg
#
#            # Check that points are arranged in counter clock-wise order
#            v0 = [x1-x0, y1-y0]
#            v1 = [x2-x1, y2-y1]
#            v2 = [x0-x2, y0-y2]
#            a0 = anglediff(v1, v0)
#            a1 = anglediff(v2, v1)
#            a2 = anglediff(v0, v2)
#
#            msg = '''Vertices (%s,%s), (%s,%s), (%s,%s) are not arranged
#            in counter clockwise order''' %(x0, y0, x1, y1, x2, y2)
#            assert a0 < pi and a1 < pi and a2 < pi, msg
#
#            # Check that normals are orthogonal to edge vectors
#            # Note that normal[k] lies opposite vertex k
#
#            normal0 = self.normals[i, 0:2]
#            normal1 = self.normals[i, 2:4]
#            normal2 = self.normals[i, 4:6]
#
#            for u, v in [ (v0, normal2), (v1, normal0), (v2, normal1) ]:
#
#                # Normalise
#                l_u = num.sqrt(u[0]*u[0] + u[1]*u[1])
#                l_v = num.sqrt(v[0]*v[0] + v[1]*v[1])
#
#                msg = 'Normal vector in triangle %d does not have unit length' %i
#                assert num.allclose(l_v, 1), msg
#
#                x = (u[0]*v[0] + u[1]*v[1])/l_u # Inner product
#
#                msg = 'Normal vector (%f,%f) is not perpendicular to' %tuple(v)
#                msg += ' edge (%f,%f) in triangle %d.' %(tuple(u) + (i,))
#                msg += ' Inner product is %e.' %x
#                assert x < epsilon, msg


        # let's try numpy constructs

        x0 = V[0::3, 0]
        y0 = V[0::3, 1]
        x1 = V[1::3, 0]
        y1 = V[1::3, 1]
        x2 = V[2::3, 0]
        y2 = V[2::3, 1]


        #print 'check areas'
        area = self.areas

        ref = -((x1*y0-x0*y1)+(x2*y1-x1*y2)+(x0*y2-x2*y0))/2


        assert num.sum(num.abs((area - ref)/area)) < epsilon, 'Error in areas'

        assert num.all(area > 0.0), 'A negative area'


        tx0 = x2 - x1
        ty0 = y2 - y1
        a0  = num.sqrt(tx0**2 + ty0**2)


        tx0 = tx0/a0
        ty0 = ty0/a0


        tx1 = x0 - x2
        ty1 = y0 - y2
        a1  = num.sqrt(tx1**2 + ty1**2)
        tx1 = tx1/a1
        ty1 = ty1/a1

        tx2 = x1 - x0
        ty2 = y1 - y0
        a2  = num.sqrt(tx2**2 + ty2**2)
        tx2 = tx2/a2
        ty2 = ty2/a2

        nx0 = self.normals[:,0]
        ny0 = self.normals[:,1]
        nx1 = self.normals[:,2]
        ny1 = self.normals[:,3]
        nx2 = self.normals[:,4]
        ny2 = self.normals[:,5]


        assert num.all(tx0*nx0 + ty0*ny0 < epsilon), 'Normal not perpendicular to edge'
        assert num.all(tx1*nx1 + ty1*ny1 < epsilon), 'Normal not perpendicular to edge'
        assert num.all(tx2*nx2 + ty2*ny2 < epsilon), 'Normal not perpendicular to edge'


        #print 'check normals are unit length'
        assert num.all(num.abs(nx0**2 + ny0**2 - 1) < epsilon), 'Normal are not normalised'
        assert num.all(num.abs(nx1**2 + ny1**2 - 1) < epsilon), 'Normal are not normalised'
        assert num.all(num.abs(nx2**2 + ny2**2 - 1) < epsilon), 'Normal are not normalised'




        # check that neighbour of neighbour is self

        # 0 neighbours
        neighs = self.neighbours
        ids = num.arange(len(neighs))

        # 0 neighbours
        nid = neighs[:,0]
        eid = self.neighbour_edges[:,0]
        nnid = num.argwhere(nid>-1).reshape(-1,)
        nid = nid[nnid]
        eid = eid[nnid]
        id  = ids[nnid]

        assert num.all(neighs[nid,eid] == id)

        # 1 neighbours
        nid = neighs[:,1]
        eid = self.neighbour_edges[:,1]
        nnid = num.argwhere(nid>-1).reshape(-1,)
        nid = nid[nnid]
        eid = eid[nnid]
        id  = ids[nnid]

        assert num.all(neighs[nid,eid] == id)

        # 2 neighbours
        nid = neighs[:,2]
        eid = self.neighbour_edges[:,2]
        nnid = num.argwhere(nid>-1).reshape(-1,)
        nid = nid[nnid]
        eid = eid[nnid]
        id  = ids[nnid]

        assert num.all(neighs[nid,eid] == id)




#        # Check neighbour structure
#        for i in xrange(N):
#            # For each triangle
#
#            for k, neighbour_id in enumerate(self.neighbours[i,:]):
#
#                #Assert that my neighbour's neighbour is me
#                #Boundaries need not fulfill this
#                if neighbour_id >= 0:
#                    edge = self.neighbour_edges[i, k]
#                    msg = 'Triangle %d has neighbour %d but it does not point back. \n' %(i,neighbour_id)
#                    msg += 'Only points to (%s)' %(self.neighbours[neighbour_id,:])
#                    assert self.neighbours[neighbour_id, edge] == i ,msg



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
        V = num.sort(V)
        assert num.allclose(V, range(3*N))

        assert num.sum(self.number_of_triangles_per_node) ==\
                       len(self.vertex_value_indices)

        
        # Check number of triangles per node
#        count = [0]*self.number_of_nodes
#        for triangle in self.triangles:
#            for i in triangle:
#                count[i] += 1

        count  = num.bincount(self.triangles.flat)


        ncount = len(count)
        #print len(count)
        #print len(self.number_of_triangles_per_node)


        number_of_lone_nodes = self.number_of_nodes - len(self.number_of_triangles_per_node)


        assert num.allclose(count, self.number_of_triangles_per_node[:ncount])


        from neighbour_mesh_ext import check_integrity_c


        #print self.vertex_value_indices.shape
        #print self.triangles.shape
        #print self.node_index.shape
        #print self.number_of_triangles_per_node.shape

        check_integrity_c(self.vertex_value_indices,
                          self.triangles,
                          self.node_index,
                          self.number_of_triangles_per_node)



#        # Check integrity of vertex_value_indices
#        current_node = 0
#        k = 0 # Track triangles touching on node
#        for index in self.vertex_value_indices:
#
#            if self.number_of_triangles_per_node[current_node] == 0:
#                # Node is lone - i.e. not part of the mesh
#                continue
#
#            k += 1
#
#            volume_id = index / 3
#            vertex_id = index % 3
#
#            msg = 'Triangle %d, vertex %d points to %d. Should have been %d'\
#                  %(volume_id, vertex_id, self.triangles[volume_id, vertex_id], current_node)
#            assert self.triangles[volume_id, vertex_id] == current_node, msg
#
#            if self.number_of_triangles_per_node[current_node] == k:
#                # Move on to next node
#                k = 0
#                current_node += 1


    def get_lone_vertices(self):
        """Return a list of vertices that are not connected to any triangles.

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



    def statistics(self, nbins=10):
        """Output statistics about mesh
        """

        from anuga.utilities.numerical_tools import histogram, create_bins

        vertex_coordinates = self.vertex_coordinates # Relative coordinates
        areas = self.areas
        x = vertex_coordinates[:,0]
        y = vertex_coordinates[:,1]


        #Setup 10 bins for area histogram
        #print "nbins",nbins
        bins = create_bins(areas, 10)
        #print "size bins",bins
        #m = max(areas)
        #bins = arange(0., m, m/10)
        hist = histogram(areas, bins)

        str =  '------------------------------------------------\n'
        str += 'Mesh statistics:\n'
        str += '  Number of triangles = %d\n' %len(self)
        str += '  Extent [m]:\n'
        str += '    x in [%8.5e, %8.5e]\n' %(num.amin(x), num.amax(x))
        str += '    y in [%8.5e, %8.5e]\n' % (num.amin(y), num.amax(y))
        str += '  Areas [m^2]:\n'
        str += '    A in [%8.5e, %8.5e]\n' %(num.amin(areas), num.amax(areas))
        str += '    number of distinct areas: %d\n' %(len(areas))
        str += '    Histogram:\n'

        hi = bins[0]
        for i, count in enumerate(hist):
            lo = hi
            if i+1 < len(bins):
                #Open upper interval
                hi = bins[i+1]
                str += '      [%8.5e, %8.5e[: %d\n' %(lo, hi, count)
            else:
                #Closed upper interval
                hi = num.max(areas)
                str += '      [%8.5e, %8.5e]: %d\n' %(lo, hi, count)

        N = len(areas)
        if N > 10:
            str += '    Percentiles (%g percent):\n' % (100/nbins)
            areas = areas.tolist()
            areas.sort()

            k = 0
            lower = num.min(areas)
            for i, a in enumerate(areas):
                if i % (N/10) == 0 and i != 0: #For every 10% of the sorted areas
                    str += '      %d triangles in [%8.5e, %8.5e]\n' %(i-k, lower, a)
                    lower = a
                    k = i

            str += '      %d triangles in [%8.5e, %8.5e]\n'\
                   %(N-k, lower, max(areas))


        str += '  Boundary:\n'
        str += '    Number of boundary segments == %d\n' %(len(self.boundary))
        str += '    Boundary tags == %s\n' %self.get_boundary_tags()
        str += '------------------------------------------------\n'


        return str


    def get_triangle_containing_point(self, point):
        """Return triangle id for triangle containing specified point (x,y)

        If point isn't within mesh, raise exception

        """

        # FIXME(Ole): This function is currently brute force
        # because I needed it for diagnostics.
        # We should make it fast - probably based on the
        # quad tree structure.
        from anuga.geometry.polygon import is_outside_polygon,\
             is_inside_polygon

        polygon = self.get_boundary_polygon()

        if is_outside_polygon(point, polygon):
            msg = 'Point %s is outside mesh' %str(point)
            raise Exception(msg)


        V = self.get_vertex_coordinates(absolute=True)

        # FIXME: Horrible brute force
        for i, triangle in enumerate(self.triangles):
            poly = V[3*i:3*i+3]

            if is_inside_polygon(point, poly, closed=True):
                return i

        msg = 'Point %s not found within a triangle' %str(point)
        raise Exception(msg)




    def get_intersecting_segments(self, polyline,
                                  use_cache=False,
                                  verbose=False):
        """Find edges intersected by polyline

        Input:
            polyline - list of points forming a segmented line
            use_cache
            verbose

        Output:
            list of instances of class Triangle_intersection

        The polyline may break inside any triangle causing multiple
        segments per triangle - consequently the same triangle may
        appear in several entries.

        If a polyline segment coincides with a triangle edge,
        the the entire shared segment will be used.
        Onle one of the triangles thus intersected will be used and that
        is the first one encountered.

        Intersections with single vertices are ignored.

        Resulting segments are unsorted
        """

        V = self.get_vertex_coordinates()
        N = len(self)

        # Adjust polyline to mesh spatial origin
        polyline = self.geo_reference.get_relative(polyline)

        if use_cache is True:
            segments = cache(get_intersecting_segments,
                             (V, N, polyline),
                             {'verbose': verbose},
                             verbose=verbose)
        else:
            segments = get_intersecting_segments(V, N, polyline,
                                                 verbose=verbose)


        return segments



    def get_triangle_neighbours(self, tri_id):
        """ Given a triangle id, Return an array of the
        3 neighbour triangle id's.

        Negative returned triangle id's represent a boundary as a neighbour.

        If the given triangle id is bad, return an empty list.
        """

        try:
            return self.neighbours[tri_id,:]
        except IndexError:
            return []


    def get_interpolation_object(self):
        """Get object I that will allow linear interpolation using this mesh

        This is a time consuming process but it needs only to be
        once for the mesh.

        Interpolation can then be done using

        result = I.interpolate_block(vertex_values, interpolation_points)

        where vertex values have been obtained from a quantity using
        vertex_values, triangles = self.get_vertex_values()
        """

        if hasattr(self, 'interpolation_object'):
            I = self.interpolation_object
        else:
            from anuga.fit_interpolate.interpolate import Interpolate

            # Get discontinuous mesh - this will match internal
            # representation of vertex values
            triangles = self.get_disconnected_triangles()
            vertex_coordinates = self.get_vertex_coordinates()

            I = Interpolate(vertex_coordinates, triangles)
            self.interpolation_object = I

        return I


class Triangle_intersection:
    """Store information about line segments intersecting a triangle

    Attributes are

        segment: Line segment intersecting triangle [[x0,y0], [x1, y1]]
        normal: [a,b] right hand normal to segment
        length: Length of intersecting segment
        triangle_id: id (in mesh) of triangle being intersected

    """


    def __init__(self,
                 segment=None,
                 normal=None,
                 length=None,
                 triangle_id=None):
        self.segment = segment
        self.normal = normal
        self.length = length
        self.triangle_id = triangle_id


    def __repr__(self):
        s = 'Triangle_intersection('
        s += 'segment=%s, normal=%s, length=%s, triangle_id=%s)'\
             %(self.segment,
               self.normal,
               self.length,
               self.triangle_id)

        return s



def _get_intersecting_segments(V, N, line,
                               verbose=False):
    """Find edges intersected by line

    Input:
        V: Vertex coordinates as obtained by mesh.get_vertex_coordinates()
        N: Number of triangles in mesh
        line - list of two points forming a segmented line
        verbose
    Output:
        list of instances of class Triangle_intersection

    This method is used by the public method
    get_intersecting_segments(self, polyline) which also contains
    more documentation.
    """

    from anuga.geometry.polygon import intersection
    from anuga.geometry.polygon import is_inside_polygon

    msg = 'Line segment must contain exactly two points'
    assert len(line) == 2, msg

    # Origin of intersecting line to be used for
    # establishing direction
    xi0 = line[0][0]
    eta0 = line[0][1]


    # Check intersection with edge segments for all triangles
    # FIXME (Ole): This should be implemented in C
    triangle_intersections={} # Keep track of segments already done
    for i in range(N):
        # Get nodes and edge segments for each triangle
        x0, y0 = V[3*i, :]
        x1, y1 = V[3*i+1, :]
        x2, y2 = V[3*i+2, :]


        edge_segments = [[[x0,y0], [x1, y1]],
                          [[x1,y1], [x2, y2]],
                          [[x2,y2], [x0, y0]]]

        # Find segments that are intersected by line

        intersections = {} # Use dictionary to record points only once
        for edge in edge_segments:

            status, value = intersection(line, edge)
            #if value is not None: log.critical('Triangle %d, status=%s, '
            #                                   'value=%s'
            #                                   % (i, str(status), str(value)))

            if status == 1:
                # Normal intersection of one edge or vertex
                intersections[tuple(value)] = i

                # Exclude singular intersections with vertices
                #if not(allclose(value, edge[0]) or\
                #       allclose(value, edge[1])):
                #    intersections.append(value)

            if status == 2:
                # Edge is sharing a segment with line

                # This is usually covered by the two
                # vertices that would have been picked up
                # under status == 1.
                # However, if coinciding line stops partway
                # along this edge, it will be recorded here.
                intersections[tuple(value[0,:])] = i
                intersections[tuple(value[1,:])] = i


        if len(intersections) == 1:
            # Check if either line end point lies fully within this triangle
            # If this is the case accept that as one end of the intersecting
            # segment

            poly = V[3*i:3*i+3]
            if is_inside_polygon(line[1], poly, closed=False):
                intersections[tuple(line[1])] = i
            elif is_inside_polygon(line[0], poly, closed=False):
                intersections[tuple(line[0])] = i
            else:
                # Ignore situations where one vertex is touch, for instance
                continue


        msg = 'There can be only two or no intersections'
        assert len(intersections) in [0,2], msg


        if len(intersections) == 2:

            # Calculate attributes for this segment


            # End points of intersecting segment
            points = intersections.keys()
            x0, y0 = points[0]
            x1, y1 = points[1]


            # Determine which end point is closer to the origin of the line
            # This is necessary for determining the direction of
            # the line and the normals

            # Distances from line origin to the two intersections
            z0 = num.array([x0 - xi0, y0 - eta0], num.float)
            z1 = num.array([x1 - xi0, y1 - eta0], num.float)
            d0 = num.sqrt(num.sum(z0**2))
            d1 = num.sqrt(num.sum(z1**2))

            if d1 < d0:
                # Swap
                xi, eta = x0, y0
                x0, y0 = x1, y1
                x1, y1 = xi, eta

            # (x0,y0) is now the origin of the intersecting segment


            # Normal direction:
            # Right hand side relative to line direction
            vector = num.array([x1 - x0, y1 - y0], num.float) # Segment vector
            length = num.sqrt(num.sum(vector**2))      # Segment length
            normal = num.array([vector[1], -vector[0]], num.float)/length


            segment = ((x0,y0), (x1, y1))
            T = Triangle_intersection(segment=segment,
                                      normal=normal,
                                      length=length,
                                      triangle_id=i)


            # Add segment unless it was done earlier
            if not triangle_intersections.has_key(segment):
                triangle_intersections[segment] = T


    # Return segments as a list
    return triangle_intersections.values()


def get_intersecting_segments(V, N, polyline,
                              verbose=False):
    """Internal function to find edges intersected by Polyline

    Input:
        V: Vertex coordinates as obtained by mesh.get_vertex_coordinates()
        N: Number of triangles in mesh
        polyline - list of points forming a segmented line
        verbose
    Output:
        list of instances of class Triangle_intersection

    This method is used by the public method
    get_intersecting_segments(self, polyline) which also contains
    more documentation.
    """

    msg = 'Polyline must contain at least two points'
    assert len(polyline) >= 2, msg


    # For all segments in polyline
    triangle_intersections = []
    for i, point0 in enumerate(polyline[:-1]):

        point1 = polyline[i+1]
        if verbose:
            log.critical('Extracting mesh intersections from line:')
            log.critical('(%.2f, %.2f) - (%.2f, %.2f)'
                         % (point0[0], point0[1], point1[0], point1[1]))

        line = [point0, point1]
        triangle_intersections += _get_intersecting_segments(V, N, line,
                                                             verbose=verbose)


    msg = 'No segments found'
    assert len(triangle_intersections) > 0, msg


    return triangle_intersections





def segment_midpoints(segments):
    """Calculate midpoints of all segments

    Inputs:
       segments: List of instances of class Segment

    Ouputs:
       midpoints: List of points
    """

    midpoints = []
    msg = 'Elements of input list to segment_midpoints must be of class Triangle_intersection'
    for segment in segments:
        assert isinstance(segment, Triangle_intersection), msg

        midpoint = num.sum(num.array(segment.segment, num.float), axis=0)/2
        midpoints.append(midpoint)

    return midpoints



