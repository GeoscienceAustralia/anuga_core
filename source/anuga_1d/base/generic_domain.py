"""Class Domain - 1D domains for finite-volume computations of
   the shallow water wave equation


   Copyright 2004
   Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
   Geoscience Australia
"""

import numpy
from generic_boundary_conditions import *


class Generic_domain:

    def __init__(self,
                 coordinates,
                 boundary = None,
                 conserved_quantities = [],
                 evolved_quantities = [],
                 other_quantities = [],
                 tagged_elements = None):
        """
        Build 1D elements from x coordinates
        """

        from config import timestepping_method
        from config import CFL
        
        #Store Points
        self.coordinates = numpy.array(coordinates, numpy.float)


        #Register number of Elements
        self.number_of_elements = N = len(self.coordinates)-1

        self.set_CFL(CFL)
        self.set_timestepping_method(timestepping_method)


        # work array for special vertices (like wet or dry cells
        self.special_vertices = numpy.zeros((N,2), numpy.int) # should this be here

        #Allocate space for neighbour and boundary structures
        self.neighbours = numpy.zeros((N, 2), numpy.int)
        #self.neighbour_edges = numpy.zeros((N, 2), numpy.int)
        self.neighbour_vertices   = numpy.zeros((N, 2), numpy.int)
        self.number_of_boundaries = numpy.zeros(N, numpy.int)
        self.surrogate_neighbours = numpy.zeros((N, 2), numpy.int)
        
        #Allocate space for geometric quantities
        self.vertices  = numpy.zeros((N, 2), numpy.float)
        self.centroids = numpy.zeros(N, numpy.float)
        self.areas     = numpy.zeros(N, numpy.float)

        self.max_speed_array = numpy.zeros(N, numpy.float)
      

        self.normals = numpy.zeros((N, 2), numpy.float)


        xl = self.coordinates[0:-1]
        xr = self.coordinates[1: ]

        self.vertices[:,0] = xl
        self.vertices[:,1] = xr

        self.centroids[:] = (xl+xr)/2.0

        msg = 'Coordinates should be ordered, smallest to largest'
        assert (xr>xl).all(), msg

        #The normal vectors
        # - point outward from each edge
        # - are orthogonal to the edge
        # - have unit length
        # - Are enumerated by left vertex then right vertex normals


        self.normals[:,0] = -1.0
        self.normals[:,1] =  1.0

        self.areas[:] = (xr-xl)

        #Initialise Neighbours (-1 means that it is a boundary neighbour)

        self.neighbours[:, :] = -1

        #Initialise edge ids of neighbours
        self.neighbour_vertices[:, :] = -1


#        for i in range(N):
#            xl = self.coordinates[i]
#            xr = self.coordinates[i+1]
#            self.vertices[i,0] = xl
#            self.vertices[i,1] = xr
#
#            centroid = (xl+xr)/2.0
#            self.centroids[i] = centroid
#
#            msg = 'Coordinates should be ordered, smallest to largest'
#            assert xr>xl, msg
#
#            #The normal vectors
#            # - point outward from each edge
#            # - are orthogonal to the edge
#            # - have unit length
#            # - Are enumerated by left vertex then right vertex normals
#
#            nl = -1.0
#            nr =  1.0
#            self.normals[i,:] = [nl, nr]
#
#            self.areas[i] = (xr-xl)
#
## #         print 'N', N
## #         print 'Centroid', self.centroids
## #         print 'Areas', self.areas
## #         print 'Vertex_Coordinates', self.vertices
#
#            #Initialise Neighbours (-1 means that it is a boundary neighbour)
#            self.neighbours[i, :] = [-1, -1]
#            #Initialise edge ids of neighbours
#            #Initialise vertex ids of neighbours
#            #In case of boundaries this slot is not used
#            #self.neighbour_edges[i, :] = [-1, -1]
#            self.neighbour_vertices[i, :] = [-1, -1]

        self.build_vertexlist()

        #Build neighbour structure
        self.build_neighbour_structure()

        #Build surrogate neighbour structure
        self.build_surrogate_neighbour_structure()          

        #Build boundary dictionary mapping (id, vertex) to symbolic tags
        self.build_boundary_dictionary(boundary)

        #Build tagged element  dictionary mapping (tag) to array of elements
        self.build_tagged_elements_dictionary(tagged_elements)
        
        from quantity import Quantity
        #from quantity_domain import Quantity, Conserved_quantity
        
        #List of quantity names entering
        #the conservation equations
        #(Must be a subset of quantities)
        if conserved_quantities is None:
            self.conserved_quantities = []
        else:
            self.conserved_quantities = conserved_quantities

        if evolved_quantities is None:
            self.evolved_quantities = self.conserved_quantities
        else:
            self.evolved_quantities = evolved_quantities

        if other_quantities is None:
            self.other_quantities = []
        else:
            self.other_quantities = other_quantities


        #Build dictionary of Quantity instances keyed by quantity names
        self.quantities = {}

        #print self.conserved_quantities
        #print self.evolved_quantities
        

        #FIXME: remove later - maybe OK, though....
        for name in self.evolved_quantities:
            self.quantities[name] = Quantity(self)
        for name in self.other_quantities:
            self.quantities[name] = Quantity(self)

        #Create an empty list for explicit forcing terms
        self.forcing_terms = []

        #Defaults
        from config import max_smallsteps, beta_w, beta_h, epsilon, CFL
        self.beta_w = beta_w
        self.beta_h = beta_h
        self.epsilon = epsilon

        #FIXME: Maybe have separate orders for h-limiter and w-limiter?
        #Or maybe get rid of order altogether and use beta_w and beta_h
        self.default_order = 1
        self.order = self.default_order

        self.default_time_order = 1
        self.time_order = self.default_time_order
        
        self.smallsteps = 0
        self.max_smallsteps = max_smallsteps
        self.number_of_steps = 0
        self.number_of_first_order_steps = 0

        #Model time
        self.time = 0.0
        self.finaltime = None
        self.min_timestep = self.max_timestep = 0.0
        self.starttime = 0 #Physical starttime if any (0 is 1 Jan 1970 00:00:00)
        #Checkpointing and storage
        from config import default_datadir
        self.set_datadir(default_datadir)
        self.filename = 'domain'
        self.checkpoint = False

    def __len__(self):
        return self.number_of_elements

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

            #a = self.triangles[i, 0]
            #b = self.triangles[i, 1]
            #c = self.triangles[i, 2]
            a = i
            b = i + 1

            #Register the vertices v as lists of
            #(triangle_id, vertex_id) tuples associated with them
            #This is used for smoothing
            #for vertex_id, v in enumerate([a,b,c]):
            for vertex_id, v in enumerate([a,b]):
                if vertexlist[v] is None:
                    vertexlist[v] = []

                vertexlist[v].append( (i, vertex_id) )

        self.vertexlist = vertexlist

        
    def build_neighbour_structure(self):
        """Update all registered triangles to point to their neighbours.

        Also, keep a tally of the number of boundaries for each triangle

        Postconditions:
          neighbours and neighbour_edges is populated
          neighbours and neighbour_vertices is populated
          number_of_boundaries integer array is defined.
        """

        #Step 1:
        #Build dictionary mapping from segments (2-tuple of points)
        #to left hand side edge (facing neighbouring triangle)

        N = self.number_of_elements
        neighbourdict = {}
        #l_edge = 0
        #r_edge = 1
        l_vertex = 0
        r_vertex = 1
        for i in range(N):

            #Register all segments as keys mapping to current triangle
            #and segment id
            #a = self.triangles[i, 0]
            #b = self.triangles[i, 1]
            #c = self.triangles[i, 2]
            a = self.vertices[i,0]
            b = self.vertices[i,1]
            
            """
            if neighbourdict.has_key((a,b)):
                    msg = "Edge 2 of triangle %d is duplicating edge %d of triangle %d.\n" %(i,neighbourdict[a,b][1],neighbourdict[a,b][0])
                    raise msg
            if neighbourdict.has_key((b,c)):
                    msg = "Edge 0 of triangle %d is duplicating edge %d of triangle %d.\n" %(i,neighbourdict[b,c][1],neighbourdict[b,c][0])
                    raise msg
            if neighbourdict.has_key((c,a)):
                    msg = "Edge 1 of triangle %d is duplicating edge %d of triangle %d.\n" %(i,neighbourdict[c,a][1],neighbourdict[c,a][0])
                    raise msg
                    """
            #neighbourdict[a,b] = (i, 2) #(id, edge)
            #neighbourdict[b,c] = (i, 0) #(id, edge)
            #neighbourdict[c,a] = (i, 1) #(id, edge)
            #neighbourdict[a,b] = (i, 1) #(id, edge)
            #neighbourdict[b,a] = (i, 0) #(id, edge)
            #neighbourdict[a,l_edge] = (i, 0) #(id, edge)
            #neighbourdict[b,r_edge] = (i, 1) #(id, edge)
            neighbourdict[a,l_vertex] = (i, 0) #(id, vertex)
            neighbourdict[b,r_vertex] = (i, 1) #(id, vertex)


        #Step 2:
        #Go through triangles again, but this time
        #reverse direction of segments and lookup neighbours.
        for i in range(N):
            #a = self.triangles[i, 0]
            #b = self.triangles[i, 1]
            #c = self.triangles[i, 2]
            
            a = self.vertices[i,0]
            b = self.vertices[i,1]
            
            #self.number_of_boundaries[i] = 3
            self.number_of_boundaries[i] = 2
            #if neighbourdict.has_key((b,l_edge)):
            if neighbourdict.has_key((b,l_vertex)):
                #self.neighbours[i, 1] = neighbourdict[b,l_edge][0]
                #self.neighbour_edges[i, 1] = neighbourdict[b,l_edge][1]
                self.neighbours[i, 1] = neighbourdict[b,l_vertex][0]
                self.neighbour_vertices[i, 1] = neighbourdict[b,l_vertex][1]
                self.number_of_boundaries[i] -= 1

            #if neighbourdict.has_key((a,r_edge)):
            if neighbourdict.has_key((a,r_vertex)):
                #self.neighbours[i, 0] = neighbourdict[a,r_edge][0]
                #self.neighbour_edges[i, 0] = neighbourdict[a,r_edge][1]
                self.neighbours[i, 0] = neighbourdict[a,r_vertex][0]
                self.neighbour_vertices[i, 0] = neighbourdict[a,r_vertex][1]
                self.number_of_boundaries[i] -= 1
                
            #if neighbourdict.has_key((b,a)):
            #    self.neighbours[i, 1] = neighbourdict[b,a][0]
            #    self.neighbour_edges[i, 1] = neighbourdict[b,a][1]
            #    self.number_of_boundaries[i] -= 1

            #if neighbourdict.has_key((c,b)):
            #    self.neighbours[i, 0] = neighbourdict[c,b][0]
            #    self.neighbour_edges[i, 0] = neighbourdict[c,b][1]
            #    self.number_of_boundaries[i] -= 1

            #if neighbourdict.has_key((a,b)):
            #    self.neighbours[i, 0] = neighbourdict[a,b][0]
            #    self.neighbour_edges[i, 0] = neighbourdict[a,b][1]
            #    self.number_of_boundaries[i] -= 1

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

        N = self.number_of_elements
        for i in range(N):
            #Find all neighbouring volumes that are not boundaries
            #for k in range(3):
            for k in range(2):    
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

        from config import default_boundary_tag

        if boundary is None:
            boundary = {}
            for vol_id in range(self.number_of_elements):
                #for edge_id in range(0, 3):
                #for edge_id in range(0, 2):
                for vertex_id in range(0, 2):    
                    #if self.neighbours[vol_id, edge_id] < 0:
                    if self.neighbours[vol_id, vertex_id] < 0:   
                        #boundary[(vol_id, edge_id)] = default_boundary_tag
                        boundary[(vol_id, vertex_id)] = default_boundary_tag
        else:
            #Check that all keys in given boundary exist
            #for vol_id, edge_id in boundary.keys():
            for vol_id, vertex_id in boundary.keys():    
                #msg = 'Segment (%d, %d) does not exist' %(vol_id, edge_id)
                msg = 'Segment (%d, %d) does not exist' %(vol_id, vertex_id)
                a, b = self.neighbours.shape
                #assert vol_id < a and edge_id < b, msg
                assert vol_id < a and vertex_id < b, msg

                #FIXME: This assert violates internal boundaries (delete it)
                #msg = 'Segment (%d, %d) is not a boundary' %(vol_id, edge_id)
                #assert self.neighbours[vol_id, edge_id] < 0, msg

            #Check that all boundary segments are assigned a tag
            for vol_id in range(self.number_of_elements):
                #for edge_id in range(0, 3):
                #for edge_id in range(0, 2):
                for vertex_id in range(0, 2):    
                    #if self.neighbours[vol_id, edge_id] < 0:
                    if self.neighbours[vol_id, vertex_id] < 0:    
                        #if not boundary.has_key( (vol_id, edge_id) ):
                        if not boundary.has_key( (vol_id, vertex_id) ):    
                            msg = 'WARNING: Given boundary does not contain '
                            #msg += 'tags for edge (%d, %d). '\
                            #       %(vol_id, edge_id)
                            msg += 'tags for vertex (%d, %d). '\
                                   %(vol_id, vertex_id)
                            msg += 'Assigning default tag (%s).'\
                                   %default_boundary_tag

                            #FIXME: Print only as per verbosity
                            #print msg

                            #FIXME: Make this situation an error in the future
                            #and make another function which will
                            #enable default boundary-tags where
                            #tags a not specified
                            #boundary[ (vol_id, edge_id) ] =\
                            boundary[ (vol_id, vertex_id) ] =\
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

        if tagged_elements is None:
            tagged_elements = {}
        else:
            #Check that all keys in given boundary exist
            for tag in tagged_elements.keys():
                tagged_elements[tag] = array(tagged_elements[tag]).astype(numpy.int)

                msg = 'Not all elements exist. '
                assert max(tagged_elements[tag]) < self.number_of_elements, msg
        #print "tagged_elements", tagged_elements
        self.tagged_elements = tagged_elements


    def set_quantities_to_be_stored(self, q):
        """Specify which quantities will be stored in the sww file.

        q must be either:
          - the name of a quantity
          - a list of quantity names
          - None

        In the two first cases, the named quantities will be stored at each
        yieldstep
        (This is in addition to the quantities elevation and friction)  
        If q is None, storage will be switched off altogether.
        """


        if q is None:
            self.quantities_to_be_stored = []            
            self.store = False
            return

        if isinstance(q, basestring):
            q = [q] # Turn argument into a list

        #Check correcness    
        for quantity_name in q:
            msg = 'Quantity %s is not a valid conserved quantity' %quantity_name
            assert quantity_name in self.conserved_quantities, msg 
        
        self.quantities_to_be_stored = q
        




    def get_boundary_tags(self):
        """Return list of available boundary tags
        """

        tags = {}
        for v in self.boundary.values():
            tags[v] = 1

        return tags.keys()
        
    def get_vertex_coordinates(self, obj = False):
        """Return all vertex coordinates.
        Return all vertex coordinates for all triangles as an Nx6 array
        (ordered as x0, y0, x1, y1, x2, y2 for each triangle)

        if obj is True, the x/y pairs are returned in a 3*N x 2 array.
        FIXME, we might make that the default.
        FIXME Maybe use keyword: continuous = False for this condition?

        
        """

        if obj is True:
        
            #V = self.vertex_coordinates
            V = self.vertices
            #return concatenate( (V[:,0:2], V[:,2:4], V[:,4:6]), axis=0)

            N = V.shape[0]
            #return reshape(V, (3*N, 2))
            return numpy.reshape(V, (N, 2))
        else:    
            #return self.vertex_coordinates
            return self.vertices
    
    def get_conserved_quantities(self, vol_id, vertex=None):#, edge=None):
        """Get conserved quantities at volume vol_id

        If vertex is specified use it as index for vertex values
        If edge is specified use it as index for edge values
        If neither are specified use centroid values
        If both are specified an exeception is raised

        Return value: Vector of length == number_of_conserved quantities

        """

        #if not (vertex is None):# or edge is None):
        #    msg = 'Values for both vertex and edge was specified.'
        #    msg += 'Only one (or none) is allowed.'
       #     raise msg

        q = numpy.zeros( len(self.conserved_quantities), numpy.float)

        for i, name in enumerate(self.conserved_quantities):
            Q = self.quantities[name]
            if vertex is not None:
                q[i] = Q.vertex_values[vol_id, vertex]
            #elif edge is not None:
            #    q[i] = Q.edge_values[vol_id, edge]
            else:
                q[i] = Q.centroid_values[vol_id]

        return q


    def get_evolved_quantities(self, vol_id, vertex=None):#, edge=None):
        """Get evolved quantities at volume vol_id

        If vertex is specified use it as index for vertex values
        If edge is specified use it as index for edge values
        If neither are specified use centroid values
        If both are specified an exeception is raised

        Return value: Vector of length == number_of_evolved quantities

        """



        #if not (vertex is None):# or edge is None):
        #    msg = 'Values for both vertex and edge was specified.'
        #    msg += 'Only one (or none) is allowed.'
       #     raise msg

        q = numpy.zeros( len(self.evolved_quantities), numpy.float)

        for i, name in enumerate(self.evolved_quantities):
            Q = self.quantities[name]
            if vertex is not None:
                q[i] = Q.vertex_values[vol_id, vertex]
            #elif edge is not None:
            #    q[i] = Q.edge_values[vol_id, edge]
            else:
                q[i] = Q.centroid_values[vol_id]

        return q
        

    def get_centroids(self):
        """Return all coordinates of centroids
        Return x coordinate of centroid for each element as a N array
        """

        return self.centroids

    def get_vertices(self):
        """Return all coordinates of centroids
        Return x coordinate of centroid for each element as a N array
        """

        return self.vertices

    def get_coordinate(self, elem_id, vertex=None):
        """Return coordinate of centroid,
        or left or right vertex.
        Left vertex (vertex=0). Right vertex (vertex=1)
        """

        if vertex is None:
            return self.centroids[elem_id]
        else:
            return self.vertices[elem_id,vertex]

    def get_area(self, elem_id=None):
        """Return area of element id
        """

        if elem_id is None:
            return sum(self.areas)
        else:
            return self.areas[elem_id]


    def get_quantity(self, name, location='vertices', indices = None):
        """Get values for named quantity

        name: Name of quantity

        In case of location == 'centroids' the dimension values must
        be a list of a numpyal array of length N, N being the number
        of elements. Otherwise it must be of dimension Nx3.

        Indices is the set of element ids that the operation applies to.

        The values will be stored in elements following their
        internal ordering.
        """

        return self.quantities[name].get_values( location, indices = indices)

    def get_centroid_coordinates(self):
        """Return all centroid coordinates.
        Return all centroid coordinates for all triangles as an Nx2 array
        (ordered as x0, y0 for each triangle)
        """
        return self.centroids


    def get_timestepping_method(self):
        return self.timestepping_method

    def set_timestepping_method(self,timestepping_method):
        
        if timestepping_method in ['euler', 'rk2', 'rk3']:
            self.timestepping_method = timestepping_method
            return

        if timestepping_method in [1, 2, 3]:
            self.timestepping_method = ['euler', 'rk2', 'rk3'][timestepping_method-1]
            return

        msg = '%s is an incorrect timestepping type'% timestepping_method
        raise Exception, msg


    def set_quantity(self, name, *args, **kwargs):
        """Set values for named quantity


        One keyword argument is documented here:
        expression = None, # Arbitrary expression

        expression:
          Arbitrary expression involving quantity names

        See Quantity.set_values for further documentation.
        """

        #FIXME (Ole): Allow new quantities here
        #from quantity import Quantity, Conserved_quantity
        #Create appropriate quantity object
        # #if name in self.conserved_quantities:
        # #    self.quantities[name] = Conserved_quantity(self)
        # #else:
        # #    self.quantities[name] = Quantity(self)


        #Do the expression stuff
        if kwargs.has_key('expression'):
            expression = kwargs['expression']
            del kwargs['expression']

            Q = self.create_quantity_from_expression(expression)
            kwargs['quantity'] = Q
        
        #Assign values
        self.quantities[name].set_values(*args, **kwargs)

    def set_boundary(self, boundary_map):
        """Associate boundary objects with tagged boundary segments.

        Input boundary_map is a dictionary of boundary objects keyed
        by symbolic tags to matched against tags in the internal dictionary
        self.boundary.

        As result one pointer to a boundary object is stored for each vertex
        in the list self.boundary_objects.
        More entries may point to the same boundary object

        Schematically the mapping is from two dictionaries to one list
        where the index is used as pointer to the boundary_values arrays
        within each quantity.

        self.boundary:          (vol_id, edge_id): tag
        boundary_map (input):   tag: boundary_object
        ----------------------------------------------
        self.boundary_objects:  ((vol_id, edge_id), boundary_object)


        Pre-condition:
          self.boundary has been built.

        Post-condition:
          self.boundary_objects is built

        If a tag from the domain doesn't appear in the input dictionary an
        exception is raised.
        However, if a tag is not used to the domain, no error is thrown.
        FIXME: This would lead to implementation of a
        default boundary condition

        Note: If a segment is listed in the boundary dictionary and if it is
        not None, it *will* become a boundary -
        even if there is a neighbouring triangle.
        This would be the case for internal boundaries

        Boundary objects that are None will be skipped.

        FIXME: If set_boundary is called multiple times and if Boundary
        object is changed into None, the neighbour structure will not be
        restored!!!
        """

        self.boundary_objects = []
        self.boundary_map = boundary_map  #Store for use with eg. boundary_stats.

        #FIXME: Try to remove the sorting and fix test_mesh.py
        x = self.boundary.keys()
        x.sort()

        #Loop through edges that lie on the boundary and associate them with
        #callable boundary objects depending on their tags
        #for k, (vol_id, edge_id) in enumerate(x):
        for k, (vol_id, vertex_id) in enumerate(x):
            #tag = self.boundary[ (vol_id, edge_id) ]
            tag = self.boundary[ (vol_id, vertex_id) ]

            if boundary_map.has_key(tag):
                B = boundary_map[tag]  #Get callable boundary object

                if B is not None:
                    #self.boundary_objects.append( ((vol_id, edge_id), B) )
                    #self.neighbours[vol_id, edge_id] = -len(self.boundary_objects)
                    self.boundary_objects.append( ((vol_id, vertex_id), B) )
                    self.neighbours[vol_id, vertex_id] = -len(self.boundary_objects)
                else:
                    pass
                    #FIXME: Check and perhaps fix neighbour structure

            else:
                msg = 'ERROR (domain.py): Tag "%s" has not been ' %tag
                msg += 'bound to a boundary object.\n'
                msg += 'All boundary tags defined in domain must appear '
                msg += 'in the supplied dictionary.\n'
                msg += 'The tags are: %s' %self.get_boundary_tags()
                raise Exception, msg



    def check_integrity(self):
        #Mesh.check_integrity(self)

        #print self.quantities
        #print self.conserved_quantities
        
        for quantity in self.conserved_quantities:
            msg = 'Conserved quantities must be a subset of all quantities'
            assert quantity in self.quantities, msg

        for quantity in self.evolved_quantities:
            msg = 'Evolved quantities must be a subset of all quantities'
            assert quantity in self.quantities, msg            

        # #assert hasattr(self, 'boundary_objects')

    def write_time(self):
        print self.timestepping_statistics()

    def timestepping_statistics(self):
        """Return string with time stepping statistics for printing or logging
        """

        msg = ''
        if self.min_timestep == self.max_timestep:
            msg += 'Time = %.4f, delta t = %.8f, steps=%d (%d)'\
                   %(self.time, self.min_timestep, self.number_of_steps,
                     self.number_of_first_order_steps)
        elif self.min_timestep > self.max_timestep:
            msg += 'Time = %.4f, steps=%d (%d)'\
                   %(self.time, self.number_of_steps,
                     self.number_of_first_order_steps)
        else:
            msg += 'Time = %.4f, delta t in [%.8f, %.8f], steps=%d (%d)'\
                   %(self.time, self.min_timestep,
                     self.max_timestep, self.number_of_steps,
                     self.number_of_first_order_steps)

        return msg

    def get_name(self):
        return self.filename

    def set_name(self, name):
        self.filename = name

    def get_datadir(self):
        return self.datadir

    def set_datadir(self, name):
        self.datadir = name

    def set_CFL(self, cfl):
        if cfl > 1.0:
            print 'WARNING: Setting CFL condition to %f which is greater than 1' % cfl
        self.CFL = cfl

    def get_CFL(self):
        return self.CFL
    
    def set_filename(self, name):
        self.filename = name

    def get_filename(self):
        return self.filename


    def get_spatial_order(self):
        return self.order

    def set_spatial_order(self, order):

        possible_orders = [1,2]

        if order in possible_orders:
            self.order = order
            return

        msg = '%s is an incorrect spatial order.\n'% limiter
        msg += 'Possible orders are: '+ ", ".join(["%s" % el for el in possible_orders])
        raise Exception, msg
    

    def get_beta(self):

        warn('limiter parameter beta associated with quantity not domain')

    def set_beta(self,beta):
        """Set the same limiter beta parameter to all evolving quantities
        """

        for name in self.evolved_quantities:
            Q = self.quantities[name]
            Q.set_beta(beta)


    def get_limiter(self):

        warn('limiter associated with quantity not domain')

    def set_limiter(self,limiter):

        for name in self.evolved_quantities:
            Q = self.quantities[name]
            Q.set_limiter(limiter)
        

        
    #--------------------------
    # Main components of evolve
    #--------------------------    

    def evolve(self, yieldstep = None,
               finaltime = None,
               duration = None,
               skip_initial_step = False):
        """Evolve model through time starting from self.starttime.


        yieldstep: numpy.interval between yields where results are stored,
                   statistics written and domain inspected or
                   possibly modified. If omitted the internal predefined
                   max timestep is used.
                   numpy.internally, smaller timesteps may be taken.

        duration: Duration of simulation

        finaltime: Time where simulation should end. This is currently
        relative time.  So it's the same as duration.

        If both duration and finaltime are given an exception is thrown.


        skip_initial_step: Boolean flag that decides whether the first
        yield step is skipped or not. This is useful for example to avoid
        duplicate steps when multiple evolve processes are dove tailed.


        Evolve is implemented as a generator and is to be called as such, e.g.

        for t in domain.evolve(yieldstep, finaltime):
            <Do something with domain and t>


        All times are given in seconds

        """

        from config import min_timestep, max_timestep, epsilon

        # FIXME: Maybe lump into a larger check prior to evolving
        msg = 'Boundary tags must be bound to boundary objects before '
        msg += 'evolving system, '
        msg += 'e.g. using the method set_boundary.\n'
        msg += 'This system has the boundary tags %s '\
               %self.get_boundary_tags()
        assert hasattr(self, 'boundary_objects'), msg


        if yieldstep is None:
            yieldstep = max_timestep
        else:
            yieldstep = numpy.float(yieldstep)

        self._order_ = self.default_order


        if finaltime is not None and duration is not None:
            # print 'F', finaltime, duration
            msg = 'Only one of finaltime and duration may be specified'
            raise msg
        else:
            if finaltime is not None:
                self.finaltime = numpy.float(finaltime)
            if duration is not None:
                self.finaltime = self.starttime + numpy.float(duration)



        N = len(self) # Number of triangles
        self.yieldtime = 0.0 # Track time between 'yields'

        # Initialise interval of timestep sizes (for reporting only)
        self.min_timestep = max_timestep
        self.max_timestep = min_timestep
        self.number_of_steps = 0
        self.number_of_first_order_steps = 0


        # Update ghosts
        self.update_ghosts()

        # Initial update of vertex and edge values
        self.distribute_to_vertices_and_edges()

        # Update extrema if necessary (for reporting)
        self.update_extrema()
        
        # Initial update boundary values
        self.update_boundary()

        # Or maybe restore from latest checkpoint
        if self.checkpoint is True:
            self.goto_latest_checkpoint()

        if skip_initial_step is False:
            yield(self.time)  # Yield initial values

        while True:

            # Evolve One Step, using appropriate timestepping method
            if self.get_timestepping_method() == 'euler':
                self.evolve_one_euler_step(yieldstep,finaltime)
                
            elif self.get_timestepping_method() == 'rk2':
                self.evolve_one_rk2_step(yieldstep,finaltime)

            elif self.get_timestepping_method() == 'rk3':
                self.evolve_one_rk3_step(yieldstep,finaltime)               

            
            # Update extrema if necessary (for reporting)
            self.update_extrema()            


            
            self.yieldtime += self.timestep
            self.number_of_steps += 1
            if self._order_ == 1:
                self.number_of_first_order_steps += 1


            # Yield results
            if finaltime is not None and self.time >= finaltime-epsilon:

                if self.time > finaltime:
                    # FIXME (Ole, 30 April 2006): Do we need this check?
                    # Probably not (Ole, 18 September 2008). Now changed to
                    # Exception
                    msg = 'WARNING (domain.py): time overshot finaltime. '
                    msg += 'Contact Ole.Nielsen@ga.gov.au'
                    raise Exception, msg
                    

                # Yield final time and stop
                self.time = finaltime
                yield(self.time)
                break

            if self.yieldtime >= yieldstep:
                # Yield (intermediate) time and allow inspection of domain

                if self.checkpoint is True:
                    self.store_checkpoint()
                    self.delete_old_checkpoints()

                # Pass control on to outer loop for more specific actions

                yield(self.time)

                # Reinitialise
                self.yieldtime = 0.0
                self.min_timestep = max_timestep
                self.max_timestep = min_timestep
                self.number_of_steps = 0
                self.number_of_first_order_steps = 0
                #self.max_speed_array = 0.0


    def evolve_one_euler_step(self, yieldstep, finaltime):
        """
        One Euler Time Step
        Q^{n+1} = E(h) Q^n
        """

        # Compute fluxes across each element edge
        self.compute_fluxes()


        # Update timestep to fit yieldstep and finaltime
        self.update_timestep(yieldstep, finaltime)     

        # Update conserved quantities
        self.update_conserved_quantities()

        # Update ghosts
        self.update_ghosts()

        # Update vertex and edge values
        self.distribute_to_vertices_and_edges()

        # Update boundary values
        self.update_boundary()

        # Update time
        self.time += self.timestep

        


    def evolve_one_rk2_step(self, yieldstep, finaltime):
        """
        One 2nd order RK timestep
        Q^{n+1} = 0.5 Q^n + 0.5 E(h)^2 Q^n
        """


        # Save initial conserved quantities values
        self.backup_conserved_quantities()            

        #--------------------------------------
        # First euler step
        #--------------------------------------

        # Compute fluxes across each element edge
        self.compute_fluxes()

        # Update timestep to fit yieldstep and finaltime
        self.update_timestep(yieldstep, finaltime)

        # Update conserved quantities
        self.update_conserved_quantities()

        # Update ghosts
        self.update_ghosts()

        # Update vertex and edge values
        self.distribute_to_vertices_and_edges()

        # Update boundary values
        self.update_boundary()

        # Update time
        self.time += self.timestep

        #------------------------------------
        # Second Euler step
        #------------------------------------
            
        # Compute fluxes across each element edge
        self.compute_fluxes()

        # Update conserved quantities
        self.update_conserved_quantities()

        #------------------------------------
        # Combine initial and final values
        # of conserved quantities and cleanup
        #------------------------------------
        
        # Combine steps
        self.saxpy_conserved_quantities(0.5, 0.5)

        #-----------------------------------
        # clean up vertex values
        #-----------------------------------
 
        # Update ghosts
        self.update_ghosts()

        # Update vertex and edge values
        self.distribute_to_vertices_and_edges()

        # Update boundary values
        self.update_boundary()



    def evolve_one_rk3_step(self, yieldstep, finaltime):
        """
        One 3rd order RK timestep
        Q^(1) = 3/4 Q^n + 1/4 E(h)^2 Q^n  (at time t^n + h/2)
        Q^{n+1} = 1/3 Q^n + 2/3 E(h) Q^(1) (at time t^{n+1})
        """

        # Save initial initial conserved quantities values
        self.backup_conserved_quantities()            

        initial_time = self.time
        
        #--------------------------------------
        # First euler step
        #--------------------------------------

        # Compute fluxes across each element edge
        self.compute_fluxes()

        # Update timestep to fit yieldstep and finaltime
        self.update_timestep(yieldstep, finaltime)

        # Update conserved quantities
        self.update_conserved_quantities()

        # Update ghosts
        self.update_ghosts()

        # Update vertex and edge values
        self.distribute_to_vertices_and_edges()

        # Update boundary values
        self.update_boundary()

        # Update time
        self.time += self.timestep

        #------------------------------------
        # Second Euler step
        #------------------------------------
            
        # Compute fluxes across each element edge
        self.compute_fluxes()

        # Update conserved quantities
        self.update_conserved_quantities()

        #------------------------------------
        #Combine steps to obtain intermediate
        #solution at time t^n + 0.5 h
        #------------------------------------

        # Combine steps
        self.saxpy_conserved_quantities(0.25, 0.75)
 
        # Update ghosts
        self.update_ghosts()

        # Update vertex and edge values
        self.distribute_to_vertices_and_edges()

        # Update boundary values
        self.update_boundary()

        # Set substep time
        self.time = initial_time + self.timestep*0.5

        #------------------------------------
        # Third Euler step
        #------------------------------------
            
        # Compute fluxes across each element edge
        self.compute_fluxes()

        # Update conserved quantities
        self.update_conserved_quantities()

        #------------------------------------
        # Combine final and initial values
        # and cleanup
        #------------------------------------
        
        # Combine steps
        self.saxpy_conserved_quantities(2.0/3.0, 1.0/3.0)
 
        # Update ghosts
        self.update_ghosts()

        # Update vertex and edge values
        self.distribute_to_vertices_and_edges()

        # Update boundary values
        self.update_boundary()

        # Set new time
        self.time = initial_time + self.timestep       
        

    def backup_conserved_quantities(self):
        N = len(self) # Number_of_triangles

        # Backup conserved_quantities centroid values
        for name in self.conserved_quantities:
            Q = self.quantities[name]
            Q.backup_centroid_values()        

    def saxpy_conserved_quantities(self,a,b):
        N = len(self) #number_of_triangles

        # Backup conserved_quantities centroid values
        for name in self.conserved_quantities:
            Q = self.quantities[name]
            Q.saxpy_centroid_values(a,b)        



        
    def distribute_to_vertices_and_edges(self):
        """Extrapolate conserved quantities from centroid to
        vertices and edge-midpoints for each volume

        Default implementation is straight first order,
        i.e. constant values throughout each element and
        no reference to non-conserved quantities.
        """

        for name in self.conserved_quantities:
            Q = self.quantities[name]
            if self.order == 1:
                Q.extrapolate_first_order()
            elif self.order == 2:
                Q.extrapolate_second_order()
                #Q.limit()
            else:
                raise 'Unknown order'
            #Q.interpolate_from_vertices_to_edges()


    def update_ghosts(self):
        pass
    
    def update_boundary(self):
        """Go through list of boundary objects and update boundary values
        for all conserved quantities on boundary.
        """

        #FIXME: Update only those that change (if that can be worked out)
        #FIXME: Boundary objects should not include ghost nodes.
        #for i, ((vol_id, edge_id), B) in enumerate(self.boundary_objects):
        #    q = B.evaluate(vol_id, edge_id)
        for i, ((vol_id, vertex_id), B) in enumerate(self.boundary_objects):
            q = B.evaluate(vol_id, vertex_id)
            #print 'q ',q
            for j, name in enumerate(self.evolved_quantities):
                #print 'name %s j = %f \n'%(name,j)
                Q = self.quantities[name]

                Q.boundary_values[i] = q[j]
                #print 'Q=',Q
    
    def update_timestep(self, yieldstep, finaltime):

        from config import min_timestep, max_timestep

        # self.timestep is calculated from speed of characteristics
        # Apply CFL condition here
        timestep = min(self.CFL*self.flux_timestep, max_timestep)

        #Record maximal and minimal values of timestep for reporting
        self.max_timestep = max(timestep, self.max_timestep)
        self.min_timestep = min(timestep, self.min_timestep)

        #Protect against degenerate time steps
        if timestep < min_timestep:

            #Number of consecutive small steps taken b4 taking action
            self.smallsteps += 1

            if self.smallsteps > self.max_smallsteps:
                self.smallsteps = 0 #Reset

                if self.order == 1:
                    msg = 'WARNING: Too small timestep %.16f reached '\
                          %timestep
                    msg += 'even after %d steps of 1 order scheme'\
                           %self.max_smallsteps
                    print msg
                    timestep = min_timestep  #Try enforcing min_step

                    #raise msg
                else:
                    #Try to overcome situation by switching to 1 order
                    print "changing Order 1"
                    self.order = 1

        else:
            self.smallsteps = 0
            if self.order == 1 and self.default_order == 2:
                self.order = 2


        #Ensure that final time is not exceeded
        if finaltime is not None and self.time + timestep > finaltime:
            timestep = finaltime-self.time

        #Ensure that model time is aligned with yieldsteps
        if self.yieldtime + timestep > yieldstep:
            timestep = yieldstep-self.yieldtime

        self.timestep = timestep

    def update_extrema(self):
        pass

    def compute_forcing_terms(self):
        """If there are any forcing functions driving the system
        they should be defined in Domain subclass and appended to
        the list self.forcing_terms
        """
        #Clears explicit_update needed for second order method
        if self.time_order == 2:
            for name in self.conserved_quantities:
                Q = self.quantities[name]
                Q.explicit_update[:] = 0.0
                
        #for name in self.conserved_quantities:
        #    Q = self.quantities[name]
        #    Q.semi_implicit_update[:] = 0.0

        for f in self.forcing_terms:
            f(self)
            

    def update_derived_quantites(self):
        pass
    
    #def update_conserved_quantities(self):
    def update_conserved_quantities(self):
        """Update vectors of conserved quantities using previously
        computed fluxes specified forcing functions.
        """


#        Stage      = self.quantities['stage']
#        Xmom       = self.quantities['xmomentum']
#        Bed        = self.quantities['elevation']
#        Height     = self.quantities['height']
#        Velocity   = self.quantities['velocity']
#
#        #Arrays
#        w_C   = Stage.centroid_values
#        uh_C  = Xmom.centroid_values
#        z_C   = Bed.centroid_values
#        h_C   = Height.centroid_values
#        u_C   = Velocity.centroid_values
#
#        import numpy as np
#        print '== before forcing ======='
#        print 'w_C', np.any(np.isnan(w_C))
#        print 'uh_C', np.any(np.isnan(uh_C))
#        print 'z_C', np.any(np.isnan(z_C))
#        print 'h_C', np.any(np.isnan(h_C))
#        print 'u_C', np.any(np.isnan(u_C))
#        print 'w_ex_update', np.any(np.isnan(Stage.explicit_update))

        timestep = self.timestep



        #Compute forcing terms
        self.compute_forcing_terms()

#        print '==after forcing ======='
#        print 'w_C', np.any(np.isnan(w_C))
#        print 'uh_C', np.any(np.isnan(uh_C))
#        print 'z_C', np.any(np.isnan(z_C))
#        print 'h_C', np.any(np.isnan(h_C))
#        print 'u_C', np.any(np.isnan(u_C))
#        print 'w_ex_update', np.any(np.isnan(Stage.explicit_update))

        #Update conserved_quantities
        for name in self.conserved_quantities:
            Q = self.quantities[name]
            Q.update(timestep)

#        print '==after quantity update  ======='
#        print 'w_C', np.any(np.isnan(w_C))
#        print 'uh_C', np.any(np.isnan(uh_C))
#        print 'z_C', np.any(np.isnan(z_C))
#        print 'h_C', np.any(np.isnan(h_C))
#        print 'u_C', np.any(np.isnan(u_C))
#        print 'w_ex_update', np.any(np.isnan(Stage.explicit_update))

if __name__ == "__main__":

    points1 = [0.0, 1.0, 2.0, 3.0]
    D1 = Domain(points1)

    print D1.get_coordinate(0)
    print D1.get_coordinate(0,1)
    print 'Number of Elements = ',D1.number_of_elements

    try:
        print D1.get_coordinate(3)
    except:
        pass
    else:
        msg =  'Should have raised an out of bounds exception'
        raise msg

    #points2 = [0.0, 1.0, 2.0, 3.0, 2.5]
    #D2 = Domain(points2)
