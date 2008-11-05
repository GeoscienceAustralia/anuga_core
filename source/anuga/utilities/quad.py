
"""quad.py - quad tree data structure for fast indexing of points in the plane


"""

from treenode import TreeNode
import string, types, sys

#FIXME verts are added one at a time. 
#FIXME add max min x y in general_mesh

class Cell(TreeNode):
    """class Cell

    One cell in the plane delimited by southern, northern,
    western, eastern boundaries.

    Public Methods:
        prune()
        insert(point)
        search(x, y)
        collapse()
        split()
        store()
        retrieve()
        count()
    """
  
    def __init__(self, southern, northern, western, eastern, mesh,
		 name = 'cell',
    	         max_points_per_cell = 4):
  
        # Initialise base classes
        TreeNode.__init__(self, string.lower(name))
	
        # Initialise cell
        self.southern = round(southern,5)    
        self.northern = round(northern,5)
        self.western = round(western,5)    
        self.eastern = round(eastern,5)
        self.mesh = mesh

        # The points in this cell     
        self.points = []
	
	self.max_points_per_cell = max_points_per_cell
        
	
    def __repr__(self):
        return self.name  


    def spawn(self):
        """Create four child cells unless they already exist
        """

        if self.children:
            return
        else:
            self.children = []

        # convenience variables
        cs = self.southern    
        cn = self.northern
        cw = self.western    
        ce = self.eastern   
        mesh = self.mesh

        # create 4 child cells
        self.AddChild(Cell((cn+cs)/2,cn,cw,(cw+ce)/2,mesh,self.name+'_nw'))
        self.AddChild(Cell((cn+cs)/2,cn,(cw+ce)/2,ce,mesh,self.name+'_ne'))
        self.AddChild(Cell(cs,(cn+cs)/2,(cw+ce)/2,ce,mesh,self.name+'_se'))
        self.AddChild(Cell(cs,(cn+cs)/2,cw,(cw+ce)/2,mesh,self.name+'_sw'))
        
 
    def search(self, x, y, get_vertices=False):
        """Find all point indices sharing the same cell as point (x, y)
        """
        branch = []
        points = []
        if self.children:
            for child in self:
                if child.contains(x,y):
                    brothers = list(self.children)
                    brothers.remove(child)
                    branch.append(brothers)
                    points, branch = child.search_branch(x,y, branch,
                                                  get_vertices=get_vertices)
        else:
            # Leaf node: Get actual waypoints
            points = self.retrieve(get_vertices=get_vertices)
        self.branch = branch   
        return points


    def search_branch(self, x, y, branch, get_vertices=False):
        """Find all point indices sharing the same cell as point (x, y)
        """
        points = []
        if self.children:
            for child in self:
                if child.contains(x,y):
                    brothers = list(self.children)
                    brothers.remove(child)
                    branch.append(brothers)
                    points, branch = child.search_branch(x,y, branch,
                                                  get_vertices=get_vertices)
                    
        else:
            # Leaf node: Get actual waypoints
            points = self.retrieve(get_vertices=get_vertices)      
        return points, branch


    def expand_search(self, get_vertices=False):
        """Find all point indices 'up' one cell from the last search
        """
        
        points = []
        if self.branch == []:
            points = []
        else:
            three_cells = self.branch.pop()
            for cell in three_cells:
                #print "cell ", cell.show() 
                points += cell.retrieve(get_vertices=get_vertices)
        return points, self.branch


    def contains(*args):    
        """True only if P's coordinates lie within cell boundaries
        This methods has two forms:
        
        cell.contains(index)
        #True if cell contains indexed point
        cell.contains(x, y)
        #True if cell contains point (x,y)	
        """
        self = args[0]
        if len(args) == 2:
            point_id = int(args[1])
            x, y = self.mesh.get_node(point_id, absolute=True)
        elif len(args) == 3:
            x = float(args[1])
            y = float(args[2])
        else:
            msg = 'Number of arguments to method must be two or three'
            raise msg    		      
	
        if y <  self.southern: return False
        if y >= self.northern: return False
        if x <  self.western:  return False
        if x >= self.eastern:  return False
        return True
    
    
    def insert(self, points, split = False):
        """insert point(s) in existing tree structure below self
           and split if requested
        """

        # Call insert for each element of a list of points
        if type(points) == types.ListType:
            for point in points:
                self.insert(point, split)
        else:
            #Only one point given as argument	
            point = points
	
            # Find appropriate cell
            if self.children is not None:
                for child in self:
                    if child.contains(point):
                        child.insert(point, split)
                        break
            else:
                # self is a leaf cell: insert point into cell
                if self.contains(point):
                    self.store(point)
                else:
                    # Have to take into account of georef.
                    #x = self.mesh.coordinates[point][0]
                    #y = self.mesh.coordinates[point][1]
                    node = self.mesh.get_node(point, absolute=True)
                    print "node", node
                    print "(" + str(node[0]) + "," + str(node[1]) + ")"
                    raise 'point not in region: %s' %str(point)
		
		
        #Split datastructure if requested        
	if split is True:
            self.split()
		


    def store(self,objects):
        
        if type(objects) not in [types.ListType,types.TupleType]:
            self.points.append(objects)
        else:
            self.points.extend(objects)


    def retrieve_triangles(self):
        """return a list of lists. For the inner lists,
        The first element is the triangle index,
        the second element is a list.for this list
           the first element is a list of three (x, y) vertices,
           the following elements are the three triangle normals.

        This info is used in searching for a triangle that a point is in.

        Post condition
        No more points can be added to the quad tree, since the
        points data structure is removed.
        """
        # FIXME Tidy up the structure that is returned.
        # if the triangles att has been made
        # return it.
        if not hasattr(self,'triangles'):
            # use a dictionary to remove duplicates
            triangles = {}
            verts = self.retrieve_vertices()
            for vert in verts:
                triangle_list = self.mesh.get_triangles_and_vertices_per_node(vert)
                for k, _ in triangle_list:
                    if not triangles.has_key(k):
                        # print 'k',k
                        tri = self.mesh.get_vertex_coordinates(k, absolute=True)
                        n0 = self.mesh.get_normal(k, 0)
                        n1 = self.mesh.get_normal(k, 1)
                        n2 = self.mesh.get_normal(k, 2) 
                        triangles[k]=(tri, (n0, n1, n2))
            self.triangles = triangles.items()
            # Delete the old cell data structure to save memory
            del self.points 
        return self.triangles
            
    def retrieve_vertices(self):
         return self.points


    def retrieve(self, get_vertices=True):
         objects = []
         if self.children is None:
             if get_vertices is True:
                 objects = self.retrieve_vertices()
             else:
                 objects =  self.retrieve_triangles()
         else:  
             for child in self:
                 objects += child.retrieve(get_vertices=get_vertices)
         return objects
        

    def count(self, keywords=None):
        """retrieve number of stored objects beneath this node inclusive
        """
        
        num_waypoint = 0
        if self.children:
            for child in self:
                num_waypoint = num_waypoint + child.count()
        else:
            num_waypoint = len(self.points)
        return num_waypoint
  

    def clear(self):
        self.Prune()   # TreeNode method


    def clear_leaf_node(self):
        """Clears storage in leaf node.
	Called from Treenod.
	Must exist.	
	"""
        self.points = []
	
	
    def clear_internal_node(self):
        """Called from Treenode.    
	Must exist.
	"""
        pass



    def split(self, threshold=None):
        """
        Partition cell when number of contained waypoints exceeds 
        threshold.  All waypoints are then moved into correct 
        child cell.
        """
        if threshold == None:
           threshold = self.max_points_per_cell
	    
        #FIXME, mincellsize removed.  base it on side length, if needed
        
        #Protect against silly thresholds such as -1
        if threshold < 1:
            return
	
        if not self.children:               # Leaf cell
            if self.count() > threshold :   
                #Split is needed
                points = self.retrieve()    # Get points from leaf cell
                self.clear()                # and remove them from storage
                    
                self.spawn()                # Spawn child cells and move
                for p in points:            # points to appropriate child
                    for child in self:
                        if child.contains(p):
                            child.insert(p) 
                            break
                        
        if self.children:                   # Parent cell
            for child in self:              # split (possibly newly created) 
                child.split(threshold)      # child cells recursively
                


    def collapse(self,threshold=None):
        """
        collapse child cells into immediate parent if total number of contained waypoints
        in subtree below is less than or equal to threshold.
        All waypoints are then moved into parent cell and
        children are removed. If self is a leaf node initially, do nothing.
        """
        
        if threshold is None:
            threshold = self.max_points_per_cell	


        if self.children:                   # Parent cell    
            if self.count() <= threshold:   # collapse
                points = self.retrieve()    # Get all points from child cells
                self.clear()                # Remove children, self is now a leaf node
                self.insert(points)         # Insert all points in local storage
            else:                          
                for child in self:          # Check if any sub tree can be collapsed
                    child.collapse(threshold)


    def Get_tree(self,depth=0):
        """Traverse tree below self
           Print for each node the name and
           if it is a leaf the number of objects
        """
        s = ''
        if depth == 0:
            s = '\n'
            
        s += "%s%s:" % ('  '*depth, self.name)
        if self.children:
            s += '\n'
            for child in self.children:
                s += child.Get_tree(depth+1)
        else:
            s += '(#wp=%d)\n' %(self.count())

        return s

	
    def show(self, depth=0):
        """Traverse tree below self
           Print for each node the name and
           if it is a leaf the number of objects
        """
        if depth == 0:
            print 
        print "%s%s" % ('  '*depth, self.name),
        if self.children:
            print
            for child in self.children:
                child.show(depth+1)
        else:
            print '(xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f): [%d]'\
		  %(self.western, self.eastern,
		    self.southern, self.northern,
		    self.count()) 


    def show_all(self,depth=0):
        """Traverse tree below self
           Print for each node the name and if it is a leaf all its objects
        """
        if depth == 0:
            print 
        print "%s%s:" % ('  '*depth, self.name),
        if self.children:
            print
            for child in self.children:
                child.show_all(depth+1)
        else:
            print '%s' %self.retrieve()


    def stats(self,depth=0,min_rad=sys.maxint,max_depth=0,max_points=0):
        """Traverse tree below self and find minimal cell radius,
           maximumtree depth and maximum number of waypoints per leaf.
        """

        if self.children:
            for child in self.children:
                min_rad, max_depth, max_points =\
                         child.Stats(depth+1,min_rad,max_depth,max_points)
        else:
            #FIXME remvoe radius stuff
            #min_rad = sys.maxint
            #if self.radius < min_rad:   min_rad = self.radius
            if depth > max_depth: max_depth = depth
            num_points = self.count()
            if num_points > max_points: max_points = num_points

        #return min_rad, max_depth, max_points    
	return max_depth, max_points    
	

    #Class initialisation method
    # this is bad.  It adds a huge memory structure to the class.
    # When the instance is deleted the mesh hangs round (leaks).
    #def initialise(cls, mesh):
    #    cls.mesh = mesh

    #initialise = classmethod(initialise)

def build_quadtree(mesh, max_points_per_cell = 4):
    """Build quad tree for mesh.

    All vertices in mesh are stored in quadtree and a reference
    to the root is returned.
    """

    #Make root cell
    #print mesh.coordinates

    xmin, xmax, ymin, ymax = mesh.get_extent(absolute=True)
    
    # Ensure boundary points are fully contained in region
    # It is a property of the cell structure that
    # points on xmax or ymax of any given cell
    # belong to the neighbouring cell.
    # Hence, the root cell needs to be expanded slightly
    ymax += (ymax-ymin)/10
    xmax += (xmax-xmin)/10

    # To avoid round off error
    ymin -= (ymax-ymin)/10
    xmin -= (xmax-xmin)/10   

    #print "xmin", xmin 
    #print "xmax", xmax
    #print "ymin", ymin 
    #print "ymax", ymax
    
    #FIXME: Use mesh.filename if it exists
    # why?
    root = Cell(ymin, ymax, xmin, xmax,mesh,
                max_points_per_cell = max_points_per_cell)

    #root.show()
    
    #Insert indices of all vertices
    root.insert( range(mesh.number_of_nodes) )

    #Build quad tree and return
    root.split()

    return root
