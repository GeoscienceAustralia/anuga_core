"""quad.py - quad tree data structure for fast indexing of regions in the plane.

This is a generic structure that can be used to store any geometry in a quadtree.


"""

from treenode import TreeNode
import string, types, sys
import anuga.utilities.log as log

# Allow children to be slightly bigger than their parents to prevent straddling of a boundary
SPLIT_BORDER_RATIO    = 0.55

class AABB:
    """Axially-aligned bounding box class.
    """
    
    def __init__(self, xmin, xmax, ymin, ymax):
        self.xmin = xmin    
        self.xmax = xmax
        self.ymin = ymin    
        self.ymax = ymax

    def __repr__(self):
        return '(xmin:%f, xmax:%f, ymin:%f, ymax:%f)' \
               % (round(self.xmin,1), round(self.xmax,1), round(self.ymin,1), round(self.ymax, 1)) 
        
    def size(self):
        """return size as (w,h)"""
        return self.xmax - self.xmin, self.ymax - self.ymin
        
    def split(self, border=SPLIT_BORDER_RATIO):
        """Split along shorter axis.
           return 2 subdivided AABBs.
        """
        
        width, height = self.size()
        assert width >= 0 and height >= 0
        
        if (width > height):
            # split vertically
            return AABB(self.xmin, self.xmin+width*border, self.ymin, self.ymax), \
                   AABB(self.xmax-width*border, self.xmax, self.ymin, self.ymax)
        else:
            # split horizontally       
            return AABB(self.xmin, self.xmax, self.ymin, self.ymin+height*border), \
                   AABB(self.xmin, self.xmax, self.ymax-height*border, self.ymax)    
    
    def is_trivial_in(self, test):
        if (test.xmin < self.xmin) or (test.xmax > self.xmax):
            return False        
        if (test.ymin < self.ymin) or (test.ymax > self.ymax):
            return False        
        return True
 
    def contains(self, x, y):
        return (self.xmin <= x <= self.xmax) and (self.ymin <= y <= self.ymax)
            
class Cell(TreeNode):
    """class Cell

    One cell in the plane delimited by southern, northern,
    western, eastern boundaries.

    Public Methods:
        insert(point)
        search(x, y)
        split()
        store()
        retrieve()
        count()
    """
  
    def __init__(self, extents,
         name = 'cell'):
  
        # Initialise base classes
        TreeNode.__init__(self, string.lower(name))
    
        self.extents = extents
        
        # The points in this cell     
        self.leaves = []
        self.children = None
        
    
    def __repr__(self):
        str = '%s: leaves: %d' \
               % (self.name , len(self.leaves))	
        if self.children:
            str += ', children: %d' % (len(self.children))
        return str

   

    def clear(self):
        self.Prune()   # TreeNode method


    def clear_leaf_node(self):
        """Clears storage in leaf node.
    Called from Treenode.
    Must exist.    
    """
        self.leaves = []
    
    
    def clear_internal_node(self):
        """Called from Treenode.    
    Must exist.
    """
        self.leaves = []


    def insert(self, new_leaf):
        # process list items sequentially
        if type(new_leaf)==type(list()):
            ret_val = []
            for leaf in new_leaf:
                self._insert(leaf)
        else:
            self._insert(new_leaf)


    def _insert(self, new_leaf):   
        new_region, data = new_leaf
        
        # recurse down to any children until we get an intersection
        if self.children:
            for child in self.children:
                if child.extents.is_trivial_in(new_region):
                    child._insert(new_leaf)
                    return
        else:            
            # try splitting this cell and see if we get a trivial in
            subregion1, subregion2 = self.extents.split()
            if subregion1.is_trivial_in(new_region):
                self.children = [Cell(subregion1), Cell(subregion2)]    
                self.children[0]._insert(new_leaf)
                return
            elif subregion2.is_trivial_in(new_region):
                self.children = [Cell(subregion1), Cell(subregion2)]    
                self.children[1]._insert(new_leaf)
                return                
    
        # recursion ended without finding a fit, so attach it as a leaf
        self.leaves.append(new_leaf)
        
     
    def retrieve(self):
        """Get all leaves from this tree. """
        
        leaves_found = list(self.leaves)
        
        if not self.children:
            return leaves_found

        for child in self.children:
            leaves_found.extend(child.retrieve())
            
        return leaves_found

    def count(self):
        """Count all leaves from this tree. """
        
        leaves_found = len(self.leaves)
        
        if not self.children:
            return leaves_found

        for child in self.children:
            leaves_found += child.count()
            
        return leaves_found        

    def show(self, depth=0):
        """Traverse tree below self
        """
        if depth == 0:
            log.critical() 
        print '%s%s' % ('  '*depth, self.name), self.extents,' [', self.leaves, ']'
        if self.children:
            log.critical()
            for child in self.children:
                child.show(depth+1)
 

    def search(self, x, y, get_vertices = False):
        """return a list of possible intersections with geometry"""
        
        intersecting_regions = []
        
        # test all leaves to see if they intersect the point
        for leaf in self.leaves:
            aabb, data = leaf
            if aabb.contains(x, y):
                if get_vertices:
                    intersecting_regions.append(leaf)
                else:
                    intersecting_regions.append(data)
        
        # recurse down into nodes that the point passes through
        if self.children:
            for child in self.children:    
                if child.extents.contains(x, y):
                    intersecting_regions.extend(child.search(x, y, get_vertices))
             
        return intersecting_regions
        
#from anuga.pmesh.mesh import Mesh
    
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
    
    root = Cell(AABB(xmin, xmax, ymin, ymax))
    
    N = len(mesh)

    # Get x,y coordinates for all vertices for all triangles
    V = mesh.get_vertex_coordinates(absolute=True)
	
    # Check each triangle
    for i in range(N):
        x0, y0 = V[3*i, :]
        x1, y1 = V[3*i+1, :]
        x2, y2 = V[3*i+2, :]

        # insert a tuple with an AABB, and the triangle index as data
        root._insert((AABB(min([x0, x1, x2]), max([x0, x1, x2]), \
                         min([y0, y1, y2]), max([y0, y1, y2])), \
                         i))

    return root
