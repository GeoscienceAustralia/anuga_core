"""quad.py - quad tree data structure for fast indexing of regions in the plane.

This generic structure can be used to store any geometry in a quadtree.
It is naive, and does not exploit any coherency - it merely tests a point
against all bounding boxes in its heirarchy.

It returns a list of bounding boxes which intersect with the test point, which
may then be iterated over with a proper intersection test to detect actual
geometry intersections.

As of June 2010 this module has a pylint quality rating of 10/10.

"""


from builtins import object
import anuga.utilities.log as log

            
class Cell:
    """ One cell in the plane.
        A cell is defined by an AABB, and can have smaller AABB children.
        The children can be rapidly searched for intersections in log(n) time.
    """
  
    def __init__(self, extents, parent,
         name = 'cell'):
        """ Construct a new cell.
            extents is an AABB defining a region on the plane.
            parent is the node above this one, or None if it is root.
        """
    
        self.extents = extents
        self.parent = parent
        
        # The points in this cell     
        self.leaves = []
        self.children = None
        
    
    def __repr__(self):
        """ String representation of the quadtree. """
        ret_str = '%s: leaves: %d' \
               % (self.name , len(self.leaves))    
        if self.children:
            ret_str += ', children: %d' % (len(self.children))
        return ret_str


    def insert(self, new_leaf):
        """ Insert a leaf into the quadtree.
            new_leaf is a tuple of (AABB extents, data), where data can
                     be any user data (geometry, triangle index, etc.).
        """
        if type(new_leaf)==type(list()):
            for leaf in new_leaf:
                self.insert_item(leaf)
        else:
            self.insert_item(new_leaf)


    def insert_item(self, new_leaf):   
        """ Internal recursive insert a single item.
            new_leaf is a tuple of (AABB extents, data), where data can
                     be any user data (geometry, triangle index, etc.).       
        """
        new_region, _ = new_leaf
        
        # recurse down to any children until we get an intersection
        if self.children:
            for child in self.children:
                if child.extents.is_trivial_in(new_region):
                    child.insert_item(new_leaf)
                    return
        else:            
            # try splitting this cell and see if we get a trivial in
            subregion1, subregion2 = self.extents.split()
            
            # option 1 - try splitting 4 ways
            #subregion11, subregion12 = subregion1.split()    
            #subregion21, subregion22 = subregion2.split()
            #regions = [subregion11, subregion12, subregion21, subregion22]
            #for region in regions:
                #if region.is_trivial_in(new_region):
                    #self.children = [Cell(x, parent=self) for x in regions]
                    #self.insert_item(new_leaf)
                    #return               

            # option 2 - try splitting 2 ways - no performance difference
            # noticed in practise between this and the above option.
            if subregion1.is_trivial_in(new_region):
                self.children = [Cell(subregion1, self), \
                                 Cell(subregion2, self)]    
                self.children[0].insert_item(new_leaf)
                return
            elif subregion2.is_trivial_in(new_region):
                self.children = [Cell(subregion1, self), \
                                 Cell(subregion2, self)]    
                self.children[1].insert_item(new_leaf)
                return                
    
        # recursion ended without finding a fit, so attach it as a leaf
        self.leaves.append(new_leaf)
        
     
    def retrieve(self):
        """Get all leaves from this tree.
           return a traversal of the entire tree.
        """
        
        leaves_found = list(self.leaves)
        
        if not self.children:
            return leaves_found

        for child in self.children:
            leaves_found.extend(child.retrieve())
            
        return leaves_found

    def count(self):
        """Count all leaves from this tree.
           return num of leaves in the tree.
        """
        
        leaves_found = len(self.leaves)
        
        if not self.children:
            return leaves_found

        for child in self.children:
            leaves_found += child.count()
            
        return leaves_found        

    def show(self, depth=0):
        """Traverse tree below self, dumping all information.
        """
        if depth == 0:
            log.critical() 
        print('%s%s'  % ('  '*depth, self.name), self.extents, ' [', \
            self.leaves, ']')
        if self.children:
            log.critical()
            for child in self.children:
                child.show(depth+1)

    def search(self, point):
        """
            Search the tree for intersection with leaves
            point is a test point.
            return a list of possible intersections with geometry.
        """
        intersecting_regions = self.test_leaves(point)
        
        # recurse down into nodes that the point passes through
        if self.children:
            for child in self.children:    
                if child.extents.contains(point):
                    intersecting_regions.extend(child.search(point))
             
        return intersecting_regions
 
 
    def test_leaves(self, point):
        """ Test all leaves on this node to see if they intersect x.
            Does not recurse into children.
            x is a point to test
            return a list of leaves that intersect x
        """
        intersecting_regions = []
        
        # test all leaves to see if they intersect the point
        for leaf in self.leaves:
            aabb, data = leaf
            if aabb.contains(point):
                intersecting_regions.append([data, self])
                
        return intersecting_regions                
 
 
    def get_siblings(self):
        """ return siblings of this node. If there is no parent, it
                   returns an empty list.
        """
        if not self.parent:
            return []
         
        siblings = list(self.parent.children)
        siblings.remove(self)
        return siblings
                
