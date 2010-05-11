"""quad.py - quad tree data structure for fast indexing of regions in the plane.

This is a generic structure that can be used to store any geometry in a quadtree.


"""

from anuga.utilities.treenode import TreeNode
import string, types, sys
import anuga.utilities.log as log
from aabb import AABB

            
class Cell(TreeNode):
    """class Cell

    One cell in the plane delimited by southern, northern,
    western, eastern boundaries.
    """
  
    def __init__(self, extents, parent,
         name = 'cell'):
  
        # Initialise base classes
        TreeNode.__init__(self, string.lower(name))
    
        self.extents = extents
        self.parent = parent
        
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
            
            # option 1 - try splitting 4 ways
            subregion11, subregion12 = subregion1.split()                                
            subregion21, subregion22 = subregion2.split()
            regions = [subregion11, subregion12, subregion21, subregion22]
            for region in regions:
                if region.is_trivial_in(new_region):
                    self.children = [Cell(x, parent=self) for x in regions]
                    self._insert(new_leaf)
                    return               

            # option 2 - try splitting 2 ways
            #if subregion1.is_trivial_in(new_region):
                #self.children = [Cell(subregion1, self), Cell(subregion2, self)]    
                #self.children[0]._insert(new_leaf)
                #return
            #elif subregion2.is_trivial_in(new_region):
                #self.children = [Cell(subregion1, self), Cell(subregion2, self)]    
                #self.children[1]._insert(new_leaf)
                #return                
    
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
        print '%s%s'  % ('  '*depth, self.name), self.extents,' [', self.leaves, ']'
        if self.children:
            log.critical()
            for child in self.children:
                child.show(depth+1)

    def search(self, x):
        """return a list of possible intersections with geometry"""
      
        intersecting_regions = self.test_leaves(x)
        
        # recurse down into nodes that the point passes through
        if self.children:
            for child in self.children:    
                if child.extents.contains(x):
                    intersecting_regions.extend(child.search(x))
             
        return intersecting_regions
 
 
    def test_leaves(self, x):
        intersecting_regions = []
        
        # test all leaves to see if they intersect the point
        for leaf in self.leaves:
            aabb, data = leaf
            if aabb.contains(x):
                intersecting_regions.append([data, self])
                
        return intersecting_regions                
 
 
    def get_siblings(self):
        """ return parent and siblings of this node.
        """
        if not self.parent:
            return []
         
        siblings = list(self.parent.children)
        siblings.remove(self)
        return siblings
                
