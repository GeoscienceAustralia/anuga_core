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
        
