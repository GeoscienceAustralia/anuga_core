
import string
import exceptions
import types


class TreeNode:
    """class TreeNode

    Base class for nodes in a tree-like structure.
    Public Methods:
        GetChildren()
        GetParent()
        AddChild(child)
        GetLeaves()
        DeleteChild(child)
        Prune(children)
        Show()
    """


    def __init__(self,name=None):
        self.children = None
        self.parent = None

        # Subclasses can implement these attributes as functions, called
        # when a node (either leaf or internal) is deleted
        
        # FIXME (Ole): The hasattr statements below were commented out on 12 July 2009 by Ole Nielsen due to 
        # errors appearing when using Python2.6. Tests pass, so I think this was just superfluous anyway.
        #
        # Excerpt from What's New in Python 2.6
        #
        #The hasattr() function was catching and ignoring all errors,
        #under the assumption that they meant a __getattr__() method
        #was failing somehow and the return value of hasattr() would
        #therefore be False. This logic shouldn't be applied to
        #KeyboardInterrupt and SystemExit, however; Python 2.6 will no
        #longer discard such exceptions when hasattr() encounters
        #them. (Fixed by Benjamin Peterson; issue 2196.)
        
        #if not hasattr(self,'ClearLeafNode'):
        #    self.ClearLeafNode = None
        #if not hasattr(self,'ClearInternalNode'):
        #    self.ClearInternalNode = None
        
        # # When not already set by derived class, set unique instance name        
        #if not hasattr(self,'name'):
        #    if name:
        #        self.name = name
        #    else:
        #        self.name = id(self)        

        # When not already set by derived class, set unique instance name
        if name:
            self.name = name
        else:
            self.name = id(self)        
        



    # Make treenode elements appear as sequences such that one
    # can iterate over them             
    def __iter__(self):
        self.index = -1   # incremented on first call to next()
        return self

    def next(self):
        try:
            self.index += 1
            return self.children[self.index]
        except IndexError:
            self.index = -1
            raise StopIteration



    #---- public methods ------------------------------------

    def GetChildren(self):
        return self.children


    def GetParent(self):
        return self.parent


    def AddChild(self,child):
        child.parent = self
        if hasattr(self,'storage'):
            child.storage = self.storage
        if self.children is None:
            self.children = []
        self.children.append(child)


    def DeleteChild(self,child):
        try:
            index = self.children.index(child)
            child.Prune()
            del self.children[index]
            if self.children == []:
                self.children = None
        except ValueError:
            raise TreeNodeError


    def Prune(self):
        if self.children is None:                  # Leaf node
            if callable(self.clear_leaf_node):
                self.clear_leaf_node()
	    else:
	        msg = 'Node must have a method named "clear_leaf_node"' 
		raise msg	
        else:
            if callable(self.clear_internal_node):   # Internal node
                self.clear_internal_node()
	    else:
	        msg = 'Node must have a method named "clear_internal_node"' 
		raise msg	
				
            for child in self.children:
                child.Prune()
            self.children = None


    #FIXME: I think this one is redundant and should could be removed (OMN)
    # Really?  Where is the equivalent method? (DEE)
    # OK - I reliased that we can use it for testing and statistics (OMN)
    def GetLeaves(self):
        if self.children is None:
            myleaves = [self]
        else:
            myleaves = []
            for child in self.children:
                myleaves += child.GetLeaves()

        return myleaves


    def CountLeaves(self):
        if self.children is None:
            count = 1
        else:
            count = 0
            for child in self.children:
                count += child.CountLeaves()

        return count


    def Show(self,depth=0):
        """Traverse tree and print node names.
        """        
        print "%s%s" % ('  '*depth, self.name)
        if self.children:
            for child in self.children:
                child.Show(depth+1)
        else:
            print ': %d' %self.Count()

         
#----------------------------------------------------------------
class TreeNodeError(exceptions.Exception):
    pass
 
