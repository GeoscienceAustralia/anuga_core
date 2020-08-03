"""
    Axially aligned bounding box.
    Contains a class describing a bounding box. It contains methods
    to split itself and return child boxes.
    
    As of June 2010 this module has a pylint quality rating of 10/10.
"""

# Allow children to be slightly bigger than their parents to prevent
# straddling of a boundary
from builtins import object
SPLIT_BORDER_RATIO    = 0.55

class AABB(object):
    """ Axially-aligned bounding box class.
        An axially aligned bounding box (or AABB) defines a box-shaped region
        of the plane which contains any other arbitrary geometry.
        It is useful because intersections can be tested against it very
        rapidly. Once all misses are trivially rejected, a more expensive
        collision detection can be done against the remaining geometry.
        
        This class defines a box which can check for intersections with other
        boxes, or points. It can also subdivide itself, returning
        smaller AABBs. This can be used as the basis of a recursive tree
        structure for optimising geometry searches.
    """
    
    def __init__(self, xmin, xmax=None, ymin=None, ymax=None):
        """ Define axially-algned bounding box.
            xmin is minimum x
            xmax is maximum x (absolute coord, ie, not size)
            ymin is minimum y
            ymax is maximum y (absolute coord, ie, not size)
        """
        if xmax is None:
            # try treating first arg as a list of points
            try:
                xmin[0][0]
            except:
                raise Exception('Single parameter to AABB must be point list.')
                
            self.xmin, self.ymin = self.xmax, self.ymax = xmin[0]
            self.include(xmin[1:])
        else:
            self.xmin = xmin    
            self.xmax = xmax
            self.ymin = ymin    
            self.ymax = ymax


    def __repr__(self):
        return 'AABB(xmin:%f, xmax:%f, ymin:%f, ymax:%f)' \
               % (round(self.xmin,1), round(self.xmax,1), \
                  round(self.ymin,1), round(self.ymax, 1)) 


    def grow(self, amount):
        """ Expand region by given amount.
            amount is a multiplier, ie 1.1 will expand border by 10%.
        """
        self.ymax += (self.ymax-self.ymin)*amount
        self.xmax += (self.xmax-self.xmin)*amount
        self.ymin -= (self.ymax-self.ymin)*amount
        self.xmin -= (self.xmax-self.xmin)*amount    

        
    def size(self):
        """return size as (w,h)"""
        return self.xmax - self.xmin, self.ymax - self.ymin

        
    def split(self, border=SPLIT_BORDER_RATIO):
        """Split box along shorter axis.
           return 2 subdivided AABBs. This is useful for recursive
           algorithms.
           
           border is the overlap between the 2 regions - if set to 0.5
                  it will subdivide them exactly, > 0.5 will create a
                  shared overlapping area. 
        """
        
        width, height = self.size()
        assert width >= 0 and height >= 0
        
        if (width > height):
            # split vertically
            split1 = self.xmin+width*border
            split2 = self.xmax-width*border
            return AABB(self.xmin, split1, self.ymin, self.ymax), \
                   AABB(split2, self.xmax, self.ymin, self.ymax)
        else:
            # split horizontally       
            split1 = self.ymin+height*border
            split2 = self.ymax-height*border
            return AABB(self.xmin, self.xmax, self.ymin, split1), \
                   AABB(self.xmin, self.xmax, split2, self.ymax)    

    
    def is_trivial_in(self, test):
        """ Is trivial in.
            test a box to test against the bounding box
            return True if the test box falls fully within the
                   bounding box without intersecting it.
        """
        if (test.xmin < self.xmin) or (test.xmax > self.xmax):
            return False        
        if (test.ymin < self.ymin) or (test.ymax > self.ymax):
            return False        
        return True
 
    def contains(self, point):
        """ is point within box
            point is a test point
            return True if the point is contained within the box.
        """
        return (self.xmin <= point[0] <= self.xmax) \
                and (self.ymin <= point[1] <= self.ymax)
        
    def include(self, point_list):
        """ Include points in AABB.
            Bounding box will be expanded to include passed points
        
            point_list is a list of points.
        """
        for point in point_list:
            pt_x, pt_y = point
            if pt_x < self.xmin:
                self.xmin = pt_x
            if pt_x > self.xmax:
                self.xmax = pt_x
            if pt_y < self.ymin:
                self.ymin = pt_y                
            if pt_y > self.ymax:
                self.ymax = pt_y                
