# Allow children to be slightly bigger than their parents to prevent straddling of a boundary
SPLIT_BORDER_RATIO    = 0.55

class AABB:
    """Axially-aligned bounding box class.
       Defines a box which can check for intersections with other boxes,
       or if a point lies within it.
    """
    
    def __init__(self, xmin, xmax, ymin, ymax):
        """ Define axially-algned bounding box.
            xmin is minimum x
            xmax is maximum x (absolute coord, ie, not size)
            ymin is minimum y
            ymax is maximum y (absolute coord, ie, not size)
        """
        self.xmin = xmin    
        self.xmax = xmax
        self.ymin = ymin    
        self.ymax = ymax


    def __repr__(self):
        return 'AABB(xmin:%f, xmax:%f, ymin:%f, ymax:%f)' \
               % (round(self.xmin,1), round(self.xmax,1), round(self.ymin,1), round(self.ymax, 1)) 


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
            return AABB(self.xmin, self.xmin+width*border, self.ymin, self.ymax), \
                   AABB(self.xmax-width*border, self.xmax, self.ymin, self.ymax)
        else:
            # split horizontally       
            return AABB(self.xmin, self.xmax, self.ymin, self.ymin+height*border), \
                   AABB(self.xmin, self.xmax, self.ymax-height*border, self.ymax)    

    
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
 
    def contains(self, x):
        """ is point within box
            x is a test point
            return True if the point is contained within the box.
        """
        return (self.xmin <= x[0] <= self.xmax) and (self.ymin <= x[1] <= self.ymax)
        
