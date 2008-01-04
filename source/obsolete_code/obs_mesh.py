
SET_COLOUR='red'

class Triangle(MeshObject):
    """
    A triangle element, defined by 3 vertices.
    Attributes based on the Triangle program.
    """

    def __init__(self, vertex1, vertex2, vertex3, attribute = None, neighbors = None ):
        """
        Vertices, the initial arguments, are listed in counterclockwise order.
        """
        self.vertices= [vertex1,vertex2, vertex3 ]
        
        if attribute is None:
            self.attribute =""
        else:
            self.attribute = attribute #this is a string
            
        if neighbors is None:
            self.neighbors=[]
        else:
            self.neighbors=neighbors

    def replace(self,new_triangle):
        self = new_triangle

    def longestSideID(self):
        ax = self.vertices[0].x
        ay = self.vertices[0].y
        
        bx = self.vertices[1].x
        by = self.vertices[1].y
        
        cx = self.vertices[2].x
        cy = self.vertices[2].y

        lenA = ((cx-bx)**2+(cy-by)**2)**0.5
        lenB = ((ax-cx)**2+(ay-cy)**2)**0.5
        lenC = ((bx-ax)**2+(by-ay)**2)**0.5
 
        len = [lenA,lenB,lenC]
        return len.index(max(len))

    def rotate(self,offset):
        """
        permute the order of the sides of the triangle
        offset must be 0,1 or 2
        """

        if offset == 0:
            pass
        else:
            if offset == 1:
                self.vertices = [self.vertices[1],self.vertices[2],self.vertices[0]]
                self.neighbors = [self.neighbors[1],self.neighbors[2],self.neighbors[0]]
            if offset == 2:
                self.vertices = [self.vertices[2],self.vertices[0],self.vertices[1]]
                self.neighbors = [self.neighbors[2],self.neighbors[0],self.neighbors[1]]

    def rotate_longest_side(self):
        self.rotate(self.longestSideID())

    def getVertices(self):
        return self.vertices
    
    def get_vertices(self):
        """
        Return a list of the vertices.  The x and y values will be relative
        Easting and Northings for the zone of the current geo_ref.
        """
        return self.vertices
    
    def calcArea(self):
        ax = self.vertices[0].x
        ay = self.vertices[0].y
        
        bx = self.vertices[1].x
        by = self.vertices[1].y
        
        cx = self.vertices[2].x
        cy = self.vertices[2].y
        
        return abs((bx*ay-ax*by)+(cx*by-bx*cy)+(ax*cy-cx*ay))/2
    
    def calcP(self):
        #calculate the perimeter
        ax = self.vertices[0].x
        ay = self.vertices[0].y
        
        bx = self.vertices[1].x
        by = self.vertices[1].y
        
        cx = self.vertices[2].x
        cy = self.vertices[2].y

        a =  ((cx-bx)**2+(cy-by)**2)**0.5 
        b =  ((ax-cx)**2+(ay-cy)**2)**0.5 
        c =  ((bx-ax)**2+(by-ay)**2)**0.5 

        return a+b+c
            
    def setNeighbors(self,neighbor1 = None, neighbor2 = None, neighbor3 = None):
        """
        neighbor1 is the triangle opposite vertex1 and so on.
        Null represents no neighbor
        """
        self.neighbors = [neighbor1, neighbor2, neighbor3]

    def setAttribute(self,attribute):
        """
        neighbor1 is the triangle opposite vertex1 and so on.
        Null represents no neighbor
        """
        self.attribute = attribute #this is a string
        
    def __repr__(self):
        return "[%s,%s]" % (self.vertices,self.attribute)
        

    def draw(self, canvas, tags, scale=1, xoffset = 0, yoffset =0, colour = "green"):
        """
        Draw a triangle, returning the objectID
        """
        return canvas.create_polygon(scale*(self.vertices[1].x + xoffset),
                                     scale*-1*(self.vertices[1].y + yoffset),
                                     scale*(self.vertices[0].x + xoffset),
                                     scale*-1*(self.vertices[0].y + yoffset),
                                     scale*(self.vertices[2].x + xoffset),
                                     scale*-1*(self.vertices[2].y + yoffset),
                                     tags = tags,
                                     outline = colour,fill = '')


        self.setID={}
        #a dictionary of names.
        #multiple sets are allowed, but the gui does not yet
        #support this
        
        self.setID['None']=0
        #contains the names of the sets pointing to the indexes
        #in the list. 
        
	self.sets=[[]]
        #Contains the lists of triangles (triangle sets)

      ##################

        
    def refineSet(self,setName):
        Triangles = self.sets[self.setID[setName]]
        Refine(self,Triangles)

    def selectAllTriangles(self):
        A=[]
        A.extend(self.meshTriangles)
        if not('All' in self.setID.keys()):
            self.setID['All']=len(self.sets)
            self.sets.append(A)
        else:
            self.sets[self.setID['All']]=A
        return 'All'
        # and objectIDs


    def clearSelection(self):
        A = []
        if not('None' in self.setID.keys()):
            self.setID['None']=len(self.sets)
            self.sets.append(A)
        return 'None'

    def drawSet(self,canvas,setName,SCALE,colour=SET_COLOUR):
    #FIXME Draws over previous triangles - may bloat canvas
        Triangles = self.sets[self.setID[setName]]
        for triangle in Triangles:
            triangle.draw(canvas,1,
                          scale = SCALE,
                          colour = colour)
	    
    def undrawSet(self,canvas,setName,SCALE,colour='green'):
    #FIXME Draws over previous lines - may bloat canvas
        Triangles = self.sets[self.setID[setName]]
        for triangle in Triangles:
            triangle.draw(canvas,1,
                          scale = SCALE,
                          colour = colour)

    def weed(self,Vertices,Segments):
        #Depreciated
        #weed out existing duplicates
        print 'len(self.getUserSegments())'
        print len(self.getUserSegments())
        print 'len(self.getUserVertices())'
        print len(self.getUserVertices())

        point_keys = {}
        for vertex in Vertices:
            point = (vertex.x,vertex.y)
            point_keys[point]=vertex
        #inlined would looks very ugly

        line_keys = {}
        for segment in Segments:
            vertex1 = segment.vertices[0]
            vertex2 = segment.vertices[1]
            point1 = (vertex1.x,vertex1.y)
            point2 = (vertex2.x,vertex2.y)
            segment.vertices[0]=point_keys[point1]
            segment.vertices[1]=point_keys[point2]
            vertex1 = segment.vertices[0]
            vertex2 = segment.vertices[1]
            point1 = (vertex1.x,vertex1.y)
            point2 = (vertex2.x,vertex2.y)
            line1 = (point1,point2)
            line2 = (point2,point1)
            if not (line_keys.has_key(line1) \
                 or line_keys.has_key(line2)):
                 line_keys[line1]=segment
        Vertices=point_keys.values()
        Segments=line_keys.values()
        return Vertices,Segments

    def segs_to_dict(self,segments):
        dict={}
        for segment in segments:
            vertex1 = segment.vertices[0]
            vertex2 = segment.vertices[1]
            point1 = (vertex1.x,vertex1.y)
            point2 = (vertex2.x,vertex2.y)
            line = (point1,point2)
            dict[line]=segment
        return dict

    def seg2line(self,s):
        return ((s.vertices[0].x,s.vertices[0].y,)\
                (s.vertices[1].x,s.vertices[1].y))

    def line2seg(self,line,tag=None):
        point0 = self.point2ver(line[0])
        point1 = self.point2ver(line[1])
        return Segment(point0,point1,tag=tag)

    def ver2point(self,vertex):
        return (vertex.x,vertex.y)

    def point2ver(self,point):
        return Vertex(point[0],point[1])

    def smooth_polySet(self,min_radius=0.05):
        #for all pairs of connecting segments:
        #    propose a new segment that replaces the 2

        #    If the difference between the new segment
        #    and the old lines is small: replace the
        #    old lines.

        seg2line = self.seg2line
        ver2point= self.ver2point
        line2seg = self.line2seg
        point2ver= self.point2ver

        #create dictionaries of lines -> segments
        userSegments = self.segs_to_dict(self.userSegments)
        alphaSegments = self.segs_to_dict(self.alphaUserSegments)

        #lump user and alpha segments
        for key in alphaSegments.keys():
            userSegments[key]=alphaSegments[key]

        #point_keys = tuple -> vertex
        #userVertices = vertex -> [line,line] - lines from that node
        point_keys = {}
        userVertices={}
        for vertex in self.getUserVertices():
            point = ver2point(vertex)
            if not point_keys.has_key(point):
                point_keys[point]=vertex
                userVertices[vertex]=[]
        for key in userSegments.keys():
            line = key
            point_0 = key[0]
            point_1 = key[1]
            userVertices[point_keys[point_0]].append(line)
            userVertices[point_keys[point_1]].append(line)

        for point in point_keys.keys():
            try:
            #removed keys can cause keyerrors
                vertex = point_keys[point]
		lines = userVertices[vertex]
    
                #if there are 2 lines on the node
		if len(lines)==2:
		    line_0 = lines[0]
		    line_1 = lines[1]
    
                    #if the tags are the the same on the 2 lines
		    if userSegments[line_0].tag == userSegments[line_1].tag:
			tag = userSegments[line_0].tag 
    
                        #point_a is one of the next nodes, point_b is the other
			if point==line_0[0]:
			    point_a = line_0[1]
			if point==line_0[1]:
			    point_a = line_0[0]
			if point==line_1[0]:
			    point_b = line_1[1]
			if point==line_1[1]:
			    point_b = line_1[0]
    
    
                        #line_2 is proposed
			line_2 = (point_a,point_b)

			#calculate the area of the triangle between
			#the two existing segments and the proposed
			#new segment
			ax = point_a[0]
			ay = point_a[1]
			bx = point_b[0]
			by = point_b[1]
			cx = point[0]
			cy = point[1]
			area=abs((bx*ay-ax*by)+(cx*by-bx*cy)+(ax*cy-cx*ay))/2

			#calculate the perimeter
			len_a =  ((cx-bx)**2+(cy-by)**2)**0.5 
			len_b =  ((ax-cx)**2+(ay-cy)**2)**0.5 
			len_c =  ((bx-ax)**2+(by-ay)**2)**0.5 
			perimeter = len_a+len_b+len_c

			#calculate the radius
			r = area/(2*perimeter)

			#if the radius is small: then replace the existing
			#segments with the new one
			if r < min_radius:
			    if len_c < min_radius: append = False
			    else: append = True
			    #if the new seg is also time, don't add it
			    if append:
				segment = self.line2seg(line_2,tag=tag)

			    list_a=userVertices[point_keys[point_a]]
			    list_b=userVertices[point_keys[point_b]]

			    if line_0 in list_a:
				list_a.remove(line_0)
			    else:
				list_a.remove(line_1)

			    if line_0 in list_b:
				list_b.remove(line_0)
			    else:
				list_b.remove(line_1)

			    if append:
				list_a.append(line_2)
				list_b.append(line_2)
			    else:
				if len(list_a)==0:
				    userVertices.pop(point_keys[point_a])
				    point_keys.pop(point_a)
				if len(list_b)==0:
				    userVertices.pop(point_keys[point_b])
				    point_keys.pop(point_b)

			    userVertices.pop(point_keys[point])
			    point_keys.pop(point)
			    userSegments.pop(line_0)
			    userSegments.pop(line_1)

			    if append:
				userSegments[line_2]=segment
            except:
                pass

        #self.userVerticies = userVertices.keys()
        #self.userSegments = []
        #for key in userSegments.keys():
        #    self.userSegments.append(userSegments[key])
        #self.alphaUserSegments = []

        self.userVerticies = []
        self.userSegments = []
        self.alphaUserSegments = []

        return userVertices,userSegments,alphaSegments

    def triangles_to_polySet(self,setName):
        #self.smooth_polySet()

        seg2line = self.seg2line
        ver2point= self.ver2point
        line2seg = self.line2seg
        point2ver= self.point2ver

        from Numeric import array,allclose
        #turn the triangles into a set
        Triangles = self.sets[self.setID[setName]]
        Triangles_dict = {}
        for triangle in Triangles:
            Triangles_dict[triangle]=None
 

        #create a dict of points to vertexes (tuple -> object)
        #also create a set of vertexes (object -> True)
        point_keys = {}
        userVertices={}
        for vertex in self.getUserVertices():
            point = ver2point(vertex)
            if not point_keys.has_key(point):
                point_keys[point]=vertex
                userVertices[vertex]=True

        #create a dict of lines to segments (tuple -> object)
        userSegments = self.segs_to_dict(self.userSegments)
        #append the userlines in an affine linespace
        affine_lines = Affine_Linespace()
        for line in userSegments.keys():
            affine_lines.append(line)
        alphaSegments = self.segs_to_dict(self.alphaUserSegments)
        for line in alphaSegments.keys():
            affine_lines.append(line)
        
        for triangle in Triangles:
            for i in (0,1,2):
                #for every triangles neighbour:
                if not Triangles_dict.has_key(triangle.neighbors[i]):
                #if the neighbour is not in the set:
                    a = triangle.vertices[i-1]
                    b = triangle.vertices[i-2]
                    #Get possible matches:
                    point_a = ver2point(a)
                    point_b = ver2point(b)
                    midpoint = ((a.x+b.x)/2,(a.y+b.y)/2)
                    line = (point_a,point_b)
                    tag = None


                    #this bit checks for matching lines
                    possible_lines = affine_lines[line] 
                    possible_lines = unique(possible_lines)
                    found = 0                            
                    for user_line in possible_lines:
                        if self.point_on_line(midpoint,user_line):
                            found+=1
                            assert found<2
                            if userSegments.has_key(user_line):
                                parent_segment = userSegments.pop(user_line)
                            if alphaSegments.has_key(user_line):
                                parent_segment = alphaSegments.pop(user_line)
                            tag = parent_segment.tag
                            offspring = [line]
                            offspring.extend(self.subtract_line(user_line,
                                                                line))
                            affine_lines.remove(user_line)
                            for newline in offspring:
                                line_vertices = []
                                for point in newline:
                                    if point_keys.has_key(point):
                                        vert = point_keys[point]
                                    else:
                                        vert = Vertex(point[0],point[1])
                                        userVertices[vert]=True
                                        point_keys[point]=vert
                                    line_vertices.append(vert)
                                segment = Segment(line_vertices[0],
                                                  line_vertices[1],tag)
                                userSegments[newline]=segment
                                affine_lines.append(newline)
                            #break
                    assert found<2



                    #if no matching lines
                    if not found:
                        line_vertices = []
		        for point in line:
			    if point_keys.has_key(point):
			        vert = point_keys[point]
			    else:
				vert = Vertex(point[0],point[1])
				userVertices[vert]=True
				point_keys[point]=vert
			    line_vertices.append(vert)
			segment = Segment(line_vertices[0],
                                          line_vertices[1],tag)
                        userSegments[line]=segment
                        affine_lines.append(line)
        
        self.userVerticies = []
        self.userSegments = []
        self.alphaUserSegments = []

        return userVertices,userSegments,alphaSegments

    def subtract_line(self,parent,child):
        #Subtracts child from parent
        #Requires that the child is a 
        #subline of parent to work.

        from Numeric import allclose,dot,array
        A= parent[0]
        B= parent[1]
        a = child[0]
        b = child[1]

        A_array = array(parent[0])
        B_array = array(parent[1])
        a_array   = array(child[0])
        b_array   = array(child[1])

        assert not A == B
        assert not a == b

        answer = []

        #if the new line does not share a 
        #vertex with the old one
        if not (allclose(A_array,a_array)\
             or allclose(B_array,b_array)\
             or allclose(A_array,b_array)\
             or allclose(a_array,B_array)):
            if dot(A_array-a_array,A_array-a_array) \
            < dot(A_array-b_array,A_array-b_array):
                sibling1 = (A,a)
                sibling2 = (B,b)
                return [sibling1,sibling2]
            else:
                sibling1 = (A,b)
                sibling2 = (B,a)
                return [sibling1,sibling2]

        elif allclose(A_array,a_array):
            if allclose(B_array,b_array):
                return []
	    else:
                sibling = (b,B)
                return [sibling]
    	elif allclose(B_array,b_array):
            sibling = (a,A)
            return [sibling]

        elif allclose(A_array,b_array):
            if allclose(B,a):
                return []
	    else:
                sibling = (a,B)
                return [sibling]
        elif allclose(a_array,B_array):
            sibling = (b,A)
            return [sibling]

    def point_on_line(self,point,line):
        #returns true within a tolerance of 3 degrees
        x=point[0]
        y=point[1]
        x0=line[0][0]
        x1=line[1][0]
        y0=line[0][1]
        y1=line[1][1]
        from Numeric import array, dot, allclose
        from math import sqrt
        tol = 3. #DEGREES
        tol = tol*3.1415/180

        a = array([x - x0, y - y0]) 
        a_normal = array([a[1], -a[0]])      
        len_a_normal = sqrt(sum(a_normal**2)) 

        b = array([x1 - x0, y1 - y0])  	       
        len_b = sqrt(sum(b**2)) 
    
        if abs(dot(a_normal, b)/(len_b*len_a_normal))< tol:
            #Point is somewhere on the infinite extension of the line

            len_a = sqrt(sum(a**2))     	       	                
            if dot(a, b) >= 0 and len_a <= len_b:
               return True
            else:   
               return False
        else:
          return False

    def line_length(self,line):
        x0=line[0][0]
        x1=line[1][0]
        y0=line[0][1]
        y1=line[1][1]
        return ((x1-x0)**2-(y1-y0)**2)**0.5     

    def threshold(self,setName,min=None,max=None,attribute_name='elevation'):
        """
        threshold using  d
        """
        triangles = self.sets[self.setID[setName]]
        A = []

        if attribute_name in self.attributeTitles:
            i = self.attributeTitles.index(attribute_name)
        else: i = -1#no attribute
        if not max == None:
            for t in triangles:
                if (min<self.av_att(t,i)<max):
                    A.append(t)
        else:
            for t in triangles:
                if (min<self.av_att(t,i)):
                    A.append(t)
        self.sets[self.setID[setName]] = A

    def general_threshold(self,setName,min=None,max=None\
              ,attribute_name = 'elevation',function=None):
        """
        Thresholds the triangles 
        """
        from visual.graph import arange,ghistogram,color as colour
        triangles = self.sets[self.setID[setName]]
        A = []
        data=[]
        #data is for the graph

        if attribute_name in self.attributeTitles:
            i = self.attributeTitles.index(attribute_name)
        else: i = -1
        if not max == None:
            for t in triangles:
                value=function(t,i)
                if (min<value<max):
                    A.append(t)
                data.append(value)
        else:
            for t in triangles:
                value=function(t,i)
                if (min<value):
                    A.append(t)
                data.append(value)
        self.sets[self.setID[setName]] = A

        if self.visualise_graph:
            if len(data)>0:
                max=data[0]
                min=data[0]
                for value in data:
                    if value > max:
                        max = value
                    if value < min:
                        min = value

                inc = (max-min)/100

                histogram = ghistogram(bins=arange(min,max,inc),\
                             color = colour.red)
                histogram.plot(data=data)
        
    def av_att(self,triangle,i):
        if i==-1: return 1
        else:
            #evaluates the average attribute of the vertices of a triangle.
            V = triangle.getVertices()
            a0 = (V[0].attributes[i])
            a1 = (V[1].attributes[i])
            a2 = (V[2].attributes[i])
            return (a0+a1+a2)/3

    def Courant_ratio(self,triangle,index):
        """
        Uses the courant threshold
        """
        e = self.av_att(triangle,index)
        A = triangle.calcArea()
        P = triangle.calcP()
        r = A/(2*P)
        e = max(0.1,abs(e))
        return r/e**0.5

    def Gradient(self,triangle,index):
        V = triangle.vertices
        x0, y0, x1, y1, x2, y2, q0, q1, q2 = V[0].x,V[0].y,V[1].x,V[1].y,V[2].x,V[2].y,V[0].attributes[index],V[1].attributes[index],V[2].attributes[index]
        grad_x,grad_y = gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2)
        if ((grad_x**2)+(grad_y**2))**(0.5)<0:
            print ((grad_x**2)+(grad_y**2))**(0.5)
        return ((grad_x**2)+(grad_y**2))**(0.5)
    

    def append_triangle(self,triangle):
        self.meshTriangles.append(triangle)

    def replace_triangle(self,triangle,replacement):
        i = self.meshTriangles.index(triangle)
        self.meshTriangles[i]=replacement
        assert replacement in self.meshTriangles

 
"""Refines triangles

   Implements the #triangular bisection?# algorithm.
 
   
"""

def Refine(mesh, triangles):
    """
    Given a general_mesh, and a triangle number, split
    that triangle in the mesh in half. Then to prevent
    vertices and edges from meeting, keep refining 
    neighbouring triangles until the mesh is clean.
    """
    state = BisectionState(mesh)
    for triangle in triangles:
        if not state.refined_triangles.has_key(triangle):
            triangle.rotate_longest_side()
            state.start(triangle)
            Refine_mesh(mesh, state)

def Refine_mesh(mesh, state):
    """
    """
    state.getState(mesh)
    refine_triangle(mesh,state)
    state.evolve()
    if not state.end:
        Refine_mesh(mesh,state)

def refine_triangle(mesh,state):
    split(mesh,state.current_triangle,state.new_point)
    if state.case == 'one':
        state.r[3]=state.current_triangle#triangle 2

        new_triangle_id = len(mesh.meshTriangles)-1
        new_triangle = mesh.meshTriangles[new_triangle_id]

        split(mesh,new_triangle,state.old_point)
        state.r[2]=new_triangle#triangle 1.2
        state.r[4]=mesh.meshTriangles[len(mesh.meshTriangles)-1]#triangle 1.1
        r = state.r
        state.repairCaseOne()

    if state.case == 'two':
        state.r[2]=mesh.meshTriangles[len(mesh.meshTriangles)-1]#triangle 1

        new_triangle = state.current_triangle

        split(mesh,new_triangle,state.old_point)

        state.r[3]=mesh.meshTriangles[len(mesh.meshTriangles)-1]#triangle 2.1
        state.r[4]=new_triangle#triangle 2.2
        r = state.r

        state.repairCaseTwo()

    if state.case == 'vertex':
        state.r[2]=state.current_triangle#triangle 2
        state.r[3]=mesh.meshTriangles[len(mesh.meshTriangles)-1]#triangle 1
        r = state.r
        state.repairCaseVertex()
        
    if state.case == 'start':
        state.r[2]=mesh.meshTriangles[len(mesh.meshTriangles)-1]#triangle 1
        state.r[3]=state.current_triangle#triangle 2

    if state.next_case == 'boundary':
        state.repairCaseBoundary()


def split(mesh, triangle, new_point):
	"""
	Given a mesh, triangle_id and a new point,
	split the corrosponding triangle into two
	new triangles and update the mesh.
	"""

	new_triangle1 = Triangle(new_point,triangle.vertices[0],
                                 triangle.vertices[1],
                                 attribute = triangle.attribute,
                                 neighbors = None)
	new_triangle2 = Triangle(new_point,triangle.vertices[2],
                                 triangle.vertices[0],
                                 attribute = triangle.attribute,
                                 neighbors = None)

        new_triangle1.setNeighbors(triangle.neighbors[2],None,new_triangle2)
        new_triangle2.setNeighbors(triangle.neighbors[1],new_triangle1,None)

        mesh.meshTriangles.append(new_triangle1)

        triangle.vertices = new_triangle2.vertices
        triangle.neighbors = new_triangle2.neighbors


class State:

    def __init__(self):
        pass

class BisectionState(State):


    def __init__(self,mesh):
        self.len = len(mesh.meshTriangles)
        self.refined_triangles = {}
        self.mesh = mesh
        self.current_triangle = None
        self.case = 'start'
        self.end = False
        self.r = [None,None,None,None,None]

    def start(self, triangle):
        self.current_triangle = triangle
        self.case = 'start'
        self.end = False
        self.r = [None,None,None,None,None]

    def getState(self,mesh):
        if not self.case == 'vertex':
            self.new_point=self.getNewVertex(mesh, self.current_triangle)
            #self.neighbour=self.getNeighbour(mesh, self.current_triangle)
            self.neighbour = self.current_triangle.neighbors[0]
            if not self.neighbour is None:
                self.neighbour.rotate_longest_side()
            self.next_case = self.get_next_case(mesh,self.neighbour,
                                                self.current_triangle)
        if self.case == 'vertex':
            self.new_point=self.old_point


    def evolve(self):
        if self.case == 'vertex':
            self.end = True

        self.last_case = self.case
        self.case = self.next_case

        self.old_point = self.new_point
        self.current_triangle = self.neighbour

        if self.case == 'boundary':
            self.end = True
        self.refined_triangles[self.r[2]]=1
        self.refined_triangles[self.r[3]]=1
        if not self.r[4] is None:
            self.refined_triangles[self.r[4]]=1
        self.r[0]=self.r[2]
        self.r[1]=self.r[3]


    def getNewVertex(self,mesh,triangle):
	coordinate1 = triangle.vertices[1]
	coordinate2 = triangle.vertices[2]
	a = ([coordinate1.x*1.,coordinate1.y*1.])
	b = ([coordinate2.x*1.,coordinate2.y*1.])
        attributes = []
        for i in range(len(coordinate1.attributes)):
            att = (coordinate1.attributes[i]+coordinate2.attributes[i])/2
            attributes.append(att)
	new_coordinate = [((a[0]-b[0])/2+b[0]),((a[1]-b[1])/2+b[1])]
	newVertex = Vertex(new_coordinate[0],new_coordinate[1],
                           attributes = attributes)
        mesh.maxVertexIndex+=1
        newVertex.index = mesh.maxVertexIndex
        mesh.meshVertices.append(newVertex)
	return newVertex

    def get_next_case(self, mesh,neighbour,triangle):
	"""
	Given the locations of two neighbouring triangles, 
	examine the interior indices of their vertices (i.e. 
	0,1 or 2) to determine what how the neighbour needs
	to be refined.
	"""
	if (neighbour is None):
		next_case = 'boundary'
	else:
		if triangle.vertices[1].x==neighbour.vertices[2].x:
		    if triangle.vertices[1].y==neighbour.vertices[2].y:
			next_case = 'vertex'
                if triangle.vertices[1].x==neighbour.vertices[0].x:
		    if triangle.vertices[1].y==neighbour.vertices[0].y:
			next_case = 'two'
		if triangle.vertices[1].x==neighbour.vertices[1].x:
		    if triangle.vertices[1].y==neighbour.vertices[1].y:
			next_case = 'one'
	return next_case



    def repairCaseVertex(self):

        r = self.r


        self.link(r[0],r[2])
        self.repair(r[0])

        self.link(r[1],r[3])
        self.repair(r[1])

        self.repair(r[2])

        self.repair(r[3])


    def repairCaseOne(self):
        r = self.rkey


        self.link(r[0],r[2])
        self.repair(r[0])

        self.link(r[1],r[4])
        self.repair(r[1])

        self.repair(r[4])

    def repairCaseTwo(self):
        r = self.r

        self.link(r[0],r[4])
        self.repair(r[0])

        self.link(r[1],r[3])
        self.repair(r[1])

        self.repair(r[4])

    def repairCaseBoundary(self):
        r = self.r
        self.repair(r[2])
        self.repair(r[3])



    def repair(self,triangle):
        """
        Given a triangle that knows its neighbours, this will
        force the neighbours to comply.

        However, it needs to compare the vertices of triangles
        for this implementation 

        But it doesn't work for invalid neighbour structures
        """
        n=triangle.neighbors
        for i in (0,1,2):
            if not n[i] is None:
                for j in (0,1,2):#to find which side of the list is broken
                    if not (n[i].vertices[j] in triangle.vertices):
                    #ie if j is the side of n that needs fixing
                        k = j
                n[i].neighbors[k]=triangle



    def link(self,triangle1,triangle2):
        """
        make triangle1 neighbors point to t
                #count = 0riangle2
        """
        count = 0
        for i in (0,1,2):#to find which side of the list is broken
            if not (triangle1.vertices[i] in triangle2.vertices):
                j = i
                count+=1
        assert count == 1
        triangle1.neighbors[j]=triangle2

class Discretised_Tuple_Set:
    """
    if a={(0.0):[(0.01),(0.02)],(0.2):[(0.17)]}
    a[(0.01)]=a[(0.0)]=[(0.01),(0.02)]
    a[(10000)]=[] #NOT KEYERROR

    a.append[(0.01)]
    => {0.0:[(0.01),(0.02),(0.01)],0.2:[(0.17)]}

    #NOT IMPLIMENTED
    a.remove[(0.01)]
    => {(0.0):[(0.02),(0.01)],0.2:[(0.17)]}

    a.remove[(0.17)]
    => {(0.0):[(0.02),(0.01)],0.2:[]}
    #NOT IMPLIMENTED
    at a.precision = 2:
        a.round_up_rel[0.0]=
        a.round_flat[0.0]=
        a.round_down_rel[0.0]=

        a.up((0.1,2.04))=

    If t_rel = 0, nothing gets rounded into
    two bins. If t_rel = 0.5, everything does.

    Ideally, precision can be set high, so that multiple
    entries are rarely in the same bin. And t_rel should
    be low (<0.1 for 1 dimension!,<(0.1/n) for small n!!)
    so that it is rare to put items in mutiple bins.

    Ex bins per entry = product(a,b,c...,n)
    a = 1 or 2 s.t. Ex(a) = 1+2*t_rel
    b = 1 or 2 ... 

    BUT!!! to avoid missing the right bin:
    (-10)**(precision+1)*t_rel must be greater than the 
    greatest possible variation that an identical element
    can display.


    Note that if tol = 0.5 (the max allowed) 0.6 will round to .7 and .5
    but not .6 - this looks wrong, but note that *everything* will round,
    so .6 wont be missed as everything close to it will check in .7 and .5.
    """
    def __init__(self,p_rel = 6,t_rel = 0.01):
        self.__p_rel__ = p_rel
        self.__t_rel__ = t_rel

        self.__p_abs__ = p_rel+1
        self.__t_abs__ = t_rel

        assert t_rel <= 0.5
        self.__items__ = {}
        from math import frexp
        self.frexp = frexp
        roundings = [self.round_up_rel,\
        self.round_down_rel,self.round_flat_rel,\
        self.round_down_abs,self.round_up_abs,\
        self.round_flat_abs]#

        self.roundings = roundings

    def __repr__(self):
        return '%s'%self.__items__

    def rounded_keys(self,key):
        key = tuple(key)
        keys = [key]
        keys = self.__rounded_keys__(key)
        return (keys)

    def __rounded_keys__(self,key):
        keys = []
        rounded_key=list(key)
        rounded_values=list(key)

        roundings = list(self.roundings)

        #initialise rounded_values
        round = roundings.pop(0)
        for i in range(len(rounded_values)):
            rounded_key[i]=round(key[i])
            rounded_values[i]={}
            rounded_values[i][rounded_key[i]]=None
        keys.append(tuple(rounded_key))
        
        for round in roundings:
            for i in range(len(rounded_key)):
                rounded_value=round(key[i])
                if not rounded_values[i].has_key(rounded_value):
                    #ie unless round_up_rel = round_down_rel
                    #so the keys stay unique
                    for j in range(len(keys)):
                        rounded_key = list(keys[j])
                        rounded_key[i]=rounded_value
                        keys.append(tuple(rounded_key))
        return keys

    def append(self,item):
        keys = self.rounded_keys(item)
        for key in keys:
            if self.__items__.has_key(key):
                self.__items__[key].append(item)
            else:
                self.__items__[key]=[item]

    def __getitem__(self,key):
        answer = []
        keys = self.rounded_keys(key)
        for key in keys:
            if self.__items__.has_key(key):
                answer.extend(self.__items__[key])
        #if len(answer)==0:
        #    raise KeyError#FIXME or return KeyError
        #                  #FIXME or just return []?
        else:
            return answer #FIXME or unique(answer)?

    def __delete__(self,item):
        keys = self.rounded_keys(item)
        answer = False
        #if any of the possible keys contains
        #a list, return true
        for key in keys:        
            if self.__items__.has_key(key):
                if item in self.__items__[key]:
                    self.__items__[key].remove(item)

    def remove(self,item):
        self.__delete__(item)

    def __contains__(self,item):

        keys = self.rounded_keys(item)
        answer = False
        #if any of the possible keys contains
        #a list, return true
        for key in keys:        
            if self.__items__.has_key(key):
                if item in self.__items__[key]:
                    answer = True
        return answer


    def has_item(self,item):
        return self.__contains__(item)

    def round_up_rel2(self,value):
         t_rel=self.__t_rel__
         #Rounding up the value
         m,e = self.frexp(value)
         m = m/2
         e = e + 1
         #m is the mantissa, e the exponent
         # 0.5 < |m| < 1.0
         m = m+t_rel*(10**-(self.__p_rel__))
         #bump m up
         m = round(m,self.__p_rel__)
         return m*(2.**e)

    def round_down_rel2(self,value):
         t_rel=self.__t_rel__
         #Rounding down the value
         m,e = self.frexp(value)
         m = m/2
         e = e + 1
         #m is the mantissa, e the exponent
         # 0.5 < m < 1.0
         m = m-t_rel*(10**-(self.__p_rel__))
         #bump the |m| down, by 5% or whatever
         #self.p_rel dictates
         m = round(m,self.__p_rel__)
         return m*(2.**e)

    def round_flat_rel2(self,value):
    #redundant
         m,e = self.frexp(value)
         m = m/2
         e = e + 1
         m = round(m,self.__p_rel__)
         return m*(2.**e)

    def round_up_rel(self,value):
         t_rel=self.__t_rel__
         #Rounding up the value
         m,e = self.frexp(value)
         #m is the mantissa, e the exponent
         # 0.5 < |m| < 1.0
         m = m+t_rel*(10**-(self.__p_rel__))
         #bump m up
         m = round(m,self.__p_rel__)
         return m*(2.**e)

    def round_down_rel(self,value):
         t_rel=self.__t_rel__
         #Rounding down the value
         m,e = self.frexp(value)
         #m is the mantissa, e the exponent
         # 0.5 < m < 1.0
         m = m-t_rel*(10**-(self.__p_rel__))
         #bump the |m| down, by 5% or whatever
         #self.p_rel dictates
         m = round(m,self.__p_rel__)
         return m*(2.**e)

    def round_flat_rel(self,value):
    #redundant
         m,e = self.frexp(value)
         m = round(m,self.__p_rel__)
         return m*(2.**e)

    def round_up_abs(self,value):
         t_abs=self.__t_abs__
         #Rounding up the value
         m = value+t_abs*(10**-(self.__p_abs__))
         #bump m up
         m = round(m,self.__p_abs__)
         return m

    def round_down_abs(self,value):
         t_abs=self.__t_abs__
         #Rounding down the value
         m = value-t_abs*(10**-(self.__p_abs__))
         #bump the |m| down, by 5% or whatever
         #self.p_rel dictates
         m = round(m,self.__p_abs__)
         return m

    def round_flat_abs(self,value):
    #redundant?
         m = round(value,self.__p_abs__)
         return m

    def keys(self):
        return self.__items__.keys()


class Mapped_Discretised_Tuple_Set(Discretised_Tuple_Set):
    """
    This is a discretised tuple set, but 
    based on a mapping. The mapping MUST
    return a sequence.

    example: 
    def weight(animal):
        return [animal.weight]
    
    a = Mapped_Discretised_Tuple_Set(weight)
    a.append[cow]
    a.append[fox]
    a.append[horse]

    a[horse]    -> [cow,horse]
    a[dog]      -> [fox]
    a[elephant] -> []
    """
    def __init__(self,mapping,p_rel = 6, t_rel=0.01):
        Discretised_Tuple_Set.__init__\
         (self,p_rel,t_rel = t_rel)
        self.mapping = mapping

    def rounded_keys(self,key):
        mapped_key = tuple(self.mapping(key))
        keys = self.__rounded_keys__(mapped_key)
        return keys

class Affine_Linespace(Mapped_Discretised_Tuple_Set):
    """
    The affine linespace creates a way to record and compare lines.
    Precision is a bit of a hack, but it creates a way to avoid 
    misses caused by round offs (between two lines of different
    lenghts, the short one gets rounded off more). 
    I am starting to think that a quadratic search would be faster.
    Nearly.
    """
    def __init__(self,p_rel=4,t_rel=0.2):
        Mapped_Discretised_Tuple_Set.__init__\
            (self,self.affine_line,\
            p_rel=p_rel,t_rel=t_rel)

        roundings = \
        [self.round_down_rel,self.round_up_rel,self.round_flat_rel]
        self.roundings = roundings
        #roundings = \
        #[self.round_down_abs,self.round_up_abs,self.round_flat_abs]
        #self.roundings = roundings

    def affine_line(self,line):
        point_1 = line[0]
        point_2 = line[1]
        #returns the equation of a line
        #between two points, in the from
        #(a,b,-c), as in ax+by-c=0
        #or line *dot* (x,y,1) = (0,0,0)

        #Note that it normalises the line
        #(a,b,-c) so Line*Line = 1.
        #This does not change the mathematical
        #properties, but it makes comparism
        #easier.

        #There are probably better algorithms.
        x1 = point_1[0]
        x2 = point_2[0]
        y1 = point_1[1]
        y2 = point_2[1]
        dif_x = x1-x2
        dif_y = y1-y2

        if dif_x == dif_y == 0:
            msg = 'points are the same'
            raise msg
        elif abs(dif_x)>=abs(dif_y):
            alpha = (-dif_y)/dif_x
            #a = alpha * b
            b = -1.
            c = (x1*alpha+x2*alpha+y1+y2)/2.
            a = alpha*b
        else:
            beta = dif_x/(-dif_y)
            #b = beta * a
            a = 1.
            c = (x1+x2+y1*beta+y2*beta)/2.
            b = beta*a
        mag = abs(a)+abs(b)
        #This does not change the mathematical
        #properties, but it makes comparism possible.

        #note that the gradient is b/a, or (a/b)**-1.
        #so 

        #if a == 0:
        #    sign_a = 1.
        #else:
        #    sign_a = a/((a**2)**0.5)
        #if b == 0:
        #    sign_b = 1.
        #else:
        #    sign_b = b/((b**2)**0.5)
        #if c == 0:
        #    sign_c = 1.
        #else:
        #    sign_c = c/((c**2)**0.5)
        #a = a/mag*sign_a
        #b = b/mag*sign_b
        #c = c/mag*sign_c
        a = a/mag
        b = b/mag
        c = c/mag
        return a,b,c
