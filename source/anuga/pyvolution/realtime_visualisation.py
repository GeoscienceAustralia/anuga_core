from visual import *

def normal(v1,v2):
    """Safe computation of normalised cross product
    """
    try:
        n = norm(cross(v1, v2))
    except ZeroDivisionError:
        print v1, v2
        raise Exception
        #n = vector(0,0,0)

    return n


class Triangle:
    def __init__(self,v0,v1,v2,color=(1,1,1), frame=None):
        """Create one flat triangular panel with two sides
        """

        #Outward normal

        self.normal=normal(v1-v0, v2-v0)

        self.panel = faces(frame = frame)
        self.backpanel = faces(frame = frame)
        for v in (v0,v1,v2):
            self.panel.append(pos=v, normal = self.normal, color=color)

        for v in (v2,v1,v0):
            self.backpanel.append(pos=v, normal = -self.normal, color=color)


    def move(self,v):
        """Move panel in direction given by vector v
        """

        for vertex in self.panel.pos:
            vertex += v

    def update_height(self,d,fix_baseline=False):
        """Change height of triangle by displacing the third vertex by
        d in the direction perpendicular to the baseline (v1-v0) and the
        outward normal vector
        """

        v0 = self.panel.pos[0]
        v1 = self.panel.pos[1]
        v2 = self.panel.pos[2]

        n = normal(v1-v0, self.normal)

        if fix_baseline:
            self.panel.pos[2] -= d*n
        else:
            self.panel.pos[:2] += d*n

    def update_color(self, c):
        """Change color to c
        """

        self.panel.color = c


    def set_vertexheights(self, heights, floor_heights = None):
        from Numeric import zeros, Float

        from anuga.config import minimum_allowed_height as hmin
        if floor_heights is None:
            floor_heights = zeros(heights.shape, Float)

        all_vertices_below_threshold = True
        for k in range(3):
            w = heights[k]
            z = floor_heights[k]

            if w-z >= hmin:
                all_vertices_below_threshold = False

            vertex = self.panel.pos[k]
            vertex[2] = w

            vertex = self.backpanel.pos[2-k]
            vertex[2] = w


        #Update color to visualise dry areas
        if all_vertices_below_threshold:
            self.panel.color = self.bottom_color
            self.backpanel.color = self.bottom_color
        else:
            self.panel.color = self.top_color
            self.backpanel.color = self.top_color


        #update normal
        v0 = self.panel.pos[0]
        v1 = self.panel.pos[1]
        v2 = self.panel.pos[2]


        n = normal(v1-v0, v2-v0)
        self.panel.normal=n
        self.backpanel.normal=-n




def create_surface(domain):
    """Link surface of domains to their visual representations
    """

    fr = frame() #Default frame to contain all objects
    s=0

    Q = domain.quantities['stage']
    try:
        Z = domain.quantities['elevation']
        bed = True
    except:
        bed = False

    N = Q.vertex_values.shape[0]

    domain.visuals = []
    for i in range(N):

        z0 = Q.vertex_values[i, 0]
        z1 = Q.vertex_values[i, 1]
        z2 = Q.vertex_values[i, 2]

        x0, y0 = domain.get_vertex_coordinate(i,0)
        x1, y1 = domain.get_vertex_coordinate(i,1)
        x2, y2 = domain.get_vertex_coordinate(i,2)

        V0 = vector(x0, y0, z0)
        V1 = vector(x1, y1, z1)
        V2 = vector(x2, y2, z2)


        #Top surface
        c0 = 0.1
        c1 = 0.4
        c2 = 0.99
        if s:
            c2 = 0.7*c2   #To show triangles in slightly different shades

        s = 1-s

        col = (c0,c1,c2)

        visual_top = Triangle(V0,V1,V2,color=col,frame=fr)
        visual_top.top_color = col


        #Bottom surface
        #v0, v1, v2 = volume.vertices
        #v0 = vector(v0.x, v0.y, v0.field_values[index])
        #v1 = vector(v1.x, v1.y, v1.field_values[index])
        #v2 = vector(v2.x, v2.y, v2.field_values[index])
        if bed:
            z0 = Z.vertex_values[i, 0]
            z1 = Z.vertex_values[i, 1]
            z2 = Z.vertex_values[i, 2]
        else:
            z0 = z1 = z2 = 0.0

        V0 = vector(x0, y0, z0)
        V1 = vector(x1, y1, z1)
        V2 = vector(x2, y2, z2)

        c0 = 0.3
        c1 = 0.3
        c2 = 0.3
        if s:
            c2 = 0.4*c2   #To show triangles in slightly different shades

        col = (c0,c1,c2)
        visual_bottom = Triangle(V0, V1, V2, color=col,frame=fr)
        visual_top.bottom_color=col

        domain.visuals.append( (visual_top, visual_bottom) )

    update(domain)
    #print 'Scale', scene.scale


def update(domain):
    """Update vertex heights.
    The argument index refers to which conserved quantity to visualise.
    If domain.smooth is set True, vertex values will be averaged
    yielding a smoother surface.
    """

    from Numeric import array


    Q = domain.quantities['stage']
    N = Q.vertex_values.shape[0]

    #print scene.forward
    #FIXME: Use smoother from anuga.pyvolution instead
    if domain.smooth:
        #Get all average point values
        vertex_heights = {}

        for k in range(N):
            for i in range(3):
                vertex = domain.triangles[k, i]
                if vertex_heights.has_key(vertex):
                    vertex_heights[vertex].append(
                        Q.vertex_values[k, i])
                else:
                    vertex_heights[vertex] = []
                    vertex_heights[vertex].append(
                        Q.vertex_values[k, i])



    for k in range(N):
        if domain.smooth:
            #The averages
            x = zeros(3, Float)
            for i in range(3):
                vertex = domain.triangles[k, i]
                A = array(vertex_heights[vertex])
                x[i] = sum(A)/len(A)
        else:
            #The true values
            x = [Q.vertex_values[k, 0],
                 Q.vertex_values[k, 1],
                 Q.vertex_values[k, 2]]


        #Do it
        floor_heights = array([pos[2] for pos in domain.visuals[k][1].panel.pos])

        domain.visuals[k][0].set_vertexheights(x, floor_heights)




scene.width = 1000
scene.height = 800

#Original
scene.center = (0.5,0.5,0)
scene.forward = vector(0.0, 0.5, -0.5)

#Temporary (for bedslope)
#scene.forward = vector(0.0006, 0.7, -0.03)


#Temporary for hackett - begin
#scene.autoscale = 0
#scene.scale = (0.002, 0.002, 0.01) #Scale z so that countours stand out more
#scene.center = (300.0,500.0,-10)
#Temporary for hackett - end


scene.ambient = 0.4
scene.lights = [(0.6, 0.3, 0.2), (0.1, -0.5, 0.4), (-0.1, 0.1, -0.4),
                (-0.2, 0.2, 0.1)]




