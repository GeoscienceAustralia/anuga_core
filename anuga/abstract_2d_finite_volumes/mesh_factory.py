"""Library of standard meshes and facilities for reading various
mesh file formats
"""

from builtins import object
import anuga.utilities.log as log
import numpy as num


def rectangular_old(m, n, len1=1.0, len2=1.0, origin = (0.0, 0.0)):

    """Setup a rectangular grid of triangles
    with m+1 by n+1 grid points
    and side lengths len1, len2. If side lengths are omitted
    the mesh defaults to the unit square.

    len1: x direction (left to right)
    len2: y direction (bottom to top)

    Return to lists: points and elements suitable for creating a Mesh or
    FVMesh object, e.g. Mesh(points, elements)
    """

    from anuga.config import epsilon

    deltax = float(len1)/m
    deltay = float(len2)/n

    #Dictionary of vertex objects
    vertices = {}
    points = []

    for i in range(m+1):
        for j in range(n+1):
            vertices[i,j] = len(points)
            points.append([i*delta1 + origin[0], j*delta2 + origin[1]])


    #Construct 2 triangles per rectangular element and assign tags to boundary
    elements = []
    boundary = {}
    for i in range(m):
        for j in range(n):
            v1 = vertices[i,j+1]
            v2 = vertices[i,j]
            v3 = vertices[i+1,j+1]
            v4 = vertices[i+1,j]

            #Update boundary dictionary and create elements
            if i == m-1:
                boundary[(len(elements), 2)] = 'right'
            if j == 0:
                boundary[(len(elements), 1)] = 'bottom'
            elements.append([v4,v3,v2]) #Lower element

            if i == 0:
                boundary[(len(elements), 2)] = 'left'
            if j == n-1:
                boundary[(len(elements), 1)] = 'top'
            elements.append([v1,v2,v3]) #Upper element

    return points, elements, boundary

def rectangular(m, n, len1=1.0, len2=1.0, origin = (0.0, 0.0)):

    """Setup a rectangular grid of triangles
    with m+1 by n+1 grid points
    and side lengths len1, len2. If side lengths are omitted
    the mesh defaults to the unit square.

    len1: x direction (left to right)
    len2: y direction (bottom to top)

    Return to lists: points and elements suitable for creating a Mesh or
    FVMesh object, e.g. Mesh(points, elements)
    """

    from anuga.config import epsilon

    delta1 = float(len1)/m
    delta2 = float(len2)/n

    #Calculate number of points
    Np = (m+1)*(n+1)

    class Index(object):

        def __init__(self, n,m):
            self.n = n
            self.m = m

        def __call__(self, i,j):
            return j+i*(self.n+1)


    index = Index(n,m)

    points = num.zeros((Np, 2), float)

    for i in range(m+1):
        for j in range(n+1):

            points[index(i,j),:] = [i*delta1 + origin[0], j*delta2 + origin[1]]

    #Construct 2 triangles per rectangular element and assign tags to boundary
    #Calculate number of triangles
    Nt = 2*m*n


    elements = num.zeros((Nt, 3), int)
    boundary = {}
    nt = -1
    for i in range(m):
        for j in range(n):
            nt = nt + 1
            i1 = index(i,j+1)
            i2 = index(i,j)
            i3 = index(i+1,j+1)
            i4 = index(i+1,j)


            #Update boundary dictionary and create elements
            if i == m-1:
                boundary[nt, 2] = 'right'
            if j == 0:
                boundary[nt, 1] = 'bottom'
            elements[nt,:] = [i4,i3,i2] #Lower element
            nt = nt + 1

            if i == 0:
                boundary[nt, 2] = 'left'
            if j == n-1:
                boundary[nt, 1] = 'top'
            elements[nt,:] = [i1,i2,i3] #Upper element

    return points, elements, boundary

def rectangular_cross(m, n, len1=1.0, len2=1.0, origin = (0.0, 0.0)):
    """Setup a rectangular grid of triangles
    with m+1 by n+1 grid points
    and side lengths len1, len2. If side lengths are omitted
    the mesh defaults to the unit square.

    len1: x direction (left to right)
    len2: y direction (bottom to top)

    Return to lists: points and elements suitable for creating a Mesh or
    Domain object, e.g. Mesh(points, elements)
    """

    len1 = float(len1)
    len2 = float(len2)

    m = int(m)
    n = int(n)

    params = []
    params.append(m)
    params.append(n)
    params.append(len1)
    params.append(len2)

    arrParams = num.array(params, dtype=float)
    arrOrigin = num.array(origin, dtype=float)
    
    points = num.empty([(m+1)*(n+1)+m*n,2], dtype=float)
    elements = num.empty([4*m*n,3], dtype=int)

    from .mesh_factory_ext import rectangular_cross_construct
    boundary = rectangular_cross_construct(arrParams, arrOrigin, points, elements)

    #points = list(arrPoints)
    #elements = list(arrElements)

    return points, elements, boundary

def rectangular_cross_python(m, n, len1=1.0, len2=1.0, origin = (0.0, 0.0)):

    """Setup a rectangular grid of triangles
    with m+1 by n+1 grid points
    and side lengths len1, len2. If side lengths are omitted
    the mesh defaults to the unit square.

    len1: x direction (left to right)
    len2: y direction (bottom to top)

    Return to lists: points and elements suitable for creating a Mesh or
    Domain object, e.g. Mesh(points, elements)
    """

    from anuga.config import epsilon

    delta1 = float(len1)/m
    delta2 = float(len2)/n

    #Dictionary of vertex objects
    vertices = {}
    points = []

    for i in range(m+1):
        for j in range(n+1):
            vertices[i,j] = len(points)
            points.append([delta1*i + origin[0], delta2*j  + origin[1]])

    # Construct 4 triangles per element
    elements = []
    boundary = {}
    for i in range(m):
        for j in range(n):
            v1 = vertices[i,j+1]
            v2 = vertices[i,j]
            v3 = vertices[i+1,j+1]
            v4 = vertices[i+1,j]
            x = (points[v1][0]+points[v2][0]+points[v3][0]+points[v4][0])*0.25
            y = (points[v1][1]+points[v2][1]+points[v3][1]+points[v4][1])*0.25

            # Create centre point
            v5 = len(points)
            points.append([x, y])

            #Create left triangle
            if i == 0:
                boundary[(len(elements), 1)] = 'left'
            elements.append([v2,v5,v1])

            #Create bottom triangle
            if j == 0:
                boundary[(len(elements), 1)] = 'bottom'
            elements.append([v4,v5,v2])

            #Create right triangle
            if i == m-1:
                boundary[(len(elements), 1)] = 'right'
            elements.append([v3,v5,v4])

            #Create top triangle
            if j == n-1:
                boundary[(len(elements), 1)] = 'top'
            elements.append([v1,v5,v3])


    return points, elements, boundary





def rectangular_cross_slit(m, n, len1=1.0, len2=1.0, origin = (0.0, 0.0)):

    """Setup a rectangular grid of triangles
    with m+1 by n+1 grid points
    and side lengths len1, len2. If side lengths are omitted
    the mesh defaults to the unit square.

    len1: x direction (left to right)
    len2: y direction (bottom to top)

    Return to lists: points and elements suitable for creating a Mesh or
    FVMesh object, e.g. Mesh(points, elements)
    """

    from anuga.config import epsilon

    delta1 = float(len1)/m
    delta2 = float(len2)/n

    #Dictionary of vertex objects
    vertices = {}
    points = []

    for i in range(m+1):
        for j in range(n+1):
            vertices[i,j] = len(points)
            points.append([delta1*i + origin[0], delta2*j  + origin[1]])

    # Construct 4 triangles per element
    elements = []
    boundary = {}
    for i in range(m):
        for j in range(n):
            v1 = vertices[i,j+1]
            v2 = vertices[i,j]
            v3 = vertices[i+1,j+1]
            v4 = vertices[i+1,j]
            x = (points[v1][0]+points[v2][0]+points[v3][0]+points[v4][0])*0.25
            y = (points[v1][1]+points[v2][1]+points[v3][1]+points[v4][1])*0.25

            # Create centre point
            v5 = len(points)
            points.append([x, y])

            #Create left triangle
            if i == 0:
                boundary[(len(elements), 1)] = 'left'
            elements.append([v2,v5,v1])

            #Create bottom triangle
            if j == 0:
                boundary[(len(elements), 1)] = 'bottom'
            elements.append([v4,v5,v2])

            #Create right triangle
            if i == m-1:
                boundary[(len(elements), 1)] = 'right'
            elements.append([v3,v5,v4])

            #Create top triangle
            if j == n-1:
                boundary[(len(elements), 1)] = 'top'
            elements.append([v1,v5,v3])


    return points, elements, boundary



def rectangular_periodic(m_g, n_g, len1_g=1.0, len2_g=1.0, origin_g = (0.0, 0.0)):


    """Setup a rectangular grid of triangles
    with m+1 by n+1 grid points
    and side lengths len1, len2. If side lengths are omitted
    the mesh defaults to the unit square, divided between all the
    processors

    len1: x direction (left to right)
    len2: y direction (bottom to top)

    """

    processor = 0
    numproc   = 1


    n = n_g
    m_low  = -1
    m_high = m_g +1

    m = m_high - m_low

    delta1 = float(len1_g)/m_g
    delta2 = float(len2_g)/n_g

    len1 = len1_g*float(m)/float(m_g)
    len2 = len2_g
    origin = ( origin_g[0]+float(m_low)/float(m_g)*len1_g, origin_g[1] )

    #Calculate number of points
    Np = (m+1)*(n+1)

    class VIndex(object):

        def __init__(self, n,m):
            self.n = n
            self.m = m

        def __call__(self, i,j):
            return j+i*(self.n+1)

    class EIndex(object):

        def __init__(self, n,m):
            self.n = n
            self.m = m

        def __call__(self, i,j):
            return 2*(j+i*self.n)


    I = VIndex(n,m)
    E = EIndex(n,m)

    points = num.zeros( (Np,2), float)

    for i in range(m+1):
        for j in range(n+1):

            points[I(i,j),:] = [i*delta1 + origin[0], j*delta2 + origin[1]]

    #Construct 2 triangles per rectangular element and assign tags to boundary
    #Calculate number of triangles
    Nt = 2*m*n


    elements = num.zeros( (Nt,3), int)
    boundary = {}
    Idgl = []
    Idfl = []
    Idgr = []
    Idfr = []

    full_send_dict = {}
    ghost_recv_dict = {}
    nt = -1
    for i in range(m):
        for j in range(n):

            i1 = I(i,j+1)
            i2 = I(i,j)
            i3 = I(i+1,j+1)
            i4 = I(i+1,j)

            #Lower Element
            nt = E(i,j)
            if i == 0:
                Idgl.append(nt)

            if i == 1:
                Idfl.append(nt)

            if i == m-2:
                Idfr.append(nt)

            if i == m-1:
                Idgr.append(nt)

            if i == m-1:
                if processor == numproc-1:
                    boundary[nt, 2] = 'right'
                else:
                    boundary[nt, 2] = 'ghost'

            if j == 0:
                boundary[nt, 1] = 'bottom'
            elements[nt,:] = [i4,i3,i2]

            #Upper Element
            nt = E(i,j)+1
            if i == 0:
                Idgl.append(nt)

            if i == 1:
                Idfl.append(nt)

            if i == m-2:
                Idfr.append(nt)

            if i == m-1:
                Idgr.append(nt)

            if i == 0:
                if processor == 0:
                    boundary[nt, 2] = 'left'
                else:
                    boundary[nt, 2] = 'ghost'
            if j == n-1:
                boundary[nt, 1] = 'top'
            elements[nt,:] = [i1,i2,i3]

    Idfl.extend(Idfr)
    Idgr.extend(Idgl)

    Idfl = num.array(Idfl, int)
    Idgr = num.array(Idgr, int)

    full_send_dict[processor]  = [Idfl, Idfl]
    ghost_recv_dict[processor] = [Idgr, Idgr]


    return  points, elements, boundary, full_send_dict, ghost_recv_dict


def oblique(m, n, lenx = 1.0, leny = 1.0, theta = 8.95, origin = (0.0, 0.0)):
    """Setup a oblique grid of triangles
    with m segments in the x-direction and n segments in the y-direction

    """

    # FIXME (Ole): Someone wrote this but didn't add a test for it. Anything could happen here
    
    import math

    from anuga.config import epsilon


    deltax = lenx/float(m)
    deltay = leny/float(n)
    a = 0.75*lenx*math.tan(theta/180.*math.pi)
    x1 = lenx
    y1 = 0
    x2 = lenx
    y2 = leny
    x3 = 0.25*lenx
    y3 = leny
    x4 = x3
    y4 = 0
    a2 = a/(x1-x4)
    a1 = -a2*x4
    a4 = ((a1 + a2*x3)/y3-(a1 + a2*x2)/y2)/(x2-x3)
    a3 = 1. - (a1 + a2*x3)/y3 - a4*x3

    # Dictionary of vertex objects
    vertices = {}
    points = []

    for i in range(m+1):
        x = deltax*i
        for j in range(n+1):
            y = deltay*j
            if x > 0.25*lenx:
                y = a1 + a2*x + a3*y + a4*x*y

            vertices[i,j] = len(points)
            points.append([x + origin[0], y + origin[1]])

    # Construct 2 triangles per element
    elements = []
    boundary = {}
    for i in range(m):
        for j in range(n):
            v1 = vertices[i,j+1]
            v2 = vertices[i,j]
            v3 = vertices[i+1,j+1]
            v4 = vertices[i+1,j]

            #Update boundary dictionary and create elements
            if i == m-1:
                boundary[(len(elements), 2)] = 'right'
            if j == 0:
                boundary[(len(elements), 1)] = 'bottom'
            elements.append([v4,v3,v2]) #Lower

            if i == 0:
                boundary[(len(elements), 2)] = 'left'
            if j == n-1:
                boundary[(len(elements), 1)] = 'top'

            elements.append([v1,v2,v3]) #Upper

    return points, elements, boundary


def circular(m, n, radius=1.0, center = (0.0, 0.0)):
    """Setup a circular grid of triangles with m concentric circles and
    with n radial segments. If radius is are omitted the mesh defaults to
    the unit circle radius.

    radius: radius of circle

    #FIXME: The triangles become degenerate for large values of m or n.
    """



    from math import pi, cos, sin

    radius = float(radius) #Ensure floating point format

    #Dictionary of vertex objects and list of points
    vertices = {}
    points = [[0.0, 0.0]] #Center point
    vertices[0, 0] = 0

    for i in range(n):
        theta = 2*i*pi/n
        x = cos(theta)
        y = sin(theta)
        for j in range(1,m+1):
            delta = j*radius/m
            vertices[i,j] = len(points)
            points.append([delta*x, delta*y])

    #Construct 2 triangles per element
    elements = []
    for i in range(n):
        for j in range(1,m):

            i1 = (i + 1) % n  #Wrap around

            v1 = vertices[i,j+1]
            v2 = vertices[i,j]
            v3 = vertices[i1,j+1]
            v4 = vertices[i1,j]

            elements.append([v4,v2,v3]) #Lower
            elements.append([v1,v3,v2]) #Upper


    #Do the center
    v1 = vertices[0,0]
    for i in range(n):
        i1 = (i + 1) % n  #Wrap around
        v2 = vertices[i,1]
        v3 = vertices[i1,1]

        elements.append([v1,v2,v3]) #center

    return points, elements


def from_polyfile(name):
    """Read mesh from .poly file, an obj like file format
    listing first vertex coordinates and then connectivity
    """

    from anuga.utilities.numerical_tools import anglediff
    from math import pi
    import os.path
    root, ext = os.path.splitext(name)

    if ext == 'poly':
        filename = name
    else:
        filename = name + '.poly'


    fid = open(filename)

    points = []    #x, y
    values = []    #z
    ##vertex_values = []    #Repeated z
    triangles = [] #v0, v1, v2

    lines = fid.readlines()

    keyword = lines[0].strip()
    msg = 'First line in .poly file must contain the keyword: POINTS'
    assert keyword == 'POINTS', msg

    offending = 0
    i = 1
    while keyword == 'POINTS':
        line = lines[i].strip()
        i += 1

        if line == 'POLYS':
            keyword = line
            break

        fields = line.split(':')
        assert int(fields[0]) == i-1, 'Point indices not consecutive'

        #Split the three floats
        xyz = fields[1].split()

        x = float(xyz[0])
        y = float(xyz[1])
        z = float(xyz[2])

        points.append([x, y])
        values.append(z)


    k = i
    while keyword == 'POLYS':
        line = lines[i].strip()
        i += 1

        if line == 'END':
            keyword = line
            break


        fields = line.split(':')
        assert int(fields[0]) == i-k, 'Poly indices not consecutive'

        #Split the three indices
        vvv = fields[1].split()

        i0 = int(vvv[0])-1
        i1 = int(vvv[1])-1
        i2 = int(vvv[2])-1

        #Check for and exclude degenerate areas
        x0 = points[i0][0]
        y0 = points[i0][1]
        x1 = points[i1][0]
        y1 = points[i1][1]
        x2 = points[i2][0]
        y2 = points[i2][1]

        area = abs((x1*y0-x0*y1)+(x2*y1-x1*y2)+(x0*y2-x2*y0))/2
        if area > 0:

            #Ensure that points are arranged in counter clock-wise order
            v0 = [x1-x0, y1-y0]
            v1 = [x2-x1, y2-y1]
            v2 = [x0-x2, y0-y2]

            a0 = anglediff(v1, v0)
            a1 = anglediff(v2, v1)
            a2 = anglediff(v0, v2)


            if a0 < pi and a1 < pi and a2 < pi:
                #all is well
                j0 = i0
                j1 = i1
                j2 = i2
            else:
                #Swap two vertices
                j0 = i1
                j1 = i0
                j2 = i2

            triangles.append([j0, j1, j2])
            ##vertex_values.append([values[j0], values[j1], values[j2]])
        else:
            offending +=1

    log.critical('Removed %d offending triangles out of %d'
                 % (offending, len(lines)))
    return points, triangles, values



def strang_mesh(filename):
    """Read Strang generated mesh.
    """

    from math import pi
    from anuga.utilities.numerical_tools import anglediff


    fid = open(filename)
    points = []    # List of x, y coordinates
    triangles = [] # List of vertex ids as listed in the file

    for line in fid.readlines():
        fields = line.split()
        if len(fields) == 2:
            # we are reading vertex coordinates
            points.append([float(fields[0]), float(fields[1])])
        elif len(fields) == 3:
            # we are reading triangle point id's (format ae+b)
            triangles.append([int(float(fields[0]))-1,
                              int(float(fields[1]))-1,
                              int(float(fields[2]))-1])
        else:
            raise Excetion('wrong format in %s' % filename)

    elements = [] #Final list of elements

    for t in triangles:
        #Get vertex coordinates
        v0 = t[0]
        v1 = t[1]
        v2 = t[2]

        x0 = points[v0][0]
        y0 = points[v0][1]
        x1 = points[v1][0]
        y1 = points[v1][1]
        x2 = points[v2][0]
        y2 = points[v2][1]

        #Check that points are arranged in counter clock-wise order
        vec0 = [x1-x0, y1-y0]
        vec1 = [x2-x1, y2-y1]
        vec2 = [x0-x2, y0-y2]

        a0 = anglediff(vec1, vec0)
        a1 = anglediff(vec2, vec1)
        a2 = anglediff(vec0, vec2)

        if a0 < pi and a1 < pi and a2 < pi:
            elements.append([v0, v1, v2])
        else:
            elements.append([v0, v2, v1])

    return points, elements

# #Map from edge number to indices of associated vertices
# edge_map = ((1,2), (0,2), (0,1))

def contracting_channel(m, n, W_upstream = 1., W_downstream = 0.75,
                            L_1 = 5.0, L_2 = 2.0, L_3 = 10, origin = (0.0, 0.0)):
    """Setup a contracting channel grid of triangles
    with m segments in the x-direction and n segments in the y-direction

    """

    import math

    from anuga.config import epsilon


    lenx = L_1 + L_2 + L_3
    leny = W_upstream
    deltax = lenx/float(m)
    deltay = leny/float(n)

    x1 = 0
    y1 = 0
    x2 = L_1
    y2 = 0
    x3 = L_1 + L_2
    y3 = (W_upstream - W_downstream)/2
    x4 = L_1 + L_2 + L_3
    y4 = y3
    x5 = x4
    y5 = y4 + W_downstream
    x6 = L_1 + L_2
    y6 = y5
    x7 = L_1
    y7 = W_upstream
    x8 = 0
    y8 = W_upstream
    a1 = 0
    a2 = (W_upstream - W_downstream)/(2*L_2)
    a3 = 1
    a4 = (W_downstream - W_upstream)/(L_2*W_upstream)

    # Dictionary of vertex objects
    vertices = {}
    points = []

    for i in range(m+1):
        x = deltax*i
        for j in range(n+1):
            y = deltay*j
            if x > L_1 and x <= (L_1 + L_2):
                y = a1 + a2*(x - L_1) + a3*y + a4*(x - L_1)*y
            elif x > L_1 + L_2:
                y = (W_upstream - W_downstream)/2 + deltay*j*W_downstream/W_upstream

            vertices[i,j] = len(points)
            points.append([x + origin[0], y + origin[1]])

    # Construct 2 triangles per element
    elements = []
    boundary = {}
    for i in range(m):
        for j in range(n):
            v1 = vertices[i,j+1]
            v2 = vertices[i,j]
            v3 = vertices[i+1,j+1]
            v4 = vertices[i+1,j]

            #Update boundary dictionary and create elements
            if i == m-1:
                boundary[(len(elements), 2)] = 'right'
            if j == 0:
                boundary[(len(elements), 1)] = 'bottom'
            elements.append([v4,v3,v2]) #Lower

            if i == 0:
                boundary[(len(elements), 2)] = 'left'
            if j == n-1:
                boundary[(len(elements), 1)] = 'top'

            elements.append([v1,v2,v3]) #Upper

    return points, elements, boundary


def contracting_channel_cross(m, n, W_upstream = 1., W_downstream = 0.75,
                              L_1 = 5.0, L_2 = 2.0, L_3 = 10, origin = (0.0, 0.0)):
    """Setup a contracting channel grid of triangles
    with m segments in the x-direction and n segments in the y-direction

    """

    import math

    from anuga.config import epsilon


    lenx = L_1 + L_2 + L_3
    leny = W_upstream
    deltax = lenx/float(m)
    deltay = leny/float(n)

    x1 = 0
    y1 = 0
    x2 = L_1
    y2 = 0
    x3 = L_1 + L_2
    y3 = (W_upstream - W_downstream)/2
    x4 = L_1 + L_2 + L_3
    y4 = y3
    x5 = x4
    y5 = y4 + W_downstream
    x6 = L_1 + L_2
    y6 = y5
    x7 = L_1
    y7 = W_upstream
    x8 = 0
    y8 = W_upstream
    a1 = 0
    a2 = (W_upstream - W_downstream)/(2*L_2)
    a3 = 1
    a4 = (W_downstream - W_upstream)/(L_2*W_upstream)

    # Dictionary of vertex objects
    vertices = {}
    points = []

    for i in range(m+1):
        x = deltax*i
        for j in range(n+1):
            y = deltay*j
            if x > L_1 and x <= (L_1 + L_2):
                y = a1 + a2*(x - L_1) + a3*y + a4*(x - L_1)*y
            elif x > L_1 + L_2:
                y = (W_upstream - W_downstream)/2 + deltay*j*W_downstream/W_upstream

            vertices[i,j] = len(points)
            points.append([x + origin[0], y + origin[1]])

    # Construct 4 triangles per element
    elements = []
    boundary = {}
    for i in range(m):
        for j in range(n):
            v1 = vertices[i,j+1]
            v2 = vertices[i,j]
            v3 = vertices[i+1,j+1]
            v4 = vertices[i+1,j]
            x = (points[v1][0]+points[v2][0]+points[v3][0]+points[v4][0])*0.25
            y = (points[v1][1]+points[v2][1]+points[v3][1]+points[v4][1])*0.25
            v5 = len(points)
            points.append([x, y])

            #Create left triangle
            if i == 0:
                boundary[(len(elements), 1)] = 'left'
            elements.append([v2,v5,v1])

            #Create bottom triangle
            if j == 0:
                boundary[(len(elements), 1)] = 'bottom'
            elements.append([v4,v5,v2])

            #Create right triangle
            if i == m-1:
                boundary[(len(elements), 1)] = 'right'
            elements.append([v3,v5,v4])

            #Create top triangle
            if j == n-1:
                boundary[(len(elements), 1)] = 'top'
            elements.append([v1,v5,v3])


    return points, elements, boundary




def oblique_cross(m, n, lenx = 1.0, leny = 1.0, theta = 8.95, origin = (0.0, 0.0)):
    """Setup a oblique grid of triangles
    with m segments in the x-direction and n segments in the y-direction

    """

    import math

    from anuga.config import epsilon


    deltax = lenx/float(m)
    deltay = leny/float(n)
    a = 0.75*lenx*math.tan(theta/180.*math.pi)
    x1 = lenx
    y1 = 0
    x2 = lenx
    y2 = leny
    x3 = 0.25*lenx
    y3 = leny
    x4 = x3
    y4 = 0
    a2 = a/(x1-x4)
    a1 = -a2*x4
    a4 = ((a1 + a2*x3)/y3-(a1 + a2*x2)/y2)/(x2-x3)
    a3 = 1. - (a1 + a2*x3)/y3 - a4*x3

    # Dictionary of vertex objects
    vertices = {}
    points = []

    for i in range(m+1):
        x = deltax*i
        for j in range(n+1):
            y = deltay*j
            if x > 0.25*lenx:
                y = a1 + a2*x + a3*y + a4*x*y

            vertices[i,j] = len(points)
            points.append([x + origin[0], y + origin[1]])

    # Construct 4 triangles per element
    elements = []
    boundary = {}
    for i in range(m):
        for j in range(n):
            v1 = vertices[i,j+1]
            v2 = vertices[i,j]
            v3 = vertices[i+1,j+1]
            v4 = vertices[i+1,j]
            x = (points[v1][0]+points[v2][0]+points[v3][0]+points[v4][0])*0.25
            y = (points[v1][1]+points[v2][1]+points[v3][1]+points[v4][1])*0.25
            v5 = len(points)
            points.append([x, y])

            #Update boundary dictionary and create elements
                        #Create left triangle
            if i == 0:
                boundary[(len(elements), 1)] = 'left'
            elements.append([v2,v5,v1])

            #Create bottom triangle
            if j == 0:
                boundary[(len(elements), 1)] = 'bottom'
            elements.append([v4,v5,v2])

            #Create right triangle
            if i == m-1:
                boundary[(len(elements), 1)] = 'right'
            elements.append([v3,v5,v4])

            #Create top triangle
            if j == n-1:
                boundary[(len(elements), 1)] = 'top'
            elements.append([v1,v5,v3])


    return points, elements, boundary
