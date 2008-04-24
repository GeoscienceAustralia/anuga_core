"""parallel-meshes -
2D triangular domains for parallel finite-volume computations of
the advection equation, with extra structures to define the
sending and receiving communications define in dictionaries
full_send_dict and ghost_recv_dict


Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
Geoscience Australia, 2005

Modified by Linda Stals, March 2006, to include ghost boundaries

"""


import sys
from Numeric import array, zeros, Float, Int, ones, sum

import pypar

from anuga.config import epsilon





def parallel_rectangle(m_g, n_g, len1_g=1.0, len2_g=1.0, origin_g = (0.0, 0.0)):


    """Setup a rectangular grid of triangles
    with m+1 by n+1 grid points
    and side lengths len1, len2. If side lengths are omitted
    the mesh defaults to the unit square, divided between all the
    processors

    len1: x direction (left to right)
    len2: y direction (bottom to top)

    """

    processor = pypar.rank()
    numproc   = pypar.size()

    m_low, m_high = pypar.balance(m_g, numproc, processor)
    
    n = n_g
    m_low  = m_low-1
    m_high = m_high+1
    m = m_high - m_low

    delta1 = float(len1_g)/m_g
    delta2 = float(len2_g)/n_g

    len1 = len1_g*float(m)/float(m_g)
    len2 = len2_g
    origin = ( origin_g[0]+float(m_low)/float(m_g)*len1_g, origin_g[1] )

    #Calculate number of points
    Np = (m+1)*(n+1)

    class VIndex:

        def __init__(self, n,m):
            self.n = n
            self.m = m

        def __call__(self, i,j):
            return j+i*(self.n+1)

    class EIndex:

        def __init__(self, n,m):
            self.n = n
            self.m = m

        def __call__(self, i,j):
            return 2*(j+i*self.n)


    I = VIndex(n,m)
    E = EIndex(n,m)

    points = zeros( (Np,2), Float)

    for i in range(m+1):
        for j in range(n+1):

            points[I(i,j),:] = [i*delta1 + origin[0], j*delta2 + origin[1]]

    #Construct 2 triangles per rectangular element and assign tags to boundary
    #Calculate number of triangles
    Nt = 2*m*n


    elements = zeros( (Nt,3), Int)
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

    if numproc==1:
        Idfl.extend(Idfr)
        Idgr.extend(Idgl)
        Idfl = array(Idfl,Int)
        Idgr = array(Idgr,Int)
        full_send_dict[processor]  = [Idfl, Idfl]
        ghost_recv_dict[processor] = [Idgr, Idgr]
    elif numproc == 2:
        Idfl.extend(Idfr)
        Idgr.extend(Idgl)
        Idfl = array(Idfl,Int)
        Idgr = array(Idgr,Int)
        full_send_dict[(processor-1)%numproc]  = [Idfl, Idfl]
        ghost_recv_dict[(processor-1)%numproc] = [Idgr, Idgr]
    else:
        Idfl = array(Idfl,Int)
        Idgl = array(Idgl,Int)

        Idfr = array(Idfr,Int)
        Idgr = array(Idgr,Int)

        full_send_dict[(processor-1)%numproc]  = [Idfl, Idfl]
        ghost_recv_dict[(processor-1)%numproc] = [Idgl, Idgl]
        full_send_dict[(processor+1)%numproc]  = [Idfr, Idfr]
        ghost_recv_dict[(processor+1)%numproc] = [Idgr, Idgr]



        
    
    return  points, elements, boundary, full_send_dict, ghost_recv_dict



def rectangular_periodic(m, n, len1=1.0, len2=1.0, origin = (0.0, 0.0)):


    """Setup a rectangular grid of triangles
    with m+1 by n+1 grid points
    and side lengths len1, len2. If side lengths are omitted
    the mesh defaults to the unit square.

    len1: x direction (left to right)
    len2: y direction (bottom to top)

    Return to lists: points and elements suitable for creating a Mesh or
    FVMesh object, e.g. Mesh(points, elements)
    """

    delta1 = float(len1)/m
    delta2 = float(len2)/n

    #Calculate number of points
    Np = (m+1)*(n+1)

    class VIndex:

        def __init__(self, n,m):
            self.n = n
            self.m = m

        def __call__(self, i,j):
            return j+i*(self.n+1)

    class EIndex:

        def __init__(self, n,m):
            self.n = n
            self.m = m

        def __call__(self, i,j):
            return 2*(j+i*self.n)


    I = VIndex(n,m)
    E = EIndex(n,m)

    points = zeros( (Np,2), Float)

    for i in range(m+1):
        for j in range(n+1):

            points[I(i,j),:] = [i*delta1 + origin[0], j*delta2 + origin[1]]

    #Construct 2 triangles per rectangular element and assign tags to boundary
    #Calculate number of triangles
    Nt = 2*m*n


    elements = zeros( (Nt,3), Int)
    boundary = {}
    ghosts = {}
    nt = -1
    for i in range(m):
        for j in range(n):

            i1 = I(i,j+1)
            i2 = I(i,j)
            i3 = I(i+1,j+1)
            i4 = I(i+1,j)

            #Lower Element
            nt = E(i,j)
            if i == m-1:
                ghosts[nt] = E(1,j)
            if i == 0:
                ghosts[nt] = E(m-2,j)

            if j == n-1:
                ghosts[nt] = E(i,1)

            if j == 0:
                ghosts[nt] = E(i,n-2)

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
            if i == m-1:
                ghosts[nt] = E(1,j)+1
            if i == 0:
                ghosts[nt] = E(m-2,j)+1

            if j == n-1:
                ghosts[nt] = E(i,1)+1

            if j == 0:
                ghosts[nt] = E(i,n-2)+1

            if i == 0:
                if processor == 0:
                    boundary[nt, 2] = 'left'
                else:
                    boundary[nt, 2] = 'ghost'
            if j == n-1:
                boundary[nt, 1] = 'top'
            elements[nt,:] = [i1,i2,i3]

    #bottom left
    nt = E(0,0)
    nf = E(m-2,n-2)
    ghosts[nt]   = nf
    ghosts[nt+1] = nf+1

    #bottom right
    nt = E(m-1,0)
    nf = E(1,n-2)
    ghosts[nt]   = nf
    ghosts[nt+1] = nf+1

    #top left
    nt = E(0,n-1)
    nf = E(m-2,1)
    ghosts[nt]   = nf
    ghosts[nt+1] = nf+1

    #top right
    nt = E(m-1,n-1)
    nf = E(1,1)
    ghosts[nt]   = nf
    ghosts[nt+1] = nf+1

    return points, elements, boundary, ghosts

def rectangular_periodic_lr(m, n, len1=1.0, len2=1.0, origin = (0.0, 0.0)):


    """Setup a rectangular grid of triangles
    with m+1 by n+1 grid points
    and side lengths len1, len2. If side lengths are omitted
    the mesh defaults to the unit square.

    len1: x direction (left to right)
    len2: y direction (bottom to top)

    Return to lists: points and elements suitable for creating a Mesh or
    Domain object, e.g. Mesh(points, elements)
    """

    delta1 = float(len1)/m
    delta2 = float(len2)/n

    #Calculate number of points
    Np = (m+1)*(n+1)

    class VIndex:

        def __init__(self, n,m):
            self.n = n
            self.m = m

        def __call__(self, i,j):
            return j+i*(self.n+1)

    class EIndex:

        def __init__(self, n,m):
            self.n = n
            self.m = m

        def __call__(self, i,j):
            return 2*(j+i*self.n)


    I = VIndex(n,m)
    E = EIndex(n,m)

    points = zeros( (Np,2), Float)

    for i in range(m+1):
        for j in range(n+1):

            points[I(i,j),:] = [i*delta1 + origin[0], j*delta2 + origin[1]]

    #Construct 2 triangles per rectangular element and assign tags to boundary
    #Calculate number of triangles
    Nt = 2*m*n


    elements = zeros( (Nt,3), Int)
    boundary = {}
    ghosts = {}
    nt = -1
    for i in range(m):
        for j in range(n):

            i1 = I(i,j+1)
            i2 = I(i,j)
            i3 = I(i+1,j+1)
            i4 = I(i+1,j)

            #Lower Element
            nt = E(i,j)
            if i == m-1:
                ghosts[nt] = E(1,j)
            if i == 0:
                ghosts[nt] = E(m-2,j)

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
            if i == m-1:
                ghosts[nt] = E(1,j)+1
            if i == 0:
                ghosts[nt] = E(m-2,j)+1

            if i == 0:
                if processor == 0:
                    boundary[nt, 2] = 'left'
                else:
                    boundary[nt, 2] = 'ghost'
            if j == n-1:
                boundary[nt, 1] = 'top'
            elements[nt,:] = [i1,i2,i3]


    return points, elements, boundary, ghosts
