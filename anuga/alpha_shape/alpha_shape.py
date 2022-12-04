"""Alpha shape
Determine the shape of a set of points.

From website by Kaspar Fischer:
As mentionned in Edelsbrunner's and Muecke's paper, one can
intuitively think of an alpha-shape as the following:

Imagine a huge mass of ice-cream making up the space and containing
the points S as ``hard'' chocolate pieces. Using one of these
sphere-formed ice-cream spoons we carve out all parts of the ice-cream
block we can reach without bumping into chocolate pieces, even
carving out holes in the inside (eg. parts not reachable by simply
moving the spoon from the outside). We will eventually end up with a
(not necessarily convex) object bounded by caps, arcs and points. If
we now straighten all ``round'' faces to triangles and line segments,
we have an intuitive description of what is called the alpha-shape.

Author: Vanessa Robins, ANU
"""


from builtins import filter
from builtins import range
from builtins import object
import random

from anuga.load_mesh.loadASCII import export_boundary_file
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.utilities import log 

import numpy as num

# Python 2.7 Hack
try:
    from exceptions import Exception
except:
    pass
class AlphaError(Exception):pass
class PointError(AlphaError): pass
class FlagError(AlphaError): pass

OUTPUT_FILE_TITLE = "# The alpha shape boundary defined by point index pairs of edges"
INF = pow(10,9)
EPSILON = 1.0e-12

def alpha_shape_via_files(point_file, boundary_file, alpha= None):
    """
    Load a point file and return the alpha shape boundary as a boundary file.
    
    Inputs:
    point_file: File location of the input file, points format (.csv or .pts)
    boundary_file: File location of the generated output file
    alpha: The alpha value can be optionally specified.  If it is not specified
    the optimum alpha value will be used.
    """
    geospatial = Geospatial_data(point_file)
    points = geospatial.get_data_points(absolute=False)
    
    AS = Alpha_Shape(points, alpha)
    AS.write_boundary(boundary_file)
    

class Alpha_Shape(object):

    def __init__(self, points, alpha = None):
        """
        An Alpha_Shape requires input of a set of points. Other class routines
        return the alpha shape boundary.
        
        Inputs:       
          points: List of coordinate pairs [[x1, y1],[x2, y2]..] 
          alpha: The alpha value can be optionally specified.  If it is
          not specified the optimum alpha value will be used.
        """
        self._set_points(points)
        self._alpha_shape_algorithm(alpha)

    def _set_points(self, points):
        """
        Create self.points array, do Error checking
        Inputs:       
          points: List of coordinate pairs [[x1, y1],[x2, y2]..] 
        """
        if len (points) <= 2:
            raise PointError("Too few points to find an alpha shape")
        if len(points)==3:
            #check not in a straight line
            # FIXME check points 1,2,3 if straingt, check if points 2,3,4, ect
            x01 = points[0][0] - points[1][0]
            y01 = points[0][1] - points[1][1]
            x12 = points[1][0] - points[2][0]
            y12 = points[1][1] - points[2][1]
            crossprod = x01*y12 - x12*y01
            if crossprod==0:
                raise PointError("Three points on a straight line")
        
        #Convert input to numeric arrays
        self.points = num.array(points, float)

    
    def write_boundary(self,file_name):
        """
        Write the boundary to a file
        """
        export_boundary_file(file_name, self.get_boundary(),
                             OUTPUT_FILE_TITLE, delimiter = ',')
    
    def get_boundary(self):
        """
        Return a list of tuples.
        Each tuple represents a segment in the boundary
        by the index of its two end points.
        The list of tuples represents the alpha shape boundary.
        """
        return self.boundary

    def set_boundary_type(self,raw_boundary=True,
                          remove_holes=False,
                          smooth_indents=False,
                          expand_pinch=False,
                          boundary_points_fraction=0.2):
        """
        Use the flags to set constraints on the boundary:
        raw_boundary    Return raw boundary i.e. the regular edges of the
                        alpha shape. 
        remove_holes    filter to remove small holes
                        (small is defined by  boundary_points_fraction ) 
        smooth_indents  remove sharp triangular indents in boundary
        expand_pinch    test for pinch-off and correct 
                           i.e. a boundary vertex with more than two edges.
        """

        if raw_boundary:
            # reset alpha shape boundary 
            reg_edge = self.get_regular_edges(self.alpha)
            self.boundary = [self.edge[k] for k in reg_edge]
            self._init_boundary_triangles()
        if remove_holes:
            #remove small holes
            self.boundary = self._remove_holes(boundary_points_fraction)
        if smooth_indents:
            #remove sharp triangular indents
            self.boundary = self._smooth_indents()
        if expand_pinch:
            #deal with pinch-off
            self.boundary = self._expand_pinch()
        

    def get_delaunay(self):
        """
        """
        return self.deltri

    def get_optimum_alpha(self):
        """
        """
        return self.optimum_alpha

    def get_alpha(self):
        """
        Return current alpha value.
        """
        return self.alpha

    def set_alpha(self,alpha):
        """
        Set alpha and update alpha-boundary. 
        """
        self.alpha = alpha
        reg_edge = self.get_regular_edges(alpha)    
        self.boundary = [self.edge[k] for k in reg_edge]
        self._init_boundary_triangles()
    
            
    def _alpha_shape_algorithm(self, alpha=None):
        """
        Given a set of points (self.points) and an optional alpha value
        determines the alpha shape boundary (stored in self.boundary,
        accessed by get_boundary).
        
        Inputs:
          alpha: The alpha value can be optionally specified.  If it is
          not specified the optimum alpha value will be used.
        """

        self.alpha = alpha

        ## Build Delaunay triangulation
        from anuga.mesh_engine.mesh_engine import generate_mesh
        points = []
        seglist = []
        holelist = []
        regionlist = []
        pointattlist = []
        segattlist = []

        points = [(pt[0], pt[1]) for pt in self.points]
        pointattlist = [ [] for i in range(len(points)) ] 
        mode = "Qzcn"
        tridata = generate_mesh(points,seglist,holelist,regionlist,
                                 pointattlist,segattlist,mode)
        self.deltri = tridata['generatedtrianglelist']
        self.deltrinbr = tridata['generatedtriangleneighborlist']
        self.hulledges = tridata['generatedsegmentlist'] 

        ## Build Alpha table
        ## the following routines determine alpha thresholds for the
        ## triangles, edges, and vertices of the delaunay triangulation
        self._tri_circumradius()
        self._edge_intervals()
        self._vertex_intervals()

        if alpha==None:
            # Find optimum alpha
            # Ken Clarkson's hull program uses smallest alpha so that
            # every vertex is non-singular so... 
            self.optimum_alpha = max([iv[0] for iv in self.vertexinterval \
                                      if iv!=[] ])
            alpha = self.optimum_alpha
        self.alpha = alpha
        reg_edge = self.get_regular_edges(self.alpha)
        self.boundary = [self.edge[k] for k in reg_edge]
        self._init_boundary_triangles()
        
        return

    def _tri_circumradius(self):
        """
        Compute circumradii of the delaunay triangles
        """
        
        x = self.points[:,0]
        y = self.points[:,1]
        ind1 = [self.deltri[j][0] for j in range(len(self.deltri))]
        ind2 = [self.deltri[j][1] for j in range(len(self.deltri))]
        ind3 = [self.deltri[j][2] for j in range(len(self.deltri))]

        x1 = num.array([x[j] for j in ind1])
        y1 = num.array([y[j] for j in ind1])
        x2 = num.array([x[j] for j in ind2])
        y2 = num.array([y[j] for j in ind2])
        x3 = num.array([x[j] for j in ind3])
        y3 = num.array([y[j] for j in ind3])

        x21 = x2-x1
        x31 = x3-x1
        y21 = y2-y1
        y31 = y3-y1

        dist21 = x21*x21 + y21*y21
        dist31 = x31*x31 + y31*y31

        denom = x21*y31 - x31*y21

        # dx/2, dy/2 give circumcenter relative to x1,y1. 
        # dx = (y31*dist21 - y21*dist31)/denom
        # dy = (x21*dist31 - x31*dist21)/denom
        # first need to check for near-zero values of denom 
        delta = 0.00000001
        zeroind = [k for k in range(len(denom)) if \
                   (denom[k]< EPSILON and  denom[k] > -EPSILON)]
        # if some denom values are close to zero,
        # we perturb the associated vertices and recalculate 
        while zeroind!=[]:
            random.seed()
            log.critical("Warning: degenerate triangles found in alpha_shape.py, results may be inaccurate.")
            for d in zeroind:
                x1[d] = x1[d]+delta*(random.random()-0.5)
                x2[d] = x2[d]+delta*(random.random()-0.5)
                x3[d] = x3[d]+delta*(random.random()-0.5)
                y1[d] = y1[d]+delta*(random.random()-0.5)
                y2[d] = y2[d]+delta*(random.random()-0.5)
                y3[d] = y3[d]+delta*(random.random()-0.5)
            x21 = x2-x1
            x31 = x3-x1
            y21 = y2-y1
            y31 = y3-y1
            dist21 = x21*x21 + y21*y21
            dist31 = x31*x31 + y31*y31
            denom = x21*y31 - x31*y21
            zeroind = [k for k in range(len(denom)) if \
                       (denom[k]< EPSILON and  denom[k] > -EPSILON)]

        if num.any(denom == 0.0):
            raise AlphaError

        dx = num.divide(y31*dist21 - y21*dist31, denom)
        dy = num.divide(x21*dist31 - x31*dist21, denom)

        self.triradius = 0.5*num.sqrt(dx*dx + dy*dy)

    def _edge_intervals(self):
        """
         for each edge, find triples
         (length/2, min_adj_triradius, max_adj_triradius) if unattached
         (min_adj_triradius, min_adj_triradius, max_adj_triradius) if attached.
         An edge is attached if it is opposite an obtuse angle 
        """
        
        # It should be possible to rewrite this routine in an array-friendly
        # form like _tri_circumradius()  if we need to speed things up.
        # Hard to do though.

        edges = []
        edgenbrs = []
        edgeinterval = []
        for t in range(len(self.deltri)):
            tri = self.deltri[t]
            trinbr = self.deltrinbr[t]
            dx = num.array([self.points[tri[(i+1)%3],0] -
                            self.points[tri[(i+2)%3],0] for i in [0,1,2]])
            dy = num.array([self.points[tri[(i+1)%3],1] -
                            self.points[tri[(i+2)%3],1] for i in [0,1,2]])
            elen = num.sqrt(dx*dx+dy*dy)
            # really only need sign - not angle value:
            anglesign = num.array([(-dx[(i+1)%3]*dx[(i+2)%3]-
                                    dy[(i+1)%3]*dy[(i+2)%3]) for i in [0,1,2]])
            
            for i in [0,1,2]:
                j = (i+1)%3
                k = (i+2)%3
                if trinbr[i]==-1:
                    edges.append((tri[j], tri[k]))
                    edgenbrs.append((t, -1))
                    edgeinterval.append([0.5*elen[i], self.triradius[t], INF])
                elif (tri[j]<tri[k]):
                    edges.append((tri[j], tri[k]))
                    edgenbrs.append((t, trinbr[i]))
                    edgeinterval.append([0.5*elen[i],\
                       min(self.triradius[t],self.triradius[trinbr[i]]),\
                          max(self.triradius[t],self.triradius[trinbr[i]]) ])
                else:
                    continue
                if anglesign[i] < 0: 
                    edgeinterval[-1][0] = edgeinterval[-1][1]
                    
        self.edge = edges
        self.edgenbr = edgenbrs
        self.edgeinterval = edgeinterval

    def _vertex_intervals(self):
        """
        for each vertex find pairs
        (min_adj_triradius, max_adj_triradius)
        """
        nv = len(self.points)
        vertexnbrs = [ [] for i in range(nv)]
        vertexinterval = [ [] for i in range(nv)] 
        for t in range(len(self.deltri)):
            for j in self.deltri[t]:
                vertexnbrs[int(j)].append(t)
        for h in range(len(self.hulledges)):
            for j in self.hulledges[h]:
                vertexnbrs[int(j)].append(-1)
        
        for i in range(nv):
            radii = [ self.triradius[t] for t in vertexnbrs[i] if t>=0 ]
            try:
                vertexinterval[i] = [min(radii), max(radii)]
                if vertexnbrs[i][-1]==-1:
                    vertexinterval[i][1]=INF
            except ValueError:
                raise AlphaError

        self.vertexnbr = vertexnbrs
        self.vertexinterval = vertexinterval
        
    def get_alpha_triangles(self,alpha):
        """
        Given an alpha value,
        return indices of triangles in the alpha-shape
        """
        def tri_rad_lta(k):
            return self.triradius[k]<=alpha

        return list(filter(tri_rad_lta, list(range(len(self.triradius)))))

    def get_regular_edges(self,alpha):
        """
        Given an alpha value,
        return the indices of edges on the boundary of the alpha-shape
        """
        def reg_edge(k):
            return self.edgeinterval[k][1]<=alpha and \
                   self.edgeinterval[k][2]>alpha

        return list(filter(reg_edge, list(range(len(self.edgeinterval)))))

    def get_exposed_vertices(self,alpha):
        """
        Given an alpha value,
        return the vertices on the boundary of the alpha-shape
        """
        def exp_vert(k):
            return self.vertexinterval[k][0]<=alpha and \
                   self.vertexinterval[k][1]>alpha

        return list(filter(exp_vert, list(range(len(self.vertexinterval)))))        

    def _vertices_from_edges(self,elist):
        """
        Returns the list of unique vertex labels from edges in elist
        """

        v1 = [elist[k][0] for k in range(len(elist))]
        v2 = [elist[k][1] for k in range(len(elist))]
        v = v1+v2
        v.sort()
        vertices = [v[k] for k in range(len(v)) if v[k]!=v[k-1] ]
        return vertices

    def _init_boundary_triangles(self):
        """
        Creates the initial list of triangle indices
        exterior to and touching the boundary of the alpha shape
        """
        def tri_rad_gta(k):
            return self.triradius[k]>self.alpha

        extrind = list(filter(tri_rad_gta, list(range(len(self.triradius)))))

        bv = self._vertices_from_edges(self.boundary)
        
        btri = []
        for et in extrind:
            v0 = self.deltri[et][0]
            v1 = self.deltri[et][1]
            v2 = self.deltri[et][2]
            if v0 in bv or v1 in bv or v2 in bv:
                btri.append(et)

        self.boundarytriangle = btri
  
       
    def _remove_holes(self,small):
        """
        Given the edges in self.boundary, finds the largest components.
        The routine does this by implementing a union-find algorithm. 
        """

        bdry = self.boundary
        
        def findroot(i):
            if vptr[i] < 0:
                return i
            k = findroot(vptr[i])
            vptr[i] = k    # this produces "path compression" in the
                           # union-find tree. 
            return k



        # get a list of unique vertex labels:
        verts = self._vertices_from_edges(bdry)

        # vptr represents the union-find tree.
        # if vptr[i] = EMPTY < 0, the vertex verts[i] has not been visited yet
        # if vptr[i] = j > 0, then j verts[j] is the parent of verts[i].
        # if vptr[i] = n < 0, then verts[i] is a root vertex and
        #                       represents a connected component of n vertices.
        
        #initialise vptr to negative number outside range 
        EMPTY = -max(verts)-len(verts)
        vptr = [EMPTY for k in range(len(verts))]

        #add edges and maintain union tree
        for i in range(len(bdry)):
            vl = verts.index(bdry[i][0])
            vr = verts.index(bdry[i][1])
            rvl = findroot(vl)
            rvr = findroot(vr)
            if not(rvl==rvr):
                if (vptr[vl]==EMPTY):
                    if (vptr[vr]==EMPTY):
                        vptr[vl] = -2
                        vptr[vr] = vl
                    else:
                        vptr[vl] = rvr
                        vptr[rvr] = vptr[rvr]-1
                else:
                    if (vptr[vr]==EMPTY):
                        vptr[vr] = rvl
                        vptr[rvl] = vptr[rvl]-1
                    else:
                        if vptr[rvl] > vptr[rvr]:
                            vptr[rvr] = vptr[rvr] + vptr[rvl]
                            vptr[rvl] = rvr
                            vptr[vl] = rvr
                        else:
                            vptr[rvl] = vptr[rvl] + vptr[rvr]
                            vptr[rvr] = rvl
                            vptr[vr] = rvl 
        # end edge loop

        if vptr.count(EMPTY):
            raise FlagError("We didn't hit all the vertices in the boundary")
        
        # discard the edges in the little components
        # (i.e. those components with less than 'small' fraction of bdry points)
        cutoff = round(small*len(verts))
        largest_component = -min(vptr)
        if cutoff > largest_component:
            cutoff = round((1-small)*largest_component)

        # littleind has root indices for small components 
        littleind = [k for k in range(len(vptr)) if \
                     (vptr[k]<0 and vptr[k]>-cutoff)] 
        if littleind:
            # littlecomp has all vptr indices in the small components
            littlecomp = [k for k in range(len(vptr)) if \
                          findroot(k) in littleind]
            # vdiscard has the vertex indices corresponding to vptr indices  
            vdiscard = [verts[k] for k in littlecomp] 
            newbdry = [e for e in bdry if \
                       not((e[0] in vdiscard) and (e[1] in vdiscard))]

            newverts = self._vertices_from_edges(newbdry)
            # recompute the boundary triangles
            newbt = []
            for bt in self.boundarytriangle:
                v0 = self.deltri[bt][0]
                v1 = self.deltri[bt][1]
                v2 = self.deltri[bt][2]
                if (v0 in newverts or v1 in newverts or v2 in newverts):
                    newbt.append(bt)

            self.boundarytriangle = newbt           
        else:
            newbdry = bdry
    
        return newbdry


    def _smooth_indents(self):
        """
        Given edges in bdry, test for acute-angle triangular indents
        and remove them.
        """
        
        bdry = self.boundary
        bdrytri = self.boundarytriangle
        
        # find boundary triangles that have two edges in bdry
        # v2ind has the place index relative to the triangle deltri[ind]  
        # for the bdry vertex where the two edges meet.

        verts = self._vertices_from_edges(bdry)

        b2etri = []
        for ind in bdrytri:
            bect = 0
            v2ind = [0,1,2]
            for j in [0,1,2]:
                eda = (self.deltri[ind][(j+1)%3], self.deltri[ind][(j+2)%3])
                edb = (self.deltri[ind][(j+2)%3], self.deltri[ind][(j+1)%3])
                if eda in bdry or edb in bdry:
                    bect +=1
                    v2ind.remove(j)
            if bect==2:
                b2etri.append((ind,v2ind[0]))

        # test the bdrytri triangles for acute angles 
        acutetri = []
        for tind in b2etri:
            tri = self.deltri[tind[0]]
            
            dx = num.array([self.points[tri[(i+1)%3],0] - \
                           self.points[tri[(i+2)%3],0] for i in [0,1,2]])
            dy = num.array([self.points[tri[(i+1)%3],1] - \
                           self.points[tri[(i+2)%3],1] for i in [0,1,2]])
            anglesign = num.array([(-dx[(i+1)%3]*dx[(i+2)%3]-\
                                   dy[(i+1)%3]*dy[(i+2)%3]) for i in [0,1,2]])
            # record any triangle that has an acute angle spanned by
            #two edges along the boundary..
            if anglesign[tind[1]] > 0:
                acutetri.append(tind[0])

        # adjust the bdry edges and triangles by adding
        #in the acutetri triangles
        for pind in acutetri:
            bdrytri.remove(pind) 
            tri = self.deltri[pind]
            for i in [0,1,2]:
                bdry.append((tri[(i+1)%3], tri[(i+2)%3]))

        newbdry = []
        for ed in bdry:
            numed = bdry.count(ed)+bdry.count((ed[1],ed[0]))
            if numed%2 == 1:
                newbdry.append(ed)
        
        return newbdry

    def _expand_pinch(self):
        """
        Given edges in bdry, test for vertices with more than 2 incident edges.
        Expand by adding back in associated triangles.
        """

        bdry = self.boundary
        bdrytri = self.boundarytriangle
        
        v1 = [bdry[k][0] for k in range(len(bdry))]
        v2 = [bdry[k][1] for k in range(len(bdry))]
        v = v1+v2
        v.sort()
        probv = [v[k] for k in range(len(v)) \
                 if (v[k]!=v[k-1] and v.count(v[k])>2) ]

        # find boundary triangles that have at least one vertex in probv
        probtri = []
        for ind in bdrytri:
            v0 = self.deltri[ind][0]
            v1 = self.deltri[ind][1]
            v2 = self.deltri[ind][2]
            if v0 in probv or v1 in probv or v2 in probv:
                probtri.append(ind)


        # "add in" the problem triangles
        for pind in probtri:
            bdrytri.remove(pind) 
            tri = self.deltri[pind]
            for i in [0,1,2]:
                bdry.append((tri[(i+1)%3], tri[(i+2)%3]))

        newbdry = []
        for ed in bdry:
            numed = bdry.count(ed)+bdry.count((ed[1],ed[0]))
            if numed%2 == 1:
                newbdry.append(ed)
        
        return newbdry


#-------------------------------------------------------------
if __name__ == "__main__":
    """
    Load in a data point file.
    Determine the alpha shape boundary
    Save the boundary to a file.

    usage: alpha_shape.py point_file.csv boundary_file.bnd [alpha]
    
    The alpha value is optional.
    """
    
    import os, sys
    usage = "usage: %s point_file.csv boundary_file.bnd [alpha]"%os.path.basename(sys.argv[0])
    if len(sys.argv) < 3:
        print(usage)
    else:
        point_file = sys.argv[1]
        boundary_file = sys.argv[2]
        if len(sys.argv) > 4:
            alpha = sys.argv[3]
        else:
            alpha = None

        #print "about to call alpha shape routine \n"
        alpha_shape_via_files(point_file, boundary_file, alpha)
        
