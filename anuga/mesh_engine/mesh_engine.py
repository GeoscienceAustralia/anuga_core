#!/usr/bin/env python

import sys

from types import ListType, TupleType

import exceptions

class NoTrianglesError(exceptions.Exception): pass
import anuga.mesh_engine.mesh_engine_c_layer as triang
#import anuga.mesh_engine.list_dic as triang

import numpy as num

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.anuga_exceptions import ANUGAError
    
def generate_mesh(points=None,
                  segments=None,holes=None,regions=None,
                  pointatts=None,segatts=None,
                  mode=None, dummy_test=None):
    """
    pointatts can be a list of lists.

    generatedtriangleattributelist is used to represent tagged regions.
    #FIXME (DSG-DSG): add comments
    """
    #FIXME (DSG-DSG): Catch parameters that are lists,
    #instead of lists of lists
    # check shape[1] is 2 etc

    if points is None:
        points = []

    if segments is None:
        segments = []

    if holes is None:
        holes = []
        
    if regions is None:
        regions = []

    if dummy_test is None:
        dummy_test  = []
        
    try:
        points =  ensure_numeric(points, num.float)
    except ValueError:
        msg = 'ERROR: Inconsistent points array.'
        raise ANUGAError(msg)
    if points.shape[1] <>2:
        msg = 'ERROR: Bad shape points array.'
        raise ANUGAError(msg)

    # This is after points is numeric
    if pointatts is None or pointatts == []:
        pointatts = [[] for x in range(points.shape[0])]
        
    try:
        # If num.int is used, instead of num.int32, it fails in Linux
        segments = ensure_numeric(segments, num.int32)
        
    except ValueError:
        msg = 'ERROR: Inconsistent segments array.'
        raise ANUGAError(msg)
    
    # This is after segments is numeric
    if segatts is None or segatts == []:
        segatts = [0 for x in range(segments.shape[0])]
        
    try:
        holes = ensure_numeric(holes, num.float)
    except ValueError:
        msg = 'ERROR: Inconsistent holess array.'
        raise ANUGAError(msg)

   
    regions = add_area_tag(regions)
    try:
        regions = ensure_numeric(regions, num.float)
    except  (ValueError, TypeError):
        msg = 'ERROR: Inconsistent regions array.'
        raise ANUGAError(msg)
        
    if not regions.shape[0] == 0 and regions.shape[1] <= 2:
        msg = 'ERROR: Bad shape points array.'
        raise ANUGAError(msg)
    
    try:
        pointatts = ensure_numeric(pointatts, num.float)
    except (ValueError, TypeError):
        msg = 'ERROR: Inconsistent point attributes array.'
        raise ANUGAError(msg)

    if pointatts.shape[0] <> points.shape[0]:
        msg = """ERROR: Point attributes array not the same shape as
        point array."""
        raise ANUGAError(msg)
    if len(pointatts.shape) == 1:
        pointatts = num.reshape(pointatts,(pointatts.shape[0],1))
    
    try:
        segatts = ensure_numeric(segatts, num.int32)
    except ValueError:
        msg = 'ERROR: Inconsistent point attributes array.'
        raise ANUGAError(msg)
    if segatts.shape[0] <> segments.shape[0]:
        msg = """ERROR: Segment attributes array not the same shape as
        segment array."""
        raise ANUGAError(msg)
    
    if mode.find('n'):
        #pass
        mode = 'j' + mode
        # j- Jettisons vertices that are not part of the final
        #    triangulation from the output .node file (including duplicate
        #    input vertices and vertices ``eaten'' by holes).  - output a
        #    list of neighboring triangles
        # EG handles lone verts!
    #
    # GD (June 2014): We get segfaults in some cases with breakLines, unless
    # we remove repeated values in 'points', and adjust segments accordingly
    # 
    pts_complex=points[:,0]+1j*points[:,1] # Use to match points
    i=0 # Use this as a counter, since the length of 'points' changes as we go
    while (i < len(pts_complex)-1):
        i=i+1 # Maximum i = len(pts_complex)-1 = largest index of points
        #
        # Check if points[i,] is the same as a previous point
        if(any(pts_complex[i]==pts_complex[0:i])):
            # Get index of previous occurrence
            earlierInd=(pts_complex[i]==pts_complex[0:i]).nonzero()[0][0]
            # Remove the ith point, and adjust the segments
            for ii in range(len(segments)):
                for j in range(2):
                    if(segments[ii,j] == i):
                        # Segment will use previous occurrence of this point
                        segments[ii,j]=earlierInd
                    if(segments[ii,j]>i):
                        # Decrement the index (since we will remove point[i,])
                        segments[ii,j] = segments[ii,j]-1
            # Remove ith point
            points=num.delete(points, i, 0)
            pointatts=num.delete(pointatts,i,0)
            # Recompute the complex number points for matching
            pts_complex=points[:,0]+1j*points[:,1]
            i=i-1 # Repeat for the last value of i = next point

    trianglelist, pointlist, pointmarkerlist, pointattributelist, triangleattributelist, segmentlist, segmentmarkerlist, neighborlist = triang.genMesh(points,segments,holes,regions,
                          pointatts,segatts, mode)
    mesh_dict = {}
    # the values as arrays
    mesh_dict['generatedtrianglelist'] = trianglelist
    mesh_dict['generatedpointlist'] = pointlist
    # WARNING - generatedpointmarkerlist IS UNTESTED
    mesh_dict['generatedpointmarkerlist'] = pointmarkerlist
    mesh_dict['generatedpointattributelist'] = pointattributelist
    mesh_dict['generatedsegmentlist'] = segmentlist 
    mesh_dict['generatedsegmentmarkerlist'] =  segmentmarkerlist
    mesh_dict['generatedtriangleneighborlist'] = neighborlist
    mesh_dict['qaz'] = 1 #debugging

    #mesh_dict['triangleattributelist'] = triangleattributelist
    if True: 
        mesh_dict['generatedtriangleattributelist'] = triangleattributelist

        if mesh_dict['generatedtriangleattributelist'].shape[1] == 0:
            mesh_dict['generatedtriangleattributelist'] = None
            
        if mesh_dict['generatedpointattributelist'].shape[1] == 0:
            mesh_dict['generatedpointattributelist'] = None
            
        if mesh_dict['generatedtriangleneighborlist'].shape[1] == 0:
            mesh_dict['generatedtriangleneighborlist'] = None
            
        if trianglelist.shape[0] == 0:
            # There are no triangles.
            # this is used by urs_ungridded2sww
            raise NoTrianglesError
    a = mesh_dict['generatedtriangleattributelist']
    # the structure of generatedtriangleattributelist is an list of
    # list of integers.  It is transformed into a list of list of
    # strings later on.  This is then inputted into an triangle
    # object.  The triangle object outputs a list of strings.  Note
    # the subtle change!  How should I handle this?  For my code, when
    # the list of list of integers is transformed, transform it into a
    # list of strings, not a list of list of strings.
    
    return mesh_dict

def add_area_tag(regions):
    """
    So, what is the format?
    A list with
    [x,y,region_tag,area] OR [x,y,region_tag]
    if it's [x,y,region_tag], add a 4th element, value of 0.0.
    """
    if isinstance(regions, ListType):
        for i, region in enumerate(regions):
            if len(region) == 3:
                if isinstance(region, TupleType):
                    #FIXME: How do you convert a tuple to a list?
                    # I can do it a stupid way..
                    tuple = region[:]
                    regions[i] = []
                    for j in tuple:
                        regions[i].append(j)
                    regions[i].append(0.0)
                else:
                    regions[i].append(0.0)
                    
                # let ensure numeric catch this..    
                #len(region) <= 2:
                #msg = 'ERROR: Inconsistent regions array.'
                #raise Exception(msg)
            #elif
    return regions

if __name__ == "__main__":
    pass 
