#!/usr/bin/env python

from builtins import range
import sys

class NoTrianglesError(Exception): pass
#import anuga.mesh_engine.mesh_engine_c_layer as triang
#import anuga.mesh_engine.list_dic as triang

try:
    import meshpy.triangle as triang
    TRILIB = 'meshpy'
except:
    import triangle as triang
    TRILIB = 'triangle'



import numpy as num

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.anuga_exceptions import ANUGAError

def generate_mesh(points=None,
                  segments=None,holes=None,regions=None,
                  pointatts=None,segatts=None,
                  mode=None, dummy_test=None, verbose=False):
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
        points =  ensure_numeric(points, float)
    except ValueError:
        msg = 'ERROR: Inconsistent points array.'
        raise ANUGAError(msg)
    if points.shape[1] !=2:
        msg = 'ERROR: Bad shape points array.'
        raise ANUGAError(msg)

    # This is after points is numeric
    if pointatts is None or pointatts == []:
        pointatts = [[] for x in range(points.shape[0])]

    try:
        # If int is used, instead of num.int32, it fails in Linux
        segments = ensure_numeric(segments, num.int32)

    except ValueError:
        msg = 'ERROR: Inconsistent segments array.'
        raise ANUGAError(msg)

    # This is after segments is numeric
    if segatts is None or segatts == []:
        segatts = [0 for x in range(segments.shape[0])]

    try:
        holes = ensure_numeric(holes, float)
    except ValueError:
        msg = 'ERROR: Inconsistent holes array.'
        raise ANUGAError(msg)


    regions = add_area_tag(regions)
    try:
        regions = ensure_numeric(regions, float)
    except  (ValueError, TypeError):
        msg = 'ERROR: Inconsistent regions array.'
        raise ANUGAError(msg)

    if not regions.shape[0] == 0 and regions.shape[1] <= 2:
        msg = 'ERROR: Bad shape points array.'
        raise ANUGAError(msg)

    try:
        pointatts = ensure_numeric(pointatts, float)
    except (ValueError, TypeError):
        msg = 'ERROR: Inconsistent point attributes array.'
        raise ANUGAError(msg)

    if pointatts.shape[0] != points.shape[0]:
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
    if segatts.shape[0] != segments.shape[0]:
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

    #print(points,segments,holes,regions, pointatts,segatts)

    if TRILIB == 'triangle':
        in_tri = ({'vertices':points})
        if segments.size != 0:
            in_tri['segments'] = segments
        if holes.size != 0:
            in_tri['holes'] = holes
        if regions.size != 0:
            in_tri['regions'] = regions
        if pointatts.size != 0:
            in_tri['vertex_attributes'] = pointatts
        if segatts.size != 0:
            in_tri['segment_markers'] = segatts

        if verbose:
            print('TRIANGLE (triangle):')
            print('  Input sizes')
            print('    vertices: %g'% points.size)
            print('    segments: %g'% segments.size)
            print('    holes:    %g'% holes.size)
            print('    regions:  %g'% regions.size)

        out_tri = triang.triangulate(in_tri,mode)

        if 'vertices' in out_tri:
            pointlist = num.ascontiguousarray(out_tri['vertices'])
        else:
            pointlist = num.empty((0,2),dtype=float)
        if 'vertex_markers' in out_tri:
            pointmarkerlist = num.ascontiguousarray(out_tri['vertex_markers'].reshape(-1))
        else:
            pointmarkerlist = num.empty(pointlist.shape[0],dtype=num.int32)
        if 'triangles' in out_tri:
            trianglelist = num.ascontiguousarray(out_tri['triangles'])
        else:
            trianglelist = num.empty((0,3),dtype=num.int32)
        if 'vertex_attributes' in out_tri:
            pointattributelist = num.ascontiguousarray(out_tri['vertex_attributes'])
        else:
            pointattributelist = num.empty((pointlist.shape[0],0),dtype=float)
        if 'triangle_attributes' in out_tri:
            triangleattributelist = num.ascontiguousarray(out_tri['triangle_attributes'])
        else:
            triangleattributelist = num.empty((trianglelist.shape[0],0),dtype=float)
        if 'segments' in out_tri:
            segmentlist = num.ascontiguousarray(out_tri['segments'])
        else:
            segmentlist = num.empty((0,2),dtype=num.int32)
        if 'segment_markers' in out_tri:
            segmentmarkerlist = num.ascontiguousarray(out_tri['segment_markers'].reshape(-1))
        else:
            segmentmarkerlist = num.empty(segmentlist.shape[0],dtype=num.int32)
        if 'neighbors' in out_tri:
            neighborlist = num.ascontiguousarray(out_tri['neighbors'])
        else:
            neighborlist = num.empty((trianglelist.shape[0],3),dtype=num.int32)


    if TRILIB == 'meshpy':

        in_tri = triang.MeshInfo()

        in_tri.set_points(points)
        in_tri.set_facets(segments, facet_markers=segatts)

        if holes.size != 0:
            in_tri.set_holes(holes)

        if regions.size != 0:
            in_tri.regions.resize(len(regions))
            for i, region in enumerate(regions):
                in_tri.regions[i] = region

        #for i, att in enumerate(in_tri.point_attributes):
        #    print(i,att)
        

        if pointatts.size != 0:
            in_tri.number_of_point_attributes = pointatts.shape[1]
            in_tri.point_attributes.setup()
            for i, pointatt in enumerate(pointatts):
                # meshpy seems to act differently for scalar atts and vector atts
                if pointatts.shape[1] == 1:
                    in_tri.point_attributes[i] = pointatt[0]
                else:
                    in_tri.point_attributes[i] = pointatt

                


        if verbose:
            print('TRIANGLE (meshpy):')
            print('  Input sizes')
            print('    vertices: %g'% points.size)
            print('    segments: %g'% segments.size)
            print('    holes:    %g'% holes.size)
            print('    regions:  %g'% regions.size)


        if verbose:
            print(70*'=')
            print('in_tri')
            print(70*'=')
            in_tri.dump()

        try:
            import locale
        except ImportError:
            have_locale = False
        else:
            have_locale = True
            prev_num_locale = locale.getlocale(locale.LC_NUMERIC)
            locale.setlocale(locale.LC_NUMERIC, "C")

        try:
            out_tri = triang.MeshInfo()
            import meshpy._internals as internals
            internals.triangulate(mode, in_tri, out_tri, triang.MeshInfo(), None)
        finally:
            # restore previous locale if we've changed it
            if have_locale:
                locale.setlocale(locale.LC_NUMERIC, prev_num_locale)

        """
            "points", "point_attributes", "point_markers",
            "elements", "element_attributes", "element_volumes",
            "neighbors",
            "facets", "facet_markers",
            "holes",
            "regions",
            "faces", "face_markers",
            "normals",
        """

        if verbose:
            print(70*'=')
            print('out_tri')
            print(70*'=')
            out_tri.dump()

        if len(out_tri.points) != 0:
            pointlist = num.ascontiguousarray(out_tri.points).reshape(-1,2)
        else:
            pointlist = num.empty((0,2),dtype=float)

        if len(out_tri.point_markers) != 0:
            pointmarkerlist = num.ascontiguousarray(out_tri.point_markers).reshape(-1)
        else:
            pointmarkerlist = num.empty(pointlist.shape[0],dtype=num.int32)

        try:
            unit = out_tri.point_attributes.unit
            pointattributelist = num.ascontiguousarray(out_tri.point_attributes).reshape(-1,unit)
        except:
            pointattributelist = num.empty((pointlist.shape[0],0),dtype=float)

        if len(out_tri.elements) != 0:
            trianglelist = num.ascontiguousarray(out_tri.elements).reshape(-1,3)
        else:
            trianglelist = num.empty((0,3),dtype=num.int32)            

        try:
            unit = out_tri.element_attributes.unit
            triangleattributelist = num.ascontiguousarray(out_tri.element_attributes).reshape(-1,unit)
        except:
            triangleattributelist = num.empty((trianglelist.shape[0],0),dtype=float)

        if len(out_tri.facets) != 0:
            segmentlist = num.ascontiguousarray(out_tri.facets).reshape(-1,2)
        else:
            segmentlist = num.empty((0,2),dtype=num.int32)

        if len(out_tri.facet_markers) != 0:
            segmentmarkerlist = num.ascontiguousarray(out_tri.facet_markers).reshape(-1)
        else:
            segmentmarkerlist = num.empty(segmentlist.shape[0],dtype=num.int32)

        try:
            neighborlist = num.ascontiguousarray(out_tri.neighbors).reshape(-1,3)
        except:
            neighborlist = num.zeros((trianglelist.shape[0],3),dtype=num.int32)        

    if verbose:
        print('  Output sizes')
        print('    vertices:  %g'% pointlist.size)
        print('    triangles: %g' % trianglelist.size)
        print('    segments:  %g'% segmentlist.size)

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

    if verbose:
        from pprint import pprint
        pprint(mesh_dict)

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
    if isinstance(regions, list):
        for i, region in enumerate(regions):
            if len(region) == 3:
                if isinstance(region, tuple):
                    #FIXME: How do you convert a tuple to a list?
                    # I can do it a stupid way..
                    tuple0 = region[:]
                    regions[i] = []
                    for j in tuple0:
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
