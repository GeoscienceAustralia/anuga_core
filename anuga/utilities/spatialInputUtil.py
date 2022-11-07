"""

Routines to ease the import of spatial data to ANUGA

Key routines:
    readShp_1PolyGeo, readShp_1LineGeo -- read SIMPLE shapefile geometries into ANUGA as a list of points.
                                          The supported geometries are quite restricted, see help
    read_polygon -- Reads either a polygon or line shapefile or a csv as a polygon
    readShpPointsAndAttributes -- read a multi-point shapefile with its attributes into ANUGA

    ListPts2Wkb -- (Probably for internal use) Convert a list of points to a
                    Wkb geometry, allowing us to use GDALs geometry tools
    Wbk2ListPts -- reverse of ListPts2Wkb

    addIntersectionPtsToLines -- (Probably for internal use, see add_intersections_to_domain_features) 
                                 Given 2 intersecting lines, add their intersection points to them, or
                                 move existing points if they are within a given distance of the intersection.

    add_intersections_to_domain_features -- Add intersections to bounding_polygon, breaklines, riverwalls
                                            This sets them up for being passed to the mesh generator.
                                            Points within a given distance of an intersection can be replaced
                                            with the intersection point if desired

    rasterValuesAtPoints -- (Probably for internal use, see quantityRasterFun in quantity_setting_functions)
                            Quite efficiently get raster cell values at points in any gdal-compatible raster
                            [gridPointsInPolygon could in future be linked with this to get raster values in a region,
                             if we develop a version of ANUGA with sub-grid topography]

    readRegionPtAreas -- read a shapefile containing regionPtAreas -- xy coordinates + 1 attribute, which is
                         the mesh triangle side length (or area) limit. Can be passed as regionPtAreas in
                         the mesh generation stage.

    readListOfBreakLines -- takes a list of shapefile names, each containing a single line geometry, and reads
                            it in as a dictionary of breaklines. 
    
    readListOfRiverWalls -- takes a list of csv names, each containing a xyz polylines defining riverwalls, with
                            an optional first line defining the non-default riverwall par, and returns a dictionary
                            of the riverwalls and the riverwall_Par

    polygon_from_matching_breaklines -- given a pattern (string) matching exactly 2 breaklines in a directory,
                                        convert them to a single polygon. This is useful to e.g. make
                                        a polygon defining the extent of a channel that is defined from 2 breaklines

    
"""

import sys
import os
import os.path
import copy
import struct
import numpy
#from matplotlib import path
import anuga
from anuga.geometry.polygon import inside_polygon

try:
    import osgeo.gdal as gdal
    import osgeo.ogr as ogr
    gdal_available = True
    #import osr # Not needed here but important in general
except ImportError as err:
    gdal_available = False


#####################################  
if gdal_available:


    def readShp_1PolyGeo(shapefile, dropLast=True):
        """
            Read a "single polygon" shapefile into an ANUGA_polygon object
            (a list of lists, each containing a polygon coordinate)
    
            The attribute table is ignored, and there can only be a single geometry in the shapefile
            
            INPUTS: shapefile -- full path to shapefile name
                    dropLast -- Logical. If so, cut the last point (which closes
                                the polygon). ANUGA uses this by default for e.g. 
                                bounding polygons
        """
    
        # Get the data
        driver=ogr.GetDriverByName("ESRI Shapefile")
        dataSrc=driver.Open(shapefile, 0)
        #dataSrc=ogr.Open(shapefile)
        layer=dataSrc.GetLayer()

        # Check it is a polygon
        layerType=ogr.GeometryTypeToName(layer.GetGeomType())
        if not layerType=='Polygon':
            msg= shapefile +' is not a polygon shapefile'
            raise Exception(msg)

        # Need a single polygon 
        try:
            assert(len(layer)==1)
        except:
            print(shapefile)
        
        boundary_poly=[]
        for feature in layer:
            #geom=feature.GetGeometryReg()
            boundary=feature.GetGeometryRef().Boundary().GetPoints()
            boundary=[list(pts) for pts in boundary]
            boundary_poly.extend(boundary)
    
        if(dropLast):
            # Return a list of points, without last point [= first point]
            return boundary_poly[:-1]
        else:
            return boundary_poly
    
    ####################
    
    def readShp_1LineGeo(shapefile):
        """
            Read a "single-line" shapefile into a list of lists (each containing a point),
                resembling an ANUGA polygon object
    
            The attribute table is ignored, and there can only be a single geometry in the shapefile
            
            INPUTS: shapefile -- full path to shapefile name
        """
    
        driver=ogr.GetDriverByName("ESRI Shapefile")
        dataSrc=driver.Open(shapefile, 0)
        #dataSrc=ogr.Open(shapefile)
        layer=dataSrc.GetLayer()
     
        # Check it is a line
        layerType=ogr.GeometryTypeToName(layer.GetGeomType())
        if not layerType=='Line String':
            msg= shapefile +' is not a line shapefile'
            raise Exception(msg) 

        # Need a single line 
        try:
            assert len(layer)==1
        except:
            print(shapefile)
        
        line_all=[]
        for feature in layer:
            line=feature.GetGeometryRef().GetPoints()
            line=[list(pts) for pts in line]
            line_all.extend(line)
    
        return line_all

    ###########################################################################
    def read_csv_optional_header(filename):
        """Read a csv file of numbers, which optionally has a single header 
            row containing alphabetical characters (which is ignored if it 
            exists)
           
            INPUT:
            @param filename -- name of appropriate file with ',' delimiter
            
            OUTPUT:
            A numpy array with the numeric data
        """
           
        f=open(filename)
        firstLine=f.readline()
        f.close()
        hasLetters=any(c.isalpha() for c in firstLine)
        outPol = numpy.genfromtxt(
            filename, delimiter=',',skip_header=int(hasLetters))

        return outPol
           
    ###########################################################################
    def read_polygon(filename, close_polygon_shapefiles=False):
        """

        Read a shapefile (polygon or line) or a csv file as an anuga polygon

        In the case of a csv file, we permit either 1 or 0 header rows, and
        the file can have > 2 columns [but only the first 2 are used]

        Try to automatically do the correct thing based on the filename

        """ 
        # Check the file exists
        msg= 'Could not read '+ filename
        assert os.path.isfile(filename), msg 

        # Get the file extension type
        fileNameNoExtension , fileExtension = os.path.splitext(filename)

        if fileExtension == '.shp':
            # Read as either a polygon or line shapefile
            try:
                outPol=readShp_1PolyGeo(filename, dropLast = (not close_polygon_shapefiles))
                assert len(outPol)>1
            except:
                try:
                    outPol= readShp_1LineGeo(filename)
                    assert len(outPol)>1
                except:
                    msg= 'Could not read '+ filename +\
                         ' as either polygon or line shapefile'
                    raise Exception(msg)
        else:
            try:
                outPol = read_csv_optional_header(filename)
                # Only take the first 2 columns
                outPol = outPol[:,0:2].tolist()
            except:
                msg = 'Failed reading polygon '+ filename +\
                      ' with anuga.utilities.spatialInputUtils.read_polygon'
                raise Exception(msg)

        return outPol 

    #################### 

    def readShpPtsAndAttributes(shapefile):
        """
            Read a point shapefile with an attribute table into a list
    
            INPUT: shapefile -- name of shapefile to read
            
            OUTPUT: List with 
            [ list_of_points, list_of_attributes, names_of_attributes]
        """
        
        driver=ogr.GetDriverByName("ESRI Shapefile")
        dataSrc=driver.Open(shapefile, 0)
        #dataSrc=ogr.Open(shapefile)
        layer=dataSrc.GetLayer()
    
        pts=[]
        attribute=[]
        for i, feature in enumerate(layer):
            if(i==0):
                attributeNames=list(feature.keys())
            pt=feature.GetGeometryRef().GetPoints()
            pt=[list(p) for p in pt]
            pts.extend(pt)
            att=[feature[i] for i in attributeNames]
            attribute.extend(att)
    
        return [pts, attribute, attributeNames]

    ########################################
    def readShpPts(shapefile):
        """
            Wrapper around readShpPtsAndAttributes
            Only get the points
        """

        out=readShpPtsAndAttributes(shapefile)[0]
        return out

    ########################################
    def read_points(filename):
        """
            Reads x,y geometries from a point shapefile or csv file (anuga polygon type format),
            and returns as a list of lists
        """
        # Check the file exists
        msg= 'Could not read '+ filename
        assert os.path.isfile(filename), msg 

        # Get the file extension type
        fileNameNoExtension , fileExtension = os.path.splitext(filename)

        if fileExtension=='.shp':
            try:
                points = readShpPts(filename)
            except:
                msg = 'Could not read points from ' + filename 
                raise Exception(msg)
        else:
            # Assume txt format
            try:
                #points=numpy.genfromtxt(filename,delimiter=',',skip_header=1)
                #points=points[:,0:2].tolist()
                points=anuga.read_polygon(filename)
            except:
                msg = 'Could not read points from ' + filename +\
                      '. Make sure it has a single header row, ' +\
                      'with comma separator, and the first 2 columns are x,y'
                raise Exception(msg)
        return points
    
    ########################################
    def ListPts2Wkb( ptsIn, geometry_type='line', appendFirstOnEnd=None):
        """
            Convert list of points to a GDAl Wkb format
            Can be either points, lines, or polygon
    
            Purpose is that once data in in Wkb format, we can use GDAL's geometric operations
                (e.g. to test for intersections, etc)
            
            INPUT: ptsIn -- list of points in the format [[x0,y0], [x1, y1], ..., [xn, yn]]
                          Actually it is also ok if [x0,y0], .. is a tuple instead
                   geometry_type -- 'point' or 'line' or 'polygon'
    
                   appendFirstOnEnd -- logical. If true, add the first point to the
                        end. Probably wanted for polygons when they are unclosed in ANUGA
    
            OUTPUT:
                    The points as a Wkb Geometry. 
                    geometry_type='point' produces a MULTIPOINT Geometry
                                 ='line' produces a LINESTRING Geometry
                                 ='polygon' produces a POLYGON Geometry
    
            FIXME: This might not sensibly-use the gdal geometry types (although it
                    passes our tests) -- consider revising
        """
        # Avoid modifying ptsIn
        pts=copy.copy(ptsIn)
    
        if appendFirstOnEnd is None:
            if(geometry_type=='polygon'):
                appendFirstOnEnd=True
            else:
                appendFirstOnEnd=False
    
        if appendFirstOnEnd:
            pts.append(pts[0])    
        
        if(geometry_type=='point'):
            data=ogr.Geometry(ogr.wkbMultiPoint)
        elif(geometry_type=='line'):
            data=ogr.Geometry(ogr.wkbLineString)
        elif(geometry_type=='polygon'):
            data=ogr.Geometry(ogr.wkbLinearRing)
        else:
            msg = "Type must be either 'point' or 'line' or 'polygon'"
            raise Exception(msg)
         
        for i in range(len(pts)):
            if(len(pts[i])==2):
                if(geometry_type=='point'):
                    newPt=ogr.Geometry(ogr.wkbPoint)
                    newPt.AddPoint(pts[i][0], pts[i][1]) 
                    data.AddGeometryDirectly(newPt)
                else:
                    data.AddPoint(pts[i][0], pts[i][1])
            elif(len(pts[i])==3):
                if(geometry_type=='point'):
                    newPt=ogr.Geometry(ogr.wkbPoint)
                    newPt.AddPoint(pts[i][0], pts[i][1], pts[i][2]) 
                    data.AddGeometryDirectly(newPt)
                else:
                    data.AddPoint(pts[i][0], pts[i][1], pts[i][2])
            else:
                raise Exception('Points must be either 2 or 3 dimensional')
    
        if(geometry_type=='polygon'):   
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(data)
            data=poly
        
        return(data)
    
    ############################################################################
    def Wkb2ListPts(wkb_geo, removeLast=False, drop_third_dimension=False):
        """
            Reverse of ListPts2Wkb
        """
        if(wkb_geo.GetGeometryName()=='POLYGON'):
            X=wkb_geo.GetBoundary()
            new=[ list(X.GetPoints()[i]) for i in range(len(X.GetPoints())) ]
        elif(wkb_geo.GetGeometryName()=='MULTIPOINT'):
            new=[ [feature.GetX(), feature.GetY(),feature.GetZ()] for feature in wkb_geo]
        elif(wkb_geo.GetGeometryName()=='LINESTRING'):
            new=[ list(wkb_geo.GetPoints()[i]) for i in range(len(wkb_geo.GetPoints())) ]
        else:
            raise Exception('Geometry type not supported')
        
        if(removeLast):
            new=new[:-1]
        if(drop_third_dimension):
            new = [new[i][0:2] for i in range(len(new))]
        return new
    
    ############################################################################
    def compute_squared_distance_to_segment(pt, line):
        """
            Compute the squared distance between a point and a [finite] line segment
    
            INPUT: pt -- [x,y]
                   line -- [[x0,y0],[x1,y1]] -- 2 points defining a line segment
    
            OUTPUT: The distance^2 of pt to the line segment
    
        """
        p0 = line[0]
        p1 = line[1]
        #
        # Get unit vector along segment
        seg_unitVec_x = float(p1[0]-p0[0])
        seg_unitVec_y = float(p1[1]-p0[1])
        segLen = (seg_unitVec_x**2+seg_unitVec_y**2)**0.5
        if(segLen == 0.):
            raise Exception('Line has repeated points: Line %s Pt %s' % (str(line),str(pt)))

        seg_unitVec_x = seg_unitVec_x/segLen
        seg_unitVec_y = seg_unitVec_y/segLen

        # Get vector from pt to p0 
        pt_p0_vec_x = float(pt[0]-p0[0])
        pt_p0_vec_y = float(pt[1]-p0[1])
        pt_p0_vec_len_squared = (pt_p0_vec_x**2 + pt_p0_vec_y**2)

        # Get dot product of above vector with unit vector
        pt_dot_segUnitVec = (pt_p0_vec_x)*seg_unitVec_x + (pt_p0_vec_y)*seg_unitVec_y

        if( (pt_dot_segUnitVec < segLen) and (pt_dot_segUnitVec > 0.)):
            # The nearest point on the line is actually between p0 and p1, so we have a 'real' candidate
            # Get distance^2
            output = pt_p0_vec_len_squared - pt_dot_segUnitVec**2.
        else:
            # Distance is the min distance from p0 and p1. 
            output = min( pt_p0_vec_len_squared,  (float(pt[0]-p1[0])**2+float(pt[1]-p1[1])**2))

        if(output < -1.0e-06):
            print('Diagnostic numbers follow: ')
            print(output)
            print(pt_p0_vec_len_squared)
            print(pt_dot_segUnitVec)
            print(pt)
            print(p1)
            print(p0)
            raise Exception('round-off in compute_squared_distance_to_segment')
        if(output < 0.):
            output=0.
        return output
    
    ############################################################################
    
    def find_nearest_segment(pt, segments):
        """
            Given a point and a line, find the line segment nearest to the line
    
            NOTE: The answer can be ambiguous if one of the segment endpoints is
                the nearest point. In that case, the behaviour is determined by the behaviour
                of numpy.argmin. Won't be a problem for this application
    
            INPUT: pt -- [x,y]
                   segments -- [[x0,y0], [x1,y1], ...]
                             A list of points, consecutive points are interpreted
                             as joined and so defining line segments
            
            OUTPUT: The squared distance, and the index i of the segment
                    [x_i,y_i],[x_i+1,y_i+1] in segments which is closest to pt
        """
        ll=len(segments)
        if(ll<=1):
            raise Exception('Segments must have length > 1 in find_nearest_segment')
       
        ptDist_sq=numpy.zeros(ll-1) # Hold the squared distance from the point to the line segment
        for i in range(len(segments)-1):
            # Compute distance from segment
            ptDist_sq[i]=compute_squared_distance_to_segment(pt, [segments[i],segments[i+1]])
        
        return [ptDist_sq.min(), ptDist_sq.argmin()]
    
    ######################################################
    
    def shift_point_on_line(pt, lineIn, nearest_segment_index):
        """
            Support pt is a point, which is near to the 'nearest_segment_index'
                segment of the line 'lineIn'
    
            This routine moves the nearest end point of that segment on line to pt.
    
            INPUTS: pt -- [x,y] point
                    lineIn -- [ [x0, y0], [x1, y1], ..., [xN,yN]] 
                    nearest_segment_index = index where the distance of pt to
                        the line from [x_i,y_i] to [x_i+1,y_i+1] is minimum
    
            OUTPUT: The new line
        """
        # Avoid Changing line
        line=copy.copy(lineIn)
    
        # Replace the nearest point on L1 with the intersection point
        p0 = line[nearest_segment_index]
        p1 = line[nearest_segment_index+1]
        d_p0 = ( (pt[0]-p0[0])**2 + (pt[1]-p0[1])**2)
        d_p1 = ( (pt[0]-p1[0])**2 + (pt[1]-p1[1])**2)
        changeP1=(d_p1<d_p0)
        line[nearest_segment_index+changeP1][0] = pt[0]
        line[nearest_segment_index+changeP1][1] = pt[1]
    
        return line
    
    #################################################################################
    def insert_intersection_point(intersectionPt, line_pts, point_movement_threshold, verbose=False):
        """
            Add intersectionPt to line_pts, either by inserting it, or if a point on line_pts is
            closer than point_movement_threshold, then by moving that point to the intersection point
    
            INPUTS:
                    intersectionPt -- the intersection point [known to lie on line_pts]
                    line_pts -- ordered list of [x,y] points making a line. 
                    point_movement_threshold -- used to decide to move or add intersectionPt
    
            OUTPUT:
                new version of lint_pts with the point added
        """
        # Avoid pointer/copy issues
        L1_pts = copy.copy(line_pts)
        iP=copy.copy(intersectionPt)
    
        # Adjust L1
        tmp = find_nearest_segment(intersectionPt, L1_pts)
        # Compute the distance from the end points of the segment to the
        # intersection point. Based on this we decide to add or move the point
        p0 = L1_pts[tmp[1]]
        p1 = L1_pts[tmp[1]+1]
        endPt_Dist_Sq=min( ( (p0[0]-iP[0])**2 + (p0[1]-iP[1])**2), 
                           ( (p1[0]-iP[0])**2 + (p1[1]-iP[1])**2))
        #
        if(endPt_Dist_Sq>point_movement_threshold**2):
            # Insert the intersection point. We do this in a tricky way to
            # account for the possibility of L1_pts having > 2 coordinates
            # (riverWalls)
            if verbose:
                print('      Inserting new point')
            dummyPt=copy.copy(L1_pts[tmp[1]])
            L1_pts.insert(tmp[1]+1,dummyPt) 
            L1_pts[tmp[1]+1][0]=iP[0]
            L1_pts[tmp[1]+1][1]=iP[1]
            if(len(L1_pts[tmp[1]+1])==3):
                # Treat 3rd coordinate
                # Find distance of inserted point from neighbours, and 
                # Set 3rd coordinate as distance-weighted average of the others
                d0=((L1_pts[tmp[1]][0]-L1_pts[tmp[1]+1][0])**2.+\
                   (L1_pts[tmp[1]][1]-L1_pts[tmp[1]+1][1])**2.)**0.5
                d1=((L1_pts[tmp[1]+2][0]-L1_pts[tmp[1]+1][0])**2.+\
                   (L1_pts[tmp[1]+2][1]-L1_pts[tmp[1]+1][1])**2.)**0.5
                L1_pts[tmp[1]+1][2] = (d0*L1_pts[tmp[1]+2][2] + d1*L1_pts[tmp[1]][2])/(d0+d1)
    
        else:
            if verbose:
                print('      Shifting existing point')
            # Move a point already on L1
            L1_pts=shift_point_on_line(iP, L1_pts, tmp[1])
    
        return L1_pts
    
    #######################################################################################################
    def check_polygon_is_small(intersection, buf, tol2=100.):
        """
            
            Elsewhere in the code, we check whether lines intersect by buffering them
            to polygons with a small buffer = buf, then getting the intersection.
             [since intersection with polygons is supported by gdal, but apparently
            not directly with lines].  
           
            The logic of our code only works with point intersections,  
            and it will fails if 2 lines overlap in a line.
    
            We crudely check for this situation here, by ensuring that the intersection polygon is 'small'
    
            WARNING: The gdal geometry routines may be a bit rough (?)
                     Intersections not very precise, etc (?)
    
            INPUT: intersection -- intersection of the 2 lines [gdal geometry]
                   buf -- a length scale giving the size of the intersection extent that we expect for a point
                   tol2 -- Throw an error if the x or y extent is greater than buf*tol2. Seems this needs to be
                           large sometimes -- this might reflect the stated weaknesses of GDALs geometry routines?
            OUTPUT: True/False
                    False should suggest that the intersection is not describing a point
        """
    
        extent=intersection.GetEnvelope()
        assert(len(extent)==4) # Make sure this assumption is valid
        if( (abs(extent[0]-extent[1])>buf*tol2) or (abs(extent[2]-extent[3]) > buf*tol2)):
            return False
        else:
            return True
    
    #######################################################################################################
    
    def addIntersectionPtsToLines(L1,L2, point_movement_threshold=0.0, buf=1.0e-05, tol2 = 100,
                                  verbose=True, nameFlag=''):
        """
            Add intersection points to lines L1 and L2 if they intersect each other
    
            This is useful e.g. so that intersections can be exact (important for
                 mesh generation in ANUGA)
    
            It currently only supports point intersection of 2 lines.
                Line intersections should fail gracefully
            
            INPUTS:  L1, L2 = Wkb LineString geometries
    
                     point_movement_threshold -- if the distance from the nearest
                        point on L1 or L2 to the intersection is < point_movement_threshold, then the
                        nearest point has its coordinates replaced with the intersection point. This is
                        to prevent points being too close to each other
        
                     buf = tolerence that is used to buffer lines to find
                            intersections. Probably doesn't need modification
    
                     tol2 = [see check_polygon_is_small] Probably doesn't need to change
    
                     nameFlag = This will be printed if intersection occurs.
                            Convenient way to display intersecting filenames
            
            OUTPUTS: L1,L2 with intersection points added in the right places
        """
    
        if(L1.Intersects(L2)):
            # Get points on the lines 
            L1_pts=Wkb2ListPts(L1) 
            L2_pts=Wkb2ListPts(L2) 
            
            # Buffer lines by a small amount
            L1_buf=L1.Buffer(buf)
            L2_buf=L2.Buffer(buf)
                
            # Get intersection point[s]    
            L1_L2_intersect=L1_buf.Intersection(L2_buf)
            if(L1_L2_intersect.GetGeometryCount()==1):
                if(not check_polygon_is_small(L1_L2_intersect, buf, tol2)):
                    msg = 'line intersection is not allowed. ' + \
                          'Envelope %s '% str(L1_L2_intersect.GetEnvelope())
                    raise Exception(msg)
                # Seems to need special treatment with only 1 intersection point
                intersectionPts=[L1_L2_intersect.Centroid().GetPoint()]
            else:
                intersectionPts=[]
                for feature in L1_L2_intersect:
                    if(not check_polygon_is_small(feature, buf, tol2)):
                        print(feature.GetEnvelope())
                        raise Exception('line intersection is not allowed')
                    intersectionPts.append(feature.Centroid().GetPoint())
    
            if(verbose):
                print(nameFlag)
                print('    Treating intersections in ', len(intersectionPts) , ' locations')
                print(intersectionPts)
    
            # Insert the points into the line segments
            for i in range(len(intersectionPts)):
                L1_pts = insert_intersection_point(intersectionPts[i], L1_pts, 
                            point_movement_threshold, verbose=verbose)
                L2_pts = insert_intersection_point(intersectionPts[i], L2_pts, 
                            point_movement_threshold, verbose=verbose)
                
            # Convert to the input format
            L1_pts=ListPts2Wkb(L1_pts,geometry_type='line')
            L2_pts=ListPts2Wkb(L2_pts,geometry_type='line')
    
            return [L1_pts, L2_pts]
        else:
            return [L1, L2]
   
    ###########################################################
    def getRasterExtent(rasterFile, asPolygon=False):
        """
            Sometimes we need to know the extent of a raster
            i.e. the minimum x, maximum x, minimum y, and maximum y values
            
            INPUT:
                rasterFile -- a gdal compatible rasterfile
                asPolygon -- if False, return [xmin,xmax,ymin,ymax].
                             If True, return [ [xmin,ymin],[xmax,ymin],[xmax,ymax],[xmin,ymax]]
            OUTPUT
                The extent as defined above

        """
        raster = gdal.Open(rasterFile)
        transform=raster.GetGeoTransform()
        xOrigin = transform[0]
        yOrigin = transform[3]
        xPixels = raster.RasterXSize
        yPixels = raster.RasterYSize

        # Compute the other extreme corner
        x2 = xOrigin + xPixels * transform[1] + yPixels * transform[2]
        y2 = yOrigin + xPixels * transform[4] + yPixels * transform[5]
        
        xmin=min(xOrigin,x2) 
        xmax=max(xOrigin,x2)

        ymin=min(yOrigin,y2)
        ymax=max(yOrigin,y2)

        if(asPolygon):
            return [ [xmin,ymin], [xmax,ymin], [xmax,ymax], [xmin,ymax]]
        else:
            return [xmin,xmax,ymin,ymax]

 
    ###########################################################
    def rasterValuesAtPoints(
        xy, 
        rasterFile, 
        band=1, 
        nodata_rel_tol = 1.0e-08,
        interpolation = 'pixel'):
        """
            Get raster values at point locations.
            Can be used to e.g. set quantity values
       
            INPUT: 
            @param xy = numpy array with point locations

            @param rasterFile = Filename of the gdal-compatible raster

            @param band = band of the raster to get

            @param nodata_rel_tol = Values are treated as nodata if
                ( abs(elev - nodataval) < nodata_rel_tol*abs(nodataval) )
                This allows for truncation errors in nodata values which seem
                to be introduced by some file-type conversions

            @param interpolation 'pixel' or 'bilinear' determines how the
                    raster cell values are used to set the point value
    
            OUTPUT:
            1d numpy array with raster values at xy
    
        """
        # Raster info
        raster = gdal.Open(rasterFile)
        rasterBand = raster.GetRasterBand(band)
        rasterBandType = gdal.GetDataTypeName(rasterBand.DataType)
        nodataval = rasterBand.GetNoDataValue()
    
        # Projection info
        transform = raster.GetGeoTransform()
        xOrigin = transform[0]
        yOrigin = transform[3]
        pixelWidth = transform[1]
        pixelHeight = transform[5] # Negative
        
        # Get coordinates in pixel values
        px = (xy[:,0] - xOrigin)/pixelWidth
        py = (xy[:,1] - yOrigin)/pixelHeight
      
        # Hold elevation 
        elev = px*0. 
    
        # Get the right character for struct.unpack
        if (rasterBandType == 'Int16'):
            CtypeName='h'
        elif (rasterBandType == 'Float32'):
            CtypeName='f'
        elif (rasterBandType == 'Float64'):
            CtypeName='d'
        elif (rasterBandType == 'Byte'):
            CtypeName='B'
        elif (rasterBandType == 'Int32'):
            CtypeName='i'
        else:
            print('unrecognized DataType:', rasterBandType)
            print('You might need to edit this code to read the data type')
            raise Exception('Stopping')
  
        # Upper bounds for pixel values, so we can fail gracefully
        xMax = raster.RasterXSize
        yMax = raster.RasterYSize
        if(px.max() < xMax and px.min() >= 0 and py.max() < yMax and py.min() >= 0):
            pass
        else:
            msg = 'Trying to extract point values that exceed the raster extent'
            raise Exception(msg)

        # Get values -- seems we have to loop, but it is efficient enough
        for i in range(len(px)):

            if(interpolation == 'pixel'):
                # Pixel coordinates
                xC = int(numpy.floor(px[i]))
                yC = int(numpy.floor(py[i]))

                structval = rasterBand.ReadRaster(xC,yC,1,1,
                    buf_type=rasterBand.DataType)
                elev[i] = struct.unpack(CtypeName, structval)[0]

            elif(interpolation=='bilinear'):
                # Pixel coordinates
                xl = int(numpy.floor(px[i]))
                yl = int(numpy.floor(py[i]))

                # Find neighbours required for bilinear interpolation
                # l = lower, u = upper
                if(px[i] - xl > 0.5):
                    xu = min(xl + 1, xMax - 1)
                else:
                    # Swap xl for xu
                    xu = xl + 0
                    xl = max(xu - 1, 0)

                if(py[i] - yl > 0.5):
                    yu = min(yl + 1, yMax - 1)
                else:
                    yu = yl + 0
                    yl = max(yu - 1, 0)

                # Map x,y to unit square
                if(xu > xl):
                    x = px[i] - (xl + 0.5)
                else:
                    x = 0.

                if(yu > yl):
                    y = py[i] - (yl + 0.5)
                else:
                    y = 0.

                if not ( (x>=0.) & (x<=1.)):
                    print('x-values error: ', x, xl, xu, px[i], xMax)
                    raise Exception('x out of bounds')

                if not ( (y>=0.) & (y<=1.)):
                    print('y-values error: ', y, yl, yu, py[i])
                    raise Exception('y out of bounds')

                # Lower-left
                structval = rasterBand.ReadRaster(xl,yl,1,1,
                    buf_type=rasterBand.DataType)
                r00 = struct.unpack(CtypeName, structval)[0]
                # Upper left
                structval = rasterBand.ReadRaster(xl,yu,1,1,
                    buf_type=rasterBand.DataType)
                r01 = struct.unpack(CtypeName, structval)[0]
                # Lower-right
                structval = rasterBand.ReadRaster(xu,yl,1,1,
                    buf_type=rasterBand.DataType)
                r10 = struct.unpack(CtypeName, structval)[0]
                # Upper right
                structval = rasterBand.ReadRaster(xu,yu,1,1,
                    buf_type=rasterBand.DataType)
                r11 = struct.unpack(CtypeName, structval)[0]

                # Bilinear interpolation
                elev[i] = r00*(1.-x)*(1.-y) + r01*(1.-x)*y +\
                          r10*x*(1.-y) + r11*x*y

                # Deal with nodata. This needs to be in the loop
                # Just check if any of the pixels are nodata
                if nodataval is not None:
                    if numpy.isfinite(nodataval):
                        rij = numpy.array([r00, r01, r10, r11])
                        rel_tol = ( abs(rij - nodataval) < nodata_rel_tol*abs(nodataval) )
                        missing = (rel_tol).nonzero()[0]
                        if len(missing) > 0:
                            elev[i] = numpy.nan
            else:
                raise Exception('Unknown value of "interpolation"')            

        # Deal with nodata for pixel based interpolation [efficient treatment
        # outside of loop]
        if (interpolation == 'pixel'):
            if nodataval is not None:
                if numpy.isfinite(nodataval):
                    rel_tol = ( abs(elev - nodataval) < nodata_rel_tol*abs(nodataval) )
                    missing = (rel_tol).nonzero()[0]
                    if len(missing) > 0:
                        elev[missing] = numpy.nan

        return elev


    def gridPointsInPolygon(polygon, approx_grid_spacing=[1.,1.], eps=1.0e-06):
        """
            Get a 'grid' of points inside a polygon. Could be used with rasterValuesAtPoints to 
               get a range of raster values inside a polygon
    
            Approach: A 'trial-grid' of points is created which is 'almost'
                      covering the range of the polygon (xmin-xmax,ymin-ymax). 
    
                      (Actually it is just inside this region, to avoid polygon-boundary issues, see below)
    
                      Then we find those points which are actually inside the polygon.
    
                      The x/y point spacing on the trial-grid will be close to
                      approx_grid_spacing, but we ensure there are at least 4x4 points on the trial grid.
                      Also, we reduce the spacing so that the min_x+R and max_x-R
                        are both similarly close to the polygon extents [see definition of R below]
    
            INPUTS:
                polygon -- the polygon in ANUGA format (list of lists of ordered xy points)
    
                approx_grid_spacing -- the approximate x,y grid spacing
    
                eps -- 'trial-grid' of points has x range from min_polygon_x+R to
                        max_polygon_x - R, where R = (max_polygon_x-min_polygon_x)*eps
                        ( and similarly for y).
    
                       This makes it more likely that points are inside the
                        polygon, not on the boundaries. Points on the boundaries can confuse the
                        point-in-polygon routine
    
            OUTPUTS: A n x 2 numpy array of points in the polygon
        """
    
        # Get polygon extent
        polygonArr = numpy.array(polygon)
        poly_xmin = polygonArr[:,0].min()
        poly_xmax = polygonArr[:,0].max()
        poly_ymin = polygonArr[:,1].min()
        poly_ymax = polygonArr[:,1].max()
    
        # Make a 'grid' of points which covers the polygon
        xGridCount = max( numpy.ceil( (poly_xmax-poly_xmin)/approx_grid_spacing[0]+1. ).astype(int), 4)
        R = (poly_xmax-poly_xmin)*eps
        Xvals = numpy.linspace(poly_xmin+R,poly_xmax-R, xGridCount)
        yGridCount = max( numpy.ceil( (poly_ymax-poly_ymin)/approx_grid_spacing[1]+1. ).astype(int), 4)
        R = (poly_ymax-poly_ymin)*eps
        Yvals = numpy.linspace(poly_ymin+R,poly_ymax-R, yGridCount)
    
        xGrid, yGrid = numpy.meshgrid(Xvals,Yvals)
        Grid = numpy.vstack([xGrid.flatten(),yGrid.flatten()]).transpose()
    
        keepers = inside_polygon(Grid, polygon)
        if(len(keepers) == 0):
            raise Exception('No points extracted from polygon')
        xyInside = Grid[keepers,:]
    
        return(xyInside)
    
    #########################################################################
    # Function to search for pattern matches in a string (turns out to be useful)
    def matchInds(pattern, stringList):
        """
            Find indices in stringList which match pattern
        """
        #matches=[ (pattern in stringList[i]) for i in range(len(stringList))]
        matches = []
        for i in range(len(stringList)):
            if pattern in stringList[i]:
                matches.append(i)
        return matches
    
    
    ###########################################################################
    #
    # Less generic utilities below
    # 
    # These are more 'anuga-specific' than above, aiming to make nice interfaces
    # in ANUGA scripts
    #
    ############################################################################
    
    
    def add_intersections_to_domain_features(
        bounding_polygonIn,
        breakLinesIn={ }, 
        riverWallsIn={ }, 
        point_movement_threshold=0.,
        verbose=True):
        """
            If bounding polygon / breaklines /riverwalls intersect with each 
            other, then add intersection points.

            INPUTS:
                bounding_polygonIn -- the bounding polygon in ANUGA format
                breakLinesIn -- the breaklines dictionary
                riverWallsIn -- the riverWalls dictionary
                point_movement_threshold -- if points on lines 
                    are < this distance from intersection points, then they are
                    replaced with the intersection point. This can prevent
                    overly close points from breaking the mesh generation

            OUTPUT:
                List with bounding_polygon,breakLines,riverwalls
        """

        bounding_polygon = copy.copy(bounding_polygonIn)
        breakLines = copy.copy(breakLinesIn)
        riverWalls = copy.copy(riverWallsIn)

        # Quick exit 
        if (breakLines == {}) and (riverWalls == {}): 
            return [bounding_polygon, breakLines, riverWalls]

        # Clean intersections of breakLines with itself
        if(verbose): 
            print('Cleaning breakline intersections')
        if(len(breakLines)>0):
            kbl = list(breakLines.keys())
            for i in range(len(kbl)):
                n1 = kbl[i]
                for j in range(len(kbl)):
                    if(i >= j):
                        continue
                    n2 = kbl[j]
                    # Convert breaklines to wkb format
                    bl1 = ListPts2Wkb(breakLines[n1],geometry_type='line')
                    bl2 = ListPts2Wkb(breakLines[n2],geometry_type='line')
                    # Add intersection points
                    bl1, bl2 = addIntersectionPtsToLines(
                        bl1, bl2,
                        point_movement_threshold=point_movement_threshold,
                        verbose=verbose, nameFlag=n1+' intersects '+ n2)
                    breakLines[n1] = Wkb2ListPts(bl1)
                    breakLines[n2] = Wkb2ListPts(bl2)


        # Clean intersections of riverWalls with itself
        if(verbose): 
            print('Cleaning riverWall intersections')
        if(len(riverWalls)>0):
            krw=list(riverWalls.keys())
            for i in range(len(krw)):
                n1=krw[i]
                for j in range(len(krw)):
                    if(i>=j):
                        continue 
                    n2 = krw[j]
                    # Convert breaklines to wkb format
                    rw1 = ListPts2Wkb(riverWalls[n1],geometry_type='line')
                    rw2 = ListPts2Wkb(riverWalls[n2],geometry_type='line')
                    # Add intersection points
                    rw1, rw2 = addIntersectionPtsToLines(rw1, rw2,\
                                    point_movement_threshold=point_movement_threshold,\
                                    verbose=verbose, nameFlag=n1+' intersects '+ n2)
                    riverWalls[n1] = Wkb2ListPts(rw1)
                    riverWalls[n2] = Wkb2ListPts(rw2)
    
        # Clean intersections of breaklines with riverwalls
        if(verbose): 
            print('Cleaning breakLine-riverWall intersections')
        if( (len(riverWalls)>0) and (len(breakLines)>0)):
            krw = list(riverWalls.keys())
            kbl = list(breakLines.keys())
            for i in range(len(krw)):
                n1 = krw[i]
                for j in range(len(kbl)):
                    n2 = kbl[j]
                    # Convert breaklines to wkb format
                    rw1 = ListPts2Wkb(riverWalls[n1],geometry_type='line')
                    bw2 = ListPts2Wkb(breakLines[n2],geometry_type='line')
                    # Add intersection points
                    rw1, bw2 = addIntersectionPtsToLines(rw1, bw2,\
                                    point_movement_threshold=point_movement_threshold,\
                                    verbose=verbose, nameFlag=n1+' intersects '+ n2)
                    riverWalls[n1] = Wkb2ListPts(rw1)
                    breakLines[n2] = Wkb2ListPts(bw2)
                    
    
        # Clean intersections of bounding polygon and riverwalls
        if(verbose): 
            print('Cleaning bounding_poly-riverWall intersections')
        if( (len(riverWalls)>0)):
            krw = list(riverWalls.keys())
            for i in range(len(krw)):
                n1 = krw[i]
                # Convert breaklines to wkb format
                rw1 = ListPts2Wkb(riverWalls[n1],geometry_type='line')
                bp2 = ListPts2Wkb(bounding_polygon,geometry_type='line', appendFirstOnEnd=True)
                # Add intersection points
                rw1, bp2 = addIntersectionPtsToLines(rw1, bp2,\
                                point_movement_threshold=point_movement_threshold,\
                                verbose=verbose, nameFlag='Bounding Pol intersects '+ n1)
                riverWalls[n1] = Wkb2ListPts(rw1)
                # Since the bounding polygon is a loop, the first/last points are the same
                # If one of these was moved, the other should be moved too. Since we
                # will drop the last bounding_polygon point, we only need to worry about the first
                bounding_polygon = Wkb2ListPts(bp2,removeLast=False)
                if(bounding_polygon[-1] is not bounding_polygon[0]):
                    bounding_polygon[0] = bounding_polygon[-1]
                # Drop the last point
                bounding_polygon = bounding_polygon[:-1]
    
        # Clean intersections of bounding polygon and breaklines
        if(verbose):
            print('Cleaning bounding_poly-breaklines intersections')
        if( (len(breakLines)>0)):
            kbl = list(breakLines.keys())
            for i in range(len(kbl)):
                n1 = kbl[i]
                # Convert breaklines to wkb format
                bl1 = ListPts2Wkb(breakLines[n1],geometry_type='line')
                bp2 = ListPts2Wkb(bounding_polygon,geometry_type='line', appendFirstOnEnd=True)
                # Add intersection points
                bl1, bp2 = addIntersectionPtsToLines(bl1, bp2,\
                                point_movement_threshold=point_movement_threshold,
                                verbose=verbose, nameFlag='Bounding Pol intersects '+n1)
                breakLines[n1] = Wkb2ListPts(bl1)
                # Since the bounding polygon is a loop, the first/last points are the same
                # If one of these was moved, the other should be moved too. Since we
                # will drop the last bp2 point, we only need to worry about the first
                bounding_polygon = Wkb2ListPts(bp2,removeLast=False)
                if(bounding_polygon[-1] is not bounding_polygon[0]):
                    bounding_polygon[0] = bounding_polygon[-1]
                # Drop the last point
                bounding_polygon = bounding_polygon[:-1]

        # Remove the extra 0.0 from bounding polygon [this cannot have 3 coordinates] 
        bounding_polygon = [ bounding_polygon[i][0:2] for i in range(len(bounding_polygon))]
        # Same for breaklines [although might not matter]
        for n1 in list(breakLines.keys()):
            breakLines[n1] = [breakLines[n1][i][0:2] for i in range(len(breakLines[n1]))]

        # Check that all mesh-boundary points are inside the bounding polygon
        from anuga.geometry.polygon import outside_polygon
        for blCat in [riverWalls, breakLines]:
            for n1 in list(blCat.keys()):
                l = len(blCat[n1])
                # Test every point -- means we can strip 3rd coordinate if needed
                for j in range(l):
                    isOut = outside_polygon(blCat[n1][j][0:2], bounding_polygon)
                    if(len(isOut)>0):
                        msg = 'Breakline/riverwall point '+str(blCat[n1][j][0:2]) +' on '+ n1+\
                            ' is outside the bounding polygon.\n'+\
                            'Check that it exceeds the bounding polygon'+\
                            ' by a distance < point_movement_threshold \n'+\
                            ' so it can be moved back onto the polygon'
                        print('Polygon\n ')
                        print(bounding_polygon)
                        print('Line \n')
                        print(blCat[n1])
                        raise Exception(msg)

        return [bounding_polygon, breakLines, riverWalls]

    ###################################################################

    def readRegionPtAreas(shapefile, convert_length_to_area=False):
        """
            Read a point shapefile to define the ANUGA mesh resoutions. 

            MUST HAVE A SINGLE ATTRIBUTE REPRESENTING THE LENGTHS OF TRIANGLES IN
             REGIONS

            INPUT: shapefile -- name of the input shapefile (or a 3 column csv file with x, y, res)
                   convert_length_to_area -- if True, res values are assumed to
                          represent triangle side lengths, and are converted to areas with 0.5*res0*res0
                          Note that this might not ensure that the max triangle side length really is res0, but
                          it will be of similar magnitude
                          If False, attribute values are assumed to represent triangle areas

            OUTPUT: list of the form  [ [x0,y0,res0], [x1, y1, res1], ...]
        """

        if shapefile[-4:] == '.shp':
            ptData = readShpPtsAndAttributes(shapefile)

            # Must have only 1 attribute
            if not (len(ptData[2]) == 1):
                msg = 'Shapefile ' + shapefile + ' does not contain exactly 1 ' +\
                      'attribute, so cannot be read as a regionPointArea'
                raise Exception(msg)

        else:
            # Assume file is csv
            f = open(shapefile)
            firstLine = f.readline()
            f.close()
            # If the top line has any letters, assume it is a header
            hasLetters = any(c.isalpha() for c in firstLine)
            ptData = numpy.genfromtxt(shapefile, delimiter=",", skip_header=int(hasLetters))
            if len(ptData[0,:]) != 3:
                msg = 'Region point areas text file must have exactly 3 columns separated by commas'
                raise Exception(msg)

            ptData = [ptData[:,0:2].tolist(), ptData[:,2].tolist()]

        # Convert to the required format
        numPts = len(ptData[0])
        outData = []
        for i in range(numPts):
            if(convert_length_to_area):
                newDat = [ptData[0][i][0], ptData[0][i][1], 0.5*float(ptData[1][i])**2]
            else:
                newDat = [ptData[0][i][0], ptData[0][i][1], float(ptData[1][i])]
            outData.append(newDat)

        return outData

    #########################################
    def readListOfBreakLines(fileList):
        """
            Take a list with the names of shapefiles or anuga_polygon csv files

            They are assumed to be '2D breaklines', so we just read their
                coordinates into a dict with their names

            Read them in

            INPUT: fileList -- a list of shapefile and/or anuga_polygon csv filenames 
                               [e.g. from glob.glob('GIS/Breaklines/*.shp')]

            OUTPUT: dictionary with breaklines [filenames are keys]
        """

        allBreakLines = {}
        for shapefile in fileList:
            allBreakLines[shapefile] = read_polygon(shapefile, 
                close_polygon_shapefiles=True)

        return allBreakLines

    #########################################
    def readListOfRiverWalls(rwfileList):
        """
            Take a list with the names of riverwall input files 
            [should be comma-separated x,y,elevation files]

            The input file can OPTIONALLY have a first line defining some
            hydraulic parameters. A valid example is

            Qfactor: 1.5, s1: 0.94 
            200., 300., 0.5
            300., 400., 0.7
            ....and so on..

            Read their coordinates into a dict with their names, read for use by ANUGA

            INPUT: rwfileList -- a list of riverwall filenames 
                        [e.g. from glob.glob('GIS/RiverWalls/*.csv')]

            OUTPUT: 
                dictionary with riverwalls [filenames are keys] AND
                dictionary with hydraulic parameters [filenames are keys]
        """
        import numpy

        allRiverWalls = {}
        allRiverWallPar = {}
        for rwfile in rwfileList:
            f = open(rwfile)
            firstLine = f.readline()
            f.close()
            # If the top line has any letters, assume it is a hydraulic
            # variables line
            hasLetters = any(c.isalpha() for c in firstLine)
            if(not hasLetters):
                allRiverWalls[rwfile] = \
                    numpy.genfromtxt(rwfile,delimiter=",").tolist()
                allRiverWallPar[rwfile] = {}
            else:
                # Get the wall geometry
                allRiverWalls[rwfile] = \
                    numpy.genfromtxt(rwfile,delimiter=",",skip_header=1).tolist()
                # Get the hydraulic par
                firstLine = firstLine.replace(' ', '') # No whitespace
                wallPar = firstLine.split(',')
                allRiverWallPar[rwfile] = {}
                for wp in wallPar:
                    keyNameValue = wp.split(':')
                    allRiverWallPar[rwfile][keyNameValue[0]] = \
                        float(keyNameValue[1])
        
        return allRiverWalls, allRiverWallPar

    ############################################################################

    def combine_breakLines_and_riverWalls_for_mesh(breakLines, riverWalls):
        """
        Combine breaklines and riverwalls for input in mesh generation,
            ensuring both have 2 coordinates only
        """
        mesh_breakLines=list(riverWalls.values())+list(breakLines.values())
        for i in range(len(mesh_breakLines)):
            mesh_breakLines[i] =\
             [mesh_breakLines[i][j][0:2] for j in range(len(mesh_breakLines[i]))]
        return mesh_breakLines
    
    ############################################################################
    def polygon_from_matching_breaklines(pattern,breakLinesIn, reverse2nd=None):
        """ We sometimes have breaklines defining 2 edges of a channel,
            and wish to make a polygon from them

            Can do this with the current function

            INPUTS: pattern == character string containing pattern which 
                        is inside exactly 2 keys in breaklines

                    breakLinesIn = the breakLines dictionary

                    reverse2nd = True/False or None. Reverse the order of the 
                       2nd set of edges before making the polygon.
                       If None, then we compute the distance between the
                       first point on breakline1 and the first/last
                       points on breakline2, and reverse2nd if the
                           'distance from the first point' <
                           'distance from the last point'

            OUTPUT: Polygon made with the 2 breaklines
        """

        breakLines=copy.copy(breakLinesIn)
        bk=list(breakLines.keys())

        # They can be pathnames from glob, and sometimes / and \\ get mixed up
        # Fix that here
        pattern_norm = os.path.normpath(pattern)
        bk_norm = [ os.path.normpath(bk_i) for bk_i in bk ]

        matchers=matchInds(pattern_norm, bk_norm)

        if(len(matchers)==0):
            msg = 'Cannot match ' + pattern + ' in breaklines file names'
            raise Exception(msg)

        if(len(matchers)!=2):
            print('Need exactly 2 matches, but pattern matched these', bk[matchers])

        # There are 2 matches

        if(reverse2nd is None):
            # Automatically compute whether we should reverse the 2nd breakline
            bl1_0=breakLines[bk[matchers[0]]][0]
            bl2_0=breakLines[bk[matchers[1]]][0]
            bl2_N=breakLines[bk[matchers[1]]][-1]
            d0 = ((bl1_0[0]-bl2_0[0])**2 + (bl1_0[1]-bl2_0[1])**2)
            dN = ((bl1_0[0]-bl2_N[0])**2 + (bl1_0[1]-bl2_N[1])**2)
            if(d0<dN):
                reverse2nd = True
            else:
                reverse2nd = False

        if(reverse2nd):
            breakLines[bk[matchers[1]]].reverse()
        polyOut=breakLines[bk[matchers[0]]] +  breakLines[bk[matchers[1]]]

        # If the breakLines values have > 2 entries (i.e. like for riverwalls),
        # remove the third
        if(len(polyOut[0])>2):
            polyOut=[polyOut[i][0:2] for i in range(len(polyOut))]

        return polyOut
    ###################

else: # gdal_available == False
    msg='Failed to import gdal/ogr modules --'\
        + 'perhaps gdal python interface is not installed.'



    def readShp_1PolyGeo(shapefile, dropLast=True):
        raise ImportError(msg)
    
    def readShp_1LineGeo(shapefile):
        raise ImportError(msg)

    def read_csv_optional_header(filename):
        raise ImportError(msg)
    
    def read_polygon(filename):
        raise ImportError(msg)
    
    def readShpPtsAndAttributes(shapefile):
        raise ImportError(msg)
    
    def read_points(filename):
        raise ImportError(msg)
    
    def ListPts2Wkb( ptsIn, geometry_type='line', appendFirstOnEnd=None):
        raise ImportError(msg)
    
    def Wkb2ListPts(wkb_geo, removeLast=False, drop_third_dimension=False):
        raise ImportError(msg)
    
    def compute_squared_distance_to_segment(pt, line):
        raise ImportError(msg)
    
    def find_nearest_segment(pt, segments):
        raise ImportError(msg)
    
    def shift_point_on_line(pt, lineIn, nearest_segment_index):
        raise ImportError(msg)
    
    def insert_intersection_point(intersectionPt, line_pts, 
                                  point_movement_threshold,verbose=False):
        raise ImportError(msg)

    def check_polygon_is_small(intersection, buf, tol2=100.):
        raise ImportError(msg)
    
    def addIntersectionPtsToLines(L1,L2, point_movement_threshold=0.0, 
                                  buf=1.0e-06, tol2 = 100,
                                  verbose=True, nameFlag=''):
        raise ImportError(msg)
   
    def getRasterExtent(rasterFile, asPolygon=False): 
        raise ImportError(msg)

    def rasterValuesAtPoints(xy, rasterFile, band=1):
        raise ImportError(msg)
    
    
    def gridPointsInPolygon(polygon, approx_grid_spacing=[1.,1.], eps=1.0e-06):
        raise ImportError(msg)
    

    def matchInds(pattern, stringList):
        raise ImportError(msg)
    
    
    def add_intersections_to_domain_features(bounding_polygonIn,
                breakLinesIn={ }, riverWallsIn={ }, point_movement_threshold=0.,
                verbose=True):
        raise ImportError(msg)
    
    
    def readRegionPtAreas(shapefile, convert_length_to_area=False):
        raise ImportError(msg)
    
    def readListOfBreakLines(shapefileList):
        raise ImportError(msg)

    def combine_breakLines_and_riverwalls_for_mesh(breakLines, riverWalls):
        raise ImportError(msg)
    
    def polygon_from_matching_breaklines(pattern,breakLinesIn, reverse2nd=None):
        raise ImportError(msg)
    ###################    


