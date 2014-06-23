"""

Routines to ease the import of spatial data to ANUGA

Key routines:
    readShp_1PolyGeo, readShp_1LineGeo -- read SIMPLE shapefile geometries into ANUGA as a list of points.
                                          The supported geometries are quite restricted, see help
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

    readRegionPtAreas -- read a shapefile containin regionPtAreas -- xy coordinates + 1 attribute, which is
                         the mesh triangle side length (or area) limit. Can be passed as regionPtAreas in
                         the mesh generation stage.

    readListOfBreakLines -- takes a list of shapefile names, each containing a single line geometry, and reads
                            it in as a dictionary of breaklines. 

    polygon_from_matching_breaklines -- given a pattern (string) matching exactly 2 breaklines in a directory,
                                        convert them to a single polygon. This is useful to e.g. make
                                        a polygon defining the extent of a channel that is defined from 2 breaklines

    
"""
import sys
import os
import copy
import struct
import numpy
from matplotlib import path

try:
    import gdal
    import ogr
    gdal_available = True
    #import osr # Not needed here but important in general
except ImportError, err:
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
        layer=dataSrc.GetLayer()
       
        # Need a single polygon 
        try:
            assert(len(layer)==1)
        except:
            print shapefile
        
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
        layer=dataSrc.GetLayer()
       
        # Need a single line 
        try:
            assert(len(layer)==1)
        except:
            print shapefile
        
        line_all=[]
        for feature in layer:
            line=feature.GetGeometryRef().GetPoints()
            line=[list(pts) for pts in line]
            line_all.extend(line)
    
        return line_all
           
    ####################
    
    def readShpPtsAndAttributes(shapefile):
        """
            Read a point shapefile with an attribute table into a list
    
            INPUT: shapefile -- name of shapefile to read
            
            OUTPUT: List with [ list_of_points, list_of_attributes, names_of_attributes]
        """
        
        driver=ogr.GetDriverByName("ESRI Shapefile")
        dataSrc=driver.Open(shapefile, 0)
        layer=dataSrc.GetLayer()
    
        pts=[]
        attribute=[]
        for i, feature in enumerate(layer):
            if(i==0):
                attributeNames=feature.keys()
            pt=feature.GetGeometryRef().GetPoints()
            pt=[list(p) for p in pt]
            pts.extend(pt)
            att=[feature[i] for i in attributeNames]
            attribute.extend(att)
    
        return [pts, attribute, attributeNames]
    
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
            raise Exception, "Type must be either 'point' or 'line' or 'polygon'"
         
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
                raise Exception, 'Points must be either 2 or 3 dimensional'
    
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
            raise Exception, 'Geometry type not supported'
        
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
        seg_unitVec_x=float(p1[0]-p0[0])
        seg_unitVec_y=float(p1[1]-p0[1])
        segLen=(seg_unitVec_x**2+seg_unitVec_y**2)**0.5
        if(segLen==0.):
            #print line
            #print 'Pt'
            #print pt
            raise Exception, 'Line has repeated points: Line %s Pt %s' % (str(line),str(pt))
        seg_unitVec_x=seg_unitVec_x/segLen
        seg_unitVec_y=seg_unitVec_y/segLen
        #
        # Get vector from pt to p0 
        pt_p0_vec_x=float(pt[0]-p0[0])
        pt_p0_vec_y=float(pt[1]-p0[1])
        pt_p0_vec_len_squared=(pt_p0_vec_x**2 + pt_p0_vec_y**2)
        # Get dot product of above vector with unit vector
        pt_dot_segUnitVec=(pt_p0_vec_x)*seg_unitVec_x+(pt_p0_vec_y)*seg_unitVec_y
        #
        if( (pt_dot_segUnitVec<segLen) and (pt_dot_segUnitVec > 0.)):
            # The nearest point on the line is actually between p0 and p1, so we have a 'real' candidate
            # Get distance^2
            output = pt_p0_vec_len_squared - pt_dot_segUnitVec**2.
        else:
            # Distance is the min distance from p0 and p1. 
            output = min( pt_p0_vec_len_squared,  (float(pt[0]-p1[0])**2+float(pt[1]-p1[1])**2))
        if(output < -1.0e-06):
            raise Exception, 'round-off in compute_squared_distance_to_segment'
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
            
            OUTPUT: The squared distance, and the index i of the segment [x_i,y_i],[x_i+1,y_i+1] in segments which is closest to pt
        """
        ll=len(segments)
        if(ll<=1):
            raise Exception, 'Segments must have length > 1 in find_nearest_segment'
       
        ptDist_sq=numpy.zeros(ll-1) # Hold the squared distance from the point to the line segment
        for i in range(len(segments)-1):
            # Compute distance from segment
            ptDist_sq[i]=compute_squared_distance_to_segment(pt, [segments[i],segments[i+1]])
        
        return [ptDist_sq.min(), ptDist_sq.argmin()]
    
    ######################################################
    
    def shift_point_on_line(pt, lineIn, nearest_segment_index):
        """
            Support pt is a point, which is near to the 'nearest_segment_index' segment of the line 'lineIn'
    
            This routine moves the nearest end point of that segment on line to pt.
    
            INPUTS: pt -- [x,y] point
                    lineIn -- [ [x0, y0], [x1, y1], ..., [xN,yN]] 
                    nearest_segment_index = index where the distance of pt to the line from [x_i,y_i] to [x_i+1,y_i+1] is minimum
    
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
    def insert_intersection_point(intersectionPt, line_pts, point_movement_threshold):
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
        endPt_Dist_Sq=min( ( (p0[0]-iP[0])**2 + (p0[1]-iP[1])**2), ( (p1[0]-iP[0])**2 + (p1[1]-iP[1])**2))
        #
        if(endPt_Dist_Sq>point_movement_threshold**2):
            # Insert the intersection point. We do this in a tricky way to
            # account for the possibility of L1_pts having > 2 coordinates
            # (riverWalls)
            dummyPt=copy.copy(L1_pts[tmp[1]])
            L1_pts.insert(tmp[1]+1,dummyPt) 
            L1_pts[tmp[1]+1][0]=iP[0]
            L1_pts[tmp[1]+1][1]=iP[1]
            if(len(L1_pts[tmp[1]+1])==3):
                # Treat 3rd coordinate
                # Find distance of inserted point from neighbours, and 
                # Set 3rd coordinate as distance-weighted average of the others
                d0=((L1_pts[tmp[1]][0]-L1_pts[tmp[1]+1][0])**2.+\
                   (L1_pts[tmp[1]][0]-L1_pts[tmp[1]+1][1])**2.)**0.5
                d1=((L1_pts[tmp[1]+2][0]-L1_pts[tmp[1]+1][0])**2.+\
                   (L1_pts[tmp[1]+2][1]-L1_pts[tmp[1]+1][1])**2.)**0.5
                L1_pts[tmp[1]+1][2] = (d0*L1_pts[tmp[1]+2][2] + d1*L1_pts[tmp[1]][2])/(d0+d1)
    
        else:
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
    
    def addIntersectionPtsToLines(L1,L2, point_movement_threshold=0.0, buf=1.0e-06, tol2 = 100,
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
    
                     nameFlag = This will be printed if intersection occurs. Convenient way to display intersecting filenames
            
            OUTPUTS: L1,L2 with intersection points added in the right places
        """
    
        if(L1.Intersects(L2)):
            # Get points on the lines 
            L1_pts=Wkb2ListPts(L1) #L1.GetPoints()
            L2_pts=Wkb2ListPts(L2) #L2.GetPoints()
            
            # Buffer lines by a small amount
            L1_buf=L1.Buffer(buf)
            L2_buf=L2.Buffer(buf)
                
            # Get intersection point[s]    
            L1_L2_intersect=L1_buf.Intersection(L2_buf)
            if(L1_L2_intersect.GetGeometryCount()==1):
                if(not check_polygon_is_small(L1_L2_intersect, buf, tol2)):
                    #print L1_L2_intersect.GetEnvelope()
                    raise Exception, 'line intersection is not allowed. Envelope %s '% str(L1_L2_intersect.GetEnvelope())
                # Seems to need special treatment with only 1 intersection point
                intersectionPts=[L1_L2_intersect.Centroid().GetPoint()]
            else:
                intersectionPts=[]
                for feature in L1_L2_intersect:
                    if(not check_polygon_is_small(feature, buf, tol2)):
                        print feature.GetEnvelope()
                        raise Exception, 'line intersection is not allowed'
                    intersectionPts.append(feature.Centroid().GetPoint())
    
            if(verbose):
                print nameFlag
                print '    Treating intersections in ', len(intersectionPts) , ' locations'
                print intersectionPts
    
            # Insert the points into the line segments
            for i in range(len(intersectionPts)):
                L1_pts = insert_intersection_point(intersectionPts[i], L1_pts, point_movement_threshold)
                L2_pts = insert_intersection_point(intersectionPts[i], L2_pts, point_movement_threshold)
                
            # Convert to the input format
            L1_pts=ListPts2Wkb(L1_pts,geometry_type='line')
            L2_pts=ListPts2Wkb(L2_pts,geometry_type='line')
    
            return [L1_pts, L2_pts]
        else:
            return [L1, L2]
    
    ###########################################################
    def rasterValuesAtPoints(xy, rasterFile, band=1):
        """
            Get raster values at point locations.
            Can be used to e.g. set quantity values
       
            INPUT: 
            xy = numpy array with point locations
            rasterFile = Filename of the gdal-compatible raster
            band = band of the raster to get
    
            OUTPUT:
            1d numpy array with raster values at xy
    
        """
        # Raster info
        raster = gdal.Open(rasterFile)
        rasterBand=raster.GetRasterBand(band)
        rasterBandType=gdal.GetDataTypeName(rasterBand.DataType)
    
        # Projection info
        transform=raster.GetGeoTransform()
        xOrigin = transform[0]
        yOrigin = transform[3]
        pixelWidth = transform[1]
        pixelHeight = transform[5] # Negative
        
        # Get coordinates in pixel values
        px = ((xy[:,0] - xOrigin) / pixelWidth).astype(int) #x 
        py = ((xy[:,1] - yOrigin) / pixelHeight).astype(int) #y 
      
        # Hold elevation 
        elev=px*0. 
    
        # Get the right character for struct.unpack
        if (rasterBandType == 'Int16'):
            CtypeName='h'
        elif (rasterBandType == 'Float32'):
            CtypeName='f'
        elif (rasterBandType == 'Byte'):
            CtypeName='B'
        else:
            print 'unrecognized DataType:', gdal.GetDataTypeName(band.DataType)
            print 'You might need to edit this code to read the data type'
            raise Exception, 'Stopping'
    
        # Get values -- seems we have to loop, but it is efficient enough
        for i in range(len(px)):
            xC=int(px[i])
            yC=int(py[i])
            structval=rasterBand.ReadRaster(xC,yC,1,1,buf_type=rasterBand.DataType)
            elev[i] = struct.unpack(CtypeName, structval)[0]
    
        return elev
    
    
    def gridPointsInPolygon(polygon, approx_grid_spacing=[1.,1.], eps=1.0e-06):
        """
            Get a 'grid' of points inside a polygon. Could be used with rasterValuesAtPoints to 
               get a range of raster values inside a polygon
    
            Approach: A 'trial-grid' of points is created which is 'almost' covering the range of the polygon (xmin-xmax,ymin-ymax). 
    
                      (Actually it is just inside this region, to avoid polygon-boundary issues, see below)
    
                      Then we find those points which are actually inside the polygon.
    
                      The x/y point spacing on the trial-grid will be close to
                      approx_grid_spacing, but we ensure there are at least 3x3 points on the trial grid.
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
        polygonArr=numpy.array(polygon)
        poly_xmin=polygonArr[:,0].min()
        poly_xmax=polygonArr[:,0].max()
        poly_ymin=polygonArr[:,1].min()
        poly_ymax=polygonArr[:,1].max()
    
        # Make a 'grid' of points which covers the polygon
        xGridCount=max( numpy.ceil( (poly_xmax-poly_xmin)/approx_grid_spacing[0]+1. ).astype(int), 3)
        R=(poly_xmax-poly_xmin)*eps
        Xvals=numpy.linspace(poly_xmin+R,poly_xmax-R, xGridCount)
        yGridCount=max( numpy.ceil( (poly_ymax-poly_ymin)/approx_grid_spacing[1]+1. ).astype(int), 3)
        R=(poly_ymax-poly_ymin)*eps
        Yvals=numpy.linspace(poly_ymin+R,poly_ymax-R, yGridCount)
    
        xGrid,yGrid=numpy.meshgrid(Xvals,Yvals)
        Grid=numpy.vstack([xGrid.flatten(),yGrid.flatten()]).transpose()
    
        # Now find the xy values which are actually inside the polygon
        polpath=path.Path(polygonArr)
        keepers=[]
        for i in range(len(Grid[:,0])):
           if(polpath.contains_point(Grid[i,:])):
                keepers=keepers+[i]
        #import pdb
        #pdb.set_trace()
        xyInside=Grid[keepers,:]
    
        return(xyInside)
    
    #########################################################################
    # Function to search for pattern matches in a string (turns out to be useful)
    def matchInds(pattern, stringList):
        """
            Find indices in stringList which match pattern
        """
        #matches=[ (pattern in stringList[i]) for i in range(len(stringList))]
        matches=[]
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
    
    
    def add_intersections_to_domain_features(bounding_polygonIn,
                breakLinesIn={ }, riverWallsIn={ }, point_movement_threshold=0.,
                verbose=True):
        """
            If bounding polygon / breaklines /riverwalls intersect with each other, then
            add intersection points.
    
            INPUTS:
                bounding_polygonIn -- the bounding polygon in ANUGA format
                breakLinesIn -- the breaklines dictionary
                riverWallsIn -- the riverWalls dictionary
    
                point_movement_threshold -- if points on lines are < this distance
                        from intersection points, then they are replaced with the intersection point.
                        This can prevent overly close points from breaking the mesh generation
    
            OUTPUT:
                List with bounding_polygon,breakLines,riverwalls
        """
    
        bounding_polygon=copy.copy(bounding_polygonIn)
        breakLines=copy.copy(breakLinesIn)
        riverWalls=copy.copy(riverWallsIn)
    
        # Clean intersections of breakLines with itself
        if(verbose): 
            print 'Cleaning breakline intersections'
        if(len(breakLines)>0):
            kbl=breakLines.keys()
            for i in range(len(kbl)):
                n1=kbl[i]
                for j in range(len(kbl)):
                    if(i>=j):
                        continue
                    n2=kbl[j]
                    # Convert breaklines to wkb format
                    bl1=ListPts2Wkb(breakLines[n1],geometry_type='line')
                    bl2=ListPts2Wkb(breakLines[n2],geometry_type='line')
                    # Add intersection points
                    bl1, bl2 =addIntersectionPtsToLines(bl1, bl2,\
                                    point_movement_threshold=point_movement_threshold,
                                    verbose=verbose, nameFlag=n1+' intersects '+ n2)
                    breakLines[n1]=Wkb2ListPts(bl1)
                    breakLines[n2]=Wkb2ListPts(bl2)
                    
    
        # Clean intersections of riverWalls with itself
        if(verbose): 
            print 'Cleaning riverWall intersections'
        if(len(riverWalls)>0):
            krw=riverWalls.keys()
            for i in range(len(krw)):
                n1=krw[i]
                for j in range(len(krw)):
                    if(i>=j):
                        continue 
                    n2=krw[j]
                    # Convert breaklines to wkb format
                    rw1=ListPts2Wkb(riverWalls[n1],geometry_type='line')
                    rw2=ListPts2Wkb(riverWalls[n2],geometry_type='line')
                    # Add intersection points
                    rw1, rw2 =addIntersectionPtsToLines(rw1, rw2,\
                                    point_movement_threshold=point_movement_threshold,\
                                    verbose=verbose, nameFlag=n1+' intersects '+ n2)
                    riverWalls[n1]=Wkb2ListPts(rw1)
                    riverWalls[n2]=Wkb2ListPts(rw2)
    
        # Clean intersections of breaklines with riverwalls
        if(verbose): 
            print 'Cleaning breakLine-riverWall intersections'
        if( (len(riverWalls)>0) and (len(breakLines)>0)):
            krw=riverWalls.keys()
            kbl=breakLines.keys()
            for i in range(len(krw)):
                n1=krw[i]
                for j in range(len(kbl)):
                    n2=kbl[j]
                    # Convert breaklines to wkb format
                    rw1=ListPts2Wkb(riverWalls[n1],geometry_type='line')
                    bw2=ListPts2Wkb(breakLines[n2],geometry_type='line')
                    # Add intersection points
                    rw1, bw2 =addIntersectionPtsToLines(rw1, bw2,\
                                    point_movement_threshold=point_movement_threshold,\
                                    verbose=verbose, nameFlag=n1+' intersects '+ n2)
                    riverWalls[n1]=Wkb2ListPts(rw1)
                    breakLines[n2]=Wkb2ListPts(bw2)
                    
    
        # Clean intersections of bounding polygon and riverwalls
        if(verbose): 
            print 'Cleaning bounding_poly-riverWall intersections'
        if( (len(riverWalls)>0)):
            krw=riverWalls.keys()
            for i in range(len(krw)):
                n1=krw[i]
                # Convert breaklines to wkb format
                rw1=ListPts2Wkb(riverWalls[n1],geometry_type='line')
                bp2=ListPts2Wkb(bounding_polygon,geometry_type='line', appendFirstOnEnd=True)
                # Add intersection points
                rw1, bp2 =addIntersectionPtsToLines(rw1, bp2,\
                                point_movement_threshold=point_movement_threshold,\
                                verbose=verbose, nameFlag='Bounding Pol intersects '+ n1)
                riverWalls[n1]=Wkb2ListPts(rw1)
                bounding_polygon=Wkb2ListPts(bp2,removeLast=True)
    
        # Clean intersections of bounding polygon and breaklines
        if(verbose):
            print 'Cleaning bounding_poly-breaklines intersections'
        if( (len(breakLines)>0)):
            kbl=breakLines.keys()
            for i in range(len(kbl)):
                n1=kbl[i]
                # Convert breaklines to wkb format
                bl1=ListPts2Wkb(breakLines[n1],geometry_type='line')
                bp2=ListPts2Wkb(bounding_polygon,geometry_type='line', appendFirstOnEnd=True)
                # Add intersection points
                bl1, bp2 =addIntersectionPtsToLines(bl1, bp2,\
                                point_movement_threshold=point_movement_threshold,
                                verbose=verbose, nameFlag='Bounding Pol intersects '+n1)
                breakLines[n1]=Wkb2ListPts(bl1)
                bounding_polygon=Wkb2ListPts(bp2, removeLast=True)
      
        # Remove the extra 0.0 from bounding polygon [this cannot have 3 coordinates] 
        bounding_polygon = [ bounding_polygon[i][0:2] for i in range(len(bounding_polygon))]
        # Same for breaklines [although might not matter]
        for n1 in breakLines.keys():
            breakLines[n1] = [breakLines[n1][i][0:2] for i in range(len(breakLines[n1]))]
     
        return [bounding_polygon, breakLines, riverWalls]
    
    ###################################################################
    
    def readRegionPtAreas(shapefile, convert_length_to_area=False):
        """
            Read a point shapefile to define the ANUGA mesh resoutions. 
    
            MUST HAVE A SINGLE ATTRIBUTE REPRESENTING THE LENGTHS OF TRIANGLES IN
             REGIONS
    
            INPUT: shapefile -- name of the input shapefile
                   convert_length_to_area -- if True, res values are assumed to
                          represent triangle side lengths, and are converted to areas with 0.5*res0*res0
                          Note that this might not ensure that the max triangle side length really is res0, but
                          it will be of similar magnitude
                          If False, attribute values are assumed to represent triangle areas
    
            OUTPUT: list of the form  [ [x0,y0,res0], [x1, y1, res1], ...]
        """
    
        ptData=readShpPtsAndAttributes(shapefile)
    
        # Must have only 1 attribute
        assert len(ptData[2])==1
    
        numPts=len(ptData[0])
        outData=[]
        for i in range(numPts):
            if(convert_length_to_area):
                newDat=[ptData[0][i][0], ptData[0][i][1], 0.5*float(ptData[1][i])**2]
            else:
                newDat=[ptData[0][i][0], ptData[0][i][1], float(ptData[1][i])]
            outData.append(newDat)
    
        return outData
    
    #########################################
    def readListOfBreakLines(shapefileList):
        """
            Take a list with the names of shapefiles
        
            They are assumed to be '2D breaklines', so we just read their
                coordinates into a dict with their names
    
            Read them in
            
            INPUT: shapefileList -- a list of shapefile names [e.g. from glob.glob('GIS/Breaklines/*.shp')]
    
            OUTPUT: dictionary with breaklines [filenames are keys]
        """
    
        allBreakLines={}
        for shapefile in shapefileList:
            allBreakLines[shapefile]=readShp_1LineGeo(shapefile)
        
        return allBreakLines
    
    ############################################################################
    def polygon_from_matching_breaklines(pattern,breakLinesIn, reverse2nd=None):
        """
            We sometimes have breaklines defining 2 edges of a channel,
            and wish to make a polygon from them
    
            Can do this with the current function
    
            INPUTS: pattern == character string containing pattern which is inside exactly 2 keys in breaklines
        
                    breakLines = the breakLinesIn dictionary
                    
                    reverse2nd = True/False or None. Reverse the order of the 2nd set of edges before making the polygon.
                                 If None, then we compute the distance between the
                                 first point on breakline1 and the first/last points on breakline2, and reverse2nd if
                                 the 'distance from the first point' < 'distance from the last point'
        
            OUTPUT: Polygon made with the 2 breaklines
        """
    
        breakLines=copy.copy(breakLinesIn)
        bk=breakLines.keys()
        matchers=matchInds(pattern, bk)
    
        if(len(matchers)==0):
            msg = 'Cannot match ' + pattern + 'in breaklines file names'
            raise Exception, msg
    
        if(len(matchers)!=2):
            print 'Need exactly 2 matches, but pattern matched these', bk[matchers]
    
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
        #polyOut=polyOut+[polyOut[0]]
        return polyOut
    ###################
    
else: # gdal_available == False
    msg='Failed to import gdal/ogr modules --'\
        + 'perhaps gdal python interface is not installed.'


    
    def readShp_1PolyGeo(shapefile, dropLast=True):
        raise ImportError, msg
    
    def readShp_1LineGeo(shapefile):
        raise ImportError, msg
    
    def readShpPtsAndAttributes(shapefile):
        raise ImportError, msg
    
    def ListPts2Wkb( ptsIn, geometry_type='line', appendFirstOnEnd=None):
        raise ImportError, msg
    
    def Wkb2ListPts(wkb_geo, removeLast=False, drop_third_dimension=False):
        raise ImportError, msg
    
    def compute_squared_distance_to_segment(pt, line):
        raise ImportError, msg
    
    def find_nearest_segment(pt, segments):
        raise ImportError, msg
    
    def shift_point_on_line(pt, lineIn, nearest_segment_index):
        raise ImportError, msg
    
    def insert_intersection_point(intersectionPt, line_pts, point_movement_threshold):
        raise ImportError, msg

    def check_polygon_is_small(intersection, buf, tol2=100.):
        raise ImportError, msg
    
    def addIntersectionPtsToLines(L1,L2, point_movement_threshold=0.0, buf=1.0e-06, tol2 = 100,
                                  verbose=True, nameFlag=''):
        raise ImportError, msg
    

    def rasterValuesAtPoints(xy, rasterFile, band=1):
        raise ImportError, msg
    
    
    def gridPointsInPolygon(polygon, approx_grid_spacing=[1.,1.], eps=1.0e-06):
        raise ImportError, msg
    

    def matchInds(pattern, stringList):
        raise ImportError, msg
    
    
    def add_intersections_to_domain_features(bounding_polygonIn,
                breakLinesIn={ }, riverWallsIn={ }, point_movement_threshold=0.,
                verbose=True):
        raise ImportError, msg
    
    
    def readRegionPtAreas(shapefile, convert_length_to_area=False):
        raise ImportError, msg
    
    def readListOfBreakLines(shapefileList):
        raise ImportError, msg
    
    def polygon_from_matching_breaklines(pattern,breakLinesIn, reverse2nd=None):
        raise ImportError, msg
    ###################    


