#!/usr/bin/env python
"""

Unit tests for spatialInputUtil

FIXME: Need to extend coverage to reading shapefiles + some other operations

"""

import unittest
import anuga
import numpy
import os
import sys
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.utilities import plot_utils as util
from anuga.config import g
from anuga.utilities import spatialInputUtil as su

# this confuses nose        
# pull -v argument from command line
#verbose = anuga.get_args().verbose
verbose = False

class Test_spatialInputUtil(unittest.TestCase):
    """
	Test the spatialInput utilities
    """

    def setUp(self):
        pass

    def tearDown(self):
        for file in ['PointData_TestData.tif']:
            try:
                os.remove(file)
            except:
                pass        

    def make_me_a_tif(self):
        # We need to make a .tif to test some functions 
        # This does the job
        #
        from anuga.utilities import plot_utils as util
        #
        # Do it with Make_Geotif
        # Pick a domain that makes sense in EPSG:32756
        x=numpy.linspace(307000., 307100., 101)     
        y=numpy.linspace(6193000., 6193100., 101)
        xG,yG=numpy.meshgrid(x, y)
        xG=xG.flatten()
        yG=yG.flatten()
        # Surface is z=x+y
        fakeZ=xG-min(xG)+yG -min(yG)
        dataToGrid=numpy.vstack([xG,yG,fakeZ]).transpose()
        #    
        util.Make_Geotif(dataToGrid, output_quantities=['TestData'],
                         EPSG_CODE=32756, output_dir='.',CellSize=1.0)
   
    def test_compute_square_distance_to_segment(self):
        # Check that this is correct for a bunch of lines

        # pt = endpoint
        pt=[0., 0.]
        seg=[ [0., 0.], [1., 2.]]
        dist_sq=su.compute_squared_distance_to_segment(pt,seg)
        assert numpy.allclose(dist_sq, 0.)
        
        pt=[1., 2.]
        seg=[ [0., 0.], [1., 2.]]
        dist_sq=su.compute_squared_distance_to_segment(pt,seg)
        assert numpy.allclose(dist_sq, 0.)

        # pt on line
        pt=[0.25, 0.5]
        seg=[ [0., 0.], [1., 2.]]
        dist_sq=su.compute_squared_distance_to_segment(pt,seg)
        assert numpy.allclose(dist_sq, 0.)
       
        # pt nearest end point 
        pt=[1., 3.]
        seg=[ [0., 0.], [1., 2.]]
        dist_sq=su.compute_squared_distance_to_segment(pt,seg)
        assert numpy.allclose(dist_sq, 1.)
        
        # pt nearest end point 
        pt=[1.01, 30.]
        seg=[ [0., 0.], [1., 2.]]
        dist_sq=su.compute_squared_distance_to_segment(pt,seg)
        expectedDistSq=( (seg[1][0]-pt[0])**2 + (seg[1][1]-pt[1])**2)
        assert numpy.allclose(dist_sq, expectedDistSq)

        # pt nearest end point 
        pt=[0., -3.]
        seg=[ [0., 0.], [1., 2.]]
        dist_sq=su.compute_squared_distance_to_segment(pt,seg)
        assert numpy.allclose(dist_sq, 9.)
       
        # Should fail because seg has 2 points which are identical 
        pt=[0., -3.]
        seg=[ [0., 0.], [0., 0.]]
        def should_fail():
            dist_sq=su.compute_squared_distance_to_segment(pt,seg)
            return
        self.assertRaises(Exception, lambda: should_fail())
      
        # This one actually uses the perpendicular distance feature 
        pt=[1.0, 1.0]
        seg=[ [0., -2.], [10., 8.]] # y=x-2
        dist_sq=su.compute_squared_distance_to_segment(pt,seg)
        assert numpy.allclose(dist_sq, 2.)
        
        return

    def test_find_nearest_segment(self):

        # Example with colinear segments
        pt=[1.01, 30.]
        seg=[ [0., 0.],[0.01, 0.02], [1., 2.]]
        nearSeg=su.find_nearest_segment(pt,seg)
        expectedDistSq=( (seg[2][0]-pt[0])**2 + (seg[2][1]-pt[1])**2)
        assert numpy.allclose(nearSeg[0], expectedDistSq)
        assert nearSeg[1]==1
        
        # Another example 
        pt=[1.01, 30.]
        seg=[ [0., 0.],[-1.01, 0.02], [1., 2.]]
        nearSeg=su.find_nearest_segment(pt,seg)
        expectedDistSq=( (seg[2][0]-pt[0])**2 + (seg[2][1]-pt[1])**2)
        assert numpy.allclose(nearSeg[0], expectedDistSq)
        assert nearSeg[1]==1
        
        # This one actually uses the perpendicular distance feature 
        # And has a double-up end point
        pt=[1.0, 1.0]
        seg=[ [0., -2.], [10., 8.], [100., -100.], [0., -2.]] 
        nearSeg=su.find_nearest_segment(pt,seg)
        assert numpy.allclose(nearSeg[0], 2.)
        assert nearSeg[1]==0
        
        # This one permutes the order of the above
        pt=[1.0, 1.0]
        seg=[ [0., -2.], [100., -100.], [10., 8.], [0., -2.]] 
        nearSeg=su.find_nearest_segment(pt,seg)
        assert numpy.allclose(nearSeg[0], 2.)
        assert nearSeg[1]==2
        
        # Get a point on the line
        pt=[100., -100.]
        seg=[ [0., -2.], [100., -100.], [10., 8.], [0., -2.]] 
        nearSeg=su.find_nearest_segment(pt,seg)
        assert numpy.allclose(nearSeg[0], 0.)
        # Here either 0 or 1 is acceptable as the nearest segment start index, for us
        assert (nearSeg[1]==0 | nearSeg[1]==1)

        return

    def test_ListPts2Wbk_conversion(self):
        # Test conversion to-from Wkb

        # Only run test if python >= 2.7
        if sys.hexversion < 0x02070000: return

        seg=[ [0., -2., 1.], [100., -100., 2.], [10., 8., 3.], [0., -2., 1.]] 
        seg_Wkb=su.ListPts2Wkb(seg,geometry_type='line')    
        seg_Wkb_List=su.Wkb2ListPts(seg_Wkb)
        for i in range(len(seg)):
            assert(seg[i]==seg_Wkb_List[i])
       
        # Check we can drop 3rd dimension ok 
        seg=[ [0., -2., 1.], [100., -100., 2.], [10., 8., 3.], [0., -2., 1.]] 
        seg_Wkb=su.ListPts2Wkb(seg,geometry_type='line')    
        seg_Wkb_List=su.Wkb2ListPts(seg_Wkb, drop_third_dimension=True)
        for i in range(len(seg)):
            assert( (seg[i][0:2]==seg_Wkb_List[i][0:2]) and (len(seg_Wkb_List[i])==2))
       
        # Check we can append first point to end ok 
        seg=[ [0., -2., 1.], [100., -100., 2.], [10., 8., 3.]] 
        seg_Wkb=su.ListPts2Wkb(seg,geometry_type='line', appendFirstOnEnd=True)    
        seg_Wkb_List=su.Wkb2ListPts(seg_Wkb)
        for i in range(len(seg)):
            assert(seg[i]==seg_Wkb_List[i])
        # Check a point was added
        assert(len(seg)==(len(seg_Wkb_List)-1))
        # Check it was the first point
        assert(seg[0]==seg_Wkb_List[len(seg_Wkb_List)-1])

        # Check this happens automatically if we assume a polygon
        seg=[ [0., -2., 1.], [100., -100., 2.], [10., 8., 3.]] 
        seg_Wkb=su.ListPts2Wkb(seg,geometry_type='polygon')    
        seg_Wkb_List=su.Wkb2ListPts(seg_Wkb)
        for i in range(len(seg)):
            assert(seg[i]==seg_Wkb_List[i])
        # Check a point was added
        assert(len(seg)==(len(seg_Wkb_List)-1))
        # Check it was the first point
        assert(seg[0]==seg_Wkb_List[len(seg_Wkb_List)-1])

        # Check that all works for points as well
        seg=[ [0., -2., 1.], [100., -100., 2.], [10., 8., 3.]] 
        seg_Wkb=su.ListPts2Wkb(seg,geometry_type='point', appendFirstOnEnd=True)    
        seg_Wkb_List=su.Wkb2ListPts(seg_Wkb)
        for i in range(len(seg)):
            assert(seg[i]==seg_Wkb_List[i])
        # Check a point was added
        assert(len(seg)==(len(seg_Wkb_List)-1))
        # Check it was the first point
        assert(seg[0]==seg_Wkb_List[len(seg_Wkb_List)-1])

        return

    def test_shift_point_on_line(self):
        #
        seg=[ [0., -2., -998.], [100., -100., 2.], [10., 8., 3.]] 
        pt= [0., 1., -999.]
        newSeg=su.shift_point_on_line(pt, seg,0)
        # Should only have changed the first 2 coordinates
        assert(newSeg[0][0:2]==pt[0:2])
        assert(newSeg[0][2]==seg[0][2])
       
        # Try with a 2D pt -- should be the same 
        seg=[ [0., -2., -998.], [100., -100., 2.], [10., 8., 3.]] 
        pt= [0., 1.]
        newSeg=su.shift_point_on_line(pt, seg,0)
        # Should only have changed the first 2 coordinates
        assert(newSeg[0][0:2]==pt[0:2])
        assert(newSeg[0][2]==seg[0][2])
        
        return

    def test_addIntersectionPtsToLines(self):
        #
        # Make an intersection, check it works

        # Only run if python >= 2.7
        if sys.hexversion < 0x02070000: return
        
        
        seg1=[ [-10., 0.], [10., 0.]]
        seg2=[ [0., -10.], [0., 10.]]
        
        seg1_wkb=su.ListPts2Wkb(seg1,geometry_type='line')
        seg2_wkb=su.ListPts2Wkb(seg2,geometry_type='line')

        newSeg1,newSeg2=su.addIntersectionPtsToLines(seg1_wkb,seg2_wkb,verbose=verbose)
        newSeg1=su.Wkb2ListPts(newSeg1)
        newSeg2=su.Wkb2ListPts(newSeg2)

        assert(numpy.allclose(newSeg1[1][0], 0.0))
        assert(numpy.allclose(newSeg1[1][1], 0.0))
        assert(numpy.allclose(newSeg2[1][0], 0.0))
        assert(numpy.allclose(newSeg2[1][1], 0.0))
       
        # As above, but for geometries with a 3rd dimension
        seg1=[ [-10., 0., 3.], [10., 0., 3.]]
        seg2=[ [0., -10., 2.], [0., 10., 2.]]
        
        seg1_wkb=su.ListPts2Wkb(seg1,geometry_type='line')
        seg2_wkb=su.ListPts2Wkb(seg2,geometry_type='line')
        newSeg1,newSeg2=su.addIntersectionPtsToLines(seg1_wkb,seg2_wkb,verbose=verbose)
        newSeg1=su.Wkb2ListPts(newSeg1)
        newSeg2=su.Wkb2ListPts(newSeg2)

        assert(numpy.allclose(newSeg1[1][0], 0.0))
        assert(numpy.allclose(newSeg1[1][1], 0.0))
        assert(numpy.allclose(newSeg2[1][0], 0.0))
        assert(numpy.allclose(newSeg2[1][1], 0.0))
        # Check that they take the 3rd dimension from the nearby points
        assert(seg2[1][2]==2.)
        assert(seg1[1][2]==3.)
        
        # Make an intersection where segments have a nearby point, and check we
        # can move nearby points in the segments to the intersection
        seg1=[ [-10., 0.], [0.001, 0.001], [10., 0.]]
        seg2=[ [0., -10.], [-0.001, -0.001], [0., 10.]]
        
        seg1_wkb=su.ListPts2Wkb(seg1,geometry_type='line')
        seg2_wkb=su.ListPts2Wkb(seg2,geometry_type='line')

        newSeg1,newSeg2=su.addIntersectionPtsToLines(seg1_wkb,seg2_wkb,verbose=verbose,\
                                point_movement_threshold=0.01)
        newSeg1=su.Wkb2ListPts(newSeg1)
        newSeg2=su.Wkb2ListPts(newSeg2)
        # Now the 2nd point in both segments should be the intersection
        assert(numpy.allclose(newSeg1[1][0], newSeg2[1][0]))
        assert(numpy.allclose(newSeg1[1][1], newSeg2[1][1]))

        return

    def test_gridPointsInPolygon(self):
        #
        # Make a polygon, check the grid points exist where we need them

        ## Simple example -- no points on the trial grid are excluded
        myPoly=[[0., 10.], [10., 10.], [10., 0.], [0., 0.]]
        
        pip=su.gridPointsInPolygon(myPoly, approx_grid_spacing=[1.,1.])
        # There should be 121 points in total
        assert(pip.shape[0]==121)
        # The min/max x and y should be inside
        assert(min(pip[:,0])>0.)
        assert(max(pip[:,0])<10.)
        assert(min(pip[:,1])>0.)
        assert(max(pip[:,1])<10.)
       
        ## Example where some points on the trial grid would be excluded 
        myPoly=[[0., 10.], [10., 10.], [10., 0.]]
        pip=su.gridPointsInPolygon(myPoly, approx_grid_spacing=[1.,1.])
        # The min/max x and y should be inside
        assert(min(pip[:,0])>0.)
        assert(max(pip[:,0])<10.)
        assert(min(pip[:,1])>0.)
        assert(max(pip[:,1])<10.)
        # x+y should always be >= 10, since the line connecting [0,10] and
        # [10,0] is x+y=10
        assert(all(pip[:,0]+pip[:,1]>=10.))


    def test_getRasterExtent(self):
        self.make_me_a_tif()

        extentOut=su.getRasterExtent('PointData_TestData.tif')
        assert(numpy.allclose(extentOut[0], 307000.-0.5))
        assert(numpy.allclose(extentOut[1], 307100.+0.5))
        assert(numpy.allclose(extentOut[2], 6193000.-0.5))
        assert(numpy.allclose(extentOut[3], 6193100.+0.5))
        
        extentOut=su.getRasterExtent('PointData_TestData.tif',asPolygon=True)
        assert(numpy.allclose(extentOut[0][0], 307000.-0.5))
        assert(numpy.allclose(extentOut[3][0], 307000.-0.5))
        assert(numpy.allclose(extentOut[1][0], 307100.+0.5))
        assert(numpy.allclose(extentOut[2][0], 307100.+0.5))
        assert(numpy.allclose(extentOut[0][1], 6193000.-0.5))
        assert(numpy.allclose(extentOut[1][1], 6193000.-0.5))
        assert(numpy.allclose(extentOut[2][1], 6193100.+0.5))
        assert(numpy.allclose(extentOut[3][1], 6193100.+0.5))

    def test_rasterValuesAtPoints(self):
        # We need to make a .tif to test this function. 
        self.make_me_a_tif()

        # Get the range of the tif
        tifRange=su.getRasterExtent('PointData_TestData.tif')

        # Now try to get some point values -- note they will be rounded to the
        # nearest cell
        xA=numpy.array([0., 10.3, 50.9, 100.])+tifRange[0]+0.5
        yA=numpy.array([0., 20.1, 75.1, 100.])+tifRange[2]+0.5
        z_predicted=numpy.round(xA)+numpy.round(yA)-tifRange[0]-tifRange[2]-1.0
        InDat=numpy.vstack([xA,yA]).transpose()

        z_fitted=su.rasterValuesAtPoints(InDat, rasterFile='PointData_TestData.tif')
        try:
            assert(numpy.allclose(z_fitted,z_predicted)) 
        except:
            raise Exception, 'Error could be in rasterValuesAtPoints or in Make_Geotif'


        # Try with bilinear interpolation
        z_fitted=su.rasterValuesAtPoints(InDat, rasterFile='PointData_TestData.tif', 
            interpolation='bilinear')
        z_predicted = xA + yA - tifRange[0] - tifRange[2] - 1.0
        try:
            assert(numpy.allclose(z_fitted,z_predicted)) 
        except:
            raise Exception, 'Error could be in rasterValuesAtPoints or in Make_Geotif'


        return

    def test_add_intersections_to_domain_features(self):

        # Only run test if python >= 2.7
        if sys.hexversion < 0x02070000: return

        bounding_polygon=[ [0., 0.], [0., 10.], [10., 10.], [10., 0.]]

        breakLines={ 'bl1': [[-0.01, 5.],[10.01, 5.]],
                     'bl2': [[5., -0.01], [5., 10.01]] }
       
        # We make an erronious riverwall, which is colinear with a breakLine
        # This is not permitted 
        riverWalls={ 'rw1': [[-0.01, 8., 2.],[10.01, 4., 3.]],
                     'rw2': [[5., -0.01, 1.], [5., 10.01, 2.]] }

        def should_fail():
            # Should fail
            newBP, newBL, newRW=su.add_intersections_to_domain_features(bounding_polygon,\
                                    breakLines, riverWalls, point_movement_threshold=0.02,\
                                    verbose=verbose)
            return
        self.assertRaises(Exception, lambda: should_fail())

        #################################################################
        # Fix the riverwall, and it should work
        riverWalls={ 'rw1': [[-0.000001, 8., 2.],[10.0000001, 4., 3.]]
                      }
        # This should work
        newBP, newBL, newRW=su.add_intersections_to_domain_features(bounding_polygon,\
                                breakLines, riverWalls, point_movement_threshold=0.02,\
                                verbose=verbose)

        # There should be several new points on breakLines + riverWalls
        assert newBL['bl1'][1]==newBL['bl2'][1]
        assert newRW['rw1'][2][0:2]==newBL['bl1'][2]
        assert newRW['rw1'][1][0:2]==newBL['bl2'][2]
        # rw1 x/y coords are close to y=8-3/10*x
        #     x/z coords are close to z=2+x/10
        assert numpy.allclose(newRW['rw1'][1][2], newRW['rw1'][1][0]/10.+2)
        assert numpy.allclose(newRW['rw1'][2][2], newRW['rw1'][2][0]/10.+2)

        return

# =========================================================================
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_spatialInputUtil, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
