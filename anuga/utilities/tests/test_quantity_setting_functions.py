#!/usr/bin/env python
"""

Unit tests for quantity_setting_functions

Coverage could be improved, but note that:

    test_elevation_from_Pt_Pol_Data_and_Raster

calls a routine which itself calls several other routines
in quantity_setting_functions, so we can catch bugs there too

"""

import unittest
import anuga
import numpy
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.utilities import plot_utils as util
from anuga.config import g
from anuga.utilities import quantity_setting_functions as qs
import os
       
# Choose minX,minY consistent a site located in UTM Zone EPSG:32756 
minX=307000.
minY=6193000.

class Test_quantity_setting_functions(unittest.TestCase):
    """
	Test the quantity_setting_functions
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def create_domain(self, InitialOceanStage, InitialLandStage, flowAlg='DE0', verbose=False):
        """
         Make the domain and set the flow algorithm for a test. Produces an sww
         that we can use for testing
        """

        boundaryPolygon=[ [minX, minY], [minX, minY+100.], 
            [minX+100., minY+100.], [minX+100., minY]]
        anuga.create_mesh_from_regions(
            boundaryPolygon, 
            boundary_tags={'left': [0],
                         'top': [1],
                         'right': [2],
                         'bottom': [3]},
            maximum_triangle_area = 1.,
            minimum_triangle_angle = 28.0,
            filename = 'test_quantity_setting_functions.msh',
            interior_regions =[ ],
            verbose=False)

        domain=anuga.create_domain_from_file('test_quantity_setting_functions.msh')
        
        os.remove('test_quantity_setting_functions.msh')

        domain.set_flow_algorithm(flowAlg)
        domain.set_name('test_quantity_setting_functions')

        domain.set_store_vertices_uniquely()
       
        def topography(x,y):
            return -x/150. 

        def stagefun(x,y):
            stg=InitialOceanStage*(x>=50.) + InitialLandStage*(x<50.)
            return stg 

        domain.set_flow_algorithm(flowAlg)
        #domain.set_quantity('elevation',topography,location='centroids')     
        domain.set_quantity('elevation',topography)     
        domain.set_quantity('friction',0.03)             
        #domain.set_quantity('stage', stagefun,location='centroids')            
        domain.set_quantity('stage', stagefun)            
       
        if(verbose):
            if(domain.store_centroids):
                print '   Centroids stored'
            else:
                print '    Centroids estimated from vertices'

        # Boundary conditions
        Br=anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom':Br})

        return domain
   
    def test_make_nearestNeighbour_quantity_function(self):

        domain=self.create_domain(1.0, 0.0)

        # Make a set of points in 'physical space' ranging from x,y,x+y-min(x)-min(y)
        xR=numpy.linspace(minX, minX+100., 101)
        yR=numpy.linspace(minY, minY+100., 101)
        xR,yR=numpy.meshgrid(xR,yR)
        xR=xR.flatten()
        yR=yR.flatten()
        zR=xR+yR-minX-minY
        inPts=numpy.vstack([xR,yR,zR]).transpose()

        F = qs.make_nearestNeighbour_quantity_function(inPts, domain, 
                threshold_distance = 9.0e+100, background_value = 9.0e+100)

        # Test that F evaluated in 'ANUGA coordinates' [lower-left = 0,0] is correct
        xGet=numpy.array([0., 10., 10., 0., 100.]) + 0.1
        yGet=numpy.array([0., 0., 20.,90., 100. ]) + 0.1
        expected=numpy.floor(xGet+yGet)
        output=F(xGet,yGet)
        assert(numpy.allclose(output,expected))


        # Test with 3 nearest neighbours at points which are exactly on data points
        F = qs.make_nearestNeighbour_quantity_function(inPts, domain, 
                threshold_distance = 9.0e+100, background_value = 9.0e+100,
                k_nearest_neighbours = 3)

        xGet=numpy.array([0., 10., 10., 0., 100.])
        yGet=numpy.array([0., 0., 20.,90., 100. ])
        expected=xGet+yGet
        output=F(xGet,yGet)
        assert(numpy.allclose(output,expected))

        # Test with 4 nearest neighbours, at points which are not exactly on data points  
        F = qs.make_nearestNeighbour_quantity_function(inPts, domain, 
                threshold_distance = 9.0e+100, background_value = 9.0e+100,
                k_nearest_neighbours = 4)
        # Test at non-trivial location
        xGet=numpy.array([0., 10., 10., 0., 99.]) + 0.5
        yGet=numpy.array([0., 0., 20.,90., 99. ]) + 0.5
        expected=xGet+yGet
        output=F(xGet,yGet)
        assert(numpy.allclose(output,expected))

        
        #################################################################################
        # Test that background_value / threshold_distance work ok

        # Make a set of points in 'physical space' ranging from x,y,x+y-min(x)-min(y)
        xR=numpy.linspace(minX+10., minX+90., 81)
        yR=numpy.linspace(minY+10., minY+90., 81)
        xR,yR=numpy.meshgrid(xR,yR)
        xR=xR.flatten()
        yR=yR.flatten()
        zR=xR+yR-minX-minY
        inPts=numpy.vstack([xR,yR,zR]).transpose()

        F = qs.make_nearestNeighbour_quantity_function(inPts, domain, 
                threshold_distance = 6., background_value = 9.0e+100)

        # Test that F evaluated in 'ANUGA coordinates' [lower-left = 0,0] is correct
        #
        # Now points 1 and 2 and 5 are outside the threshold_distance
        xGet=numpy.array( [0., 10.,   10.,10., 100. ] )
        yGet=numpy.array( [0., 3.999, 20.,90., 100. ] )
        expected=xGet+yGet
        expected[0] = 9.0e+100
        expected[1] = 9.0e+100
        expected[4] = 9.0e+100
        output=F(xGet,yGet)
        assert(numpy.allclose(output,expected))

        return

    def test_composite_quantity_setting_function(self):
        # Test the composite_quantity_setting_function
        
        domain=self.create_domain(1.0, 0.0)
        
        # Make a raster from the elevation data
        from anuga.utilities import plot_utils as util
        xs=domain.centroid_coordinates[:,0]+domain.geo_reference.xllcorner
        ys=domain.centroid_coordinates[:,1]+domain.geo_reference.yllcorner
        elev=domain.quantities['elevation'].centroid_values

        allDat=numpy.vstack([xs,ys,elev]).transpose()
        util.Make_Geotif(allDat, output_quantities=['ElevTest'], EPSG_CODE=32756, 
                        output_dir='.', CellSize=1.,k_nearest_neighbours=1)

        # Make a polygon-point pair which we use to set elevation in a 'channel'
        trenchPoly = [[minX+40., minY], [minX+40., minY+100.], 
            [minX+60., minY+100.], [minX+60., minY]]

        #################################################################
 
        # This example uses a constant, and a raster, to set the quantity           
        F=qs.composite_quantity_setting_function(
            [[trenchPoly, -1000.], ['Extent', 'PointData_ElevTest.tif']],
            domain)

        # Points where we test the function
        testPts_X=numpy.array([50., 3.])
        testPts_Y=numpy.array([1., 20.])
        fitted=F(testPts_X,testPts_Y)

        # The fitted value in the trench should be -1000.
        assert(fitted[0]==-1000.)

        # Find the nearest domain point to the second test point
        # This will have been used in constructing the elevation raster
        nearest=((domain.centroid_coordinates[:,0]-3.)**2 + 
                 (domain.centroid_coordinates[:,1]-20.)**2).argmin()
        nearest_x=domain.centroid_coordinates[nearest,0]
        assert(numpy.allclose(fitted[1],-nearest_x/150.))

        #################################################################
 
        # This example uses a constant, and a raster, to set the quantity, and
        # applies the min/max bound           
        F=qs.composite_quantity_setting_function(
            [[trenchPoly, -1000.], ['Extent', 'PointData_ElevTest.tif']],
            domain,
            clip_range = [[-500., 1.0e+100], [-1.0e+100, 1.0e+100]]) 

        # Points where we test the function
        testPts_X=numpy.array([50., 3.])
        testPts_Y=numpy.array([1., 20.])
        fitted=F(testPts_X,testPts_Y)

        # The fitted value in the trench should be -500, because of clipping
        assert(fitted[0]==-500.)

        # Find the nearest domain point to the second test point
        # This will have been used in constructing the elevation raster
        nearest = ((domain.centroid_coordinates[:,0]-3.)**2 + 
                   (domain.centroid_coordinates[:,1]-20.)**2).argmin()
        nearest_x = domain.centroid_coordinates[nearest,0]

        assert(numpy.allclose(fitted[1],-nearest_x/150.))

        #########################################################################

        # This example uses a function, and a raster, to set the quantity           
        def f0(x,y):
            return x/10.
        F = qs.composite_quantity_setting_function(
            [[trenchPoly, f0], ['Extent', 'PointData_ElevTest.tif']],
            domain) 
        fitted = F(testPts_X,testPts_Y)
        # Now the fitted value in the trench should be determined by f0
        assert(numpy.allclose(fitted[0],50./10.))
        # The second test point should be as before
        nearest = ((domain.centroid_coordinates[:,0]-3.)**2 + 
                   (domain.centroid_coordinates[:,1]-20.)**2).argmin()
        nearest_x = domain.centroid_coordinates[nearest,0]
        assert(numpy.allclose(fitted[1],-nearest_x/150.))

        ##########################################################################

        # This example uses 'All' as a polygon
        F = qs.composite_quantity_setting_function(
            [['All', f0 ], [None, 'PointData_ElevTest.tif']],
            domain) 
        fitted=F(testPts_X,testPts_Y)
        # Now the fitted value in the trench should be determined by f0
        assert(numpy.allclose(fitted[0],50./10.))
        assert(numpy.allclose(fitted[1],3./10.))

        ###########################################################################
        # This example should fail
        def should_fail():
            F=qs.composite_quantity_setting_function(
                [['All', f0], ['All', 'PointData_ElevTest.tif']],
                domain)
            # Need to call it to get the error
            F(numpy.array([3.]), numpy.array([3.]))
            return
        self.assertRaises(Exception, lambda: should_fail())

        ###########################################################################
        # This example should fail (since the clip_range minimum >= maximum)
        def should_fail_1():
            F=qs.composite_quantity_setting_function(
                [[trenchPoly, f0], ['All', 'PointData_ElevTest.tif']],
                domain,
                clip_range = [[ -500., -1000.], [-1.0e+100, 1.0e+100]]) 
            return
        self.assertRaises(Exception, lambda: should_fail_1())

        ###########################################################################
        # This example features a function with some nan return values, and uses
        # the nan_interpolation_region_polygon to try to fix it
        def f0(x,y):
            output = x/10.
            output[0] = numpy.nan
            return output 

        F = qs.composite_quantity_setting_function(
            [[trenchPoly, f0], ['Extent', 'PointData_ElevTest.tif']],
            domain,
            nan_treatment = 'fall_through',
            nan_interpolation_region_polygon = [trenchPoly],
            default_k_nearest_neighbours = 3,
            default_raster_interpolation = 'bilinear',
            verbose=False) 

        # Points where we test the function. We deliberately use many points with x=50,
        # which happens to ensure that the nan value is replaced with the same
        # value it would have had anyway
        testPts_X = numpy.array([50.,50.00, 50., 50., 97., 51., 3.])
        testPts_Y = numpy.array([1.,    2., 3. , 4  , 20., 50., 60.])
        fitted = F(testPts_X,testPts_Y)
       
        # We should have no nan values
        assert(sum(fitted!=fitted) == 0)

        # Now the fitted value in the trench should be determined by f0 because
        # the re-interpolation of nan values was designed to ensure it
        assert(numpy.allclose(fitted[0],50./10.))

        ###########################################################################
        # This example features a function with some nan return values, and uses
        # the nan_interpolation_region_polygon to try to fix it

        # Make a polygon-point pair which we use to set elevation in a 'channel'
        innerTrenchPoly = [[minX+45., minY+45.], [minX+45., minY+55.], 
            [minX+55., minY+55.], [minX+55., minY+45.]]

        def f_nan(x,y):
            output = x*0 + numpy.nan
            return(output)

        F = qs.composite_quantity_setting_function(
            [[innerTrenchPoly, f_nan], [trenchPoly, f0], ['Extent', 'PointData_ElevTest.tif']],
            domain,
            nan_treatment = 'fall_through',
            nan_interpolation_region_polygon = [trenchPoly],
            default_k_nearest_neighbours = 3,
            default_raster_interpolation = 'bilinear',
            verbose=False) 

        # Points where we test the function. We deliberately use many points with x=50,
        # which happens to ensure that the nan value is replaced with the same
        # value it would have had anyway
        testPts_X = numpy.array([50.,50.00, 50., 50., 97., 51., 3.])
        testPts_Y = numpy.array([1.,    2., 3. , 4  , 20., 50., 60.])
        fitted = F(testPts_X,testPts_Y)
       
        # We should have no nan values
        assert(sum(fitted!=fitted) == 0)

        # Now the fitted value in the trench should be determined by f0 because
        # the re-interpolation of nan values was designed to ensure it
        assert(numpy.allclose(fitted[0],50./10.))

        return

    def test_quantity_from_Pt_Pol_Data_and_Raster(self):
        # 
        # 
        #
        domain = self.create_domain(1.0, 0.0)

        # Evolve the model
        #for t in domain.evolve(yieldstep=0.2, finaltime=1.0):
        #    pass

        # Make a raster from the elevation data
        from anuga.utilities import plot_utils as util
        xs = domain.centroid_coordinates[:,0]+domain.geo_reference.xllcorner
        ys = domain.centroid_coordinates[:,1]+domain.geo_reference.yllcorner
        elev = domain.quantities['elevation'].centroid_values

        allDat = numpy.vstack([xs,ys,elev]).transpose()
        util.Make_Geotif(allDat, output_quantities=['ElevTest'], 
            EPSG_CODE=32756, output_dir='.', CellSize=1.,k_nearest_neighbours=1)

        # Make a polygon-point pair which we use to set elevation in a 'channel'
        trenchPoly = [[minX+40., minY], [minX+40., minY+100.], 
            [minX+60., minY+100.], [minX+60., minY]]
        trenchPts = numpy.array([minX+50., minY+50., -1000.])
        #
        PtPolData = [[trenchPoly, trenchPts]]
        F = qs.quantity_from_Pt_Pol_Data_and_Raster(PtPolData, 
            'PointData_ElevTest.tif', domain) 

        testPts_X = numpy.array([50., 3.])
        testPts_Y = numpy.array([1., 20.])
        fitted = F(testPts_X,testPts_Y)

        # The fitted value in the trench should be -1000.
        assert(numpy.allclose(fitted[0],-1000.))

        # Find the nearest domain point to the second test point
        # This will have been used in constructing the elevation raster
        nearest = ((domain.centroid_coordinates[:,0]-3.)**2 +\
                   (domain.centroid_coordinates[:,1]-20.)**2).argmin()
        nearest_x = domain.centroid_coordinates[nearest,0]
        assert(numpy.allclose(fitted[1],-nearest_x/150.))

        # Clean up file
        os.remove('PointData_ElevTest.tif')

        return
# =========================================================================
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_quantity_setting_functions, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
