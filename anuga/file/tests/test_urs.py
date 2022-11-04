
import unittest
import numpy as num
import os

from anuga.file.urs import calculate_boundary_points, save_boundary_as_urs
from anuga.coordinate_transforms.redfearn import degminsec2decimal_degrees
from anuga.geospatial_data.geospatial_data import Geospatial_data

class Test_Urs(unittest.TestCase):
    def setUp(self):
        self.verbose = True
        pass

    def tearDown(self):
        for filename in ['swamp.tsh']:
            try:
                os.remove(filename)
            except:
                pass
        
    def test_URS_points_needed(self):
        
        ll_lat = -21.5
        ll_long = 114.5
        grid_spacing = 1./60.
        lat_amount = 30
        long_amount = 30
        zone = 50

        boundary_polygon = [[250000,7660000],[280000,7660000],
                             [280000,7630000],[250000,7630000]]
        geo=calculate_boundary_points(boundary_polygon, zone, 
                              ll_lat, ll_long, grid_spacing, 
                              lat_amount, long_amount,
                              verbose=self.verbose)
        # to test this geo, can info from it be transfered onto the boundary
        # poly?
        #Maybe see if I can fit the data to the polygon - have to represent
        # the poly as points though.
        #geo.export_points_file("results.txt", as_lat_long=True)
        results = frozenset(geo.get_data_points(as_lat_long=True))
        #print 'results',results

        # These are a set of points that have to be in results
        points = []
        for i in range(18):
            lat = -21.0 - 8./60 - grid_spacing * i
            points.append((lat,degminsec2decimal_degrees(114,35,0))) 
            points.append((lat,degminsec2decimal_degrees(114,36,0))) 
            points.append((lat,degminsec2decimal_degrees(114,52,0))) 
            points.append((lat,degminsec2decimal_degrees(114,53,0)))
        geo_answer = Geospatial_data(data_points=points,
                                     points_are_lats_longs=True) 
        #geo_answer.export_points_file("answer.txt", as_lat_long=True)  
        answer = frozenset(points)
        
        outs = answer.difference(results)
        #print "outs", outs
        # This doesn't work.  Though visualising the results shows that
        # it is correct.
        #assert answer.issubset(results)
        # this is why;
        #point (-21.133333333333333, 114.58333333333333)
        #result (-21.133333332232368, 114.58333333300342)
        
        for point in points:
            found = False
            for result in results:
                if num.allclose(point, result):
                    found = True
                    break
            if not found:
                assert False
        
    
    def in_development_URS_points_needed(self):
        ll_lat = -21.51667
        ll_long = 114.51667
        grid_spacing = 2./60.
        lat_amount = 15
        long_amount = 15

       
        boundary_polygon = [[250000,7660000],[280000,7660000],
                             [280000,7630000],[250000,7630000]]
        calculate_boundary_points('a_test_example',boundary_polygon,
                                  ll_lat, ll_long, grid_spacing, 
                                  lat_amount, long_amount,
                                  verbose=self.verbose)
        
        
    def test_URS_points_northern_hemisphere(self):

        LL_LAT = 8.0
        LL_LONG = 97.0
        GRID_SPACING = 2.0/60.0
        LAT_AMOUNT = 2
        LONG_AMOUNT = 2
        ZONE = 47

        # 
        points = []
        for i in range(2):
            for j in range(2):
                points.append((degminsec2decimal_degrees(8,1+i*2,0),
                               degminsec2decimal_degrees(97,1+i*2,0)))
        #print "points", points
        geo_poly = Geospatial_data(data_points=points,
                                     points_are_lats_longs=True)
        poly_lat_long = geo_poly.get_data_points(as_lat_long=False,
                                       isSouthHemisphere=False)
        #print "seg_lat_long",  poly_lat_long
        
      #   geo=URS_points_needed_to_file('test_example_poly3', poly_lat_long,
#                                   ZONE,
#                                   LL_LAT, LL_LONG,
#                                   GRID_SPACING,
#                                   LAT_AMOUNT, LONG_AMOUNT,
#                                   isSouthernHemisphere=False,
#                                   export_csv=True,
#                                   verbose=self.verbose)

        geo=calculate_boundary_points(poly_lat_long,
                                  ZONE,
                                  LL_LAT, LL_LONG,
                                  GRID_SPACING,
                                  LAT_AMOUNT, LONG_AMOUNT,
                                  isSouthHemisphere=False,
                                  verbose=self.verbose)

        results = frozenset(geo.get_data_points(as_lat_long=True,
                                                isSouthHemisphere=False))

        #print 'results',results

        # These are a set of points that have to be in results
        points = [] 
        for i in range(2):
            for j in range(2):
                points.append((degminsec2decimal_degrees(8,i*2,0),
                               degminsec2decimal_degrees(97,i*2,0)))
        #print "answer points", points
        answer = frozenset(points)
        
        for point in points:
            found = False
            for result in results:
                if num.allclose(point, result):
                    found = True
                    break
            if not found:
                assert False

    def in_development_URS_points_needed_poly1(self):
        # Values used for FESA 2007 results
        # domain in southern hemisphere zone 51        
        LL_LAT = -50.0
        LL_LONG = 80.0
        GRID_SPACING = 2.0/60.0
        LAT_AMOUNT = 4800
        LONG_AMOUNT = 3600
        ZONE = 51
        
        poly1 = [[296361.89, 8091928.62],
                 [429495.07,8028278.82],
                 [447230.56,8000674.05],
                 [429661.2,7982177.6],
                 [236945.9,7897453.16],
                 [183493.44,7942782.27],
                 [226583.04,8008058.96]]

        save_boundary_as_urs('test_example_poly2', poly1,
                                  ZONE,
                                  LL_LAT, LL_LONG,
                                  GRID_SPACING,
                                  LAT_AMOUNT, LONG_AMOUNT,
                                  verbose=self.verbose)
        


    def in_development_URS_points_needed_poly2(self):
        # Values used for 2004 validation work
        # domain in northern hemisphere zone 47        
        LL_LAT = 0.0
        LL_LONG = 90.0
        GRID_SPACING = 2.0/60.0
        LAT_AMOUNT = (15-LL_LAT)/GRID_SPACING
        LONG_AMOUNT = (100-LL_LONG)/GRID_SPACING
        ZONE = 47
        
        poly2 = [[419336.424,810100.845],
                 [342405.0281,711455.8026],
                 [274649.9152,723352.9603],
                 [272089.092,972393.0131],
                 [347633.3754,968551.7784],
                 [427979.2022,885965.2313],
                 [427659.0993,875721.9386],
                 [429259.6138,861317.3083],
                 [436301.8775,840830.723]]
        
        save_boundary_as_urs('test_example_poly2', poly2,
                                  ZONE,
                                  LL_LAT, LL_LONG,
                                  GRID_SPACING,
                                  LAT_AMOUNT, LONG_AMOUNT,
                                  isSouthernHemisphere=False,
                                  verbose=self.verbose) 
        


        
################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Urs,'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)
