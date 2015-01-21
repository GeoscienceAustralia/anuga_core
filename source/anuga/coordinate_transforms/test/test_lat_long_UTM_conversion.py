
#Test of redfearns formula. Tests can be verified at
#
#http://www.cellspark.com/UTM.html
#http://www.ga.gov.au/nmd/geodesy/datums/redfearn_geo_to_grid.jsp


import unittest

from anuga.coordinate_transforms.lat_long_UTM_conversion import *
from anuga.coordinate_transforms.redfearn import degminsec2decimal_degrees, decimal_degrees2degminsec
from anuga.anuga_exceptions import ANUGAError

import numpy as num


#-------------------------------------------------------------

class TestCase(unittest.TestCase):

    def test_UTM_1(self):
        #latitude:  -37 39' 10.15610" 
        #Longitude: 143 55' 35.38390" 
        #Site Name:    GDA-MGA: (UTM with GRS80 ellipsoid) 
        #Zone:   54    
        #Easting:  758173.797  Northing: 5828674.340 
        #Latitude:   -37  39 ' 10.15610 ''  Longitude: 143  55 ' 35.38390 '' 
        #Grid Convergence:  1  47 ' 19.36 ''  Point Scale: 1.00042107 

        lat = degminsec2decimal_degrees(-37,39,10.15610)
        lon = degminsec2decimal_degrees(143,55,35.38390) 
        assert num.allclose(lat, -37.65282114)
        assert num.allclose(lon, 143.9264955)


        zone, easting, northing = LLtoUTM(lat,lon)

        assert zone == 54
        assert num.allclose(easting, 758173.797)
        assert num.allclose(northing, 5828674.340)

        lat_calced, long_calced = UTMtoLL(northing, easting, zone) 
        assert num.allclose(lat,  lat_calced)
        assert num.allclose(lon, long_calced)


    def test_UTM_2(self):
        #TEST 2

        #Latitude:  -37 57 03.7203
        #Longitude: 144 25 29.5244
        #Zone:   55    
        #Easting:  273741.297  Northing: 5796489.777 
        #Latitude:   -37  57 ' 3.72030 ''  Longitude: 144  25 ' 29.52440 '' 
        #Grid Convergence:  -1  35 ' 3.65 ''  Point Scale: 1.00023056 

        lat = degminsec2decimal_degrees(-37,57,03.7203)
        lon = degminsec2decimal_degrees(144,25,29.5244) 

        zone, easting, northing = LLtoUTM(lat,lon)

        
        assert zone == 55
        assert num.allclose(easting, 273741.297)
        assert num.allclose(northing, 5796489.777)

        lat_calced, long_calced = UTMtoLL(northing, easting, zone) 
        assert num.allclose(lat,  lat_calced)
        assert num.allclose(lon, long_calced)
        
        
    def test_UTM_3(self):
        #Test 3
        lat = degminsec2decimal_degrees(-60,0,0)
        lon = degminsec2decimal_degrees(130,0,0) 

        zone, easting, northing = LLtoUTM(lat,lon)

        assert zone == 52
        assert num.allclose(easting, 555776.267)
        assert num.allclose(northing, 3348167.264)

        Lat, Long = UTMtoLL(northing, easting, zone)

    def test_UTM_4(self):
        #Test 4 (Kobenhavn, Northern hemisphere)
        lat = 55.70248
        dd,mm,ss = decimal_degrees2degminsec(lat)

        lon = 12.58364
        dd,mm,ss = decimal_degrees2degminsec(lon)

        zone, easting, northing = LLtoUTM(lat,lon)
        
        assert zone == 33
        assert num.allclose(easting, 348157.631)
        assert num.allclose(northing, 6175612.993) 

        lat_calced, long_calced = UTMtoLL(northing, easting, zone,
                                          isSouthernHemisphere=False) 
        assert num.allclose(lat,  lat_calced)
        assert num.allclose(lon, long_calced)

    def test_UTM_5(self):
        #Test 5 (Wollongong)

        lat = degminsec2decimal_degrees(-34,30,0.)
        lon = degminsec2decimal_degrees(150,55,0.) 
        
        zone, easting, northing = LLtoUTM(lat,lon)


        assert zone == 56
        assert num.allclose(easting, 308728.009)
        assert num.allclose(northing, 6180432.601)

        lat_calced, long_calced = UTMtoLL(northing, easting, zone) 
        assert num.allclose(lat,  lat_calced)
        assert num.allclose(lon, long_calced)

#-------------------------------------------------------------

if __name__ == "__main__":
    mysuite = unittest.makeSuite(TestCase,'test')
    runner = unittest.TextTestRunner()
    runner.run(mysuite)
