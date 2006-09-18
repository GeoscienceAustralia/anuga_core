
#Test of redfearns formula. Tests can be verified at
#
#http://www.cellspark.com/UTM.html
#http://www.ga.gov.au/nmd/geodesy/datums/redfearn_geo_to_grid.jsp


import unittest
from Numeric import allclose

from redfearn import *
from anuga.utilities.anuga_exceptions import ANUGAError

#-------------------------------------------------------------

class TestCase(unittest.TestCase):

    def test_decimal_degrees_conversion(Self):
        lat = degminsec2decimal_degrees(-37,39,10.15610)
        lon = degminsec2decimal_degrees(143,55,35.38390) 
        assert allclose(lat, -37.65282114)
        assert allclose(lon, 143.9264955)

        dd,mm,ss = decimal_degrees2degminsec(-37.65282114)
        assert dd==-37
        assert mm==39
        assert allclose(ss, 10.15610)

        dd,mm,ss = decimal_degrees2degminsec(143.9264955)
        assert dd==143
        assert mm==55
        assert allclose(ss, 35.38390) 


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
        assert allclose(lat, -37.65282114)
        assert allclose(lon, 143.9264955)


        zone, easting, northing = redfearn(lat,lon)

        assert zone == 54
        assert allclose(easting, 758173.797)
        assert allclose(northing, 5828674.340)


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
        #print lat, lon

        zone, easting, northing = redfearn(lat,lon)

        assert zone == 55
        assert allclose(easting, 273741.297)
        assert allclose(northing, 5796489.777)

        
    def test_UTM_3(self):
        #Test 3
        lat = degminsec2decimal_degrees(-60,0,0)
        lon = degminsec2decimal_degrees(130,0,0) 

        zone, easting, northing = redfearn(lat,lon)
        #print zone, easting, northing

        assert zone == 52
        assert allclose(easting, 555776.267)
        assert allclose(northing, 3348167.264)


    def test_UTM_4(self):
        #Test 4 (Kobenhavn, Northern hemisphere)
        lat = 55.70248
        dd,mm,ss = decimal_degrees2degminsec(lat)

        lon = 12.58364
        dd,mm,ss = decimal_degrees2degminsec(lon)

        zone, easting, northing = redfearn(lat,lon)


        assert zone == 33
        assert allclose(easting, 348157.631)
        assert allclose(northing, 6175612.993) 


    def test_UTM_5(self):
        #Test 5 (Wollongong)

        lat = degminsec2decimal_degrees(-34,30,0.)
        lon = degminsec2decimal_degrees(150,55,0.) 
        
        zone, easting, northing = redfearn(lat,lon)

        #print zone, easting, northing

        assert zone == 56
        assert allclose(easting, 308728.009)
        assert allclose(northing, 6180432.601)


    #def test_UTM_6(self):
    #    """Test 6 (Don's Wollongong file's ref point)
    #    """
    #
    #    lat = -34.490286785873 
    #    lon = 150.79712139578
    #
    #    dd,mm,ss = decimal_degrees2degminsec(lat)
    #    print dd,mm,ss
    #    dd,mm,ss = decimal_degrees2degminsec(lon)        
    #    print dd,mm,ss
    #     
    #    zone, easting, northing = redfearn(lat,lon)
    #
    #    print zone, easting, northing
    #
    #    assert zone == 56
    #    #assert allclose(easting, 297717.36468927) #out by 10m
    #    #assert allclose(northing, 6181725.1724276)

    def test_convert_lats_longs(self):

        #Site Name:    GDA-MGA: (UTM with GRS80 ellipsoid) 
        #Zone:   56    
        #Easting:  222908.705  Northing: 6233785.284 
        #Latitude:   -34  0 ' 0.00000 ''  Longitude: 150  0 ' 0.00000 '' 
        #Grid Convergence:  -1  40 ' 43.13 ''  Point Scale: 1.00054660 

        lat_gong = degminsec2decimal_degrees(-34,30,0.)
        lon_gong = degminsec2decimal_degrees(150,55,0.)
        
        lat_2 = degminsec2decimal_degrees(-34,00,0.)
        lon_2 = degminsec2decimal_degrees(150,00,0.)
        
        lats = [lat_gong, lat_2]
        longs = [lon_gong, lon_2]
        zone, points = convert_lats_longs(lats, longs)

        assert allclose(points[0][0], 308728.009)
        assert allclose(points[0][1], 6180432.601)
        assert allclose(points[1][0],  222908.705)
        assert allclose(points[1][1], 6233785.284)
        self.failUnless(zone == 56,
                        'Bad zone error!')
        
    def test_convert_lats_longs2(self):

        #Site Name:    GDA-MGA: (UTM with GRS80 ellipsoid) 
        #Zone:   56    
        #Easting:  222908.705  Northing: 6233785.284 
        #Latitude:   -34  0 ' 0.00000 ''  Longitude: 150  0 ' 0.00000 '' 
        #Grid Convergence:  -1  40 ' 43.13 ''  Point Scale: 1.00054660 

        lat_gong = degminsec2decimal_degrees(-34,30,0.)
        lon_gong = degminsec2decimal_degrees(150,55,0.)
        
        lat_2 = degminsec2decimal_degrees(34,00,0.)
        lon_2 = degminsec2decimal_degrees(100,00,0.)
        
        lats = [lat_gong, lat_2]
        longs = [lon_gong, lon_2]
        
        try:
            zone, points = convert_lats_longs(lats, longs)
        except ANUGAError:
            pass
        else:
            self.failUnless(False,
                            'Error not thrown error!')
            
    def test_convert_lats_longs3(self):

        #Site Name:    GDA-MGA: (UTM with GRS80 ellipsoid) 
        #Zone:   56    
        #Easting:  222908.705  Northing: 6233785.284 
        #Latitude:   -34  0 ' 0.00000 ''  Longitude: 150  0 ' 0.00000 '' 
        #Grid Convergence:  -1  40 ' 43.13 ''  Point Scale: 1.00054660 

        lat_gong = "-34.5"
        lon_gong = "150.916666667"
        lat_2 = degminsec2decimal_degrees(34,00,0.)
        lon_2 = degminsec2decimal_degrees(100,00,0.)
        
        lats = [lat_gong, lat_2]
        longs = [lon_gong, lon_2]
        try:
            zone, points = convert_lats_longs(lats, longs)
        except ANUGAError:
            pass
        else:
            self.failUnless(False,
                            'Error not thrown error!')

    # Similar test for alternative interface        
    def test_convert_latlon_to_UTM1(self):

        #Site Name:    GDA-MGA: (UTM with GRS80 ellipsoid) 
        #Zone:   56    
        #Easting:  222908.705  Northing: 6233785.284 
        #Latitude:   -34  0 ' 0.00000 ''  Longitude: 150  0 ' 0.00000 '' 
        #Grid Convergence:  -1  40 ' 43.13 ''  Point Scale: 1.00054660 

        lat_gong = degminsec2decimal_degrees(-34,30,0.)
        lon_gong = degminsec2decimal_degrees(150,55,0.)
        
        lat_2 = degminsec2decimal_degrees(-34,00,0.)
        lon_2 = degminsec2decimal_degrees(150,00,0.)
        
        points = [[lat_gong, lon_gong], [lat_2, lon_2]]
        points, zone = convert_points_from_latlon_to_utm(points)

        assert allclose(points[0][0], 308728.009)
        assert allclose(points[0][1], 6180432.601)
        assert allclose(points[1][0],  222908.705)
        assert allclose(points[1][1], 6233785.284)
        self.failUnless(zone == 56,
                        'Bad zone error!')

    def test_convert_latlon_to_UTM2(self):       

        #Site Name:    GDA-MGA: (UTM with GRS80 ellipsoid) 
        #Zone:   56    
        #Easting:  222908.705  Northing: 6233785.284 
        #Latitude:   -34  0 ' 0.00000 ''  Longitude: 150  0 ' 0.00000 '' 
        #Grid Convergence:  -1  40 ' 43.13 ''  Point Scale: 1.00054660

        lat_gong = degminsec2decimal_degrees(-34,30,0.)
        lon_gong = degminsec2decimal_degrees(150,55,0.)
        
        lat_2 = degminsec2decimal_degrees(34,00,0.)
        lon_2 = degminsec2decimal_degrees(100,00,0.)

        points = [[lat_gong, lon_gong], [lat_2, lon_2]]

        try:
            points, zone = convert_points_from_latlon_to_utm(points)           
        except ANUGAError:
            pass
        else:
            self.fail('Error not thrown error!')

    def test_convert_latlon_to_UTM3(self):            

        #Site Name:    GDA-MGA: (UTM with GRS80 ellipsoid) 
        #Zone:   56    
        #Easting:  222908.705  Northing: 6233785.284 
        #Latitude:   -34  0 ' 0.00000 ''  Longitude: 150  0 ' 0.00000 '' 
        #Grid Convergence:  -1  40 ' 43.13 ''  Point Scale: 1.00054660 

        lat_gong = "-34.5"
        lon_gong = "150.916666667"
        lat_2 = degminsec2decimal_degrees(34,00,0.)
        lon_2 = degminsec2decimal_degrees(100,00,0.)

        points = [[lat_gong, lon_gong], [lat_2, lon_2]]

        try:
            points, zone = convert_points_from_latlon_to_utm(points)           
        except ANUGAError:
            pass
        else:
            self.fail('Error not thrown error!')

            
#-------------------------------------------------------------
if __name__ == "__main__":

    mysuite = unittest.makeSuite(TestCase,'test')
    runner = unittest.TextTestRunner()
    runner.run(mysuite)
