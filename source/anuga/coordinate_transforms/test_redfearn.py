
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

    def test_UTM_6_nonstandard_projection(self):
        #Test 6 (Geraldton, WA)

        #First test native projection (zone 50)
        zone, easting, northing = redfearn(-29.233299999,114.05)

        assert zone == 50
        assert allclose(easting, 213251.040253)
        assert allclose(northing, 6762559.15978)

        #Testing using the native zone
        zone, easting, northing = redfearn(-29.233299999,114.05, zone=50)

        assert zone == 50
        assert allclose(easting, 213251.040253)
        assert allclose(northing, 6762559.15978)

        #Then project to zone 49
        zone, easting, northing = redfearn(-29.233299999,114.05,zone=49)

        assert zone == 49
        assert allclose(easting, 796474.020057)
        assert allclose(northing, 6762310.25162)

       

        

        #First test native projection (zone 49)
        zone, easting, northing = redfearn(-29.1333,113.9667)

        assert zone == 49
        assert allclose(easting, 788653.192779)
        assert allclose(northing, 6773605.46384)

        #Then project to zone 50
        zone, easting, northing = redfearn(-29.1333,113.9667,zone=50)

        assert zone == 50
        assert allclose(easting, 204863.606467)
        assert allclose(northing, 6773440.04726)

        #Testing point on zone boundary
        #First test native projection (zone 50)
        zone, easting, northing = redfearn(-29.1667,114)

        assert zone == 50
        assert allclose(easting, 208199.768268)
        assert allclose(northing, 6769820.01453)

        #Then project to zone 49
        zone, easting, northing = redfearn(-29.1667,114,zone=49)

        assert zone == 49
        assert allclose(easting, 791800.231817)
        assert allclose(northing, 6769820.01453)

        #Testing furthest point in Geraldton scenario)
        #First test native projection (zone 49)
        zone, easting, northing = redfearn(-28.2167,113.4167)

        assert zone == 49
        assert allclose(easting, 737178.16131)
        assert allclose(northing, 6876426.38578)

        #Then project to zone 50
        zone, easting, northing = redfearn(-28.2167,113.4167,zone=50)

        assert zone == 50
        assert allclose(easting, 148260.567427)
        assert allclose(northing, 6873587.50926)

        #Testing outside GDA zone (New Zeland)
        #First test native projection (zone 60)
        zone, easting, northing = redfearn(-44,178)

        assert zone == 60
        assert allclose(easting, 580174.259843)
        assert allclose(northing, 5127641.114461)

        #Then project to zone 59
        zone, easting, northing = redfearn(-44,178,zone=59)

        assert zone == 59
        assert allclose(easting, 1061266.922118)
        assert allclose(northing, 5104249.395469)

        #Then skip three zones 57 (native 60)
        zone, easting, northing = redfearn(-44,178,zone=57)

        assert zone == 57
        assert allclose(easting, 2023865.527497)
        assert allclose(northing, 4949253.934967)

#Note Projecting into the Northern Hemisphere Does not coincide
#redfearn or ArcMap conversions
##        #Testing outside GDA zone (Northern Hemisphere)
##        #First test native projection (zone 57)
##        zone, easting, northing = redfearn(44,156)
##
##        assert zone == 57
##        assert allclose(easting, 259473.678944)
##        assert allclose(northing, 14876249.1268)
##
##        #Then project to zone 59
##        zone, easting, northing = redfearn(44,156,zone=56)
##
##        assert zone == 56
##        assert allclose(easting, 740526.321055)
##        assert allclose(northing, 14876249.1268)


        


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
        points, zone = convert_from_latlon_to_utm(latitudes=lats, longitudes=longs)

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
            points, zone = convert_from_latlon_to_utm(latitudes=lats, longitudes=longs)
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
            points, zone  = convert_from_latlon_to_utm(latitudes=lats, longitudes=longs)
        except ANUGAError:
            pass
        else:
            self.failUnless(False,
                            'Error not thrown error!')

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
        points, zone = convert_from_latlon_to_utm(points=points)
        #print "points",points 
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
            points, zone = convert_from_latlon_to_utm(points=points)           
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
            points, zone = convert_from_latlon_to_utm(points=points)           
        except ANUGAError:
            pass
        else:
            self.fail('Error not thrown error!')

            
#-------------------------------------------------------------
if __name__ == "__main__":

    #mysuite = unittest.makeSuite(TestCase,'test_convert_latlon_to_UTM1')
    #mysuite = unittest.makeSuite(TestCase,'test_UTM_6_nonstandard_projection')
    mysuite = unittest.makeSuite(TestCase,'test')
    runner = unittest.TextTestRunner()
    runner.run(mysuite)
