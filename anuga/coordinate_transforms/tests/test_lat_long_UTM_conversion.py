
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
    '''Test UTM and Lat Long conversions. UTMtoLL and LLtoUTM can now take a single or an array of points as input'''
    def test_UTM_1(self):
        '''test the conversion from Lat Long coordinates to utm coordinates.
	test the conversion of single point and array of points'''
	# test a single point input in the form of LLtoUTM(lat, lon)
	result0 = LLtoUTM(-25.689, 150.571)
	zone0, easting0, northing0 = result0
	utmpoint0 = num.array([easting0, northing0])
	assert zone0 == 56
	assert num.allclose(num.array([  256228.95844341,  7156515.89696371]), utmpoint0)
	
	# test a single point input in the form of a list [lat, lon]
	result1 = LLtoUTM([-25.689, 150.571])
	zone1, easting1, northing1 = result1
	utmpoint1 = num.array([easting1, northing1])
	assert zone1 == 56
	assert num.allclose(num.array([  256228.95844341,  7156515.89696371]), utmpoint1)

	# test a single point input in the form of a tuple (lat, lon)
	result2 = LLtoUTM((-25.689, 150.571))
	zone2, easting2, northing2 = result2
	utmpoint2 = num.array([easting2, northing2])
	assert zone2 == 56
	assert num.allclose(num.array([  256228.95844341,  7156515.89696371]), utmpoint2)

	# test a single point input in the form of a numpy array array([lat, lon])
	result3 = LLtoUTM(num.array([-25.689, 150.571]))
	zone3, easting3, northing3 = result3
	utmpoint3 = num.array([easting3,northing3])
	assert zone3 == 56
	assert num.allclose(num.array([  256228.95844341,  7156515.89696371]), utmpoint3)

	# test input as a list of points in the form of [[lat1, lon1], ..., [latn, lonn]]
	# the result is an array of zone numbers and an array of utm points
	result4 = LLtoUTM([[-25.689, 150.571], [-30.211, 110.532], [-120.901, 50.256]], ReferenceEllipsoid=23)
	zone4, utmpoints4 = result4
	assert zone4[0] == 56.
	assert zone4[1] == 49.
	assert zone4[2] == 39.
	assert num.allclose(num.array([[  256228.95844341,  7156515.89696371],
	       			       [  454958.02747461,  6657741.19780555],
				       [  542621.18844459, -3444616.62785417]]), utmpoints4)

	# test input as numpy array of points in the form of array([[lat1, lon1], ..., [latn, lonn]])
	# the result is an array of zone numbers and an array of utm points
	result5 = LLtoUTM(num.array([[-25.689, 150.571], [-30.211, 110.532], [-120.901, 50.256]]))
	zone5, utmpoints5 = result5
	assert zone4[0] == 56.
	assert zone4[1] == 49.
	assert zone4[2] == 39.
	assert num.allclose(num.array([[  256228.95844341,  7156515.89696371],
	       			       [  454958.02747461,  6657741.19780555],
				       [  542621.18844459, -3444616.62785417]]), utmpoints5)

	# test input as a tuple of points in the form of [(lat1, lon1), ..., (latn, lonn)]
	# the result is an array of zone numbers and an array of utm points
	result6 = LLtoUTM([(-25.689, 150.571), (-30.211, 110.532), (-120.901, 50.256)])	
	zone6, utmpoints6 = result6
	assert zone5[0] == 56.
	assert zone5[1] == 49.
	assert zone5[2] == 39.
	assert num.allclose(num.array([[  256228.95844341,  7156515.89696371],
	       			       [  454958.02747461,  6657741.19780555],
				       [  542621.18844459, -3444616.62785417]]), utmpoints6)

    def test_UTM_2(self):
	''' test the conversion from utm to Lat Long coordinates.
        test the conversion of single point and array of points'''
	# test utm to Lat Long conversion for a single point in the form of: northing, easting, zone
	northing = 7620000
	easting = 250000
	zone = 56
	latlong0 = UTMtoLL(northing, easting, zone)
	assert num.allclose(num.array([ -21.5052727,   150.58681618]), num.array(latlong0))

	# test utm to Lat Long conversion for a single point in the form of: (northing, easting, zone)
	latlong1 = UTMtoLL((northing, easting, zone))
	assert num.allclose(num.array([ -21.5052727,   150.58681618]), num.array(latlong1))

	# test utm to Lat Long conversion for a single point in the form of: [northing, easting, zone]
	latlong2 = UTMtoLL([northing, easting, zone])
	assert num.allclose(num.array([ -21.5052727,   150.58681618]), num.array(latlong2))

	# test utm to Lat Long conversion for a single point in the form of: array[northing, easting, zone]
	latlong3 = UTMtoLL(num.array([northing, easting, zone]))
	assert num.allclose(num.array([ -21.5052727,   150.58681618]), num.array(latlong3))

	# test utm to Lat Long conversion for a list of points in the form of:
	# [[northing1, easting1, zone1], ... ]
	utm_points4 = [[northing, easting, zone], [northing-1200, easting+6000, zone], [northing+12000, easting+100, zone]]
	latlong4 = UTMtoLL(utm_points4)
	assert num.allclose(num.array([[ -21.5052727,   150.58681618],
				       [ -21.51693237,  150.64452035],
			      	       [ -21.39696114,  150.5895628 ]]), latlong4)

	# test utm to Lat Long conversion for a list of points in the form of:
	# [(northing1, easting1, zone1), ... ]
	utm_points5 = [(northing, easting, zone), (northing-1200, easting+6000, zone), (northing+12000, easting+100, zone)]
	latlong5 = UTMtoLL(utm_points5, ReferenceEllipsoid=23)
	assert num.allclose(num.array([[ -21.5052727,   150.58681618],
				       [ -21.51693237,  150.64452035],
			      	       [ -21.39696114,  150.5895628 ]]), latlong5)

	# test utm to Lat Long conversion for a list of points in the form of:
	# array([[northing1, easting1, zone1], ... ])
	utm_points6 = num.array([[northing, easting, zone], [northing-1200, easting+6000, zone], [northing+12000, easting+100, zone]])
	latlong6 = UTMtoLL(utm_points6)
	assert num.allclose(num.array([[ -21.5052727,   150.58681618],
				       [ -21.51693237,  150.64452035],
			      	       [ -21.39696114,  150.5895628 ]]), latlong6)

	# test utm to Lat Long conversion for a list of points in the form of:
	# [[northing1, easting1], ... ], [zone1, ... ]
	utm_points7 = [[  7156515.89696371, 256228.95844341], [  6657741.19780555, 454958.02747461], [  -3444616.62785417, 542621.18844459]]
	zone7 = [56, 49, 39]
	latlong7 = UTMtoLL(utm_points7, zone7, isSouthernHemisphere=True, ReferenceEllipsoid=23)
	assert num.allclose(num.array([[-25.689, 150.571], 
				       [-30.211, 110.532], 
				       [-120.901, 50.256]]), latlong7)

	# test utm to Lat Long conversion for a list of points in the form of:
	# [(northing1, easting1), ... ], (zone1, ... )
	utm_points8 = [(  7156515.89696371, 256228.95844341), (  6657741.19780555, 454958.02747461), (  -3444616.62785417, 542621.18844459)]
	zone8 = (56, 49, 39)
	latlong8 = UTMtoLL(utm_points8, zone8)
	assert num.allclose(num.array([[-25.689, 150.571], 
				       [-30.211, 110.532], 
				       [-120.901, 50.256]]), latlong8)

	# test utm to Lat Long conversion for a list of points in the form of:
	# array([[northing1, easting1], ... ]), array([zone1, ... ])
	utm_points9 = num.array([[  7156515.89696371, 256228.95844341],
	       			 [  6657741.19780555, 454958.02747461],
				 [ -3444616.62785417, 542621.18844459]])
	zone9 = num.array([56, 49, 39])				
	latlong9 = UTMtoLL(utm_points9, zone9)
	assert num.allclose(num.array([[-25.689, 150.571], 
				       [-30.211, 110.532], 
				       [-120.901, 50.256]]), latlong9)

	# test utm to Lat Long conversion for a list of points in the form of:
	# array([[northing1, easting1], ... ]), zone
	# all utm points have the same zone number, which is passed as an integer
	utm_points10 = num.array(zip(northing+num.arange(0, 5*6000, 6000), easting+num.arange(0, 5*6000, 6000)))
	latlong10 = UTMtoLL(utm_points10, 56)
	assert num.allclose(num.array([[ -21.5052727,   150.58681618],
	 			       [ -21.45193465,  150.64556635],
	  			       [ -21.39857172,  150.70427386],
			   	       [ -21.3451841,   150.76293883],
	    			       [ -21.29177194,  150.82156139]]), latlong10)

    ''' test_UTM_3 to test_UTM_7 are the old test examples'''
    def test_UTM_3(self):
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


    def test_UTM_4(self):
        #TEST 4

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
        
        
    def test_UTM_5(self):
        #Test 5
        lat = degminsec2decimal_degrees(-60,0,0)
        lon = degminsec2decimal_degrees(130,0,0) 

        zone, easting, northing = LLtoUTM(lat,lon)

        assert zone == 52
        assert num.allclose(easting, 555776.267)
        assert num.allclose(northing, 3348167.264)

        Lat, Long = UTMtoLL(northing, easting, zone)

    def test_UTM_6(self):
        #Test 6 (Kobenhavn, Northern hemisphere)
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
 
    def test_UTM_7(self):
        #Test 7 (Wollongong)

        lat = degminsec2decimal_degrees(-34,30,0.)
        lon = degminsec2decimal_degrees(150,55,0.) 
        
        zone, easting, northing = LLtoUTM(lat,lon)

        assert zone == 56
        assert num.allclose(easting, 308728.009)
        assert num.allclose(northing, 6180432.601)

        lat_calced, long_calced = UTMtoLL(northing, easting, zone) 
        assert num.allclose(lat,  lat_calced)
        assert num.allclose(lon, long_calced)

#------------------------------------------------------------

if __name__ == "__main__":
    mysuite = unittest.makeSuite(TestCase,'test')
    runner = unittest.TextTestRunner()
    runner.run(mysuite)
