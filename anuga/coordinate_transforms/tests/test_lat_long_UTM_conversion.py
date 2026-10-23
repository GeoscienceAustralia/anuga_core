
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

class TestCase_extra(unittest.TestCase):
    """Cover Norway/Svalbard special zones and _UTMLetterDesignator."""

    def test_norway_special_zone(self):
        """LLtoUTM Norway special case sets ZoneNumber=32 (line 82)."""
        # Oslo approx: lat ~60°N, lon ~10°E → falls in Norway band
        lat = 60.0
        lon = 10.0
        zone, e, n = LLtoUTM(lat, lon)
        self.assertEqual(zone, 32)

    def test_svalbard_zone_31(self):
        """LLtoUTM Svalbard zone 31 (line 84)."""
        lat = 75.0
        lon = 5.0  # < 9° → zone 31
        zone, e, n = LLtoUTM(lat, lon)
        self.assertEqual(zone, 31)

    def test_svalbard_zone_33(self):
        """LLtoUTM Svalbard zone 33 (line 85)."""
        lat = 75.0
        lon = 15.0  # 9°-21° → zone 33
        zone, e, n = LLtoUTM(lat, lon)
        self.assertEqual(zone, 33)

    def test_svalbard_zone_35(self):
        """LLtoUTM Svalbard zone 35 (line 86)."""
        lat = 75.0
        lon = 25.0  # 21°-33° → zone 35
        zone, e, n = LLtoUTM(lat, lon)
        self.assertEqual(zone, 35)

    def test_svalbard_zone_37(self):
        """LLtoUTM Svalbard zone 37 (line 87)."""
        lat = 75.0
        lon = 35.0  # 33°-42° → zone 37
        zone, e, n = LLtoUTM(lat, lon)
        self.assertEqual(zone, 37)

    def test_letter_designator_all_bands(self):
        """_UTMLetterDesignator covers all latitude bands (lines 100-120)."""
        from anuga.coordinate_transforms.lat_long_UTM_conversion import \
            _UTMLetterDesignator
        # Bands: X(72-84), W(64-72), V(56-64), U(48-56), T(40-48),
        #        S(32-40), R(24-32), Q(16-24), P(8-16), N(0-8),
        #        M(-8-0), L(-16--8), K(-24--16), J(-32--24), H(-40--32),
        #        G(-48--40), F(-56--48), E(-64--56), D(-72--64), C(-80--72)
        self.assertEqual(_UTMLetterDesignator(80), 'X')
        self.assertEqual(_UTMLetterDesignator(68), 'W')
        self.assertEqual(_UTMLetterDesignator(60), 'V')
        self.assertEqual(_UTMLetterDesignator(52), 'U')
        self.assertEqual(_UTMLetterDesignator(44), 'T')
        self.assertEqual(_UTMLetterDesignator(36), 'S')
        self.assertEqual(_UTMLetterDesignator(28), 'R')
        self.assertEqual(_UTMLetterDesignator(20), 'Q')
        self.assertEqual(_UTMLetterDesignator(12), 'P')
        self.assertEqual(_UTMLetterDesignator(4), 'N')
        self.assertEqual(_UTMLetterDesignator(-4), 'M')
        self.assertEqual(_UTMLetterDesignator(-12), 'L')
        self.assertEqual(_UTMLetterDesignator(-20), 'K')
        self.assertEqual(_UTMLetterDesignator(-28), 'J')
        self.assertEqual(_UTMLetterDesignator(-36), 'H')
        self.assertEqual(_UTMLetterDesignator(-44), 'G')
        self.assertEqual(_UTMLetterDesignator(-52), 'F')
        self.assertEqual(_UTMLetterDesignator(-60), 'E')
        self.assertEqual(_UTMLetterDesignator(-68), 'D')
        self.assertEqual(_UTMLetterDesignator(-76), 'C')
        self.assertEqual(_UTMLetterDesignator(-85), 'Z')  # outside UTM limits


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCase)
    runner = unittest.TextTestRunner()
    runner.run(suite)
