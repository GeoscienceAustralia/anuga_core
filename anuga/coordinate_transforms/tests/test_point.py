
import unittest
from anuga.coordinate_transforms.point import Point
from math import fabs

#-------------------------------------------------------------

class TestCase(unittest.TestCase):

    def setUp(self):
        
        self.eps = 0.001    # Accept 0.1 % relative error
        
        self.RSISE = Point(-35.27456,149.12065)
        self.Home = Point(-35.25629,149.12494)     # 28 Scrivener Street, ACT
        self.Syd = Point(-33.93479,151.16794)      # Sydney Airport
        self.Nadi = Point(-17.75330,177.45148)     # Nadi Airport 
        self.Kobenhavn = Point(55.70248, 12.58364) # Kobenhavn, Denmark


    def testBearingNorth(self):

        eps = 1.0e-12
        
        p1 = Point(0.0,0.0)
        p2 = Point(1.0,0.0)
        p3 = Point(0.0,1.0)

        # Assert that bearing is correct within double precision
        self.assertTrue((fabs(p1.BearingTo(p2)-0) < eps),\
                        'Computed northward bearing: %d, Should have been: %d'\
                         %(p1.BearingTo(p2), 0))


    def testBearingSouth(self):

        eps = 1.0e-12
        
        p1 = Point(0.0,0.0)
        p2 = Point(1.0,0.0)
        p3 = Point(0.0,1.0)

        # Assert that bearing is correct within double precision
        self.assertTrue((fabs(p2.BearingTo(p1)-180) < eps),\
                        'Computed southhward bearing: %d, Should have been: %d'\
                         %(p2.BearingTo(p1), 180))


    def testBearingWest(self):

        eps = 1.0e-12
        
        p1 = Point(0.0,0.0)
        p2 = Point(1.0,0.0)
        p3 = Point(0.0,1.0)

        # Assert that bearing is correct within double precision
        self.assertTrue((fabs(p3.BearingTo(p1)-270) < eps),\
                        'Computed westward bearing: %d, Should have been: %d'\
                         %(p3.BearingTo(p1), 270))
        
    def testBearingEast(self):

        eps = 1.0e-12
        
        p1 = Point(0.0,0.0)
        p2 = Point(1.0,0.0)
        p3 = Point(0.0,1.0)

        # Assert that bearing is correct within double precision
        self.assertTrue((fabs(p1.BearingTo(p3)-90) < eps),\
                        'Computed eastward bearing: %d, Should have been: %d'\
                         %(p1.BearingTo(p3), 90))
        


    def testRSISE2Home(self):
        D = 2068   # True Distance to Home
        B = 11     # True Bearing to Home
        self.assertTrue((fabs(self.RSISE.DistanceTo(self.Home) - D)/D < self.eps),\
                        'Dist to Home failed')
        self.assertTrue((self.RSISE.BearingTo(self.Home) - B == 0),\
                        'Computed bearing to Home: %d, Should have been: %d'\
                         %(self.RSISE.BearingTo(self.Home), B))

        

    def testRSISE2Sydney(self):
        D = 239.5 * 1000   # True Distance to Sydney Airport
        B = 52             # True Bearing to Sydney Airport        
        self.assertTrue((fabs(self.RSISE.DistanceTo(self.Syd) - D)/D < self.eps),\
                        'Dist to Sydney failed')
        self.assertTrue((self.RSISE.BearingTo(self.Syd) - B == 0),\
                        'Computed bearing to Sydney: %d, Should have been: %d'\
                         %(self.RSISE.BearingTo(self.Syd), B))


    def testRSISE2Nadi(self):
        D = 3406.1 * 1000   # True Distance to Nadi Airport
        B = 63              # True Bearing to Nadi Airport        
        self.assertTrue((fabs(self.RSISE.DistanceTo(self.Nadi) - D)/D < self.eps),\
                        'Dist to Nadi failed')
        
        self.assertTrue((self.RSISE.BearingTo(self.Nadi) - B == 0),\
                        'Computed bearing to Nadi: %d, Should have been: %d'\
                         %(self.RSISE.BearingTo(self.Nadi), B))


    def testRSISE2Kobenhavn(self):        
        D = 16025 * 1000   # True Distance to Kobenhavn
        B = 319            # True Bearing to Kobenhavn        
        self.assertTrue((fabs(self.RSISE.DistanceTo(self.Kobenhavn) - D)/D < self.eps),\
                        'Computed Distance to Kobenhavn: %d, Should have been: %d' \
                        %(self.RSISE.DistanceTo(self.Kobenhavn), D))
        self.assertTrue((self.RSISE.BearingTo(self.Kobenhavn) - B == 0),\
                        'Computed Bearing to Kobenhavn: %d, Should have been: %d' \
                        %(self.RSISE.BearingTo(self.Kobenhavn), B))

#-------------------------------------------------------------

if __name__ == "__main__":
    mysuite = unittest.makeSuite(TestCase,'test')
    runner = unittest.TextTestRunner()
    runner.run(mysuite)

