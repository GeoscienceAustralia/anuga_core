import unittest
import numpy as num
from anuga.tsunami_source.eqf import earthquake_tsunami, Okada_func


class Test_eq(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_Okada_func(self):
        x0 = 38000.0
        y0 = 12500.0
        length = 25.0
        width = 5.0
        strike = 0.0
        depth = 3500.
        slip = 10.0
        dip = 15.0
        rake = 90.0
        test = -2205.904774487

        #print 'depth', depth
        z = Okada_func(length=length, width=width, dip=dip, \
                       x0=x0, y0=y0, strike=strike, depth=depth, \
                       slip=slip, rake=rake, test=test)

        assert num.allclose(z.length, length)
        assert num.allclose(z.width, width)
        assert num.allclose(z.x0, x0)
        assert num.allclose(z.y0, y0)
        assert num.allclose(z.strike, strike)
        assert num.allclose(z.depth, depth)
        assert num.allclose(z.slip, slip)
        assert num.allclose(z.dip, dip)
        assert num.allclose(z.rake, rake)

        #print 'in test', z.test
        assert num.allclose(z.test, -2205.904774487)
        
        #print 'hello finished okada'

    def test_earthquake_tsunami(self):

        x0 = 0.0
        y0 = 0.0
        length = 25.0
        width = 5.0
        strike = 0.0
        depth = 3500.
        slip = 10.0
        dip = 15.0
        
        eq = earthquake_tsunami(length=length, width = width, depth=depth, \
                                strike = strike, dip = dip, slip = slip, \
                                x0 = x0, y0 = y0)

#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_eq,'test_Okada_func')
    runner = unittest.TextTestRunner()
    runner.run(suite)

