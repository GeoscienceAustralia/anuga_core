"""  Test timezone
"""

import anuga
import numpy as num
import pytz
import unittest
import os
from datetime import datetime


class Test_Timzone(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        for file in ['domain.sww']:
            try:
                os.remove(file)
            except:
                pass
        


    def test_default_TZ(self):

        domain = anuga.rectangular_cross_domain(10,10)

        domainTZ = domain.get_timezone()

        UTC = pytz.utc

        assert domainTZ == UTC

    def test_set_timezone_pytz(self):

        domain = anuga.rectangular_cross_domain(10,10)
        AEST = pytz.timezone('Australia/Sydney')
        domain.set_timezone(AEST)
        domainTZ = domain.get_timezone()
        assert AEST == domainTZ


    def test_set_timezone_str(self):

        domain = anuga.rectangular_cross_domain(10,10)

        domain.set_timezone('Australia/Sydney')
        domainTZ = domain.get_timezone()

        AEST = pytz.timezone('Australia/Sydney')
        assert AEST == domainTZ

    def test_starttime_with_datetime(self):

        domain = anuga.rectangular_cross_domain(10,10)
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left' :Br , 'right' : Br, 'top' : Br, 'bottom' : Br})
        
        AEST = pytz.timezone('Australia/Sydney')
        dt = AEST.localize(datetime(2021,3,21,18,30))

        domain.set_starttime(dt)

        # The domain timezone is UTC
        assert str(domain.get_datetime()) == '2021-03-21 07:30:00+00:00'
    
        for t in domain.evolve(yieldstep=30, duration=60):
            pass

        assert str(domain.get_datetime()) == '2021-03-21 07:31:00+00:00'
        
    def test_starttime_with_naive_datetime(self):

        domain = anuga.rectangular_cross_domain(10,10)
        
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left' :Br , 'right' : Br, 'top' : Br, 'bottom' : Br})

        AEST = pytz.timezone('Australia/Sydney')
        dt = AEST.localize(datetime(2021,3,21,18,30))

        domain.set_starttime(dt)

        assert str(domain.get_datetime()) == '2021-03-21 07:30:00+00:00'
    
        for t in domain.evolve(yieldstep=30, duration=60):
            pass

        assert str(domain.get_datetime()) == '2021-03-21 07:31:00+00:00'

    def test_domainTZ_starttime_datetime(self):

        domain = anuga.rectangular_cross_domain(10,10)
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left' :Br , 'right' : Br, 'top' : Br, 'bottom' : Br})
        
        AEST = pytz.timezone('Australia/Sydney')
        domain.set_timezone(AEST)

        UTC = pytz.utc
        dt = UTC.localize(datetime(2021,3,21,7,30))
        domain.set_starttime(dt)

        # The domain timezone is AEST
        assert str(domain.get_datetime()) == '2021-03-21 18:30:00+11:00'
    
        for t in domain.evolve(yieldstep=30, duration=60):
            pass

        assert str(domain.get_datetime()) == '2021-03-21 18:31:00+11:00'
        
    def test_domainTZ_starttime_naive_datetime(self):

        domain = anuga.rectangular_cross_domain(10,10)
        
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left' :Br , 'right' : Br, 'top' : Br, 'bottom' : Br})

        AEST = pytz.timezone('Australia/Sydney')

        domain.set_timezone(AEST)

        dt = datetime(2021,3,21,18,30)

        domain.set_starttime(dt)

        assert str(domain.get_datetime()) == '2021-03-21 18:30:00+11:00'
    
        for t in domain.evolve(yieldstep=30, duration=60):
            pass

        assert str(domain.get_datetime()) == '2021-03-21 18:31:00+11:00'
            
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Timzone, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
