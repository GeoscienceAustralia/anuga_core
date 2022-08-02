"""  Test timezone
"""

import anuga
import numpy as num
import unittest
import os
from datetime import datetime

try:
    from zoneinfo import ZoneInfo
except:
    from backports.zoneinfo import ZoneInfo

#from zoneinfo import ZoneInfo

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

        UTC = ZoneInfo('UTC')

        assert domainTZ == UTC

    def test_set_timezone_zoneinfo(self):

        domain = anuga.rectangular_cross_domain(10,10)
        AEST = ZoneInfo('Australia/Sydney')
        domain.set_timezone(AEST)
        domainTZ = domain.get_timezone()
        assert AEST == domainTZ


    def test_set_timezone_str(self):

        domain = anuga.rectangular_cross_domain(10,10)

        domain.set_timezone('Australia/Sydney')
        domainTZ = domain.get_timezone()

        AEST = ZoneInfo('Australia/Sydney')
        assert AEST == domainTZ

    def test_starttime_with_datetime(self):

        domain = anuga.rectangular_cross_domain(10,10)
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left' :Br , 'right' : Br, 'top' : Br, 'bottom' : Br})
        
        AEST = ZoneInfo('Australia/Sydney')
        dt = datetime(2021,3,21,18,30, tzinfo=AEST)

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

        AEST = ZoneInfo('Australia/Sydney')
        dt = datetime(2021,3,21,18,30, tzinfo=AEST)

        domain.set_starttime(dt)

        assert str(domain.get_datetime()) == '2021-03-21 07:30:00+00:00'
    
        for t in domain.evolve(yieldstep=30, duration=60):
            pass

        assert str(domain.get_datetime()) == '2021-03-21 07:31:00+00:00'

    def test_domainTZ_starttime_datetime(self):

        domain = anuga.rectangular_cross_domain(10,10)
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left' :Br , 'right' : Br, 'top' : Br, 'bottom' : Br})
        
        AEST = ZoneInfo('Australia/Sydney')
        domain.set_timezone(AEST)

        UTC = ZoneInfo('UTC')
        dt = datetime(2021,3,21,7,30, tzinfo=UTC)
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

        AEST = ZoneInfo('Australia/Sydney')

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
