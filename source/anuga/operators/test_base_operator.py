from anuga import Domain
from anuga import Dirichlet_boundary
from anuga import Operator

import numpy as num
import unittest

class Test_Operator(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_create_operator(self):
        points = num.array([[0.0,0.0],[1.0,0.0],[0.0,1.0]])
        
        elements = num.array([[0,1,2]])
        boundary_map = {}
        boundary_map[(0,0)] = 'edge0'
        boundary_map[(0,1)] = 'edge1'
        boundary_map[(0,2)] = 'edge2'

        domain = Domain(points, elements, boundary_map)

        operator = Operator(domain)

        message = operator.statistics()
        assert message == 'You need to implement operator statistics for your operator'

        message = operator.timestepping_statistics()
        assert message == 'You need to implement timestepping statistics for your operator'

        domain.timestep = 3.0

        assert operator.get_timestep() == domain.get_timestep()
        
        try:
            operator()
        except:
            pass
        else:
            raise Exception('should have raised an exception')


################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Operator, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
