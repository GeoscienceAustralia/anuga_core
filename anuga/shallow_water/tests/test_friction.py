"""  Test environmental forcing - rain, wind, etc.
"""


import unittest, os
import anuga
from anuga import Domain
from anuga import Reflective_boundary

from anuga.shallow_water.friction import *

import numpy as np
import warnings



class Test_Friction(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        for file in ['domain.sww']:
            try:
                os.remove(file)
            except:
                pass
        


    def test_manning_friction_implicit(self):
        """Test the manning friction implicit forcing term
        """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 1.00)
        domain.set_quantity('xmomentum', 1.00)
        domain.set_quantity('ymomentum', 2.00)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        #Test friction forcing term
        domain.compute_forcing_terms()


        import pprint
        pprint.pprint(domain.quantities['xmomentum'].explicit_update)
        pprint.pprint(domain.quantities['ymomentum'].explicit_update)

        pprint.pprint(domain.quantities['xmomentum'].semi_implicit_update)
        pprint.pprint(domain.quantities['ymomentum'].semi_implicit_update)

        xmon_semi_implicit_update = np.array([-21.91346618, -21.91346618, -21.91346618, -21.91346618])
        ymon_semi_implicit_update = np.array([-43.82693236, -43.82693236, -43.82693236, -43.82693236])

        assert num.allclose(domain.quantities['stage'].explicit_update, 0.0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0.0)

        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update, xmon_semi_implicit_update)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update, ymon_semi_implicit_update)


            
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Friction, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
