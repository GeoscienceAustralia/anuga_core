#!/usr/bin/env python

import tempfile
import unittest

from model import Model
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross



class modelTestCase(unittest.TestCase):
    def setUp(self):
        # construct and name model
        self.model = Model('test_model') 
        pass

    def tearDown(self):
        pass

    def test_construction(self):
        """ Test initial setup of model. """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #              bac,     bce,     ecf,     dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        self.model.set_geometry(points, vertices)

    def test_simple_example(self):
        """ Run a new version of the channel1 example. """
        from anuga.shallow_water import Reflective_boundary
        from anuga.shallow_water import Dirichlet_boundary
        
        def topography(x, y):
            return -x/10                         

        # set quantities
        self.model.set_quantity('elevation', topography) 
        self.model.set_quantity('friction', 0.01)         
        self.model.set_quantity('stage', expression='elevation')  

        # set properties of boundaries (reflective, flow in, etc.)
        Bi = Dirichlet_boundary([0.4, 0, 0])         # Inflow
        self.model.set_boundary({'left': Bi, 'right': Bi, 'top': Bi, 'bottom': Bi})

        # set the geometry to use (may be a mesh file, or points/vertices tuple)
        self.model.set_geometry(*rectangular_cross(2, 2, len1=10.0, len2=5.0))

        # build, then run the simulation
        self.model.build()
        self.model.run(0.5, 1.0)


    def test_wrong_input_order(self):
        """The user tries to build before model is defined. """
        
        self.model.set_quantity('stage',                 # Dry bed
                            expression='elevation')          

        try:
            self.model.build()            
        except:
            pass
        else:
            msg = 'Should have raised exception for missing mesh'
            raise Exception, msg


    def test_duplicate_geometry(self):
        """The user tries to assign geometry twice. """

        self.model.set_geometry(*rectangular_cross(2, 2, len1=10.0, len2=5.0))

        try:
            self.model.set_geometry(*rectangular_cross(2, 2, len1=10.0, len2=5.0))  
        except:
            pass
        else:
            msg = 'Should have raised bad input exception'
            raise Exception, msg


    def test_input_after_build(self):
        """The user tries to change built model. """

        self.model.set_geometry(*rectangular_cross(2, 2, len1=10.0, len2=5.0))
        self.model.build()
        
        try:
            self.model.set_quantity('friction', 0.01)  
        except:
            pass
        else:
            msg = 'Should have raised exception because model already built'
            raise Exception, msg



################################################################################

if __name__ == '__main__':
    suite = unittest.makeSuite(modelTestCase,'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)

