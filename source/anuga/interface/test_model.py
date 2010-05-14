#!/usr/bin/env python

import tempfile
import unittest

from model import Model
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross



class modelTestCase(unittest.TestCase):
    def setUp(self):
        self.model = Model()     
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
        
        self.model.set_geometry(*rectangular_cross(10, 5, len1=10.0, len2=5.0))  # Create domain        
        self.model.set_name('channel1')                  # Output name

        #------------------------------------------------------------------------------
        # Setup initial conditions
        #------------------------------------------------------------------------------
        def topography(x, y):
            return -x/10                             # linear bed slope

        self.model.set_quantity('elevation', topography) # Use function for elevation
        self.model.set_quantity('friction', 0.01)        # Constant friction 
        self.model.set_quantity('stage',                 # Dry bed
                            expression='elevation')  

        self.model.build()
        
        #------------------------------------------------------------------------------
        # Setup boundary conditions
        #------------------------------------------------------------------------------
        Bi = Dirichlet_boundary([0.4, 0, 0])         # Inflow
    #    Br = Reflective_boundary(domain)             # Solid reflective wall

     #   domain.set_boundary({'left': Bi, 'right': Bi, 'top': Bi, 'bottom': Bi})

        #------------------------------------------------------------------------------
        # Evolve system through time
        #------------------------------------------------------------------------------
   #     for t in domain.evolve(yieldstep=0.2, finaltime=40.0):
   #         print domain.timestepping_statistics()        
        


################################################################################

if __name__ == '__main__':
    suite = unittest.makeSuite(modelTestCase,'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)

