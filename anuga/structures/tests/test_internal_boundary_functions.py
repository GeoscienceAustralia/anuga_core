import unittest
import os

import numpy
import anuga
from anuga.structures import internal_boundary_functions
from anuga.structures.internal_boundary_functions import hecras_internal_boundary_function
from anuga.structures.internal_boundary_functions import pumping_station_function
from anuga.utilities.system_tools import get_pathname_from_package

## Data for the domain used to test the pumping station function
boundaryPolygon = [ [0., 0.], [0., 100.], [100.0, 100.0], [100.0, 0.0]]
wallLoc = 50.
# The boundary polygon + riverwall breaks the mesh into multiple regions
# Must define the resolution in these areas with an xy point + maximum area
# Otherwise triangle.c gets confused
regionPtAreas = [ [99., 99., 20.0*20.0*0.5],
                  [1., 1., 20.0*20.0*0.5] ]

class Test_internal_boundary_functions(unittest.TestCase):

    def setUp(self):

        path = get_pathname_from_package('anuga.structures')
        self.input_hecras_file = os.path.join(path, 'tests', 'data', 'hecras_bridge_table.csv')

        return

    def tearDown(self):
        pass

        return

    def create_domain(self, wallHeight, InitialOceanStage, InitialLandStage, 
                      riverWall=None, riverWall_Par=None):
        # Riverwall = list of lists, each with a set of x,y,z (and optional QFactor) values

        if(riverWall is None):
            riverWall = { 'centralWall':
                           [ [wallLoc, 0.0, wallHeight],
                             [wallLoc, 100.0, wallHeight]] 
                        }
        if(riverWall_Par is None):
            riverWall_Par = {'centralWall':{'Qfactor':1.0}}
        # Make the domain
        anuga.create_mesh_from_regions(boundaryPolygon, 
                                 boundary_tags={'left': [0],
                                                'top': [1],
                                                'right': [2],
                                                'bottom': [3]},
                                   maximum_triangle_area = 200.,
                                   minimum_triangle_angle = 28.0,
                                   filename = 'testRiverwall.msh',
                                   interior_regions =[ ], #[ [higherResPolygon, 1.*1.*0.5],
                                                          #  [midResPolygon, 3.0*3.0*0.5]],
                                   breaklines=riverWall.values(),
                                   use_cache=False,
                                   verbose=False,
                                   regionPtArea=regionPtAreas)

        domain = anuga.create_domain_from_file('testRiverwall.msh')

        # 05/05/2014 -- riverwalls only work with DE0 and DE1
        domain.set_flow_algorithm('DE0')
        domain.set_name('test_riverwall')

        def topography(x,y):
            return -x/150. 

        def stagefun(x,y):
            stg = InitialOceanStage*(x>=wallLoc) + InitialLandStage*(x<wallLoc)
            return stg 

        # NOTE: Setting quantities at centroids is important for exactness of tests
        domain.set_quantity('elevation',topography,location='centroids')     
        domain.set_quantity('stage', stagefun,location='centroids')            
        
        domain.riverwallData.create_riverwalls(riverWall,riverWall_Par,verbose=False) 

        # Boundary conditions
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom':Br})

        return domain

    def test_basic_properties(self):

        self.hb = hecras_internal_boundary_function(self.input_hecras_file, verbose=False)
        assert self.hb.name == self.input_hecras_file
        assert len(self.hb.hw_max_given_tw) == len(self.hb.nonfree_flow_tw)
        assert len(self.hb.nonfree_flow_curves) == len(self.hb.nonfree_flow_tw)

        return

    def test_hecras_internal_boundary_function_values(self):

        self.hb = hecras_internal_boundary_function(self.input_hecras_file, verbose=False)

        # Stationary states #
        assert numpy.allclose(self.hb(-3., -3.), 0.)
        assert numpy.allclose(self.hb(-2.5, -2.5), 0.)
        assert numpy.allclose(self.hb(0., 0.), 0.)
        assert numpy.allclose(self.hb(1., 1.), 0.)
        assert numpy.allclose(self.hb(2.8, 2.8), 0.)
        assert numpy.allclose(self.hb(-2.146, -2.146), 0.)

        # Known values (looked up manually in table) #

        # Top of one nonfree-flow curve
        assert numpy.allclose(self.hb(-2.759, -2.947), 9.425)

        # Slightly above previous (i.e. free flow curve)
        assert numpy.allclose(self.hb(-2.747, -2.947), 9.736)

        # Value on most extreme nonfree flow curve
        assert numpy.allclose(self.hb(2.909, 2.894), 5.118)

        # Value on least-extreme nonfree flow curve
        assert numpy.allclose(self.hb(-3.26, -3.266), 0.27)

        # Value with hw below all curves
        assert numpy.allclose(self.hb(-4., -4.1), 0.)

        # Value below least-extreme nonfree flow curve, but hw still valid on
        # free flow curve
        assert numpy.allclose(self.hb(-2.747, -3.4), 9.736)

        # Near top of free flow curve
        assert numpy.allclose(self.hb(2.468, -1.4), 82.89)
        assert numpy.allclose(self.hb(2.468, -2.0), 82.89)
        assert numpy.allclose(self.hb(2.468, -5.0), 82.89)

        # Logical ordering
        assert self.hb(2.4, 2.0) > self.hb(2.2, 2.0)
        assert self.hb(1.5, 1.0) > self.hb(1.5, 1.01)

        # Things which should throw errors
        def tw_too_big():
            return self.hb(9.00, 2.0)

        def hw_too_big():
            return self.hb(10.00, 0.)

        self.assertRaises(Exception, lambda: tw_too_big())
        self.assertRaises(Exception, lambda: hw_too_big())


        # Check case where sign reversal is allowed
        self.hb.allow_sign_reversal=True
        Q1 = self.hb( -2.75, -2.95)
        Q2 = self.hb( -2.95, -2.75)
        assert numpy.allclose(Q1, -Q2)

        return

    def test_pumping_station_function(self):

        domain = self.create_domain(wallHeight=10., 
            InitialOceanStage=4., 
            InitialLandStage=0.)

        ps_function = pumping_station_function(
            domain=domain,
            pump_capacity=1.0,
            hw_to_start_pumping=0.0,
            hw_to_stop_pumping=-1.0,
            initial_pump_rate=0., 
            pump_rate_of_increase = 1.0, 
            pump_rate_of_decrease = 1.0, 
            verbose=False)


        # Pump should be starting at zero       
        assert(ps_function(1., 1.) == 0.)
    
        # As time increases, so will the pump rate 
        domain.time = 0.5
        assert(numpy.allclose(ps_function(1., 1.), 0.5))
        assert(numpy.allclose(ps_function(1., 1.), 0.5))
        domain.time = 0.75
        assert(numpy.allclose(ps_function(1., 1.), 0.75))

        # Pump rate maximum = 1.0
        domain.time = 1.5
        assert(numpy.allclose(ps_function(1., 1.), 1.0))

        # Force the pump rate to zero, and it should start increasing as time
        # increases
        ps_function.pump_rate = 0.
        assert (ps_function(1., 1.) == 0.)
        domain.time = 1.75
        assert (ps_function(1., 1.) == 0.25)
        domain.time = 2.5
        assert (ps_function(1., 1.) == 1.0)
        # Can't increase any more
        domain.time = 3.0 
        assert (ps_function(1., 1.) == 1.0)
       
        # Let's try to decrease the pump
        assert (ps_function(-1.5, 1.) == 1.0)
        domain.time = 3.5
        assert (ps_function(-1.5, 1.) == 0.5)
        domain.time = 3.9

        assert (numpy.allclose(ps_function(-1.5, 1.), 0.1))
        domain.time = 4.0 
        assert (numpy.allclose(ps_function(-1.5, 1.), 0.0))
        domain.time = 5.0 
        assert (numpy.allclose(ps_function(-1.5, 1.), 0.0))

        
        return 


if __name__ == "__main__":
    suite = unittest.makeSuite(Test_internal_boundary_functions, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
