import os
import unittest
import tempfile
import numpy as num

from anuga.coordinate_transforms.geo_reference import Geo_reference
from csv_file import load_csv_as_array, load_csv_as_dict
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.shallow_water.shallow_water_domain import Domain
from sww import load_sww_as_domain

# boundary functions
from anuga.shallow_water.boundaries import Reflective_boundary, \
            Field_boundary, Transmissive_momentum_set_stage_boundary, \
            Transmissive_stage_zero_momentum_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary, Dirichlet_boundary, \
            Time_boundary, File_boundary, AWI_boundary


class Test_sww(unittest.TestCase):
    def setUp(self):
        self.verbose = False
        pass

    def tearDown(self):
        pass
        
    def test_sww2domain1(self):
        ################################################
        #Create a test domain, and evolve and save it.
        ################################################
        from mesh_factory import rectangular

        #Create basic mesh

        yiel=0.01
        points, vertices, boundary = rectangular(10,10)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.geo_reference = Geo_reference(56,11,11)
        domain.smooth = False
        domain.store = True
        domain.set_name('bedslope')
        domain.default_order=2
        #Bed-slope and friction
        domain.set_quantity('elevation', lambda x,y: -x/3)
        domain.set_quantity('friction', 0.1)
        # Boundary conditions
        from math import sin, pi
        Br = Reflective_boundary(domain)
        Bt = Transmissive_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])
        Bw = Time_boundary(domain=domain,f=lambda t: [(0.1*sin(t*2*pi)), 0.0, 0.0])

        #domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})

        domain.quantities_to_be_stored['xmomentum'] = 2
        domain.quantities_to_be_stored['ymomentum'] = 2
        #Initial condition
        h = 0.05
        elevation = domain.quantities['elevation'].vertex_values
        domain.set_quantity('stage', elevation + h)

        domain.check_integrity()
        #Evolution
        #domain.tight_slope_limiters = 1
        for t in domain.evolve(yieldstep = yiel, finaltime = 0.05):
            #domain.write_time()
            pass



        filename = domain.datadir + os.sep + domain.get_name() + '.sww'
        domain2 = load_sww_as_domain(filename, None, fail_if_NaN=False,
                                        verbose=self.verbose)
        #points, vertices, boundary = rectangular(15,15)
        #domain2.boundary = boundary
        ###################
        ##NOW TEST IT!!!
        ###################

        os.remove(filename)

        bits = ['vertex_coordinates']
        for quantity in domain.quantities_to_be_stored:
            bits.append('get_quantity("%s").get_integral()' % quantity)
            bits.append('get_quantity("%s").get_values()' % quantity)

        for bit in bits:
            #print 'testing that domain.'+bit+' has been restored'
            #print bit
            #print 'done'
            assert num.allclose(eval('domain.'+bit),eval('domain2.'+bit))

        ######################################
        #Now evolve them both, just to be sure
        ######################################x
        domain.time = 0.
        from time import sleep

        final = .1
        domain.set_quantity('friction', 0.1)
        domain.store = False
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})


        for t in domain.evolve(yieldstep = yiel, finaltime = final):
            #domain.write_time()
            pass

        final = final - (domain2.starttime-domain.starttime)
        #BUT since domain1 gets time hacked back to 0:
        final = final + (domain2.starttime-domain.starttime)

        domain2.smooth = False
        domain2.store = False
        domain2.default_order=2
        domain2.set_quantity('friction', 0.1)
        #Bed-slope and friction
        # Boundary conditions
        Bd2=Dirichlet_boundary([0.2,0.,0.])
        domain2.boundary = domain.boundary
        #print 'domain2.boundary'
        #print domain2.boundary
        domain2.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})
        #domain2.set_boundary({'exterior': Bd})

        domain2.check_integrity()

        for t in domain2.evolve(yieldstep = yiel, finaltime = final):
            #domain2.write_time()
            pass

        ###################
        ##NOW TEST IT!!!
        ##################

        bits = ['vertex_coordinates']

        for quantity in ['elevation','stage', 'ymomentum','xmomentum']:
            bits.append('get_quantity("%s").get_integral()' %quantity)
            bits.append('get_quantity("%s").get_values()' %quantity)

        #print bits
        for bit in bits:
            #print bit
            #print eval('domain.'+bit)
            #print eval('domain2.'+bit)
            
            #print eval('domain.'+bit+'-domain2.'+bit)
            msg = 'Values in the two domains are different for ' + bit
            assert num.allclose(eval('domain.'+bit),eval('domain2.'+bit),
                                rtol=1.e-5, atol=3.e-8), msg



#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_sww, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
