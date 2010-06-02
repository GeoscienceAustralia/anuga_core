import os
import unittest
import tempfile
import numpy as num

from anuga.coordinate_transforms.geo_reference import Geo_reference
from csv_file import load_csv_as_array, load_csv_as_dict
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.shallow_water.shallow_water_domain import Domain
from sww import load_sww_as_domain, weed, get_mesh_and_quantities_from_file

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



    def test_get_mesh_and_quantities_from_sww_file(self):
        """test_get_mesh_and_quantities_from_sww_file(self):
        """     
        
        # Generate a test sww file with non trivial georeference
        
        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        # Setup
        from mesh_factory import rectangular

        # Create basic mesh (100m x 5m)
        width = 5
        length = 50
        t_end = 10
        points, vertices, boundary = rectangular(length, width, 50, 5)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference = Geo_reference(56,308500,6189000))

        domain.set_name('flowtest')
        swwfile = domain.get_name() + '.sww'
        domain.set_datadir('.')

        Br = Reflective_boundary(domain)    # Side walls
        Bd = Dirichlet_boundary([1, 0, 0])  # inflow

        domain.set_boundary( {'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = t_end):
            pass

        
        # Read it

        # Get mesh and quantities from sww file
        X = get_mesh_and_quantities_from_file(swwfile,
                                              quantities=['elevation',
                                                          'stage',
                                                          'xmomentum',
                                                          'ymomentum'], 
                                              verbose=False)
        mesh, quantities, time = X
        

        # Check that mesh has been recovered
        assert num.alltrue(mesh.triangles == domain.get_triangles())
        assert num.allclose(mesh.nodes, domain.get_nodes())

        # Check that time has been recovered
        assert num.allclose(time, range(t_end+1))

        # Check that quantities have been recovered
        # (sww files use single precision)
        z=domain.get_quantity('elevation').get_values(location='unique vertices')
        assert num.allclose(quantities['elevation'], z)

        for q in ['stage', 'xmomentum', 'ymomentum']:
            # Get quantity at last timestep
            q_ref=domain.get_quantity(q).get_values(location='unique vertices')

            #print q,quantities[q]
            q_sww=quantities[q][-1,:]

            msg = 'Quantity %s failed to be recovered' %q
            assert num.allclose(q_ref, q_sww, atol=1.0e-6), msg
            
        # Cleanup
        os.remove(swwfile)
        
        

    def test_weed(self):
        coordinates1 = [[0.,0.],[1.,0.],[1.,1.],[1.,0.],[2.,0.],[1.,1.]]
        volumes1 = [[0,1,2],[3,4,5]]
        boundary1= {(0,1): 'external',(1,2): 'not external',(2,0): 'external',(3,4): 'external',(4,5): 'external',(5,3): 'not external'}
        coordinates2,volumes2,boundary2=weed(coordinates1,volumes1,boundary1)

        points2 = {(0.,0.):None,(1.,0.):None,(1.,1.):None,(2.,0.):None}

        assert len(points2)==len(coordinates2)
        for i in range(len(coordinates2)):
            coordinate = tuple(coordinates2[i])
            assert points2.has_key(coordinate)
            points2[coordinate]=i

        for triangle in volumes1:
            for coordinate in triangle:
                assert coordinates2[points2[tuple(coordinates1[coordinate])]][0]==coordinates1[coordinate][0]
                assert coordinates2[points2[tuple(coordinates1[coordinate])]][1]==coordinates1[coordinate][1]

#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_sww, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
