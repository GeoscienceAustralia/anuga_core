
import unittest
import anuga
import numpy
import os

boundaryPolygon=[ [0., 0.], [0., 100.], [100.0, 100.0], [100.0, 0.0]]

verbose=False

class Test_boundary_flux_integral_operator(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('test_boundaryfluxintegral.msh')
        except:
            pass


        try:
            os.remove('test_boundaryfluxintegral.sww')
        except:
            pass

    def create_domain(self, flowalg):
        # Riverwall = list of lists, each with a set of x,y,z (and optional QFactor) values

        # Make the domain
        domain = anuga.create_domain_from_regions(boundaryPolygon, 
                                 boundary_tags={'left': [0],
                                                'top': [1],
                                                'right': [2],
                                                'bottom': [3]},
                                   mesh_filename='test_boundaryfluxintegral.msh',
                                   maximum_triangle_area = 200.,
                                   minimum_triangle_angle = 28.0,
                                   use_cache=False,
                                   verbose=verbose)


        # 05/05/2014 -- riverwalls only work with DE0 and DE1
        domain.set_flow_algorithm(flowalg)
        domain.set_name('test_boundaryfluxintegral')

        domain.set_store_vertices_uniquely()
       
        def topography(x,y):
            return -x/150. 

        # NOTE: Setting quantities at centroids is important for exactness of tests
        domain.set_quantity('elevation',topography,location='centroids')     
        domain.set_quantity('friction',0.03)             
        domain.set_quantity('stage', topography,location='centroids')            
       
        # Boundary conditions
        Br=anuga.Reflective_boundary(domain)
        Bd=anuga.Dirichlet_boundary([0., 0., 0.])
        domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom':Br})

        return domain

    def test_boundary_flux_operator_DE0(self):
        """
        A (the) boundary flux operator is instantiated when a domain is created.
        This tests the calculation for euler timestepping 
        """
  
        flowalg = 'DE0'
          
        domain=self.create_domain(flowalg)
  
        #domain.print_statistics()
        for t in domain.evolve(yieldstep=1.0,finaltime=5.0):
            if verbose: domain.print_timestepping_statistics()
            if verbose: print(domain.get_water_volume())
            pass
        # The domain was initially dry
        vol=domain.get_water_volume()
        boundaryFluxInt=domain.get_boundary_flux_integral()
  
        if verbose: print(flowalg, vol, boundaryFluxInt)        
        assert(numpy.allclose(vol,boundaryFluxInt))
         

        
    def test_boundary_flux_operator_DE1(self):
        """
        A (the) boundary flux operator is instantiated when a domain is created.
        This tests the calculation for rk2 timestepping 
        """
        flowalg = 'DE1'
                 
        domain=self.create_domain(flowalg)
        #domain.print_statistics()
        for t in domain.evolve(yieldstep=1.0,finaltime=5.0):
            if verbose: domain.print_timestepping_statistics()
            if verbose: print(domain.get_water_volume())
            pass
        # The domain was initially dry
        vol=domain.get_water_volume()
        boundaryFluxInt=domain.get_boundary_flux_integral()
         
        if verbose: print(flowalg, vol, boundaryFluxInt)
        assert(numpy.allclose(vol,boundaryFluxInt))
        
        

    def test_boundary_flux_operator_DE2(self):
        """
        A (the) boundary flux operator is instantiated when a domain is created.
        This tests the calculation for rk3 timestepping 
        """
 
        flowalg = 'DE2'
 
        domain=self.create_domain(flowalg)
        #domain.print_statistics()
        for t in domain.evolve(yieldstep=1.0,finaltime=5.0):
            if verbose: domain.print_timestepping_statistics()
            if verbose: print(domain.get_water_volume(), domain.get_boundary_flux_integral())
            pass
        # The domain was initially dry
        vol=domain.get_water_volume()
        boundaryFluxInt=domain.get_boundary_flux_integral()
 
        if verbose: print(flowalg, vol, boundaryFluxInt)
        assert(numpy.allclose(vol,boundaryFluxInt))        
         
        
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_boundary_flux_integral_operator, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

