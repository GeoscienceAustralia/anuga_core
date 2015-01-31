import unittest
import anuga
import numpy
import os

boundaryPolygon=[ [0., 0.], [0., 100.], [100.0, 100.0], [100.0, 0.0]]

verbose=False

class Test_local_extrapolation_and_flux_updating(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        for file in ['domain.sww', 'test_boundaryfluxintegral.msh', 'test_boundaryfluxintegral.sww']:
            try:
                os.remove(file)
            except:
                pass        


    def create_domain(self, flowalg):
        # Riverwall = list of lists, each with a set of x,y,z (and optional QFactor) values

        # Make the domain
        anuga.create_mesh_from_regions(boundaryPolygon, 
                                 boundary_tags={'left': [0],
                                                'top': [1],
                                                'right': [2],
                                                'bottom': [3]},
                                   maximum_triangle_area = 200.,
                                   minimum_triangle_angle = 28.0,
                                   filename = 'test_boundaryfluxintegral.msh',
                                   use_cache=False,
                                   verbose=verbose)

        domain=anuga.create_domain_from_file('test_boundaryfluxintegral.msh')

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

    def test_local_extrapolation_and_flux_updating_DE0(self):
        """

        We check that results with and without local_extrapolation_and_flux_updating are 'close enough'

        We also check that the total volume and boundary flux integral are equal for both methods
        """
        
        domain=self.create_domain('DE0')
        for t in domain.evolve(yieldstep=0.1,finaltime=20.0):
            pass
        # The domain was initially dry
        vol=domain.get_water_volume()
        boundaryFluxInt=domain.get_boundary_flux_integral()


        domain2=self.create_domain('DE0')
        domain2.set_local_extrapolation_and_flux_updating(nlevels=8)
        for t in domain2.evolve(yieldstep=0.1,finaltime=20.0):
            pass
        # The domain was initially dry
        vol2=domain2.get_water_volume()
        boundaryFluxInt2=domain2.get_boundary_flux_integral()
        
        assert(numpy.allclose(vol,vol2, rtol=0.05))
        assert(numpy.allclose(vol,boundaryFluxInt))
        assert(numpy.allclose(vol2,boundaryFluxInt2))
        assert( numpy.all(abs(domain.quantities['stage'].centroid_values-domain2.quantities['stage'].centroid_values) <0.02))
        
        return
    
    def test_local_extrapolation_and_flux_updating_DE1(self):
        """

        LEAFU should fail for DE1 since we don't have support for rk2 timestepping yet

        """
        
        domain=self.create_domain('DE1')
        Failed=True
        try:
            domain.set_local_extrapolation_and_flux_updating(nlevels=8)
            Failed=False
        except:
            pass

        assert(Failed==True)

        return

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_local_extrapolation_and_flux_updating, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

