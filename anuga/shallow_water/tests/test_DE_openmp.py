"""  Test environmental forcing - rain, wind, etc.
"""

import unittest, os

import anuga

from anuga import Reflective_boundary
from anuga import rectangular_cross_domain

from anuga import Domain

import numpy as num
import warnings
import time
import math



class Test_DE_openmp(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        for file in ['runup_openmp.sww']:
            try:
                os.remove(file)
            except:
                pass 


    def test_runup_openmp(self):
        """ Run a version of the validation test runup_sinusoid
        to ensure limiting solution has small velocity
        """

        def create_domain(name='domain'):
            domain = anuga.rectangular_cross_domain(2,2, len1=1., len2=1.)

            domain.set_flow_algorithm('DE0')
            domain.set_low_froude(0)
        
            domain.set_name(name)  
            domain.set_datadir('.')
        
            #------------------
            # Define topography
            #------------------
            scale_me=1.0

            def topography(x,y):
                return (-x/2.0 +0.05*num.sin((x+y)*50.0))*scale_me

            def stagefun(x,y):
                stage=-0.2*scale_me #+0.01*(x>0.9)
                return stage

            domain.set_quantity('elevation',topography)     # Use function for elevation
            domain.set_quantity('friction',0.03)            # Constant friction
            domain.set_quantity('stage', stagefun)          # Constant negative initial stage

            #--------------------------
            # Setup boundary conditions
            #--------------------------
            Br=anuga.Reflective_boundary(domain)                 # Solid reflective wall
            Bd=anuga.Dirichlet_boundary([-0.1*scale_me,0.,0.])   # Constant boundary values -- not used in this example

            #----------------------------------------------
            # Associate boundary tags with boundary objects
            #----------------------------------------------
            domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom':Br})

            return domain

        print('')
        print(70*'=')
        print('Test Runup')
        print(70*'=')


        domain1 = create_domain('domain_original')
        domain1.set_multiprocessor_mode(0)

        domain2 = create_domain('domain_openmp')
        domain2.set_multiprocessor_mode(0) # will change to 2 once burn in

        #------------------------------
        #Evolve the system through time
        #------------------------------
        print('Evolve domain1')
        for t in domain1.evolve(yieldstep=0.1,finaltime=0.1):
            domain1.print_timestepping_statistics()

        print('Evolve domain2')
        for t in domain2.evolve(yieldstep=0.1,finaltime=0.1):
            domain2.print_timestepping_statistics()

        #----------------------------------------
        # Now just run the openmp code on domain2
        #----------------------------------------
        domain2.set_multiprocessor_mode(2)
        timestep = 0.1

        domain1.distribute_to_vertices_and_edges()
        domain1.compute_fluxes()
        timestep1 = domain1.flux_timestep

        domain2.distribute_to_vertices_and_edges()
        domain2.compute_fluxes()
        timestep2 = domain2.flux_timestep

        # Compare update arrays and timestep


        print('domain1 timestep ', timestep1)
        print('domain2 timestep ', timestep2)

        quantities1 = domain1.quantities
        stage1 = quantities1["stage"]
        xmom1 = quantities1["xmomentum"]
        ymom1 = quantities1["ymomentum"]


        quantities2 = domain2.quantities
        stage2 = quantities2["stage"]
        xmom2 = quantities2["xmomentum"]
        ymom2 = quantities2["ymomentum"]


        print('timestep error              ', abs(timestep1-timestep2))
        print('stage explicit update error ', num.linalg.norm(stage1.explicit_update-stage2.explicit_update))
        print('xmom  explicit update error ', num.linalg.norm(xmom1.explicit_update-xmom2.explicit_update))
        print('ymom  explicit update error ', num.linalg.norm(ymom1.explicit_update-ymom2.explicit_update))
        print('edge timestep error         ', num.linalg.norm(domain1.edge_timestep-domain2.edge_timestep))
        print('pressure work error         ', num.linalg.norm(domain1.pressuregrad_work-domain2.pressuregrad_work))
        print('edge flux work error        ', num.linalg.norm(domain1.edge_flux_work-domain2.edge_flux_work))



        assert num.allclose(timestep1,timestep2)
        assert num.allclose(stage1.explicit_update,stage2.explicit_update)
        assert num.allclose(xmom1.explicit_update,xmom2.explicit_update)
        assert num.allclose(ymom1.explicit_update,ymom2.explicit_update)
        assert num.allclose(domain1.edge_timestep,domain2.edge_timestep)
        assert num.allclose(domain1.pressuregrad_work,domain2.pressuregrad_work)
        assert num.allclose(domain1.edge_flux_work,domain2.edge_flux_work)

        # ki3 = num.argmax(num.abs(domain1.edge_flux_work-domain2.edge_flux_work))

        # ki = ki3//3
        # q = ki3%3
        # k = ki//3
        # e = ki%3

        # print('edge_flux_work ki,q,k,e ', ki, q, k, e)

        # import pprint

        # #pprint.pprint(domain1.edge_flux_work)
        # edge_flux_diff = domain2.edge_flux_work- domain1.edge_flux_work
        # edge_timestep_diff =  domain2.edge_timestep- domain1.edge_timestep
        # #pprint.pprint(domain2.edge_flux_work- domain1.edge_flux_work)

        # for k in range(domain2.number_of_elements):
        #     for i in range(3):
        #         ki = 3*k+i
        #         ki3 = 3*ki
        #         print(k,i, domain2.neighbours[k,i], edge_timestep_diff[ki], edge_flux_diff[ki3],edge_flux_diff[ki3+1],edge_flux_diff[ki3+2])



    def test_riverwall_openmp(self):

        def create_domain(name='domain'):

            bounding_polygon = [[0.0, 0.0],
                    [20.0, 0.0],
                    [20.0, 10.0],
                    [0.0, 10.0]]

            boundary_tags={'bottom': [0],
               'right': [1],
               'top': [2],
               'left': [3]
              }


            riverWalls = { 'wall1': [[5.0,0.0,   0.5], [5.0,4.0,  0.5]],
               'wall2': [[15.0,0.0, -0.5], [15.0,4.0,-0.5]],
               'wall3': [[10.0,10.0, 0.0], [10.0,6.0, 0.0]]
             }


              
            domain = anuga.create_domain_from_regions(bounding_polygon, 
                                           boundary_tags,
                                           maximum_triangle_area = 0.4,
                                           breaklines = riverWalls.values())

            domain.set_name(name)

            #Initial Conditions
            domain.set_quantity('elevation', lambda x,y : -x/10, location='centroids') # Use function for elevation
            domain.set_quantity('friction', 0.01, location='centroids')                # Constant friction 
            domain.set_quantity('stage', expression='elevation', location='centroids') # Dry Bed 

            # Boundary Conditions
            Bi = anuga.Dirichlet_boundary([0.4, 0, 0])         # Inflow
            Bo = anuga.Dirichlet_boundary([-2, 0, 0])          # Inflow
            Br = anuga.Reflective_boundary(domain)            # Solid reflective wall

            domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

            # Setup RiverWall
            domain.riverwallData.create_riverwalls(riverWalls, verbose=False)

            return domain

        print('')
        print(70*'=')
        print('Test Riverwall')
        print(70*'=')


        domain1 = create_domain('domain_original')
        domain1.set_multiprocessor_mode(0)

        domain2 = create_domain('domain_openmp')
        domain2.set_multiprocessor_mode(0) # will change to 2 once burn in

        #------------------------------
        #Evolve the system through time
        #------------------------------

        print('Evolve domain1')
        for t in domain1.evolve(yieldstep=0.1,finaltime=0.1):
            domain1.print_timestepping_statistics()

        print('Evolve domain2')
        for t in domain2.evolve(yieldstep=0.1,finaltime=0.1):
            domain2.print_timestepping_statistics()

        #----------------------------------------
        # Now just run the openmp code on domain2
        #----------------------------------------
        domain2.set_multiprocessor_mode(2)
        timestep = 0.1

        domain1.distribute_to_vertices_and_edges()
        domain1.compute_fluxes()
        timestep1 = domain1.flux_timestep

        domain2.distribute_to_vertices_and_edges()
        domain2.compute_fluxes()
        timestep2 = domain2.flux_timestep

        # Compare update arrays and timestep

        print('domain1 timestep ', timestep1)
        print('domain2 timestep ', timestep2)
        

        quantities1 = domain1.quantities
        stage1 = quantities1["stage"]
        xmom1 = quantities1["xmomentum"]
        ymom1 = quantities1["ymomentum"]


        quantities2 = domain2.quantities
        stage2 = quantities2["stage"]
        xmom2 = quantities2["xmomentum"]
        ymom2 = quantities2["ymomentum"]

        print('timestep error              ', abs(timestep1-timestep2))
        print('stage explicit update error ', num.linalg.norm(stage1.explicit_update-stage2.explicit_update))
        print('xmom  explicit update error ', num.linalg.norm(xmom1.explicit_update-xmom2.explicit_update))
        print('ymom  explicit update error ', num.linalg.norm(ymom1.explicit_update-ymom2.explicit_update))
        print('edge timestep error         ', num.linalg.norm(domain1.edge_timestep-domain2.edge_timestep))
        print('pressure work error         ', num.linalg.norm(domain1.pressuregrad_work-domain2.pressuregrad_work))
        print('edge flux work error        ', num.linalg.norm(domain1.edge_flux_work-domain2.edge_flux_work))


        assert num.allclose(timestep1,timestep2)
        assert num.allclose(stage1.explicit_update,stage2.explicit_update)
        assert num.allclose(xmom1.explicit_update,xmom2.explicit_update)
        assert num.allclose(ymom1.explicit_update,ymom2.explicit_update)
        assert num.allclose(domain1.edge_timestep,domain2.edge_timestep)
        assert num.allclose(domain1.pressuregrad_work,domain2.pressuregrad_work)
        assert num.allclose(domain1.edge_flux_work,domain2.edge_flux_work)

        import pprint

        #pprint.pprint(domain1.edge_timestep)
        #pprint.pprint(domain2.edge_timestep)    

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_DE_openmp, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
