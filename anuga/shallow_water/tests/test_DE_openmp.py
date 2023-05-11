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
            domain = anuga.rectangular_cross_domain(20,20, len1=1., len2=1.)

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

        domain0 = create_domain('domain_original')
        domain0.openmp_code=0

        domain1 = create_domain('domain_parallel_loop')
        domain1.openmp_code=1

        domain2 = create_domain('domain_parallel_loop')
        domain2.openmp_code=1 # will change to 2 once burn in

        #------------------------------
        #Evolve the system through time
        #------------------------------

        print('Evolve domain0')
        for t in domain0.evolve(yieldstep=0.1,finaltime=0.1):
            domain0.print_timestepping_statistics()

        print('Evolve domain1')
        for t in domain1.evolve(yieldstep=0.1,finaltime=0.1):
            domain1.print_timestepping_statistics()

        print('Evolve domain2')
        for t in domain2.evolve(yieldstep=0.1,finaltime=0.1):
            domain2.print_timestepping_statistics()

        #----------------------------------------
        # Now just run the openmp code on domain2
        #----------------------------------------
        domain2.openmp_code = 2
        timestep = 0.1

        from anuga.shallow_water.swDE_domain_ext import extrapolate_second_order_edge_sw
        from anuga.shallow_water.swDE_domain_ext import compute_fluxes_ext_central

        extrapolate_second_order_edge_sw(domain1)
        timestep1 = compute_fluxes_ext_central(domain1, timestep)

        extrapolate_second_order_edge_sw(domain2)
        timestep2 = compute_fluxes_ext_central(domain2, timestep)

        # Compare update arrays and timestep


        print(timestep1)
        print(timestep2)

        quantities1 = domain1.quantities
        stage1 = quantities1["stage"]
        xmom1 = quantities1["xmomentum"]
        ymom1 = quantities1["ymomentum"]


        quantities2 = domain2.quantities
        stage2 = quantities2["stage"]
        xmom2 = quantities2["xmomentum"]
        ymom2 = quantities2["ymomentum"]

        print(num.allclose(stage1.explicit_update, stage2.explicit_update))
        print(num.allclose(xmom1.explicit_update, xmom2.explicit_update))
        print(num.allclose(ymom1.explicit_update, ymom2.explicit_update))

        print(num.allclose(domain1.edge_timestep, domain2.edge_timestep))

        print(num.allclose(domain1.pressuregrad_work, domain2.pressuregrad_work))


        print(num.allclose(domain1.edge_timestep, domain2.edge_timestep))

        import pprint

        pprint.pprint(domain1.edge_timestep)
        pprint.pprint(domain2.edge_timestep)



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
                                           maximum_triangle_area = 0.2,
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

        domain0 = create_domain('domain_original')
        domain0.openmp_code=0

        domain1 = create_domain('domain_parallel_loop')
        domain1.openmp_code=1

        domain2 = create_domain('domain_parallel_loop')
        domain2.openmp_code=1 # will change to 2 once burn in

        #------------------------------
        #Evolve the system through time
        #------------------------------

        print('Evolve domain0')
        for t in domain0.evolve(yieldstep=0.1,finaltime=0.1):
            domain0.print_timestepping_statistics()

        print('Evolve domain1')
        for t in domain1.evolve(yieldstep=0.1,finaltime=0.1):
            domain1.print_timestepping_statistics()

        print('Evolve domain2')
        for t in domain2.evolve(yieldstep=0.1,finaltime=0.1):
            domain2.print_timestepping_statistics()

        #----------------------------------------
        # Now just run the openmp code on domain2
        #----------------------------------------
        domain2.openmp_code = 2
        timestep = 0.1

        from anuga.shallow_water.swDE_domain_ext import extrapolate_second_order_edge_sw
        from anuga.shallow_water.swDE_domain_ext import compute_fluxes_ext_central

        extrapolate_second_order_edge_sw(domain1)
        timestep1 = compute_fluxes_ext_central(domain1, timestep)

        extrapolate_second_order_edge_sw(domain2)
        timestep2 = compute_fluxes_ext_central(domain2, timestep)

        # Compare update arrays and timestep


        print(timestep1)
        print(timestep2)

        quantities1 = domain1.quantities
        stage1 = quantities1["stage"]
        xmom1 = quantities1["xmomentum"]
        ymom1 = quantities1["ymomentum"]


        quantities2 = domain2.quantities
        stage2 = quantities2["stage"]
        xmom2 = quantities2["xmomentum"]
        ymom2 = quantities2["ymomentum"]

        print(num.allclose(stage1.explicit_update, stage2.explicit_update))
        print(num.allclose(xmom1.explicit_update, xmom2.explicit_update))
        print(num.allclose(ymom1.explicit_update, ymom2.explicit_update))

        print(num.allclose(domain1.edge_timestep, domain2.edge_timestep))

        print(num.allclose(domain1.pressuregrad_work, domain2.pressuregrad_work))


        print(num.allclose(domain1.edge_timestep, domain2.edge_timestep))

        import pprint

        pprint.pprint(domain1.edge_timestep)
        pprint.pprint(domain2.edge_timestep)    

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_DE_openmp, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
