#!/usr/bin/env python
#

# This file was reverted from changeset:5484 to changeset:5470 on 10th July 
# by Ole.

import unittest
import copy
import os
import numpy as num


import anuga.file.sww as sww
                

class Test_read_sww(unittest.TestCase):
    # Class variable
    verbose = False

    def set_verbose(self):
        Test_read_sww.verbose = True
        
    def setUp(self):
        pass

    def tearDown(self):
        for filename in ['read_sww_test0.sww', 'read_sww_test_c0.sww']:
            try:
                os.remove(filename)
            except:
                pass

    def test_read_sww(self):
        """
        Save to an sww file and then read back the info.
        Here we store the info "uniquely"
        """

        #---------------------------------------------------------------------
        # Import necessary modules
        #---------------------------------------------------------------------
        from anuga.abstract_2d_finite_volumes.mesh_factory import \
            rectangular_cross
        from anuga.shallow_water.shallow_water_domain import Domain
        from anuga import Reflective_boundary
        from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
                            import Dirichlet_boundary, Time_boundary

        #---------------------------------------------------------------------
        # Setup computational domain
        #---------------------------------------------------------------------
        length = 8.
        width = 4.
        dx = dy = 2   # Resolution: Length of subdivisions on both axes
        
        inc = 0.05 # Elevation increment

        points, vertices, boundary = rectangular_cross(int(length/dx), 
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        domain = Domain(points, vertices, boundary)
        domain.set_name('read_sww_test'+str(domain.processor))  # Output name
        domain.set_quantities_to_be_stored({'elevation': 2,
                                            'stage': 2,
                                            'xmomentum': 2,
                                            'ymomentum': 2,
                                            'friction': 1})

        domain.set_store_vertices_uniquely(True)
        
        #---------------------------------------------------------------------
        # Setup initial conditions
        #---------------------------------------------------------------------
        domain.set_quantity('elevation', 0.0)    # Flat bed initially
        domain.set_quantity('friction', 0.01)    # Constant friction
        domain.set_quantity('stage', 0.0)        # Dry initial condition

        #------------------------------------------------------------------
        # Setup boundary conditions
        #------------------------------------------------------------------
        Bi = Dirichlet_boundary([0.4, 0, 0])          # Inflow
        Br = Reflective_boundary(domain)              # Solid reflective wall
        Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow

        domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

        #-------------------------------------------------------------------
        # Evolve system through time
        #-------------------------------------------------------------------

        for t in domain.evolve(yieldstep=1, finaltime=4.0):
            pass
            
            
        # Check that quantities have been stored correctly    
        source = domain.get_name() + '.sww'


        #x = fid.variables['x'][:]
        #y = fid.variables['y'][:]
        #stage = fid.variables['stage'][:]
        #elevation = fid.variables['elevation'][:]
        #fid.close()
                   
        #assert len(stage.shape) == 2
        #assert len(elevation.shape) == 2        
        
        #M, N = stage.shape
                
        sww_file = sww.Read_sww(source)

        #print 'last frame number',sww_file.get_last_frame_number()

        assert num.allclose(sww_file.x, domain.get_vertex_coordinates()[:,0])
        assert num.allclose(sww_file.y, domain.get_vertex_coordinates()[:,1])

        
        assert num.allclose(sww_file.time, [0.0, 1.0, 2.0, 3.0, 4.0])
        
        M = domain.get_number_of_triangles()
        
        assert num.allclose(num.reshape(num.arange(3*M), (M,3)), sww_file.vertices)

        last_frame_number = sww_file.get_last_frame_number() 
        assert last_frame_number == 4

        assert num.allclose(sww_file.get_bounds(), [0.0, length, 0.0, width])

        assert 'stage'     in sww_file.quantities.keys()
        assert 'friction'  in sww_file.quantities.keys()
        assert 'elevation' in sww_file.quantities.keys()
        assert 'xmomentum' in sww_file.quantities.keys()
        assert 'ymomentum' in sww_file.quantities.keys()


        for qname, q in sww_file.read_quantities(last_frame_number).items():
            
            #print qname
            #print num.linalg.norm(num.abs((domain.get_quantity(qname).get_values()-q).flatten()), ord=1)
            
            assert num.allclose(domain.get_quantity(qname).get_values(), q)

        #-----------------------------------------
        # Start the evolution off again at frame 3
        #-----------------------------------------
        sww_file.read_quantities(last_frame_number-1)

        points, vertices, boundary = rectangular_cross(int(length/dx), 
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        new_domain = Domain(points, vertices, boundary)
        new_domain.set_quantities_to_be_stored(None)

        new_domain.set_store_vertices_uniquely(True)

        for qname, q in sww_file.read_quantities(last_frame_number-1).items():
            new_domain.set_quantity(qname, q)    

        #------------------------------------------------------------------
        # Setup boundary conditions
        #------------------------------------------------------------------
        Bi = Dirichlet_boundary([0.4, 0, 0])          # Inflow
        Br = Reflective_boundary(new_domain)          # Solid reflective wall
        Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow

        new_domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

        #-------------------------------------------------------------------
        # Evolve system through time
        #-------------------------------------------------------------------

        for t in new_domain.evolve(yieldstep=1.0, finaltime=1.0):
             pass

        # Compare  new_domain and domain quantities 
        for quantity in domain.get_quantity_names():
            dv = domain.get_quantity(quantity).get_values()
            ndv = new_domain.get_quantity(quantity).get_values()

            #print dv-ndv

            assert num.allclose( dv, ndv, rtol=5.e-2, atol=5.e-2)

        # Clean up
        #os.remove(source)
        
    def test_read_sww_with_centroids(self):
        """
        Save to an sww file and then read back the info.
        Here we store the info "uniquely"
        """

        #---------------------------------------------------------------------
        # Import necessary modules
        #---------------------------------------------------------------------
        from anuga.abstract_2d_finite_volumes.mesh_factory import \
            rectangular_cross
        from anuga.shallow_water.shallow_water_domain import Domain
        from anuga.shallow_water.boundaries import Reflective_boundary
        from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
                            import Dirichlet_boundary, Time_boundary

        #---------------------------------------------------------------------
        # Setup computational domain
        #---------------------------------------------------------------------
        length = 8.
        width = 4.
        dx = dy = 2   # Resolution: Length of subdivisions on both axes
        
        inc = 0.05 # Elevation increment

        points, vertices, boundary = rectangular_cross(int(length/dx), 
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        domain = Domain(points, vertices, boundary)
        domain.set_name('read_sww_test_c'+str(domain.processor))  # Output name
        domain.set_quantities_to_be_stored({'elevation': 2,
                                            'stage': 2,
                                            'xmomentum': 2,
                                            'ymomentum': 2,
                                            'friction': 1})

        domain.set_store_vertices_uniquely(True)
        domain.set_store_centroids(True)
        
        #---------------------------------------------------------------------
        # Setup initial conditions
        #---------------------------------------------------------------------
        domain.set_quantity('elevation', 0.0)    # Flat bed initially
        domain.set_quantity('friction', 0.01)    # Constant friction
        domain.set_quantity('stage', 0.0)        # Dry initial condition

        #------------------------------------------------------------------
        # Setup boundary conditions
        #------------------------------------------------------------------
        Bi = Dirichlet_boundary([0.4, 0, 0])          # Inflow
        Br = Reflective_boundary(domain)              # Solid reflective wall
        Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow

        domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

        #-------------------------------------------------------------------
        # Evolve system through time
        #-------------------------------------------------------------------

        for t in domain.evolve(yieldstep=1, finaltime=4.0):
            pass
            
            
        # Check that quantities have been stored correctly    
        source = domain.get_name() + '.sww'


        #x = fid.variables['x'][:]
        #y = fid.variables['y'][:]
        #stage = fid.variables['stage'][:]
        #elevation = fid.variables['elevation'][:]
        #fid.close()
                   
        #assert len(stage.shape) == 2
        #assert len(elevation.shape) == 2        
        
        #M, N = stage.shape
                
        sww_file = sww.Read_sww(source)

        #print 'last frame number',sww_file.get_last_frame_number()

        assert num.allclose(sww_file.x, domain.get_vertex_coordinates()[:,0])
        assert num.allclose(sww_file.y, domain.get_vertex_coordinates()[:,1])

        
        assert num.allclose(sww_file.time, [0.0, 1.0, 2.0, 3.0, 4.0])
        
        M = domain.get_number_of_triangles()
        
        assert num.allclose(num.reshape(num.arange(3*M), (M,3)), sww_file.vertices)

        last_frame_number = sww_file.get_last_frame_number() 
        assert last_frame_number == 4

        assert num.allclose(sww_file.get_bounds(), [0.0, length, 0.0, width])

        #print 50*"="
        #print sww_file.quantities.keys()

        assert 'stage'     in sww_file.quantities.keys()
        assert 'friction'  in sww_file.quantities.keys()
        assert 'elevation' in sww_file.quantities.keys()
        assert 'xmomentum' in sww_file.quantities.keys()
        assert 'ymomentum' in sww_file.quantities.keys()


        for qname, q in sww_file.read_quantities(last_frame_number).items():
            assert num.allclose(domain.get_quantity(qname).get_values(), q)

        #-----------------------------------------
        # Start the evolution off again at frame 3
        #-----------------------------------------
        sww_file.read_quantities(last_frame_number-1)

        points, vertices, boundary = rectangular_cross(int(length/dx), 
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        new_domain = Domain(points, vertices, boundary)
        new_domain.set_quantities_to_be_stored(None)

        new_domain.set_store_vertices_uniquely(True)

        for qname, q in sww_file.read_quantities(last_frame_number-1).items():
            new_domain.set_quantity(qname, q)    

        #------------------------------------------------------------------
        # Setup boundary conditions
        #------------------------------------------------------------------
        Bi = Dirichlet_boundary([0.4, 0, 0])          # Inflow
        Br = Reflective_boundary(new_domain)          # Solid reflective wall
        Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow

        new_domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

        #-------------------------------------------------------------------
        # Evolve system through time
        #-------------------------------------------------------------------

        for t in new_domain.evolve(yieldstep=1.0, finaltime=1.0):
             pass

        # Compare  new_domain and domain quantities 
        for quantity in domain.get_quantity_names():
            dv = domain.get_quantity(quantity).get_values()
            ndv = new_domain.get_quantity(quantity).get_values()

            #print dv-ndv

            assert num.allclose( dv, ndv, rtol=5.e-2, atol=5.e-2)

        # Clean up
        #os.remove(source)
        

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_read_sww, 'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)
    
