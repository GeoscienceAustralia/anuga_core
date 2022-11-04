
import os
import unittest
import tempfile
import numpy as num

from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.file.csv_file import load_csv_as_array, load_csv_as_dict
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.file.sww import weed, get_mesh_and_quantities_from_file, \
                Write_sww
from anuga.file.netcdf import NetCDFFile

from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_float, default_boundary_tag

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
        for filename in ['test_get_mesh_and_quantities_from_unique_vertices_sww_file.sww', \
                         'test_get_mesh_and_quantities_from_sww_file.sww']:
            try:
                os.remove(filename)
            except:
                pass

    def test_default_boundary(self):
        """Test that default boundary is correctly assigned
        """
        
        yiel = 0.01
        points, vertices, temp_boundary = rectangular(4, 4)
        
        # Deliberately remove tag 'right' to enforce the application of the default tag
        boundary = {}
        for key in temp_boundary:
            if temp_boundary[key] != 'right':
                boundary[key] = temp_boundary[key]
                
        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        
        found_default_boundary_tag = False
        for key in domain.boundary:
            
            if domain.boundary[key] == 'right':
                msg = 'Unexpected tag \'right\' found in domain.boundary'
                raise Exception(msg)
                
            if domain.boundary[key] == default_boundary_tag:
                found_default_boundary_tag = True
        
        msg = 'Did not find default boundary tag (%s) as expected' % default_boundary_tag
        assert found_default_boundary_tag, msg
                
                
        domain.geo_reference = Geo_reference(56,11,11)
        domain.smooth = False
        domain.store = True
        domain.set_name('default_bc')
        domain.default_order=2
        

        domain.set_quantity('elevation', lambda x,y: -x/3.0)
        domain.set_quantity('friction', 0.1)
        
        # Boundary conditions
        from math import sin, pi
        Br = Reflective_boundary(domain)
        Bt = Transmissive_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])
        Bw = Time_boundary(domain=domain,function=lambda t: [(0.1*sin(t*2*pi)), 0.0, 0.0])

        # Check that using the removed tag 'right' triggers an exception
        try:
            domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})
        except Exception as ex:
            # Check error message is correct
            assert 'Tag "right" provided does not exist' in str(ex) 
        else:
            msg = 'Invalid boundary tag should have failed.'        
            raise Exception(msg)
            
        # Now set boundary conditions appropriately
        domain.set_boundary({'left': Bd, default_boundary_tag: Bd, 'top': Bd, 'bottom': Bd})        

        # Check that it is as expected
        assert domain.boundary == {(0, 1): 'bottom', (1, 2): 'left', (3, 2): 'left', (5, 2): 'left', (7, 2): 'left', (7, 1): 'top', (8, 1): 'bottom', (15, 1): 'top', (16, 1): 'bottom', (23, 1): 'top', (24, 1): 'bottom', (31, 1): 'top', (24, 2): 'exterior', (26, 2): 'exterior', (28, 2): 'exterior', (30, 2): 'exterior'}
        

        # And just check that it runs for good measure
        
        domain.quantities_to_be_stored['xmomentum'] = 2
        domain.quantities_to_be_stored['ymomentum'] = 2
        # Initial condition
        h = 0.05
        elevation = domain.quantities['elevation'].vertex_values
        domain.set_quantity('stage', elevation + h)

        domain.check_integrity()
        for t in domain.evolve(yieldstep = yiel, finaltime = 0.05):
            #domain.print_timestepping_statistics()
            pass

        os.remove(domain.get_name() + '.sww') 


    def Xtest_sww2domain1(self):
    
        # FIXME (Ole): DELETE THIS TEST
    
        ################################################
        #Create a test domain, and evolve and save it.
        ################################################
        #from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh

        yiel = 0.01
        points, vertices, boundary = rectangular(10,10)
        #print('Boundary from rectangular', boundary):  {(0, 1): 'bottom', (1, 2): 'left', (3, 2): 'lef.....

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.geo_reference = Geo_reference(56,11,11)
        domain.smooth = False
        domain.store = True
        domain.set_name('sww2domain')
        domain.default_order=2
        #Bed-slope and friction
        domain.set_quantity('elevation', lambda x,y: -x/3.0)
        domain.set_quantity('friction', 0.1)
        # Boundary conditions
        from math import sin, pi
        Br = Reflective_boundary(domain)
        Bt = Transmissive_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])
        Bw = Time_boundary(domain=domain, function=lambda t: [(0.1*sin(t*2*pi)), 0.0, 0.0])

        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})

        domain.quantities_to_be_stored['xmomentum'] = 2
        domain.quantities_to_be_stored['ymomentum'] = 2
        # Initial condition
        h = 0.05
        elevation = domain.quantities['elevation'].vertex_values
        domain.set_quantity('stage', elevation + h)

        domain.check_integrity()
        # Evolution
        #domain.tight_slope_limiters = 1
        for t in domain.evolve(yieldstep = yiel, finaltime = 0.05):
            #domain.print_timestepping_statistics()
            pass

        #print boundary


        filename = domain.datadir + os.sep + domain.get_name() + '.sww'
        domain2 = load_sww_as_domain(filename, 
                                     #boundary=domain.boundary, 
                                     boundary=None,
                                     fail_if_NaN=False,
                                     verbose=self.verbose)


        # Unfortunately we lose the boundaries top, bottom, left and right,
        # they are now all lumped into "exterior"

        #print("=============== boundary domain2 =======================")
        #print(domain2.boundary)
        

        #print domain2.get_boundary_tags()
        
        #points, vertices, boundary = rectangular(15,15)
        #domain2.boundary = boundary
        ###################
        ##NOW TEST IT!!!
        ###################

        try:
            os.remove(filename)  # Clean up
        except:
            pass        

        bits = ['vertex_coordinates']
        for quantity in ['stage']:
            bits.append('get_quantity("%s").get_integral()' % quantity)
            bits.append('get_quantity("%s").get_values()' % quantity)

        for bit in bits:
            #print 'testing that domain.'+bit+' has been restored'
            #print bit
            #print 'done'
            #print eval('domain.'+bit)
            #print eval('domain2.'+bit)
            assert num.allclose(eval('domain.'+bit),eval('domain2.'+bit))

        ######################################
        #Now evolve them both, just to be sure
        ######################################x
        from time import sleep

        final = .1
        domain.set_quantity('friction', 0.1)
        domain.store = False
        
        domain.set_boundary({'left' : Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})

        for t in domain.evolve(yieldstep = yiel, finaltime = final):
            #domain.print_timestepping_statistics()
            pass


        domain2.smooth = False
        domain2.store = False
        domain2.default_order=2
        domain2.set_quantity('friction', 0.1)
        #Bed-slope and friction
        # Boundary conditions
        Bd2=Dirichlet_boundary([0.2,0.,0.])
        domain2.boundary = domain.boundary
        #print('domain2.boundary')
        #print(domain2.boundary)

        #print('------------------------------------')
        #print('Third set_boundary')
        #print('------------------------------------')                
                
        domain2.set_boundary({'left' : Bd,  'right': Bd, 'top': Bd, 'bottom': Bd})        
        #domain2.set_boundary({'exterior' : Bd})

        #print()
        #print()
        #print(domain2.boundary_map)


        domain2.check_integrity()
        
        for t in domain2.evolve(yieldstep = yiel, finaltime = final):
            #domain2.print_timestepping_statistics()
            pass

        ###################
        ##NOW TEST IT!!!
        ##################

        bits = ['vertex_coordinates']

        for quantity in ['elevation','stage', 'ymomentum','xmomentum']:
            bits.append('get_quantity("%s").get_integral()' % quantity)
            bits.append('get_quantity("%s").get_values()' % quantity)

        #print bits
        for bit in bits:
            #print bit
            #print eval('domain.'+bit)
            #print eval('domain2.'+bit)
            
            msg = 'Values in the two domains are different for ' + bit
            assert num.allclose(eval('domain.' + bit), eval('domain2.' + bit),
                                rtol=5.e-2, atol=5.e-2), msg


    def test_sww2domain_starttime(self):
        """Test that domain start time is stored correctly in the sww file 
        """

        verbose=False
        starttime = 200.0
             
        points, vertices, boundary = rectangular(10,10)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.geo_reference = Geo_reference(56,11,11)
        domain.smooth = False
        domain.store = True
        domain.set_name('sww2domain_starttime')
        domain.default_order=2
        domain.set_starttime(starttime)

        # Bed-slope and friction
        domain.set_quantity('elevation', lambda x,y: -x/3)
        domain.set_quantity('friction', 0.1)
        
        # Boundary conditions
        from math import sin, pi
        Br = Reflective_boundary(domain)
        Bt = Transmissive_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])
        Bw = Time_boundary(domain=domain,function=lambda t: [(0.1*sin(t*2*pi)), 0.0, 0.0])

        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})

        domain.quantities_to_be_stored['xmomentum'] = 2
        domain.quantities_to_be_stored['ymomentum'] = 2
        
        # Initial conditions
        h = 0.05
        elevation = domain.quantities['elevation'].vertex_values
        domain.set_quantity('stage', elevation + h)

        domain.check_integrity()

        if verbose:
            print('evolve domain')
        
        for t in domain.evolve(yieldstep = 0.01, duration = 0.05):
            #domain.print_timestepping_statistics()
            pass

        # Test that start time was stored correctly
        if verbose:
            print('read in domain')
            
        filename = domain.datadir + os.sep + domain.get_name() + '.sww'
        fid = NetCDFFile(filename, netcdf_mode_r)  # Open sww file for read        

        stored_starttime = float(fid.starttime)
        assert stored_starttime == starttime

        try:
            os.remove(filename)  # Clean up
        except:
            pass


    def test_get_mesh_and_quantities_from_1_5_sww_file(self):
        """test_get_mesh_and_quantities_from_sww_file(self):
        """     
        
        # Generate a test sww file with non trivial georeference
        
        import time, os

        # Setup
        #from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh (100m x 5m)
        width = 5
        length = 50
        t_end = 10
        points, vertices, boundary = rectangular(length, width, 50, 5)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference = Geo_reference(56,308500,6189000))

        domain.set_name('test_get_mesh_and_quantities_from_sww_file')
        swwfile = domain.get_name() + '.sww'
        domain.set_datadir('.')
        domain.set_flow_algorithm('1_5')

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
        assert num.allclose(time, list(range(t_end+1)))

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
            assert num.allclose(q_ref, q_sww, atol=1.0e-10), msg
            
        # Cleanup
        #os.remove(swwfile)
        
    def test_get_mesh_and_quantities_from_de0_sww_file(self):
        """test_get_mesh_and_quantities_from_sww_file(self):
        """     
        
        # Generate a test sww file with non trivial georeference
        
        import time, os

        # Setup
        #from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh (100m x 5m)
        width = 5
        length = 50
        t_end = 10
        points, vertices, boundary = rectangular(length, width, 50, 5)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference = Geo_reference(56,308500,6189000))

        domain.set_name('test_get_mesh_and_quantities_from_sww_file')
        swwfile = domain.get_name() + '.sww'
        domain.set_datadir('.')
        domain.set_flow_algorithm('DE0')

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
        assert num.allclose(time, list(range(t_end+1)))

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
            assert num.allclose(q_ref, q_sww, atol=1.0e-2), msg
            
        # Cleanup
        #os.remove(swwfile)   
    
    def test_get_mesh_and_quantities_from_unique_vertices_1_5_sww_file(self):
        """test_get_mesh_and_quantities_from_unique_vertices_sww_file(self):
        """     
        
        # Generate a test sww file with non trivial georeference
        
        import time, os

        # Setup
        #from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh (100m x 5m)
        width = 5
        length = 50
        t_end = 10
        points, vertices, boundary = rectangular(10, 1, length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference = Geo_reference(56,308500,6189000))

        domain.set_name('test_get_mesh_and_quantities_from_unique_vertices_sww_file')
        swwfile = domain.get_name() + '.sww'
        domain.set_datadir('.')
        domain.set_flow_algorithm('1_5')
        domain.set_store_vertices_uniquely()

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
    

        #print quantities
        #print time

        dhash = domain.get_nodes()[:,0]*10+domain.get_nodes()[:,1]
        mhash = mesh.nodes[:,0]*10+mesh.nodes[:,1]        


        #print 'd_nodes',len(dhash)
        #print 'm_nodes',len(mhash)
        di = num.argsort(dhash)
        mi = num.argsort(mhash)
        minv = num.argsort(mi)
        dinv = num.argsort(di)

        #print 'd_tri',len(domain.get_triangles())
        #print 'm_tri',len(mesh.triangles)
        
        # Check that mesh has been recovered
        # triangle order should be ok
        assert num.allclose(mesh.nodes[mi,:],domain.get_nodes()[di,:])
        assert num.alltrue(minv[mesh.triangles] == dinv[domain.get_triangles()])


        # Check that time has been recovered
        assert num.allclose(time, list(range(t_end+1)))

        z=domain.get_quantity('elevation').get_values(location='vertices').flatten()
        

        assert num.allclose(quantities['elevation'], z)

        for q in ['stage', 'xmomentum', 'ymomentum']:
            # Get quantity at last timestep
            q_ref=domain.get_quantity(q).get_values(location='vertices').flatten()

            #print q,quantities[q]
            q_sww=quantities[q][-1,:]
            
            
            msg = 'Quantity %s failed to be recovered' %q
            assert num.allclose(q_ref, q_sww, atol=1.0e-6), msg
            
        # Cleanup
        #os.remove(swwfile)
        
        
    def test_get_mesh_and_quantities_from_unique_vertices_DE0_sww_file(self):
        """test_get_mesh_and_quantities_from_unique_vertices_sww_file(self):
        """     
        
        # Generate a test sww file with non trivial georeference
        
        import time, os

        # Setup
        #from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh (100m x 5m)
        width = 5
        length = 50
        t_end = 10
        points, vertices, boundary = rectangular(10, 1, length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference = Geo_reference(56,308500,6189000))

        domain.set_name('test_get_mesh_and_quantities_from_unique_vertices_sww_file')
        swwfile = domain.get_name() + '.sww'
        domain.set_datadir('.')
        domain.set_flow_algorithm('DE0')
        domain.set_store_vertices_uniquely()

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

        #print quantities
        #print time

        dhash = domain.get_nodes()[:,0]*10+domain.get_nodes()[:,1]
        mhash = mesh.nodes[:,0]*10+mesh.nodes[:,1]        


        #print 'd_nodes',len(dhash)
        #print 'm_nodes',len(mhash)
        di = num.argsort(dhash)
        mi = num.argsort(mhash)
        minv = num.argsort(mi)
        dinv = num.argsort(di)

        #print 'd_tri',len(domain.get_triangles())
        #print 'm_tri',len(mesh.triangles)
        
        # Check that mesh has been recovered
        # triangle order should be ok
        assert num.allclose(mesh.nodes[mi,:],domain.get_nodes()[di,:])
        assert num.alltrue(minv[mesh.triangles] == dinv[domain.get_triangles()])


        # Check that time has been recovered
        assert num.allclose(time, list(range(t_end+1)))

        z=domain.get_quantity('elevation').get_values(location='vertices').flatten()
        

        
        assert num.allclose(quantities['elevation'], z)

        for q in ['stage', 'xmomentum', 'ymomentum']:
            # Get quantity at last timestep
            q_ref=domain.get_quantity(q).get_values(location='vertices').flatten()

            #print q,quantities[q]
            q_sww=quantities[q][-1,:]
            
            msg = 'Quantity %s failed to be recovered' %q
            assert num.allclose(q_ref, q_sww, atol=1.0e-6), msg
            
        # Cleanup
        #os.remove(swwfile)
        
    def test_weed(self):
        coordinates1 = [[0.,0.],[1.,0.],[1.,1.],[1.,0.],[2.,0.],[1.,1.]]
        volumes1 = [[0,1,2],[3,4,5]]
        boundary1= {(0,1): 'external',(1,2): 'not external',(2,0): 'external',(3,4): 'external',(4,5): 'external',(5,3): 'not external'}
        coordinates2,volumes2,boundary2=weed(coordinates1,volumes1,boundary1)

        points2 = {(0.,0.):None,(1.,0.):None,(1.,1.):None,(2.,0.):None}

        assert len(points2)==len(coordinates2)
        for i in range(len(coordinates2)):
            coordinate = tuple(coordinates2[i])
            assert coordinate in points2
            points2[coordinate]=i

        for triangle in volumes1:
            for coordinate in triangle:
                assert coordinates2[points2[tuple(coordinates1[coordinate])]][0]==coordinates1[coordinate][0]
                assert coordinates2[points2[tuple(coordinates1[coordinate])]][1]==coordinates1[coordinate][1]


    def test_triangulation(self):
        # 
        #  
        
        filename = tempfile.mktemp("_data_manager.sww")
        outfile = NetCDFFile(filename, netcdf_mode_w)
        points_utm = num.array([[0.,0.],[1.,1.],[0.,1.]])
        volumes = [[0,1,2]]
        elevation = [0,1,2]
        new_origin = None
        new_origin = Geo_reference(56, 0, 0)
        times = [0, 10]
        number_of_volumes = len(volumes)
        number_of_points = len(points_utm)
        sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])
        sww.store_header(outfile, times, number_of_volumes,
                         number_of_points, description='fully sick testing',
                         verbose=self.verbose,sww_precision=netcdf_float)
        sww.store_triangulation(outfile, points_utm, volumes,
                                elevation,  new_origin=new_origin,
                                verbose=self.verbose)       
        outfile.close()
        fid = NetCDFFile(filename)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        fid.close()

        assert num.allclose(num.array(list(zip(x,y))), points_utm)

        try:
            os.remove(filename)  # Clean up
        except:
            pass        


        
    def test_triangulationII(self):
        # 
        #  

        filename = tempfile.mktemp("_data_manager.sww")
        outfile = NetCDFFile(filename, netcdf_mode_w)
        points_utm = num.array([[0.,0.],[1.,1.], [0.,1.]])
        volumes = [[0,1,2]]
        elevation = [0,1,2]
        new_origin = None
        #new_origin = Geo_reference(56, 0, 0)
        times = [0, 10]
        number_of_volumes = len(volumes)
        number_of_points = len(points_utm)
        sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])        
        sww.store_header(outfile, times, number_of_volumes,
                         number_of_points, description='fully sick testing',
                         verbose=self.verbose,sww_precision=netcdf_float)
        sww.store_triangulation(outfile, points_utm, volumes,
                                new_origin=new_origin,
                                verbose=self.verbose)
        sww.store_static_quantities(outfile, elevation=elevation)                                
                                
        outfile.close()
        fid = NetCDFFile(filename)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        results_georef = Geo_reference()
        results_georef.read_NetCDF(fid)
        
        assert results_georef == Geo_reference(zone=None, xllcorner=0, yllcorner=0)
        fid.close()

        assert num.allclose(num.array(list(zip(x,y))), points_utm)


        try:
            os.remove(filename)  # Clean up
        except:
            pass        

        
    def test_triangulation_new_origin(self):
        # 
        #  
        
        filename = tempfile.mktemp('_data_manager.sww')
        outfile = NetCDFFile(filename, netcdf_mode_w)
        points_utm = num.array([[0.,0.], [1.,1.], [0.,1.]])
        volumes = [[0,1,2]]
        elevation = [0,1,2]
        new_origin = None
        new_origin = Geo_reference(56, 1, 554354)
        points_utm = new_origin.change_points_geo_ref(points_utm)
        times = [0, 10]
        number_of_volumes = len(volumes)
        number_of_points = len(points_utm)
        sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])        
        sww.store_header(outfile, times, number_of_volumes,
                         number_of_points, description='fully sick testing',
                         verbose=self.verbose,sww_precision=netcdf_float)
        sww.store_triangulation(outfile, points_utm, volumes,
                                elevation,  new_origin=new_origin,
                                verbose=self.verbose)
        outfile.close()
        fid = NetCDFFile(filename)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        results_georef = Geo_reference()
        results_georef.read_NetCDF(fid)
        assert results_georef == new_origin
        fid.close()

        absolute = Geo_reference(56, 0, 0)
        assert num.allclose(num.array(
            absolute.change_points_geo_ref(list(zip(x,y)),
                                           new_origin)),points_utm)
        
        try:
            os.remove(filename)  # Clean up
        except:
            pass
        
    def test_triangulation_points_georeference(self):
        # 
        #  
        
        filename = tempfile.mktemp("_data_manager.sww")
        outfile = NetCDFFile(filename, netcdf_mode_w)
        points_utm = num.array([[0.,0.],[1.,1.], [0.,1.]])
        volumes = [[0,1,2]]
        elevation = [0,1,2]
        new_origin = None
        points_georeference = Geo_reference(56, 1, 554354)
        points_utm = points_georeference.change_points_geo_ref(points_utm)
        times = [0, 10]
        number_of_volumes = len(volumes)
        number_of_points = len(points_utm)
        sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])        
        sww.store_header(outfile, times, number_of_volumes,
                         number_of_points, description='fully sick testing',
                         verbose=self.verbose,sww_precision=netcdf_float)
        sww.store_triangulation(outfile, points_utm, volumes,
                                elevation,  new_origin=new_origin,
                                points_georeference=points_georeference,
                                verbose=self.verbose)       
        outfile.close()
        fid = NetCDFFile(filename)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        results_georef = Geo_reference()
        results_georef.read_NetCDF(fid)
        assert results_georef == points_georeference
        fid.close()

        assert num.allclose(num.array(list(zip(x,y))), points_utm)

        try:
            os.remove(filename)  # Clean up
        except:
            pass
        
    def test_triangulation_2_geo_refs(self):
        # 
        #  
        
        filename = tempfile.mktemp("_data_manager.sww")
        outfile = NetCDFFile(filename, netcdf_mode_w)
        points_utm = num.array([[0.,0.],[1.,1.], [0.,1.]])
        volumes = [[0,1,2]]
        elevation = [0,1,2]
        new_origin = Geo_reference(56, 1, 1)
        points_georeference = Geo_reference(56, 0, 0)
        points_utm = points_georeference.change_points_geo_ref(points_utm)
        times = [0, 10]
        number_of_volumes = len(volumes)
        number_of_points = len(points_utm)
        sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])        
        sww.store_header(outfile, times, number_of_volumes,
                         number_of_points, description='fully sick testing',
                         verbose=self.verbose,sww_precision=netcdf_float)
        sww.store_triangulation(outfile, points_utm, volumes,
                                elevation,  new_origin=new_origin,
                                points_georeference=points_georeference,
                                verbose=self.verbose)       
        outfile.close()
        fid = NetCDFFile(filename)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        results_georef = Geo_reference()
        results_georef.read_NetCDF(fid)
        assert results_georef == new_origin
        fid.close()


        absolute = Geo_reference(56, 0, 0)
        assert num.allclose(num.array(
            absolute.change_points_geo_ref(list(zip(x,y)),
                                           new_origin)),points_utm)

        try:
            os.remove(filename)  # Clean up
        except:
            pass


#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_sww, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
