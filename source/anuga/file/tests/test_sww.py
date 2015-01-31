import os
import unittest
import tempfile
import numpy as num

from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.file.csv_file import load_csv_as_array, load_csv_as_dict
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.file.sww import load_sww_as_domain, weed, get_mesh_and_quantities_from_file, \
                Write_sww
from anuga.file.netcdf import NetCDFFile

from anuga.config import netcdf_mode_w, netcdf_float

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
        
    def test_sww2domain1(self):
        ################################################
        #Create a test domain, and evolve and save it.
        ################################################
        #from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh

        yiel=0.01
        points, vertices, boundary = rectangular(10,10)

        #print "=============== boundary rect ======================="
        #print boundary

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
        Bw = Time_boundary(domain=domain,function=lambda t: [(0.1*sin(t*2*pi)), 0.0, 0.0])

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

        #print boundary


        filename = domain.datadir + os.sep + domain.get_name() + '.sww'
        domain2 = load_sww_as_domain(filename, None, fail_if_NaN=False,
                                        verbose=self.verbose)

        # Unfortunately we loss the boundaries top, bottom, left and right,
        # they are now all lumped into "exterior"

        #print "=============== boundary domain2 ======================="
        #print domain2.boundary
        

        #print domain2.get_boundary_tags()
        
        #points, vertices, boundary = rectangular(15,15)
        #domain2.boundary = boundary
        ###################
        ##NOW TEST IT!!!
        ###################

        os.remove(filename)

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
        domain.set_boundary({'exterior': Bd, 'left' : Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})


        for t in domain.evolve(yieldstep = yiel, finaltime = final):
            #domain.write_time()
            pass

        #BUT since domain1 gets time hacked back to 0:
        
        final = final + (domain2.get_starttime() - domain.get_starttime())

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
        domain2.set_boundary({'exterior': Bd, 'left' : Bd,  'right': Bd, 'top': Bd, 'bottom': Bd})
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
                                rtol=5.e-2, atol=5.e-2), msg



    def test_get_mesh_and_quantities_from_sww_file(self):
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
        #os.remove(swwfile)
        
    def test_get_mesh_and_quantities_from_unique_vertices_sww_file(self):
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
        assert num.allclose(time, range(t_end+1))

        z=domain.get_quantity('elevation').get_values(location='unique vertices')
        
        assert num.allclose(quantities['elevation'], z)

        for q in ['stage', 'xmomentum', 'ymomentum']:
            # Get quantity at last timestep
            q_ref=domain.get_quantity(q).get_values(location='unique vertices')

            #print q,quantities[q]
            q_sww=quantities[q][-1,:]
            
            msg = 'Quantity %s failed to be recovered' %q
            assert num.allclose(q_ref[di], q_sww[mi], atol=1.0e-6), msg
            
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
            assert points2.has_key(coordinate)
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

        assert num.allclose(num.array(map(None, x,y)), points_utm)
        os.remove(filename)

        
    def test_triangulationII(self):
        # 
        #  

        DEFAULT_ZONE = 0 # Not documented anywhere what this should be.
        
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
        assert results_georef == Geo_reference(DEFAULT_ZONE, 0, 0)
        fid.close()

        assert num.allclose(num.array(map(None, x,y)), points_utm)
        os.remove(filename)

        
    def test_triangulation_new_origin(self):
        # 
        #  
        
        filename = tempfile.mktemp("_data_manager.sww")
        outfile = NetCDFFile(filename, netcdf_mode_w)
        points_utm = num.array([[0.,0.],[1.,1.], [0.,1.]])
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

        absolute = Geo_reference(56, 0,0)
        assert num.allclose(num.array(
            absolute.change_points_geo_ref(map(None, x,y),
                                           new_origin)),points_utm)
        
        os.remove(filename)
        
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

        assert num.allclose(num.array(map(None, x,y)), points_utm)
        os.remove(filename)
        
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


        absolute = Geo_reference(56, 0,0)
        assert num.allclose(num.array(
            absolute.change_points_geo_ref(map(None, x,y),
                                           new_origin)),points_utm)
        os.remove(filename)

#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_sww, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
