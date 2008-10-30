#!/usr/bin/env python
#
"""
Let's do some system wide tests
"""
import tempfile
import unittest
import os

from Scientific.IO.NetCDF import NetCDFFile
from Numeric import allclose

from anuga.shallow_water import Domain
from anuga.shallow_water import Dirichlet_boundary, Time_boundary
from anuga.shallow_water import File_boundary
from anuga.pmesh.mesh import Mesh
from anuga.abstract_2d_finite_volumes.pmesh2domain import pmesh_instance_to_domain_instance



class Test_system(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def create_sww_boundary(self, boundary_starttime):
        """
        This creates a boundary file with ;
        time     stage
        0        5
        10       2.15268
        20       13.9773
        """
        tide = 5
        boundary_filename = tempfile.mktemp(".sww")
        dir, base = os.path.split(boundary_filename)
        boundary_name = base[:-4]
        
        # Setup computational domain
        mesh = Mesh()
        mesh.add_region_from_polygon([[0,0], [100,0], [100,100], [0,100]])
        mesh.generate_mesh(verbose=False)
        
        domain = pmesh_instance_to_domain_instance(mesh, Domain)
        domain.set_name(boundary_name)                 
        domain.set_datadir(dir)          
        domain.set_starttime(boundary_starttime)
        
        # Setup initial conditions
        domain.set_quantity('elevation', 0.0) 
        domain.set_quantity('stage', tide)         

        # Setup boundary conditions
        Bd = Dirichlet_boundary([tide,0.,0.]) # Constant boundary values
        Bd = Time_boundary(domain=domain,     # Time dependent boundary  
                   f=lambda t: [t, 0.0, 0.0])
        domain.set_boundary({'exterior': Bd})
        for t in domain.evolve(yieldstep = 10, finaltime = 20.0):
            pass
            #print domain.boundary_statistics('stage')
            q = Bd.evaluate()
    
            # FIXME (Ole): This test would not have passed in 
            # changeset:5846.
            msg = 'Time boundary not evaluated correctly'
            assert allclose(t, q[0]), msg
            
            #print domain.get_quantity('stage').get_values()
            #domain.write_time()
            #print "domain.time", domain.time
            
        return boundary_filename
    
    def test_boundary_time(self):
        """
        test_boundary_time(self):
        test that the starttime of a boundary condition is carried thru
        to the output sww file.
        
        """
     
        boundary_starttime = 500
        boundary_filename = self.create_sww_boundary(boundary_starttime)
        filename = tempfile.mktemp(".sww")
        #print "filename",filename 
        dir, base = os.path.split(filename)
        senario_name = base[:-4]
 
        mesh = Mesh()
        ###mesh.add_region_from_polygon([[10,10], [90,10], [90,90], [10,90]])
        mesh.add_region_from_polygon([[0,0], [100,0], [100,100], [0,100]])
        mesh.generate_mesh(verbose=False)
        
        domain = pmesh_instance_to_domain_instance(mesh, Domain) 
        domain.set_name(senario_name)                 
        domain.set_datadir(dir) 

        # Setup initial conditions
        domain.set_quantity('elevation', 0.0) 
        domain.set_quantity('stage', 0.0)         
        Bf = File_boundary(boundary_filename,
                           domain,  use_cache=False, verbose=False)

        # Setup boundary conditions
        domain.set_boundary({'exterior': Bf})

        
        for t in domain.evolve(yieldstep = 5.0, finaltime = 10.0):
            pass
            #print domain.write_time()
            #print "domain.time", domain.time

        # do an assertion on the time of the produced sww file
        fid = NetCDFFile(filename, 'r')    #Open existing file for read
        times = fid.variables['time'][:]
        #print "times", times
        #print "fid.starttime", fid.starttime
        assert allclose(fid.starttime, boundary_starttime)
        fid.close()

        # clean up
        os.remove(boundary_filename)
        os.remove(filename)
        
    def test_boundary_timeII(self):
        """
        test_boundary_timeII(self):
        Test that starttime can be set in the middle of a boundary condition
        """
        
        boundary_starttime = 500
        boundary_filename = self.create_sww_boundary(boundary_starttime)
        #print "boundary_filename",boundary_filename 
        
        filename = tempfile.mktemp(".sww")
        #print "filename",filename 
        dir, base = os.path.split(filename)
        senario_name = base[:-4]
 
        mesh = Mesh()
        mesh.add_region_from_polygon([[10,10], [90,10], [90,90], [10,90]])
        mesh.generate_mesh(verbose=False)
        
        domain = pmesh_instance_to_domain_instance(mesh, Domain) 
        domain.set_name(senario_name)                 
        domain.set_datadir(dir)
        new_starttime = 510.
        domain.set_starttime(new_starttime)

        # Setup initial conditions
        domain.set_quantity('elevation', 0.0) 
        domain.set_quantity('stage', 0.0)         
        Bf = File_boundary(boundary_filename,
                           domain,  use_cache=False, verbose=False)

        # Setup boundary conditions
        domain.set_boundary({'exterior': Bf})
        for t in domain.evolve(yieldstep = 5, finaltime = 9.0):
            pass
            #print domain.boundary_statistics()
            #domain.write_time()
            #print "domain.time", domain.time

        # do an assertion on the time of the produced sww file
        fid = NetCDFFile(filename, 'r')    #Open existing file for read
        times = fid.variables['time'][:]
        stage = fid.variables['stage'][:]
        #print stage
        #print "times", times
        #print "fid.starttime", fid.starttime
        assert allclose(fid.starttime, new_starttime)
        fid.close()
        
        #print "stage[2,0]", stage[2,0]
        msg = "This test is a bit hand crafted, based on the output file. "
        msg += "Not logic. "
        msg += "It's testing that starttime is working"
        assert allclose(stage[2,0], 11.9867153168),msg
        
        

        # clean up
        os.remove(boundary_filename)
        os.remove(filename)
        
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_system,'test')
    #suite = unittest.makeSuite(Test_system,'test_boundary_timeII')
    runner = unittest.TextTestRunner()
    runner.run(suite)
