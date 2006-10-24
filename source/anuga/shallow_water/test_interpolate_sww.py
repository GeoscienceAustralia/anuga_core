#!/usr/bin/env python
#
"""
Testing interpolate_sww, based on test_data_manageer, so there maybe code
that isn't needed, eg in the setup file
"""

import unittest
from Numeric import zeros, array, allclose, Float
from anuga.utilities.numerical_tools import mean

from interpolate_sww import *
from anuga.shallow_water import Domain, Transmissive_boundary
from anuga.shallow_water.data_manager import get_dataobject


class Test_Interpolate_sww(unittest.TestCase):
    def setUp(self):

        import time
        from mesh_factory import rectangular


        #Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order=2
        domain.beta_h = 0


        #Set some field values
        domain.set_quantity('elevation', lambda x,y: -x)
        domain.set_quantity('friction', 0.03)


        ######################
        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        ######################
        #Initial condition - with jumps

        bed = domain.quantities['elevation'].vertex_values
        stage = zeros(bed.shape, Float)

        h = 0.3
        for i in range(stage.shape[0]):
            if i % 2 == 0:
                stage[i,:] = bed[i,:] + h
            else:
                stage[i,:] = bed[i,:]

        domain.set_quantity('stage', stage)

        domain.distribute_to_vertices_and_edges()


        self.domain = domain

        C = domain.get_vertex_coordinates()
        self.X = C[:,0:6:2].copy()
        self.Y = C[:,1:6:2].copy()

        self.F = bed


    def tearDown(self):
        pass

    def test_sww_DSG(self):
        """Not a test, rather a look at the sww format
        """

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate
        from Scientific.IO.NetCDF import NetCDFFile

        #FIXME (Ole): This test still passes if commented out
        #self.domain.set_name('datatest' + str(time.time()))
        
        
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = mean

        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')
        self.domain.time = 2.
        sww.store_timestep('stage')

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(sww.filename, 'r')

        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']

        volumes = fid.variables['volumes']
        time = fid.variables['time']

        # 2D
        stage = fid.variables['stage']

        X = x[:]
        Y = y[:]
        Z = z[:]
        V = volumes[:]
        T = time[:]
        S = stage[:,:]

        if False:
            print "****************************"
            print "X ",X
            print "****************************"
            print "Y ",Y
            print "****************************"
            print "Z ",Z
            print "****************************"
            print "V ",V
            print "****************************"
            print "Time ",T
            print "****************************"
            print "Stage ",S
            print "****************************"




        fid.close()

        #Cleanup
        os.remove(sww.filename)

    def test_interpolate_sww(self):
        ### Not really a unit test, rather a system test for
        ###  

        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate, \
             transpose
        from Scientific.IO.NetCDF import NetCDFFile
        import tempfile
        from load_mesh.loadASCII import  import_points_file, \
             concatinate_attributelist

        self.domain.set_name('datatest' + str(time.time()))
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = mean

        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')
        self.domain.time = 2.
        sww.store_timestep('stage')

        interp = Interpolate_sww(sww.filename, 'depth')

        assert allclose(interp.time,[0.0,2.0])

        # create an .xya file
        point_file = tempfile.mktemp(".xya")
        fd = open(point_file,'w')
        fd.write("colour, thickness \n 0.0, 0.6,2.,4 \n 0.0, 0.9,4,8 \n 0.0,0.1,4.,8 \n 0.4,1.0,4.,8 \n")
        fd.close()

        interp.interpolate_xya(point_file)

        answer = {}
        answer['0.0'] =  [ 0.08, 0.02, 0.14,0.08]
        answer['2.0'] =  [ 0.08, 0.02, 0.14,0.08]
        
        #print "answer",answer
        #print "interp.interpolated_quantity",interp.interpolated_quantity
        assert allclose(interp.interpolated_quantity['0.0'],answer['0.0'])
        assert allclose(interp.interpolated_quantity['2.0'],answer['2.0'])

        # create an output .xya file
        point_file_out = tempfile.mktemp(".xya")
        interp.write_depth_xya(point_file_out)

        #check the output file
        xya_dict = import_points_file(point_file_out)
        assert allclose(xya_dict['attributelist']['0.0'],answer['0.0'])
        assert allclose(xya_dict['attributelist']['2.0'],answer['2.0'])
        assert allclose(xya_dict['pointlist'],[[ 0.0, 0.6],[0.0, 0.9],[ 0.0,0.1],[ 0.4,1.0]])

        # Try another quantity
        interp = Interpolate_sww(sww.filename, 'stage')
        interp.interpolate_xya(point_file)

        answer['0.0'] =  [ 0.08, 0.02, 0.14,-0.32]
        answer['2.0'] =  [ 0.08, 0.02, 0.14,-0.32]        
        #print "answer",answer
        #print "interp.interpolated_quantity",interp.interpolated_quantity
        assert allclose(interp.interpolated_quantity['0.0'],answer['0.0'])
        assert allclose(interp.interpolated_quantity['2.0'],answer['2.0'])

        # look at error catching
        try:
            interp = Interpolate_sww(sww.filename, 'funky!')
        except KeyError:
            pass
        else:
            self.failUnless(0==1,
                        'bad key did not raise an error!')

        # look at error catching
        try:
            interp = Interpolate_sww(sww.filename, 'elevation')
        except KeyError:
            pass
        else:
            self.failUnless(0==1,
                        'bad key did not raise an error!')

        #Cleanup
        os.remove(sww.filename)
        os.remove(point_file_out)
        os.remove(point_file)

    def test_Interpolate_sww_errors(self):
        import tempfile
        import os
        try:
            interpolate_sww2xya('??ffd??', 'stage','yeah','yeah.x',
                                verbose = False)
        except SystemExit:  pass
        else:
            self.failUnless(0 ==1,  'Bad file did not raise error!') 
        
    def DISABLED_TEST_test_Interpolate_sww_errorsII(self):
        """
        THIS TEST HAS BEEN DISABLED, SINCE IT PRINTS TO SCREEN,
        but is should still work.  test it sometimes!
        """
        import tempfile
        import os
        import sys
       
        sww_file = tempfile.mktemp(".sww")
        fd = open(sww_file,'w')
        fd.write("unit testing a bad .sww file \n")
        fd.close()
        
        try:
            interpolate_sww2xya(sww_file, 'stage','yeah','yeah.x',
                                verbose = False)
            
        except SystemExit:  pass
        else:
            self.failUnless(0 ==1,  'Bad file did not raise error!')        
        #clean up
        os.remove(sww_file)
        
    def test_Interpolate_sww_errorsIII(self):
        import time, os
        from Numeric import array, zeros, allclose, Float, concatenate, \
             transpose
        from Scientific.IO.NetCDF import NetCDFFile
        from load_mesh.loadASCII import  import_points_file, \
             concatinate_attributelist

        self.domain.set_name('datatest' + str(time.time()))
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = mean

        sww = get_dataobject(self.domain)
        sww.store_connectivity()
        sww.store_timestep('stage')
        self.domain.time = 2.
        sww.store_timestep('stage')
        
        try:
            interpolate_sww2xya(self.domain.get_name(),
                                'stage','yeah','yeah.x',
                                verbose = False)
        except SystemExit:  pass
        else:
            self.failUnless(0 ==1,  'Bad file did not raise error!')        
        #clean up
        os.remove(sww.filename)
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Interpolate_sww,'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
