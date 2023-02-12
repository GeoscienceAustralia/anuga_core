#!/usr/bin/env python
#

from builtins import zip
from builtins import str
import unittest
import tempfile
import os

from anuga.coordinate_transforms.geo_reference import *
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a

import numpy as num


class geo_referenceTestCase(unittest.TestCase):
    def setUp(self):
        pass
        
    def tearDown(self):
        pass

    def test_get_origin(self):
        g = Geo_reference(56,1.9,1.9)
        (z,x,y) = g.get_origin()

        self.assertTrue(z == g.get_zone(), ' failed')
        self.assertTrue(x == g.get_xllcorner(), ' failed')
        self.assertTrue(y == g.get_yllcorner(), ' failed') 
        
    def test_read_write_NetCDF(self):
        from anuga.file.netcdf import NetCDFFile
        g = Geo_reference(56,1.9,1.9)
        file_name = tempfile.mktemp(".geo_referenceTest")
        
        out_file = NetCDFFile(file_name, netcdf_mode_w)
        g.write_NetCDF(out_file)
        out_file.close()
        
        in_file = NetCDFFile(file_name, netcdf_mode_r)
        new_g = Geo_reference(NetCDFObject=in_file)
        in_file.close()
        os.remove(file_name)

        self.assertTrue(g == new_g, 'test_read_write_NetCDF failed')  
        
    def test_read_NetCDFI(self):
        # test if read_NetCDF
        from anuga.file.netcdf import NetCDFFile
        
        g = Geo_reference(56,1.9,1.9)
        file_name = tempfile.mktemp(".geo_referenceTest")
        outfile = NetCDFFile(file_name, netcdf_mode_w)
        g.write_NetCDF(outfile)
        outfile.close()
        
        in_file = NetCDFFile(file_name, netcdf_mode_r)
        new_g = Geo_reference(NetCDFObject=in_file)
        in_file.close()
        os.remove(file_name)

        self.assertTrue(g == new_g, ' failed')
        
    def test_read_write_ASCII(self):
        from anuga.file.netcdf import NetCDFFile
        g = Geo_reference(56,1.9,1.9)
        file_name = tempfile.mktemp(".geo_referenceTest")
        fd = open(file_name,'w')
        g.write_ASCII(fd)
        fd.close()
        
        fd = open(file_name,'r')
        new_g = Geo_reference(ASCIIFile=fd)
        fd.close()
        os.remove(file_name)

        self.assertTrue(g == new_g, 'test_read_write_ASCII failed')  
    
    def test_read_write_ASCII2(self):
        from anuga.file.netcdf import NetCDFFile
        g = Geo_reference(56,1.9,1.9)
        file_name = tempfile.mktemp(".geo_referenceTest")
        fd = open(file_name,'w')
        g.write_ASCII(fd)
        fd.close()       
        fd = open(file_name,'r')
        line = fd.readline()
        new_g = Geo_reference(ASCIIFile=fd, read_title=line)
        fd.close()
        os.remove(file_name)

        self.assertTrue(g == new_g, 'test_read_write_ASCII failed')
        
    def test_read_write_ASCII3(self):
        from anuga.file.netcdf import NetCDFFile
        g = Geo_reference(56,1.9,1.9)
        file_name = tempfile.mktemp(".geo_referenceTest")
        fd = open(file_name,'w')
        g.write_ASCII(fd)
        fd.close()       
        fd = open(file_name,'r')
        line = fd.readline()
        line = "fail !!"
        try:
            new_g = Geo_reference(ASCIIFile=fd, read_title=line)
            fd.close()
            os.remove(file_name)
        except TitleError:
            fd.close()
            os.remove(file_name)
        else:
            self.assertTrue(0 ==1,
                        'bad text file did not raise error!')
            
    def test_change_points_geo_ref(self):
        x = 433.0
        y = 3.0
        g = Geo_reference(56,x,y)
        lofl = [[3.0,311.0], [677.0,6.0]]
        new_lofl = g.change_points_geo_ref(lofl)

        self.assertTrue(isinstance(new_lofl, list), ' failed')
        self.assertTrue(type(new_lofl) == type(lofl), ' failed')
        for point,new_point in zip(lofl,new_lofl):
            self.assertTrue(point[0]-x==new_point[0], ' failed')
            self.assertTrue(point[1]-y==new_point[1], ' failed')
         
        
    def test_change_points_geo_ref2(self):
        x = 3.0
        y = 543.0
        g = Geo_reference(56,x,y)
        lofl = [[3.0,388.0]]
        new_lofl = g.change_points_geo_ref(lofl)

        self.assertTrue(isinstance(new_lofl, list), ' failed')
        self.assertTrue(type(new_lofl) == type(lofl), ' failed')
        for point,new_point in zip(lofl,new_lofl):
            self.assertTrue(point[0]-x==new_point[0], ' failed')
            self.assertTrue(point[1]-y==new_point[1], ' failed')
        
    def test_change_points_geo_ref3(self):
        x = 3.0
        y = 443.0
        g = Geo_reference(56,x,y)
        lofl = [3.0,345.0]
        new_lofl = g.change_points_geo_ref(lofl)

        self.assertTrue(isinstance(new_lofl, list), ' failed')
        self.assertTrue(type(new_lofl) == type(lofl), ' failed')
        for point,new_point in zip([lofl],new_lofl):
            self.assertTrue(point[0]-x==new_point[0], ' failed')
            self.assertTrue(point[1]-y==new_point[1], ' failed')
        
    
    def test_change_points_geo_ref4(self):
        x = 3.0
        y = 443.0
        g = Geo_reference(56,x,y)
        lofl = num.array([[3.0,323.0], [6.0,645.0]])
        new_lofl = g.change_points_geo_ref(lofl)

        self.assertTrue(isinstance(new_lofl, num.ndarray), ' failed')
        self.assertTrue(type(new_lofl) == type(lofl), ' failed')
        lofl[:,0] -= x
        lofl[:,1] -= y
        assert num.allclose(lofl,new_lofl)
        
    def test_change_points_geo_ref5(self):
        x = 103.0
        y = 3.0
        g = Geo_reference(56,x,y)
        lofl = num.array([[3.0,323.0]])

        new_lofl = g.change_points_geo_ref(lofl.copy())

        self.assertTrue(isinstance(new_lofl, num.ndarray), ' failed')
        self.assertTrue(type(new_lofl) == type(lofl), ' failed')


        for point,new_point in zip(lofl,new_lofl):
            self.assertTrue(point[0]-x==new_point[0], ' failed')
            self.assertTrue(point[1]-y==new_point[1], ' failed')
        
    def test_change_points_geo_ref6(self):
        x = 53.0
        y = 3.0
        g = Geo_reference(56,x,y)
        lofl = num.array([355.0,3.0])
        new_lofl = g.change_points_geo_ref(lofl.copy())        

        self.assertTrue(isinstance(new_lofl, num.ndarray), ' failed')
        self.assertTrue(type(new_lofl) == type(lofl), ' failed')
        for point,new_point in zip([lofl],new_lofl):
            self.assertTrue(point[0]-x==new_point[0], ' failed')
            self.assertTrue(point[1]-y==new_point[1], ' failed')
     
    def test_change_points_geo_ref7(self):
        x = 23.0
        y = 3.0
        point_x = 9.0
        point_y = -60.0
        g = Geo_reference(56,x,y)
        points_geo_ref = Geo_reference(56,point_x,point_y)
        lofl = [[3.0,30.0], [67.0,6.0]]
        new_lofl = g.change_points_geo_ref(lofl,points_geo_ref=points_geo_ref)

        self.assertTrue(isinstance(new_lofl, list), ' failed')
        self.assertTrue(type(new_lofl) == type(lofl), ' failed')
        for point,new_point in zip(lofl,new_lofl):
            self.assertTrue(point[0]+point_x-x==new_point[0], ' failed')
            self.assertTrue(point[1]+point_y-y==new_point[1], ' failed')

    def test_get_absolute_list(self):
        # test with supplied offsets
        x = 7.0
        y = 3.0
        
        g = Geo_reference(56, x, y)
        points = [[3.0,34.0], [64.0,6.0]]
        new_points = g.get_absolute(points)

        self.assertTrue(isinstance(new_points, list), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        for point, new_point in zip(points, new_points):
            self.assertTrue(point[0]+x == new_point[0], 'failed')
            self.assertTrue(point[1]+y == new_point[1], 'failed')

        # test with no supplied offsets
        g = Geo_reference()
        points = [[3.0,34.0], [64.0,6.0]]
        new_points = g.get_absolute(points)

        self.assertTrue(isinstance(new_points, list), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        for point, new_point in zip(points, new_points):
            self.assertTrue(point[0] == new_point[0], 'failed')
            self.assertTrue(point[1] == new_point[1], 'failed')
            
        # test that calling get_absolute twice does the right thing
        # first call
        dx = 10.0
        dy = 12.0
        g = Geo_reference(56, dx, dy)
        points = [[3.0,34.0], [64.0,6.0]]
        expected_new_points = [[3.0+dx,34.0+dy], [64.0+dx,6.0+dy]]
        new_points = g.get_absolute(points)

        self.assertTrue(isinstance(new_points, list), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        self.assertTrue(new_points == expected_new_points, 'failed')

        # and repeat from 'new_points = g.get_absolute(points)' above
        # to see if second call with same input gives same results.
        new_points = g.get_absolute(points)

        self.assertTrue(isinstance(new_points, list), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        self.assertTrue(new_points == expected_new_points, 'failed')

    def test_get_absolute_array(self):
        '''Same test as test_get_absolute_list(), but with numeric arrays.'''

        # test with supplied offsets
        x = 7.0
        y = 3.0
        
        g = Geo_reference(56, x, y)
        points = num.array([[3.0,34.0], [64.0,6.0]])
        new_points = g.get_absolute(points)

        self.assertTrue(isinstance(new_points, num.ndarray), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        msg = 'points=\n%s\nnew_points=\n%s' % (str(points), str(new_points))
        for point, new_point in zip(points, new_points):
            self.assertTrue(point[0]+x == new_point[0], msg)
            self.assertTrue(point[1]+y == new_point[1], msg)

        # test with no supplied offsets
        g = Geo_reference()
        points = num.array([[3.0,34.0], [64.0,6.0]])
        new_points = g.get_absolute(points)

        self.assertTrue(isinstance(new_points, num.ndarray), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        self.assertTrue(num.alltrue(points == new_points), 'failed')

        # test that calling get_absolute twice does the right thing
        # first call
        dx = 11.0
        dy = 13.0
        g = Geo_reference(56, dx, dy)
        points = num.array([[3.0,34.0], [64.0,6.0]])
        expected_new_points = num.array([[3.0+dx,34.0+dy], [64.0+dx,6.0+dy]])
        new_points = g.get_absolute(points)

        self.assertTrue(isinstance(new_points, num.ndarray), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        msg = ('First call of .get_absolute() returned %s\nexpected %s'
               % (str(new_points), str(expected_new_points)))
        self.assertTrue(num.alltrue(expected_new_points == new_points), msg)

        # and repeat from 'new_points = g.get_absolute(points)' above
        # to see if second call with same input gives same results.
        new_points = g.get_absolute(points)

        self.assertTrue(isinstance(new_points, num.ndarray), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        msg = ('Second call of .get_absolute() returned\n%s\nexpected\n%s'
               % (str(new_points), str(expected_new_points)))
        self.assertTrue(num.alltrue(expected_new_points == new_points), msg)

        # and repeat again to see if *third* call with same input
        # gives same results.
        new_points = g.get_absolute(points)

        self.assertTrue(isinstance(new_points, num.ndarray), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        msg = ('Third call of .get_absolute() returned %s\nexpected %s'
               % (str(new_points), str(expected_new_points)))
        self.assertTrue(num.alltrue(expected_new_points == new_points), msg)

    def test_get_relative_list(self):
        # test with supplied offsets
        x = 7.0
        y = 3.0
        
        g = Geo_reference(56, x, y)
        points = [[3.0,34.0], [64.0,6.0]]
        new_points = g.get_relative(points)

        self.assertTrue(isinstance(new_points, list), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        for point, new_point in zip(points, new_points):
            self.assertTrue(point[0]-x == new_point[0], 'failed')
            self.assertTrue(point[1]-y == new_point[1], 'failed')

        # test with no supplied offsets
        g = Geo_reference()
        points = [[3.0,34.0], [64.0,6.0]]
        new_points = g.get_relative(points)

        self.assertTrue(isinstance(new_points, list), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        for point, new_point in zip(points, new_points):
            self.assertTrue(point[0] == new_point[0], 'failed')
            self.assertTrue(point[1] == new_point[1], 'failed')
            
        # test that calling get_absolute twice does the right thing
        # first call
        dx = 10.0
        dy = 12.0
        g = Geo_reference(56, dx, dy)
        points = [[3.0,34.0], [64.0,6.0]]
        expected_new_points = [[3.0-dx,34.0-dy], [64.0-dx,6.0-dy]]
        new_points = g.get_relative(points)

        self.assertTrue(isinstance(new_points, list), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        self.assertTrue(new_points == expected_new_points, 'failed')

        # and repeat from 'new_points = g.get_absolute(points)' above
        # to see if second call with same input gives same results.
        new_points = g.get_relative(points)

        self.assertTrue(isinstance(new_points, list), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        self.assertTrue(new_points == expected_new_points, 'failed')

    def test_get_relative_array(self):
        '''Same test as test_get_relative_list(), but with numeric arrays.'''

        # test with supplied offsets
        x = 7.0
        y = 3.0
        
        g = Geo_reference(56, x, y)
        points = num.array([[3.0,34.0], [64.0,6.0]])
        new_points = g.get_relative(points)

        self.assertTrue(isinstance(new_points, num.ndarray), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        msg = 'points=\n%s\nnew_points=\n%s' % (str(points), str(new_points))
        for point, new_point in zip(points, new_points):
            self.assertTrue(point[0]-x == new_point[0], msg)
            self.assertTrue(point[1]-y == new_point[1], msg)

        # test with no supplied offsets
        g = Geo_reference()
        points = num.array([[3.0,34.0], [64.0,6.0]])
        new_points = g.get_relative(points)

        self.assertTrue(isinstance(new_points, num.ndarray), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        self.assertTrue(num.alltrue(points == new_points), 'failed')

        # test that calling get_relative twice does the right thing
        # first call
        dx = 11.0
        dy = 13.0
        g = Geo_reference(56, dx, dy)
        points = num.array([[3.0,34.0], [64.0,6.0]])
        expected_new_points = num.array([[3.0-dx,34.0-dy], [64.0-dx,6.0-dy]])
        new_points = g.get_relative(points)

        self.assertTrue(isinstance(new_points, num.ndarray), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        msg = ('First call of .get_relative() returned %s\nexpected %s'
               % (str(new_points), str(expected_new_points)))
        self.assertTrue(num.alltrue(expected_new_points == new_points), msg)

        # and repeat from 'new_points = g.get_relative(points)' above
        # to see if second call with same input gives same results.
        new_points = g.get_relative(points)

        self.assertTrue(isinstance(new_points, num.ndarray), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        msg = ('Second call of .get_relative() returned\n%s\nexpected\n%s'
               % (str(new_points), str(expected_new_points)))
        self.assertTrue(num.alltrue(expected_new_points == new_points), msg)

        # and repeat again to see if *third* call with same input
        # gives same results.
        new_points = g.get_relative(points)

        self.assertTrue(isinstance(new_points, num.ndarray), 'failed')
        self.assertTrue(type(new_points) == type(points), 'failed')
        msg = ('Third call of .get_relative() returned %s\nexpected %s'
               % (str(new_points), str(expected_new_points)))
        self.assertTrue(num.alltrue(expected_new_points == new_points), msg)

    def test_is_absolute(self):
        
        g = Geo_reference(34,0,0)
        points = [[3.0,34.0], [64.0,6.0]]

        assert g.is_absolute()

        g = Geo_reference(34,7,-6)
        assert not g.is_absolute()        

                        
    def test___cmp__(self):
        g = Geo_reference(56,1.9,1.9,)
        new_g = Geo_reference(56,1.9,1.9)
     
        self.assertTrue(g == new_g, 'test___cmp__ failed')   


    def test_reconcile(self):
        g1 = Geo_reference(56,2,5)
        g2 = Geo_reference(50,4,5)
        g3 = Geo_reference(50,66,6)        
        g_default = Geo_reference()                


        g2.reconcile_zones(g3)
        assert g2.get_zone() == g3.get_zone()

        g_default.reconcile_zones(g3)
        assert g_default.get_zone() == g3.get_zone()

        g_default = Geo_reference()                
        g3.reconcile_zones(g_default)
        assert g_default.get_zone() == g3.get_zone()                

        try:
            g1.reconcile_zones(g2)
        except:
            pass
        else:
            msg = 'Should have raised an exception'
            raise Exception(msg)  


    def test_set_hemisphere(self):
        g1 = Geo_reference(56,2,5)

        assert g1.hemisphere == 'undefined'

        g1.set_hemisphere('southern')
        assert g1.hemisphere == 'southern'

        # Generate exception with invalid hemisphere value
        try:
            g1.set_hemisphere('bogus')
        except:
            pass
        else:
            msg = 'Should have raised an exception' 

    def test_get_hemisphere(self):

        g1 = Geo_reference(56,2,5)

        assert g1.get_hemisphere() == 'undefined'

        g2 = Geo_reference()
        g2.set_hemisphere('southern')

        assert g2.get_hemisphere() == 'southern'

        assert g2.get_zone() == -1   

    def test_set_zone(self):
        g1 = Geo_reference(56,2,5)

        assert g1.zone == 56

        g1.set_zone('55')
        assert g1.zone == 55

        g1.set_zone(-1)
        assert g1.zone == -1

        # Generate exception with invalid zone value
        try:
            g1.set_zone(0)
        except:
            pass
        else:
            msg = 'Should have raised an exception'  

    def test_get_zone(self):

        g1 = Geo_reference(56,2,5)

        assert g1.get_zone() == 56

        g2 = Geo_reference()

        assert g2.get_zone() == -1

           
  
    def test_bad_ASCII_title(self):      
        # create an text file
        point_file = tempfile.mktemp(".xxx")
        fd = open(point_file,'w')
        fd.write("# hey! \n")
        fd.close()
        
        fd = open(point_file,'r')
        # 
        #new_g = Geo_reference(ASCIIFile=fd)
        try:
            new_g = Geo_reference(ASCIIFile=fd)
            fd.close()
            os.remove(point_file)
        except TitleError:
            fd.close()
            os.remove(point_file)
        else:
            self.assertTrue(0 ==1,
                        'bad text file did not raise error!')
            os.remove(point_file)

    def test_read_write_ASCII_test_and_fail(self):
        from anuga.file.netcdf import NetCDFFile

        # This is to test a fail
        g = Geo_reference(56,1.9,1.9)
        file_name = tempfile.mktemp(".geo_referenceTest")
        fd = open(file_name,'w')
        g.write_ASCII(fd)
        fd.close()       
        fd = open(file_name,'r')
        line = fd.readline()
        line = " #Geo"
        try:
            new_g = Geo_reference(ASCIIFile=fd, read_title=line)
            fd.close()
            os.remove(file_name)
        except TitleError:
            fd.close()
            os.remove(file_name)
        else:
            self.assertTrue(0 ==1,
                        'bad text file did not raise error!')

        # this tests a pass
        g = Geo_reference(56,1.9,1.9)
        file_name = tempfile.mktemp(".geo_referenceTest")
        fd = open(file_name,'w')
        g.write_ASCII(fd)
        fd.close()
        
        fd = open(file_name,'r')
        line = fd.readline()
        line = "#geo_yeah"
        new_g = Geo_reference(ASCIIFile=fd, read_title=line)
        fd.close()
        os.remove(file_name)

        self.assertTrue(g == new_g, 'test_read_write_ASCII failed')
        
        # this tests a pass
        g = Geo_reference(56,1.9,1.9)
        file_name = tempfile.mktemp(".geo_referenceTest")
        fd = open(file_name,'w')
        g.write_ASCII(fd)
        fd.close()
        
        fd = open(file_name,'r')
        line = fd.readline()
        line = "#geo crap"
        new_g = Geo_reference(ASCIIFile=fd, read_title=line)
        fd.close()
        os.remove(file_name)

        self.assertTrue(g == new_g, 'test_read_write_ASCII failed')
        
    def test_good_title(self):      
 # create an .xxx file
        point_file = tempfile.mktemp(".xxx")
        fd = open(point_file,'w')
        fd.write("#Geo crap \n 56\n ")
        fd.close()
        
        fd = open(point_file,'r')
        # 
        #new_g = Geo_reference(ASCIIFile=fd)
        try:
            new_g = Geo_reference(ASCIIFile=fd)
            fd.close()
            os.remove(point_file)
        except ValueError:
            fd.close()
            os.remove(point_file)
        else:
            self.assertTrue(0 ==1,
                        'bad text file did not raise error!')
            os.remove(point_file)

    def test_error_message_ShapeError(self):
        
        new_g = Geo_reference()
        try:
            new_g.get_absolute((8.9, 7.8, 9.0)) 
        except ShapeError:
            pass
        else:
            self.assertTrue(0 ==1,
                        'bad shape did not raise error!')
            os.remove(point_file)
            
        new_g = Geo_reference()
        try:
            new_g.get_absolute(((8.9, 7.8, 9.0))) 
        except ShapeError:
            pass
        else:
            self.assertTrue(0 ==1,
                        'bad shape did not raise error!')
            os.remove(point_file)

    def test_functionality_get_absolute(self):
        x0 = 1000.0
        y0 = 2000.0
        geo = Geo_reference(56, x0, y0)

        # iterable points (*not* num.array())
        points = ((2,3), (3,1), (5,2))
        abs_points = geo.get_absolute(points)
        # check we haven't changed 'points' itself
        self.assertFalse(num.alltrue(abs_points == points))
        new_points = abs_points.copy()
        new_points[:,0] -= x0
        new_points[:,1] -= y0
        self.assertTrue(num.alltrue(new_points == points))

        # points in num.array()
        points = num.array(((2,3), (3,1), (5,2)), float)
        abs_points = geo.get_absolute(points)
        # check we haven't changed 'points' itself
        self.assertFalse(num.alltrue(abs_points == points))
        new_points = abs_points.copy()
        new_points[:,0] -= x0
        new_points[:,1] -= y0
        self.assertTrue(num.alltrue(new_points == points))

    def test_georef_types(self):
        '''Ensure that attributes of a georeference are of correct type.
        
        zone            int
        false_easting   int
        false_northing  int
        xllcorner       float
        yllcorner       float
        '''

        from anuga.file.netcdf import NetCDFFile

        # ensure that basic instance attributes are correct
        g = Geo_reference(56, 1.8, 1.8)
        self.assertTrue(isinstance(g.zone, int),
                        "geo_ref .zone should be 'int' type, "
                        "was '%s' type" % type(g.zone))  
        self.assertTrue(isinstance(g.false_easting, int),
                        "geo_ref .false_easting should be int type, "
                        "was '%s' type" % type(g.false_easting))  
        self.assertTrue(isinstance(g.false_northing, int),
                        "geo_ref .false_northing should be int type, "
                        "was '%s' type" % type(g.false_northing))
        self.assertTrue(isinstance(g.xllcorner, float),
                        "geo_ref .xllcorner should be float type, "
                        "was '%s' type" % type(g.xllcorner))
        self.assertTrue(isinstance(g.yllcorner, float),
                        "geo_ref .yllcorner should be float type, "
                        "was '%s' type" % type(g.yllcorner))

        # now write fikle, read back and check types again
        file_name = tempfile.mktemp(".geo_referenceTest")
        
        out_file = NetCDFFile(file_name, netcdf_mode_w)
        g.write_NetCDF(out_file)
        out_file.close()
        
        in_file = NetCDFFile(file_name, netcdf_mode_r)
        new_g = Geo_reference(NetCDFObject=in_file)
        in_file.close()
        os.remove(file_name)

        self.assertTrue(isinstance(new_g.zone, int),
                        "geo_ref .zone should be 'int' type, "
                        "was '%s' type" % type(new_g.zone))  
        self.assertTrue(isinstance(new_g.false_easting, int),
                        "geo_ref .false_easting should be int type, "
                        "was '%s' type" % type(new_g.false_easting))  
        self.assertTrue(isinstance(new_g.false_northing, int),
                        "geo_ref .false_northing should be int type, "
                        "was '%s' type" % type(new_g.false_northing))
        self.assertTrue(isinstance(new_g.xllcorner, float),
                        "geo_ref .xllcorner should be float type, "
                        "was '%s' type" % type(new_g.xllcorner))
        self.assertTrue(isinstance(new_g.yllcorner, float),
                        "geo_ref .yllcorner should be float type, "
                        "was '%s' type" % type(new_g.yllcorner))
        
    def test_georef_types_coerceable(self):
        '''Ensure that attributes of a georeference are of correct type.
        
        zone            int
        false_easting   int
        false_northing  int
        xllcorner       float
        yllcorner       float
        '''

        # now provide wrong types but coerceable
        g = Geo_reference(56.0, '1.8', '1.8')
        self.assertTrue(isinstance(g.zone, int),
                        "geo_ref .zone should be 'int' type, "
                        "was '%s' type" % type(g.zone))  
        self.assertTrue(isinstance(g.false_easting, int),
                        "geo_ref .false_easting should be int type, "
                        "was '%s' type" % type(g.false_easting))  
        self.assertTrue(isinstance(g.false_northing, int),
                        "geo_ref .false_northing should be int type, "
                        "was '%s' type" % type(g.false_northing))
        self.assertTrue(isinstance(g.xllcorner, float),
                        "geo_ref .xllcorner should be float type, "
                        "was '%s' type" % type(g.xllcorner))
        self.assertTrue(isinstance(g.yllcorner, float),
                        "geo_ref .yllcorner should be float type, "
                        "was '%s' type" % type(g.yllcorner))

        
#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(geo_referenceTestCase, 'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)
    
