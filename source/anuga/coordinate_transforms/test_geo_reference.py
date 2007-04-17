#!/usr/bin/env python
#

import unittest
import tempfile
import os

from geo_reference import *
from Numeric import allclose,array

# Ignore these warnings, since we still want to test .xya code.
import warnings
warnings.filterwarnings(action = 'ignore',
                        message='.xya format is deprecated.  Please use .txt.',
                        category=DeprecationWarning)

warnings.filterwarnings(action = 'ignore',
                        message='Text file format is moving to comma se',
                        category=DeprecationWarning)


class geo_referenceTestCase(unittest.TestCase):
    def setUp(self):
        pass
        
    def tearDown(self):
        pass

    
    
    def test_get_origin(self):
        g = Geo_reference(56,1.9,1.9)
        (z,x,y) = g.get_origin()

        self.failUnless(z == g.get_zone(), ' failed')
        self.failUnless(x == g.get_xllcorner(), ' failed')
        self.failUnless(y == g.get_yllcorner(), ' failed') 
        
    def test_read_write_NetCDF(self):
        from Scientific.IO.NetCDF import NetCDFFile
        g = Geo_reference(56,1.9,1.9)
        file_name = tempfile.mktemp(".geo_referenceTest")
        
        out_file = NetCDFFile(file_name, 'w')
        g.write_NetCDF(out_file)
        out_file.close()
        
        in_file = NetCDFFile(file_name, 'r')
        new_g = Geo_reference(NetCDFObject=in_file)
        in_file.close()
        os.remove(file_name)

        self.failUnless(g == new_g, 'test_read_write_NetCDF failed')  
        
    def test_read_write_NetCDFII(self):
        from Scientific.IO.NetCDF import NetCDFFile
        g = Geo_reference(56,1.9,1.9)
        file_name = tempfile.mktemp(".geo_referenceTest")
        
        outfile = NetCDFFile(file_name, 'w')
        outfile.xllcorner = g.get_xllcorner() 
        outfile.yllcorner =  g.get_yllcorner() 
        outfile.zone = g.get_zone()
        outfile.close()
        
        in_file = NetCDFFile(file_name, 'r')
        new_g = Geo_reference(NetCDFObject=in_file)
        in_file.close()
        os.remove(file_name)

        self.failUnless(g == new_g, ' failed')
        
    def test_read_write_ASCII(self):
        from Scientific.IO.NetCDF import NetCDFFile
        g = Geo_reference(56,1.9,1.9)
        file_name = tempfile.mktemp(".geo_referenceTest")
        fd = open(file_name,'w')
        g.write_ASCII(fd)
        fd.close()
        
        fd = open(file_name,'r')
        new_g = Geo_reference(ASCIIFile=fd)
        fd.close()
        os.remove(file_name)

        self.failUnless(g == new_g, 'test_read_write_ASCII failed')  
    
    def test_read_write_ASCII2(self):
        from Scientific.IO.NetCDF import NetCDFFile
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

        self.failUnless(g == new_g, 'test_read_write_ASCII failed')
        
    def test_read_write_ASCII3(self):
        from Scientific.IO.NetCDF import NetCDFFile
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
            self.failUnless(0 ==1,
                        'bad text file did not raise error!')
            
    def test_change_points_geo_ref(self):
        x = 433.0
        y = 3.0
        g = Geo_reference(56,x,y)
        lofl = [[3.0,311.0], [677.0,6.0]]
        new_lofl = g.change_points_geo_ref(lofl)
        #print "lofl",lofl 
        #print "new_lofl",new_lofl

        self.failUnless(type(new_lofl) == types.ListType, ' failed')
        self.failUnless(type(new_lofl) == type(lofl), ' failed')
        for point,new_point in map(None,lofl,new_lofl):
            self.failUnless(point[0]-x==new_point[0], ' failed')
            self.failUnless(point[1]-y==new_point[1], ' failed')
         
        
    def test_change_points_geo_ref2(self):
        x = 3.0
        y = 543.0
        g = Geo_reference(56,x,y)
        lofl = [[3.0,388.0]]
        new_lofl = g.change_points_geo_ref(lofl)
        #print "lofl",lofl 
        #print "new_lofl",new_lofl

        self.failUnless(type(new_lofl) == types.ListType, ' failed')
        self.failUnless(type(new_lofl) == type(lofl), ' failed')
        for point,new_point in map(None,lofl,new_lofl):
            self.failUnless(point[0]-x==new_point[0], ' failed')
            self.failUnless(point[1]-y==new_point[1], ' failed')
        
    def test_change_points_geo_ref3(self):
        x = 3.0
        y = 443.0
        g = Geo_reference(56,x,y)
        lofl = [3.0,345.0]
        new_lofl = g.change_points_geo_ref(lofl)
        #print "lofl",lofl 
        #print "new_lofl",new_lofl

        self.failUnless(type(new_lofl) == types.ListType, ' failed')
        self.failUnless(type(new_lofl) == type(lofl), ' failed')
        for point,new_point in map(None,[lofl],new_lofl):
            self.failUnless(point[0]-x==new_point[0], ' failed')
            self.failUnless(point[1]-y==new_point[1], ' failed')
        
    
    def test_change_points_geo_ref4(self):
        x = 3.0
        y = 443.0
        g = Geo_reference(56,x,y)
        lofl = array([[3.0,323.0], [6.0,645.0]])
        new_lofl = g.change_points_geo_ref(lofl)
        #print "4 lofl",lofl 
        #print "4 new_lofl",new_lofl

        self.failUnless(type(new_lofl) == ArrayType, ' failed')
        self.failUnless(type(new_lofl) == type(lofl), ' failed')
        lofl[:,0] -= x
        lofl[:,1] -= y
        assert allclose(lofl,new_lofl)
        
    def test_change_points_geo_ref5(self):
        x = 103.0
        y = 3.0
        g = Geo_reference(56,x,y)
        lofl = array([[3.0,323.0]])

        
        #print "5 lofl before",lofl         
        new_lofl = g.change_points_geo_ref(lofl.copy())
        #print "5 lofl",lofl 
        #print "5 new_lofl",new_lofl

        self.failUnless(type(new_lofl) == ArrayType, ' failed')
        self.failUnless(type(new_lofl) == type(lofl), ' failed')


        for point,new_point in map(None,lofl,new_lofl):
            self.failUnless(point[0]-x==new_point[0], ' failed')
            self.failUnless(point[1]-y==new_point[1], ' failed')
        
    def test_change_points_geo_ref6(self):
        x = 53.0
        y = 3.0
        g = Geo_reference(56,x,y)
        lofl = array([355.0,3.0])
        new_lofl = g.change_points_geo_ref(lofl.copy())        
        #print "lofl",lofl 
        #print "new_lofl",new_lofl

        self.failUnless(type(new_lofl) == ArrayType, ' failed')
        self.failUnless(type(new_lofl) == type(lofl), ' failed')
        for point,new_point in map(None,[lofl],new_lofl):
            self.failUnless(point[0]-x==new_point[0], ' failed')
            self.failUnless(point[1]-y==new_point[1], ' failed')
     
    def test_change_points_geo_ref7(self):
        x = 23.0
        y = 3.0
        point_x = 9.0
        point_y = -60.0
        g = Geo_reference(56,x,y)
        points_geo_ref = Geo_reference(56,point_x,point_y)
        lofl = [[3.0,30.0], [67.0,6.0]]
        new_lofl = g.change_points_geo_ref(lofl,points_geo_ref=points_geo_ref)
        #print "lofl",lofl 
        #print "new_lofl",new_lofl

        self.failUnless(type(new_lofl) == types.ListType, ' failed')
        self.failUnless(type(new_lofl) == type(lofl), ' failed')
        for point,new_point in map(None,lofl,new_lofl):
            self.failUnless(point[0]+point_x-x==new_point[0], ' failed')
            self.failUnless(point[1]+point_y-y==new_point[1], ' failed')
      

    def test_get_absolute(self):
        x = 7.0
        y = 3.0
        
        g = Geo_reference(56,x,y)
        lofl = [[3.0,34.0], [64.0,6.0]]
        new_lofl = g.get_absolute(lofl)
        #print "lofl",lofl 
        #print "new_lofl",new_lofl

        self.failUnless(type(new_lofl) == types.ListType, ' failed')
        self.failUnless(type(new_lofl) == type(lofl), ' failed')
        for point,new_point in map(None,lofl,new_lofl):
            self.failUnless(point[0]+x==new_point[0], ' failed')
            self.failUnless(point[1]+y==new_point[1], ' failed')

            
        g = Geo_reference()
        lofl = [[3.0,34.0], [64.0,6.0]]
        new_lofl = g.get_absolute(lofl)
        #print "lofl",lofl 
        #print "new_lofl",new_lofl

        self.failUnless(type(new_lofl) == types.ListType, ' failed')
        self.failUnless(type(new_lofl) == type(lofl), ' failed')
        for point,new_point in map(None,lofl,new_lofl):
            self.failUnless(point[0]==new_point[0], ' failed')
            self.failUnless(point[1]==new_point[1], ' failed')
            
    def test_is_absolute(self):
        
        g = Geo_reference(34,0,0)
        points = [[3.0,34.0], [64.0,6.0]]

        assert g.is_absolute()

        g = Geo_reference(34,7,-6)
        assert not g.is_absolute()        

                        
    def test___cmp__(self):
        g = Geo_reference(56,1.9,1.9,)
        new_g = Geo_reference(56,1.9,1.9)
     
        self.failUnless(g == new_g, 'test___cmp__ failed')   


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
            raise msg
  
    def test_bad_ASCII_title(self):      
 # create an .xya file
        point_file = tempfile.mktemp(".xya")
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
            self.failUnless(0 ==1,
                        'bad text file did not raise error!')
            os.remove(point_file)

    def test_read_write_ASCII_test_and_fail(self):
        from Scientific.IO.NetCDF import NetCDFFile

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
            self.failUnless(0 ==1,
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

        self.failUnless(g == new_g, 'test_read_write_ASCII failed')
        
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

        self.failUnless(g == new_g, 'test_read_write_ASCII failed')
        
    def xxtest_good_title(self):      
 # create an .xya file
        point_file = tempfile.mktemp(".xya")
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
        except TitleError:
            fd.close()
            os.remove(point_file)
        else:
            self.failUnless(0 ==1,
                        'bad text file did not raise error!')
            os.remove(point_file)

    def test_error_message_ShapeError(self):
        
        new_g = Geo_reference()
        try:
            new_g.get_absolute((8.9, 7.8, 9.0)) 
        except ShapeError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad shape did not raise error!')
            os.remove(point_file)
            
        new_g = Geo_reference()
        try:
            new_g.get_absolute(((8.9, 7.8, 9.0))) 
        except ShapeError:
            pass
        else:
            self.failUnless(0 ==1,
                        'bad shape did not raise error!')
            os.remove(point_file)
        
#-------------------------------------------------------------
if __name__ == "__main__":

    suite = unittest.makeSuite(geo_referenceTestCase,'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)
    
