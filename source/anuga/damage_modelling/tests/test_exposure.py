import csv
import unittest, os
import tempfile
import numpy as num
import sys


from anuga.damage_modelling.exposure import Exposure

from anuga.anuga_exceptions import TitleValueError, \
                                    DataMissingValuesError

class Test_Exposure(unittest.TestCase):
    def test_exposure_csv_loading(self):
        file_name = tempfile.mktemp(".csv")
        file = open(file_name,"w")
        file.write("LATITUDE, LONGITUDE ,sound  , speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 40.0\n")
        file.close()
        exposure = Exposure(file_name, title_check_list = ['speed','sound'])
        exposure.get_column("sound")
       
        self.assertTrue(exposure._attribute_dic['sound'][2]==' bang',
                        'FAILED!')
        self.assertTrue(exposure._attribute_dic['speed'][2]==' 40.0',
                        'FAILED!')
        
        os.remove(file_name)
        
    def test_exposure_csv_loadingII(self):
        

        file_name = tempfile.mktemp(".txt")
        file = open(file_name,"w")
        file.write("LATITUDE, LONGITUDE ,sound  , speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 40.0\n")
        file.close()
        exposure = Exposure(file_name)
        exposure.get_column("sound")
       
        self.assertTrue(exposure._attribute_dic['sound'][2]==' bang',
                        'FAILED!')
        self.assertTrue(exposure._attribute_dic['speed'][2]==' 40.0',
                        'FAILED!')
        
        os.remove(file_name)
        
    def test_exposure_csv_loading_title_check_list(self):

        # I can't get cvs.reader to close the exposure file
        # The hacks below are to get around this.        
        if sys.platform == 'win32':
            file_name = tempfile.gettempdir() + \
                    "test_exposure_csv_loading_title_check_list.csv"
        else:
            file_name = tempfile.mktemp(".csv")
        file = open(file_name,"w")
        file.write("LATITUDE, LONGITUDE ,sound  , speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 40.0\n")
        file.close()
        try:
            exposure = Exposure(file_name, title_check_list = ['SOUND'])
        except IOError:
            pass
        else:
            self.assertTrue(0 ==1,  'Assertion not thrown error!')
            
        if not sys.platform == 'win32':
            os.remove(file_name)
        
    def test_exposure_csv_cmp(self):
        file_name = tempfile.mktemp(".csv")
        file = open(file_name,"w")
        file.write("LATITUDE, LONGITUDE ,sound  , speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 40.0\n")
        file.close()
        
        e1 = Exposure(file_name)
        e2 = Exposure(file_name)
        os.remove(file_name)

        self.assertTrue(cmp(e1,e2)==0,
                        'FAILED!')
        
        self.assertTrue(cmp(e1,"hey")==1,
                        'FAILED!')
        
        file_name = tempfile.mktemp(".csv")
        file = open(file_name,"w")
        # Note, this has less spaces in the title,
        # the instances will be the same.
        file.write("LATITUDE,LONGITUDE ,sound, speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 40.0\n")
        file.close()
        e3 = Exposure(file_name)
        os.remove(file_name)

        self.assertTrue(cmp(e3,e2)==0,
                        'FAILED!')
        
        file_name = tempfile.mktemp(".csv")
        file = open(file_name,"w")
        # Note, 40 changed to 44 .
        file.write("LATITUDE,LONGITUDE ,sound, speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 44.0\n")
        file.close()
        e4 = Exposure(file_name)
        os.remove(file_name)
        #print "e4",e4._attribute_dic 
        #print "e2",e2._attribute_dic 
        self.assertTrue(cmp(e4,e2)<>0,
                        'FAILED!')
        
        file_name = tempfile.mktemp(".csv")
        file = open(file_name,"w")
        # Note, the first two columns are swapped.
        file.write("LONGITUDE,LATITUDE ,sound, speed \n\
 -21.0,115.0, splat, 0.0\n\
 -21.7,114.0, pow, 10.0\n\
 -21.4,114.5, bang, 40.0\n")
        file.close()
        e5 = Exposure(file_name)
        os.remove(file_name)

        self.assertTrue(cmp(e3,e5)<>0,
                        'FAILED!')
        
    def test_exposure_csv_saving(self):
        

        file_name = tempfile.mktemp(".csv")
        file = open(file_name,"w")
        file.write("LATITUDE, LONGITUDE ,sound  , speed \n\
115.0, -21.0, splat, 0.0\n\
114.0, -21.7, pow, 10.0\n\
114.5, -21.4, bang, 40.0\n")
        file.close()
        e1 = Exposure(file_name)
        
        file_name2 = tempfile.mktemp(".csv")
        e1.save(file_name = file_name2)
        e2 = Exposure(file_name2)
       
        self.assertTrue(cmp(e1,e2)==0,
                        'FAILED!')
        os.remove(file_name)
        os.remove(file_name2)

    def test_exposure_csv_get_location(self):
        file_name = tempfile.mktemp(".csv")
        file = open(file_name,"w")
        file.write("LONGITUDE , LATITUDE, sound  , speed \n\
150.916666667, -34.5, splat, 0.0\n\
150.0, -34.0, pow, 10.0\n")
        file.close()
        e1 = Exposure(file_name)

        gsd = e1.get_location()
        
        points = gsd.get_data_points(absolute=True)
        
        assert num.allclose(points[0][0], 308728.009)
        assert num.allclose(points[0][1], 6180432.601)
        assert num.allclose(points[1][0],  222908.705)
        assert num.allclose(points[1][1], 6233785.284)
        self.assertTrue(gsd.get_geo_reference().get_zone() == 56,
                        'Bad zone error!')

        os.remove(file_name)
        
    def test_exposure_csv_set_column_get_column(self):
        file_name = tempfile.mktemp(".csv")
        file = open(file_name,"w")
        file.write("LONGITUDE , LATITUDE, sound  , speed \n\
150.916666667, -34.5, splat, 0.0\n\
150.0, -34.0, pow, 10.0\n")
        file.close()
        e1 = Exposure(file_name)      
        os.remove(file_name)

        new_title = "feast"
        new_values = ["chicken","soup"]
        e1.set_column(new_title, new_values)
        returned_values = e1.get_column(new_title)
        self.assertTrue(returned_values == new_values,
                        ' Error!')
        
        file_name2 = tempfile.mktemp(".csv")
        e1.save(file_name = file_name2)
        e2 = Exposure(file_name2)
        returned_values = e2.get_column(new_title)
        self.assertTrue(returned_values == new_values,
                        ' Error!')       
        os.remove(file_name2)

    def test_exposure_csv_set_column_get_column_error_checking(self):
        file_name = tempfile.mktemp(".csv")
        file = open(file_name,"w")
        file.write("LONGITUDE , LATITUDE, sound  , speed \n\
150.916666667, -34.5, splat, 0.0\n\
150.0, -34.0, pow, 10.0\n")
        file.close()
        e1 = Exposure(file_name)      
        os.remove(file_name)

        new_title = "sound"
        new_values = [12.5,7.6]
        try:
            e1.set_column(new_title, new_values)
        except TitleValueError:
            pass
        else:
            self.assertTrue(0 ==1,  'Error not thrown error!')
            
        e1.set_column(new_title, new_values, overwrite=True)
        returned_values = e1.get_column(new_title)
        self.assertTrue(returned_values == new_values,
                        ' Error!')       
        
        new2_title = "short list"
        new2_values = [12.5]
        try:
            e1.set_column(new2_title, new2_values)
        except DataMissingValuesError:
            pass
        else:
            self.assertTrue(0 ==1,  'Error not thrown error!')
            
        new2_title = "long list"
        new2_values = [12.5, 7,8]
        try:
            e1.set_column(new2_title, new2_values)
        except DataMissingValuesError:
            pass
        else:
            self.assertTrue(0 ==1,  'Error not thrown error!')
        file_name2 = tempfile.mktemp(".csv")
        e1.save(file_name = file_name2)
        e2 = Exposure(file_name2)
        returned_values = e2.get_column(new_title)
        for returned, new in map(None, returned_values, new_values):
            self.assertTrue(returned == str(new), ' Error!')
        #self.assertTrue(returned_values == new_values, ' Error!')       
        os.remove(file_name2)
        
        try:
            e1.get_column("toe jam")
        except TitleValueError:
            pass
        else:
            self.assertTrue(0 ==1,  'Error not thrown error!')
            
    def test_exposure_csv_loading_x_y(self):
        

        file_name = tempfile.mktemp(".csv")
        file = open(file_name,"w")
        file.write("x, y ,sound  , speed \n\
115.0, 7, splat, 0.0\n\
114.0, 8.0, pow, 10.0\n\
114.5, 9., bang, 40.0\n")
        file.close()
        e1 = Exposure(file_name, is_x_y_locations=True)
        gsd = e1.get_location()
        
        points = gsd.get_data_points(absolute=True)
        
        assert num.allclose(points[0][0], 115)
        assert num.allclose(points[0][1], 7)
        assert num.allclose(points[1][0], 114)
        assert num.allclose(points[1][1], 8)
        assert num.allclose(points[2][0], 114.5)
        assert num.allclose(points[2][1], 9)
        self.assertTrue(gsd.get_geo_reference().get_zone() == -1,
                        'Bad zone error!')

        os.remove(file_name)

           
    def test_exposure_csv_loading_x_y2(self):
        
        csv_file = tempfile.mktemp(".csv")
        fd = open(csv_file,'wb')
        writer = csv.writer(fd)
        writer.writerow(['x','y','STR_VALUE','C_VALUE','ROOF_TYPE','WALLS', 'SHORE_DIST'])
        writer.writerow([5.5,0.5,'199770','130000','Metal','Timber',20])
        writer.writerow([4.5,1.0,'150000','76000','Metal','Double Brick',20])
        writer.writerow([4.5,1.5,'150000','76000','Metal','Brick Veneer',20])
        fd.close()

        e1 = Exposure(csv_file)
        gsd = e1.get_location()
        
        points = gsd.get_data_points(absolute=True)
        assert num.allclose(points[0][0], 5.5)
        assert num.allclose(points[0][1], 0.5)
        assert num.allclose(points[1][0], 4.5)
        assert num.allclose(points[1][1], 1.0)
        assert num.allclose(points[2][0], 4.5)
        assert num.allclose(points[2][1], 1.5)
        self.assertTrue(gsd.get_geo_reference().get_zone() == -1,
                        'Bad zone error!')

        os.remove(csv_file)


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Exposure,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

