import os
import unittest
import tempfile
import numpy as num

from anuga.utilities.system_tools import get_pathname_from_package

from anuga.file.csv_file import load_csv_as_array, load_csv_as_dict, store_parameters, \
                        load_csv_as_matrix, load_csv_as_polygons, \
                        load_csv_as_building_polygons

class Test_csv(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def _create_csv_file(self):
        """
            Create a dummy csv file.
            Return its filename.
        """      
        filename = tempfile.mktemp(".txt")
        file = open(filename,"w")
        file.write("elevation, stage\n\
1.0, 3  \n\
0.0, 4 \n\
4.0, 3 \n\
1.0, 6 \n")
        file.close()
        return filename
        
    def test_get_data_from_file1(self):
        filename = self._create_csv_file()

        x = load_csv_as_array(filename)  
        
        os.remove(filename)

        assert num.allclose(x['elevation'], [1.0, 0.0, 4.0, 1.0])
        assert num.allclose(x['stage'], [3.0, 4.0, 3.0, 6.0])        


    def test_get_data_from_file(self):
        filename = self._create_csv_file()
        
        header,x = load_csv_as_matrix(filename)
        os.remove(filename)
        
        assert num.allclose(x[:,0], [1.0, 0.0,4.0, 1.0])

        
    def test_store_parameters(self):
        """tests store temporary file
        """
        
        from os import sep, getenv
        
        output_dir=''
        file_name='details.csv'
        
        kwargs = {'file_name':'new2.txt',
                  'output_dir':output_dir,
                  'file_name':file_name,
                  'who':'me',
                  'what':'detail',
                  'how':2,
                  'why':241,
#                  'completed':345
                  }
        store_parameters(verbose=False,**kwargs)

        temp='detail_temp.csv'
        fid = open(temp)
        file_header = fid.readline()
        file_line = fid.readline()
        fid.close()
        
        
        keys = kwargs.keys()
        keys.sort()
        line=''
        header=''
        count=0
        #used the sorted keys to create the header and line data
        for k in keys:
#            print "%s = %s" %(k, kwargs[k]) 
            header = header+str(k)
            line = line+str(kwargs[k])
            count+=1
            if count <len(kwargs):
                header = header+','
                line = line+','
        header+='\n'
        line+='\n'
        
        
        #file exists
        assert os.access(temp, os.F_OK)
        assert header == file_header
        assert line == file_line
        
        os.remove(temp)
        
    def test_store_parameters1(self):
        """tests store in temporary file and other file 
        """
        
        from os import sep, getenv
        
        output_dir=''
        file_name='details.csv'
        
        kwargs = {'file_name':'new2.txt',
                  'output_dir':output_dir,
                  'file_name':file_name,
                  'who':'me',
                  'what':'detail',
                  'how':2,
                  'why':241,
#                  'completed':345
                  }
        store_parameters(verbose=False,**kwargs)
        
        kwargs['how']=55
        kwargs['completed']=345

        keys = kwargs.keys()
        keys.sort()
        line=''
        header=''
        count=0
        #used the sorted keys to create the header and line data
        for k in keys:
#            print "%s = %s" %(k, kwargs[k]) 
            header = header+str(k)
            line = line+str(kwargs[k])
            count+=1
            if count <len(kwargs):
                header = header+','
                line = line+','
        header+='\n'
        line+='\n'
        
        kwargs['how']=55
        kwargs['completed']=345
        
        store_parameters(verbose=False,**kwargs)
        
#        temp='detail_temp.csv'
        fid = open(file_name)
        file_header = fid.readline()
        file_line1 = fid.readline()
        file_line2 = fid.readline()
        fid.close()
        
        
        #file exists
#        print 'header',header,'line',line
#        print 'file_header',file_header,'file_line1',file_line1,'file_line2',file_line2
        assert os.access(file_name, os.F_OK)
        assert header == file_header
        assert line == file_line1
        
        temp='detail_temp.csv'
        os.remove(temp)
        os.remove(file_name)        
        
    def test_store_parameters2(self):
        """tests appending the data to the end of an existing file
        """
        
        from os import sep, getenv
        
        output_dir=''
        file_name='details.csv'
        
        kwargs = {'file_name':'new2.txt',
                  'output_dir':output_dir,
                  'file_name':file_name,
                  'who':'me',
                  'what':'detail',
                  'how':2,
                  'why':241,
                  'completed':345
                  }
        store_parameters(verbose=False,**kwargs)
        
        kwargs['how']=55
        kwargs['completed']=23.54532
        
        store_parameters(verbose=False,**kwargs)
        
        keys = kwargs.keys()
        keys.sort()
        line=''
        header=''
        count=0
        #used the sorted keys to create the header and line data
        for k in keys:
#            print "%s = %s" %(k, kwargs[k]) 
            header = header+str(k)
            line = line+str(kwargs[k])
            count+=1
            if count <len(kwargs):
                header = header+','
                line = line+','
        header+='\n'
        line+='\n'
        
        fid = open(file_name)
        file_header = fid.readline()
        file_line1 = fid.readline()
        file_line2 = fid.readline()
        fid.close()
        
        assert os.access(file_name, os.F_OK)
        assert header == file_header
        assert line == file_line2
        
        os.remove(file_name)        
        

    
    def test_csv2polygons(self):
        """test loading of a csv polygon file.
        """
        
        path = get_pathname_from_package('anuga.shallow_water')                
        testfile = os.path.join(path, 'tests', 'data', 'polygon_values_example.csv')                

        polygons, values = load_csv_as_polygons(testfile, 
                                        value_name='floors')

        assert len(polygons) == 7, 'Must have seven polygons'
        assert len(values) == 7, 'Must have seven values'
            
        # Known floor values
        floors = {'1': 2, '2': 0, '3': 1, '4': 2, '5': 0, '8': 1, '9': 1}
        
        # Known polygon values
        known_polys = {'1': [[422681.61,871117.55],
                             [422691.02,871117.60],
                             [422690.87,871084.23],
                             [422649.36,871081.85],
                             [422649.36,871080.39],
                             [422631.86,871079.50],
                             [422631.72,871086.75],
                             [422636.75,871087.20],
                             [422636.75,871091.50],
                             [422649.66,871092.09],
                             [422649.83,871084.91],
                             [422652.94,871084.90],
                             [422652.84,871092.39],
                             [422681.83,871093.73],
                             [422681.61,871117.55]],
                       '2': [[422664.22,870785.46],
                             [422672.48,870780.14],
                             [422668.17,870772.62],
                             [422660.35,870777.17],
                             [422664.22,870785.46]],
                       '3': [[422661.30,871215.06],
                             [422667.50,871215.70],
                             [422668.30,871204.86],
                             [422662.21,871204.33],
                             [422661.30,871215.06]],
                       '4': [[422473.44,871191.22],
                             [422478.33,871192.26],
                             [422479.52,871186.03],
                             [422474.78,871185.14],
                             [422473.44,871191.22]],
                       '5': [[422369.69,871049.29],
                             [422378.63,871053.58],
                             [422383.91,871044.51],
                             [422374.97,871040.32],
                             [422369.69,871049.29]],
                       '8': [[422730.56,871203.13],
                             [422734.10,871204.90],
                             [422735.26,871202.18],
                             [422731.87,871200.58],
                             [422730.56,871203.13]],
                       '9': [[422659.85,871213.80],
                             [422660.91,871210.97],
                             [422655.42,871208.85],
                             [422654.36,871211.68],
                             [422659.85,871213.80]]
                       }        
        

        
                
        for id in ['1', '2', '3', '4', '5' ,'8' ,'9']:
            assert id in polygons.keys()
            assert id in values.keys()            

            assert int(values[id]) == int(floors[id])
            assert len(polygons[id]) == len(known_polys[id])
            assert num.allclose(polygons[id], known_polys[id])


    def test_csv2polygons_with_clipping(self):
        """test_csv2polygons with optional clipping
        """
        #FIXME(Ole): Not Done!!
        
        path = get_pathname_from_package('anuga.shallow_water')                
        testfile = os.path.join(path, 'tests', 'data', 'polygon_values_example.csv')                

        polygons, values = load_csv_as_polygons(testfile, 
                                        value_name='floors',
                                        clipping_polygons=None)

        assert len(polygons) == 7, 'Must have seven polygons'
        assert len(values) == 7, 'Must have seven values'
            
        # Known floor values
        floors = {'1': 2, '2': 0, '3': 1, '4': 2, '5': 0, '8': 1, '9': 1}
        
        # Known polygon values
        known_polys = {'1': [[422681.61,871117.55],
                             [422691.02,871117.60],
                             [422690.87,871084.23],
                             [422649.36,871081.85],
                             [422649.36,871080.39],
                             [422631.86,871079.50],
                             [422631.72,871086.75],
                             [422636.75,871087.20],
                             [422636.75,871091.50],
                             [422649.66,871092.09],
                             [422649.83,871084.91],
                             [422652.94,871084.90],
                             [422652.84,871092.39],
                             [422681.83,871093.73],
                             [422681.61,871117.55]],
                       '2': [[422664.22,870785.46],
                             [422672.48,870780.14],
                             [422668.17,870772.62],
                             [422660.35,870777.17],
                             [422664.22,870785.46]],
                       '3': [[422661.30,871215.06],
                             [422667.50,871215.70],
                             [422668.30,871204.86],
                             [422662.21,871204.33],
                             [422661.30,871215.06]],
                       '4': [[422473.44,871191.22],
                             [422478.33,871192.26],
                             [422479.52,871186.03],
                             [422474.78,871185.14],
                             [422473.44,871191.22]],
                       '5': [[422369.69,871049.29],
                             [422378.63,871053.58],
                             [422383.91,871044.51],
                             [422374.97,871040.32],
                             [422369.69,871049.29]],
                       '8': [[422730.56,871203.13],
                             [422734.10,871204.90],
                             [422735.26,871202.18],
                             [422731.87,871200.58],
                             [422730.56,871203.13]],
                       '9': [[422659.85,871213.80],
                             [422660.91,871210.97],
                             [422655.42,871208.85],
                             [422654.36,871211.68],
                             [422659.85,871213.80]]
                       }        
        

        
                
        for id in ['1', '2', '3', '4', '5' ,'8' ,'9']:
            assert id in polygons.keys()
            assert id in values.keys()            

            assert int(values[id]) == int(floors[id])
            assert len(polygons[id]) == len(known_polys[id])
            assert num.allclose(polygons[id], known_polys[id])




    
    def test_csv2building_polygons(self):
        """test_csv2building_polygons
        """
        
        path = get_pathname_from_package('anuga.shallow_water')                
        testfile = os.path.join(path, 'tests', 'data', 'polygon_values_example.csv')                

        polygons, values = load_csv_as_building_polygons(testfile, 
                                                 floor_height=3)

        assert len(polygons) == 7, 'Must have seven polygons'
        assert len(values) == 7, 'Must have seven values'
            
        # Known floor values
        floors = {'1': 6, '2': 0, '3': 3, '4': 6, '5': 0, '8': 3, '9': 3}
        
                
        for id in ['1', '2', '3', '4', '5' ,'8' ,'9']:
            assert id in polygons.keys()
            assert id in values.keys()            

            assert float(values[id]) == float(floors[id])



#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_csv, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
