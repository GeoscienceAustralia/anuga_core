import os
import unittest
import tempfile
import numpy as num

from csv_file import load_csv_as_array, load_csv_as_dict, store_parameters, \
                        load_csv_as_matrix

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
        



#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_csv, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
