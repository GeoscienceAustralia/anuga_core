import unittest
import tempfile
import os
import shutil

from anuga.utilities.file_utils import copy_code_files


class Test_FileUtils(unittest.TestCase):
                
    def test_copy_code_files(self):
        '''test that the copy_code_files() function is sane.'''

        def create_file(f):
            fd = open(f, 'w')
            fd.write('%s\n' % f)
            fd.close()

        # create working directories and test files
        work_dir = tempfile.mkdtemp()
        dst_dir = tempfile.mkdtemp(dir=work_dir)
        src_dir = tempfile.mkdtemp(dir=work_dir)

        f1 = 'file1'        
        filename1 = os.path.join(src_dir, f1)
        create_file(filename1)
        f2 = 'file2'        
        filename2 = os.path.join(src_dir, f2)
        create_file(filename2)
        f3 = 'file3'        
        filename3 = os.path.join(src_dir, f3)
        create_file(filename3)
        f4 = 'file4'        
        filename4 = os.path.join(src_dir, f4)
        create_file(filename4)
        f5 = 'file5'        
        filename5 = os.path.join(src_dir, f5)
        create_file(filename5)

        # exercise the copy function
        copy_code_files(dst_dir, filename1)
        copy_code_files(dst_dir, filename1, filename2)
        copy_code_files(dst_dir, (filename4, filename5, filename3))

        # test that files were actually copied
        self.failUnless(os.access(os.path.join(dst_dir, f1), os.F_OK))
        self.failUnless(os.access(os.path.join(dst_dir, f2), os.F_OK))
        self.failUnless(os.access(os.path.join(dst_dir, f3), os.F_OK))
        self.failUnless(os.access(os.path.join(dst_dir, f4), os.F_OK))
        self.failUnless(os.access(os.path.join(dst_dir, f5), os.F_OK))

        # clean up
        shutil.rmtree(work_dir)
            

#-------------------------------------------------------------

if __name__ == "__main__":
    #suite = unittest.makeSuite(Test_Data_Manager, 'test_sww2domain2')
    suite = unittest.makeSuite(Test_FileUtils, 'test_sww')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)    
