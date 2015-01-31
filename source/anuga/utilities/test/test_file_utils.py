import unittest
import tempfile
import os
import shutil
import sys

from anuga.utilities.file_utils import copy_code_files, get_all_swwfiles
from anuga.utilities.file_utils import del_dir
from anuga.utilities.sww_merge import sww_merge, _sww_merge


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
        self.assertTrue(os.access(os.path.join(dst_dir, f1), os.F_OK))
        self.assertTrue(os.access(os.path.join(dst_dir, f2), os.F_OK))
        self.assertTrue(os.access(os.path.join(dst_dir, f3), os.F_OK))
        self.assertTrue(os.access(os.path.join(dst_dir, f4), os.F_OK))
        self.assertTrue(os.access(os.path.join(dst_dir, f5), os.F_OK))

        # clean up
        shutil.rmtree(work_dir)
        
    def test_get_all_swwfiles(self):
        try:
            swwfiles = get_all_swwfiles('','test.txt')  #Invalid
        except IOError:
            pass
        else:
            raise Exception('Should have raised exception')
        
    def test_get_all_swwfiles1(self):
        
        temp_dir = tempfile.mkdtemp('','sww_test')
        filename0 = tempfile.mktemp('.sww','test',temp_dir)
        filename1 = tempfile.mktemp('.sww','test',temp_dir)
        filename2 = tempfile.mktemp('.sww','test',temp_dir)
        filename3 = tempfile.mktemp('.sww','test',temp_dir)
       
        #print'filename', filename0,filename1,filename2,filename3
        
        fid0 = open(filename0, 'w')
        fid1 = open(filename1, 'w')
        fid2 = open(filename2, 'w')
        fid3 = open(filename3, 'w')

        fid0.write('hello')
        fid1.write('hello')
        fid2.write('hello')
        fid3.write('hello')
        
        fid0.close()
        fid1.close()
        fid2.close()
        fid3.close()
        
        
        dir, name0 = os.path.split(filename0)
        #print 'dir',dir,name0
        
        iterate=get_all_swwfiles(dir,'test')
        
        del_dir(temp_dir)
#        removeall(temp_dir)

        _, name0 = os.path.split(filename0) 
        #print'name0',name0[:-4],iterate[0]    
        _, name1 = os.path.split(filename1)       
        _, name2 = os.path.split(filename2)       
        _, name3 = os.path.split(filename3)       

        assert name0[:-4] in iterate
        assert name1[:-4] in iterate
        assert name2[:-4] in iterate
        assert name3[:-4] in iterate
        
        assert len(iterate)==4            

    def test_merge_swwfiles(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular, \
                                                                    rectangular_cross
        from anuga.shallow_water.shallow_water_domain import Domain
        from anuga.file.sww import SWW_file
        from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import \
            Dirichlet_boundary

        Bd = Dirichlet_boundary([0.5, 0., 0.])

        # Create shallow water domain
        domain = Domain(*rectangular_cross(2, 2))
        domain.set_name('test1')
        domain.set_quantity('elevation', 2)
        domain.set_quantity('stage', 5)
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})
        for t in domain.evolve(yieldstep=0.5, finaltime=1):
            pass
            
        domain = Domain(*rectangular(3, 3))
        domain.set_name('test2')
        domain.set_quantity('elevation', 3)
        domain.set_quantity('stage', 50)
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})
        for t in domain.evolve(yieldstep=0.5, finaltime=1):
            pass
                
        outfile = 'test_out.sww'
        _sww_merge(['test1.sww', 'test2.sww'], outfile)
        self.assertTrue(os.access(outfile, os.F_OK))  
        
        # remove temp files
        if not sys.platform == 'win32':		
			os.remove('test1.sww')
			os.remove('test2.sww')
			os.remove(outfile)      
        
        

#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_FileUtils, 'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)    
