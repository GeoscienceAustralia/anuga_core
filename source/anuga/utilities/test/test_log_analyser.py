#!/usr/bin/env python
#

# Warning, this will not update as the log format is updated.

import os
import unittest
import tempfile
from anuga.utilities.log import TimingDelimiter
from anuga.utilities.log_analyser import analyse_log

class logTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def log_lines(self, key, value):
        
        line = '2012-01-23 13:20:10,147 INFO   general_mesh:202 |'

        #FIXME this logic is in the code somewhere.  So it 
        # shouldn't be repeated here.
        end_line = TimingDelimiter + ' ' + key + ', ' + \
            value + '\n' 
        return line + end_line

    def test_log(self):
        # create a dummy directory and log file
        root_dir = tempfile.mkdtemp('test_log_analyser')
        dir1 = tempfile.mkdtemp(dir=root_dir)
        dir2 = tempfile.mkdtemp(dir=root_dir)
        
        log_file_name = 'anuga.log'
        
        # Create a fake log file
        log_path_file1 = os.path.join(dir1, log_file_name)
        handle = file(log_path_file1, 'w')
        handle.write('yeah\n yeah\n ')
        handle.write('2012-01-23 13:20:10,147 INFO   general_mesh:202 |\n')
        #made_up = {'three':'3', 'five':'5', 'seven':'7'}
        made_up1 = {'numTriangles':'1', 'startMeshTime':'2', 'finishMeshTime':'3', \
                       'finishTime':'4', 'startMemory':'5', 'finishMemory':'6'}
        for key, val in made_up1.iteritems():
            handle.write(self.log_lines(key, val))
        handle.close()
        
        # Create another fake log file
        log_path_file2 = os.path.join(dir2, log_file_name)
        handle = file(log_path_file2, 'w')
        handle.write('yeah\n yeah\n ')
        handle.write('2012-01-23 13:20:10,147 INFO   general_mesh:202 |\n')
        #made_up = {'three':'3', 'five':'5', 'seven':'7'}
        made_up2 = {'another':'1', 'startMeshTime':'2', 'finishMeshTime':'3', \
                       'finishTime':'4', 'startMemory':'5', 'finishMemory':'6'}
        for key, val in made_up2.iteritems():
            handle.write(self.log_lines(key, val))
        handle.close()

        # output file
        (handle, output_file) = tempfile.mkstemp('.csv','get_bridges_from_dic_')
        os.close(handle)
        output_file = 'yeah.csv'

        analyse_log(root_dir, output_file, log_file_name)
        
        #FIXME check these results!!
        
        
        os.remove(log_path_file1)
        os.remove(log_path_file2)
        os.rmdir(dir1)
        os.rmdir(dir2)
        os.rmdir(root_dir)
        os.remove(output_file)
        #os.remove(log_file_name)




################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(logTestCase,'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)
    
