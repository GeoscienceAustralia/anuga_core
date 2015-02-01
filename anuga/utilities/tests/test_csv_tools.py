#!/usr/bin/env python

import os
import unittest
import tempfile
import csv

import anuga.utilities.csv_tools as csv_tools


# this dictionary sets the column header string for
# column number modulo 4.
col_text_string = {0: 'col_%d',
                   1: ' col_%d',
                   2: 'col_%d ',
                   3: ' col_%d '
                  }

class Test_CSV_utils(unittest.TestCase):

    NUM_FILES = 10
    NUM_COLS = 6
    NUM_LINES = 10
    OUTPUT_FILE = 'test.csv'
    
    def setUp(self):
        # create temporary scratch directory
        self.tmp_dir = tempfile.mkdtemp()

        # create 4 test CSV files
        self.num_files = self.NUM_FILES
        self.filenames = []
        for i in range(self.NUM_FILES):
            self.filenames.append(tempfile.mktemp('.csv'))
        for (i, fn) in enumerate(self.filenames):
            fd = open(fn, 'w')
            csv_fd = csv.writer(fd)
            # write colums row
            columns = []
            for j in range(self.NUM_COLS):
                columns.append(col_text_string[j % 4] % j)
            csv_fd.writerow(columns)

            # write data rows
            for j in xrange(self.NUM_LINES):
                data = [j, j, '%d.%d' % (j, i)] + ['qwert']*(self.NUM_COLS-3)
                csv_fd.writerow(data)
            fd.close()


    def tearDown(self):
        for fn in self.filenames:
            try:
                os.remove(fn)
            except:
                pass
        try:
            os.remove(self.OUTPUT_FILE)
        except:
            pass


    def test_merge_one_file(self):
        """Test merging a single CSV file.
        
        This is the same as a two coluymn extract, with column rename.
        """

        file_title_list = [(self.filenames[0], 'test')]
        csv_tools.merge_csv_key_values(file_title_list, self.OUTPUT_FILE,
                                       key_col='col_0', data_col='col_3')

        expected = '''col_0,test
0,qwert
1,qwert
2,qwert
3,qwert
4,qwert
5,qwert
6,qwert
7,qwert
8,qwert
9,qwert
'''

        got = self.get_file_contents(self.OUTPUT_FILE)
        msg = ('Merging one file,\n'
               'expected file=\n'
               '--------------------\n'
               '%s'
               '--------------------\n'
               'got file=\n'
               '--------------------\n'
               '%s'
               '--------------------\n'
               % (expected, got))
        self.assertTrue(self.str_cmp(got, expected), msg)


    def test_merge_two_files(self):
        """Test merging two CSV files."""

        file_title_list = [(self.filenames[0], 'test0'),
                           (self.filenames[1], 'test1')]
        csv_tools.merge_csv_key_values(file_title_list, self.OUTPUT_FILE,
                                       key_col='col_0', data_col='col_3')

        expected = '''col_0,test0,test1
0,qwert,qwert
1,qwert,qwert
2,qwert,qwert
3,qwert,qwert
4,qwert,qwert
5,qwert,qwert
6,qwert,qwert
7,qwert,qwert
8,qwert,qwert
9,qwert,qwert
'''

        got = self.get_file_contents(self.OUTPUT_FILE)
        msg = ('Merging two files,\n'
               'expected file=\n'
               '--------------------\n'
               '%s'
               '--------------------\n'
               'got file=\n'
               '--------------------\n'
               '%s'
               '--------------------\n'
               % (expected, got))
        self.assertTrue(self.str_cmp(got, expected), msg)


    def test_merge_two_files2(self):
        """Test merging two CSV files."""

        file_title_list = [(self.filenames[0], 'test0'),
                           (self.filenames[1], 'test1')]
        csv_tools.merge_csv_key_values(file_title_list, self.OUTPUT_FILE,
                                       key_col='col_0', data_col='col_2')

        expected = '''col_0,test0,test1
0,0.0,0.1
1,1.0,1.1
2,2.0,2.1
3,3.0,3.1
4,4.0,4.1
5,5.0,5.1
6,6.0,6.1
7,7.0,7.1
8,8.0,8.1
9,9.0,9.1
'''

        got = self.get_file_contents(self.OUTPUT_FILE)
        msg = ('Merging two file,\n'
               'expected file=\n'
               '--------------------\n'
               '%s'
               '--------------------\n'
               'got file=\n'
               '--------------------\n'
               '%s'
               '--------------------\n'
               % (expected, got))
        self.assertTrue(self.str_cmp(got, expected), msg)


    def test_merge_four_files(self):
        """Test merging four CSV files."""

        file_title_list = [(self.filenames[0], 'test0'),
                           (self.filenames[1], 'test1'),
                           (self.filenames[2], 'test2'),
                           (self.filenames[3], 'test3')]
        csv_tools.merge_csv_key_values(file_title_list, self.OUTPUT_FILE,
                                       key_col='col_0', data_col='col_2')

        expected = '''col_0,test0,test1,test2,test3
0,0.0,0.1,0.2,0.3
1,1.0,1.1,1.2,1.3
2,2.0,2.1,2.2,2.3
3,3.0,3.1,3.2,3.3
4,4.0,4.1,4.2,4.3
5,5.0,5.1,5.2,5.3
6,6.0,6.1,6.2,6.3
7,7.0,7.1,7.2,7.3
8,8.0,8.1,8.2,8.3
9,9.0,9.1,9.2,9.3
'''

        got = self.get_file_contents(self.OUTPUT_FILE)
        msg = ('Merging four files,\n'
               'expected file=\n'
               '--------------------\n'
               '%s'
               '--------------------\n'
               'got file=\n'
               '--------------------\n'
               '%s'
               '--------------------\n'
               % (expected, got))
        self.assertTrue(self.str_cmp(got, expected), msg)


    def test_merge_ten_files(self):
        """Test merging ten CSV files."""

        file_title_list = [(self.filenames[0], 'test0'),
                           (self.filenames[1], 'test1'),
                           (self.filenames[2], 'test2'),
                           (self.filenames[3], 'test3'),
                           (self.filenames[4], 'test4'),
                           (self.filenames[5], 'test5'),
                           (self.filenames[6], 'test6'),
                           (self.filenames[7], 'test7'),
                           (self.filenames[8], 'test8'),
                           (self.filenames[9], 'test9')]
        csv_tools.merge_csv_key_values(file_title_list, self.OUTPUT_FILE,
                                       key_col='col_1', data_col='col_2')

        expected = '''col_1,test0,test1,test2,test3,test4,test5,test6,test7,test8,test9
0,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
1,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9
2,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9
3,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9
4,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9
5,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9
6,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9
7,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9
8,8.0,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9
9,9.0,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9
'''

        got = self.get_file_contents(self.OUTPUT_FILE)
        msg = ('Merging four files,\n'
               'expected file=\n'
               '--------------------\n'
               '%s'
               '--------------------\n'
               'got file=\n'
               '--------------------\n'
               '%s'
               '--------------------\n'
               % (expected, got))
        self.assertTrue(self.str_cmp(got, expected), msg)


    def test_no_key_column(self):
        """Test merging two CSV files with expected missing key column."""

        file_title_list = [(self.filenames[0], 'test0'),
                           (self.filenames[2], 'test2')]
        self.assertRaises(Exception,
                              csv_tools.merge_csv_key_values,
                              file_title_list,
                              self.OUTPUT_FILE,
                              key_col='col_A',
                              data_col='col_2'
                             )


    def test_no_input_files(self):
        """Test merging *zero* CSV files!"""

        file_title_list = []
        self.assertRaises(Exception,
                              csv_tools.merge_csv_key_values,
                              file_title_list,
                              self.OUTPUT_FILE,
                              key_col='col_1',
                              data_col='col_A'
                             )


    def test_no_data_column(self):
        """Test merging two CSV files with expected missing data column."""

        file_title_list = [(self.filenames[0], 'test0'),
                           (self.filenames[2], 'test2')]
        self.assertRaises(Exception,
                              csv_tools.merge_csv_key_values,
                              file_title_list,
                              self.OUTPUT_FILE,
                              key_col='col_1',
                              data_col='col_A'
                             )


    def test_different_num_rows(self):
        """Test merging two CSV files with different number of rows."""

        # get data from file [1]
        fd = open(self.filenames[1], 'r')
        data = fd.readlines()
        fd.close()

        # delete a row in data and write to test file
        test_filename = 'my_test.csv'
        fd = open(test_filename, 'w')
        fd.write(''.join(data[0:-1]))
        fd.close()

        file_title_list = [(self.filenames[0], 'test0'),
                           (test_filename, 'test2')]
        self.assertRaises(Exception,
                              csv_tools.merge_csv_key_values,
                              file_title_list,
                              self.OUTPUT_FILE,
                              key_col='col_1',
                              data_col='col_A'
                             )

        try:
            os.remove(test_filename)
        except:
            pass


    def test_different_key_values(self):
        """Test merging two CSV files with different key values."""

        # get data from file [1]
        fd = open(self.filenames[1], 'r')
        data = fd.readlines()
        fd.close()

        # chnage a row key value in data and write to test file
        test_filename = 'my_test.csv'
        fd = open(test_filename, 'w')
        data[3] = '1' + data[3]
        fd.write(''.join(data))
        fd.close()

        file_title_list = [(self.filenames[0], 'test0'),
                           (test_filename, 'test2')]
        self.assertRaises(Exception,
                              csv_tools.merge_csv_key_values,
                              file_title_list,
                              self.OUTPUT_FILE,
                              key_col='col_1',
                              data_col='col_A'
                             )

        try:
            os.remove(test_filename)
        except:
            pass


    def test_latex_example(self):
        """Test merging two CSV files - example from latex doc."""

        fd = open('alpha.csv', 'w')
        csv_fd = csv.writer(fd)
        csv_fd.writerow(['time', 'hours', 'stage', 'depth'])
        csv_fd.writerow(['3600', '1.00', '100.3', '10.2'])
        csv_fd.writerow(['3636', '1.01', '100.3', '10.0'])
        csv_fd.writerow(['3672', '1.02', '100.3', '9.7'])
        csv_fd.writerow(['3708', '1.03', '100.3', '8.9'])
        csv_fd.writerow(['3744', '1.04', '100.3', '7.1'])
        fd.close()

        fd = open('beta.csv', 'w')
        csv_fd = csv.writer(fd)
        csv_fd.writerow(['time', 'hours', 'stage', 'depth'])
        csv_fd.writerow(['3600', '1.00', '100.3', '11.3'])
        csv_fd.writerow(['3636', '1.01', '100.3', '10.5'])
        csv_fd.writerow(['3672', '1.02', '100.3', '10.0'])
        csv_fd.writerow(['3708', '1.03', '100.3', '9.7'])
        csv_fd.writerow(['3744', '1.04', '100.3', '8.2'])
        fd.close()

        file_title_list = [('alpha.csv', 'alpha'),
                           ('beta.csv',  'beta')]
        csv_tools.merge_csv_key_values(file_title_list,
                                       'gamma.csv',
                                       key_col='hours',
                                       data_col='depth')

        expected = '''hours,alpha,beta
1.00,10.2,11.3
1.01,10.0,10.5
1.02,9.7,10.0
1.03,8.9,9.7
1.04,7.1,8.2
'''

        got = self.get_file_contents('gamma.csv')
        msg = ('Merging two files,\n'
               'expected file=\n'
               '--------------------\n'
               '%s'
               '--------------------\n'
               'got file=\n'
               '--------------------\n'
               '%s'
               '--------------------\n'
               % (expected, got))
        self.assertTrue(self.str_cmp(got, expected), msg)

        try:
            os.remove('alpha.csv')
            os.remove('beta.csv')
            os.remove('gamma.csv')
        except:
            pass


    def str_cmp(self, str1, str2):
        '''Compare 2 strings, removing end-of-line stuff first.'''

        s1 = str1.split('\n')
        s2 = str2.split('\n')
        for (sub1, sub2) in zip(s1, s2):
            if sub1 != sub2:
                return False
        return True


    def get_file_contents(self, filename):
        '''Return file contents as a string.'''

        fd = open(filename, 'r')
        data = fd.readlines()
        fd.close()
        return ''.join(data).replace('\r', '')

################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_CSV_utils, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
