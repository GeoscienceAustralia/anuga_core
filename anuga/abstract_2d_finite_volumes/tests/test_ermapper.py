#!/usr/bin/env python

import unittest


import anuga.abstract_2d_finite_volumes.ermapper_grids as ermapper_grids
from os import remove

import numpy as num


class Test_ERMapper(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_write_grid(self):
        header_filename = 'test_write_ermapper_grid.ers'
        data_filename = 'test_write_ermapper_grid'
        
        original_grid = num.array([[0.0, 0.1, 1.0], [2.0, 3.0, 4.0]])

        # Check that the function works when passing the filename without
        # a '.ers' extension
        ermapper_grids.write_ermapper_grid(data_filename, original_grid)
        new_grid = ermapper_grids.read_ermapper_grid(data_filename)
        assert num.allclose(original_grid, new_grid)

        # Check that the function works when passing the filename with
        # a '.ers' extension
        ermapper_grids.write_ermapper_grid(header_filename, original_grid)
        new_grid = ermapper_grids.read_ermapper_grid(header_filename)
        assert num.allclose(original_grid, new_grid)

        # Clean up created files
        remove(data_filename)
        remove(header_filename)
        
    def test_basic_single_line_grid(self):
        # Setup test data
        filename = 'test_write_ermapper_grid'
        original_grid = num.array([0.0, 0.1, 1.0, 2.0, 3.0, 4.0])

        # Write test data
        ermapper_grids.write_ermapper_data(original_grid, filename, num.float64)

        # Read in the test data
        new_grid = ermapper_grids.read_ermapper_data(filename, num.float64)

        # Check that the test data that has been read in matches the original data
        assert num.allclose(original_grid, new_grid)

        # Clean up created files
        remove(filename)
        
    def test_basic_single_line_grid_default_format(self):
        # Setup test data
        filename = 'test_write_ermapper_grid'
        original_grid = num.array([0.0, 0.1, 1.0, 2.0, 3.0, 4.0])

        # Write test data
        ermapper_grids.write_ermapper_data(original_grid, filename)

        # Read in the test data
        new_grid = ermapper_grids.read_ermapper_data(filename)

        # Check that the test data that has been read in matches the original data
        assert num.allclose(original_grid, new_grid)

        # Clean up created files
        remove(filename)
        
    def test_write_default_header(self):
        data_filename = 'test_write_ermapper_grid'

        # setup test data
        original_grid = num.array([[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]])
        # Write test data
        ermapper_grids.write_ermapper_data(original_grid, data_filename)
        # Write test header using all default values
        header_filename = data_filename + '.ers'                              
        ermapper_grids.write_ermapper_header(header_filename)

        # Check that the read in values match the default values
        header = ermapper_grids.read_ermapper_header(header_filename)

        assert header['datum'] == '"GDA94"'
        assert header['projection'] == '"GEOGRAPHIC"'  
        assert header['coordinatetype'] == 'LL'
        assert header['rotation'] == '0:0:0.0'
        assert header['celltype'] == 'IEEE4ByteReal'
        assert header['nullcellvalue'] == '-99999'
        assert header['xdimension'] == '100'
        assert header['registrationcellx'] == '0'
        assert header['ydimension'] == '100'
        assert header['registrationcelly'] == '2'
        assert header['nroflines'] == '3'
        assert header['nrofcellsperline'] == '4'
        assert header['longitude'] == '0:0:0'
        assert header['latitude'] == '0:0:0'
        assert header['nrofbands'] == '1'
        assert header['value'] == '"Default_Band"'

        # Clean up created files
        remove(data_filename)
        remove(header_filename)
        
    def test_header_creation(self):
        header = {}
        # have some values that aren't defaults
        header['nroflines'] = '2'
        header['nrofcellsperline'] = '3'

        header = ermapper_grids.create_default_header(header)


        # default values
        assert header['datum'] == '"GDA94"'
        assert header['projection'] == '"GEOGRAPHIC"'  
        assert header['coordinatetype'] == 'LL'
        assert header['rotation'] == '0:0:0.0'
        assert header['celltype'] == 'IEEE4ByteReal'
        assert header['nullcellvalue'] == '-99999'
        assert header['xdimension'] == '100'
        assert header['registrationcellx'] == '0'
        assert header['ydimension'] == '100'
        assert header['longitude'] == '0:0:0'
        assert header['latitude'] == '0:0:0'
        assert header['nrofbands'] == '1'
        assert header['value'] == '"Default_Band"'

        # non-default values
        assert header['nroflines'] == '2'
        assert header['nrofcellsperline'] == '3'
        assert header['registrationcelly'] == '1'


    def test_write_non_default_header(self):
        header_filename = 'test_write_ermapper_grid.ers'

        # setup test data
        header = {}
        # have some values that aren't defaults
        header['nroflines'] = '2'
        header['nrofcellsperline'] = '3'
        
        # Write test header using non-default values                           
        ermapper_grids.write_ermapper_header(header_filename, header)

        # Check that the read in values match the default values
        header = ermapper_grids.read_ermapper_header(header_filename)

        # default values
        assert header['datum'] == '"GDA94"'
        assert header['projection'] == '"GEOGRAPHIC"'  
        assert header['coordinatetype'] == 'LL'
        assert header['rotation'] == '0:0:0.0'
        assert header['celltype'] == 'IEEE4ByteReal'
        assert header['nullcellvalue'] == '-99999'
        assert header['xdimension'] == '100'
        assert header['registrationcellx'] == '0'
        assert header['ydimension'] == '100'
        assert header['longitude'] == '0:0:0'
        assert header['latitude'] == '0:0:0'
        assert header['nrofbands'] == '1'
        assert header['value'] == '"Default_Band"'

        # non-default values
        assert header['nroflines'] == '2'
        assert header['nrofcellsperline'] == '3'
        assert header['registrationcelly'] == '1'

        # Clean up created files
        remove(header_filename)        


# def test_default_filenames
# def test_write_header
# def test_multi_band_grid

       
#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_ERMapper,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
