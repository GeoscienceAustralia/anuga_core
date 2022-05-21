""" Modules for performing conversions between file types, or for
    resampling a given file. 
    
    In general, the naming convention follows this rule:
        <file_ext_in>2<file_ext_out>('filename_in', 'filename_out')
        
    for example:
        sww2dem('northbeach.sww', 'outfile.dem')
        
    Some formats input and output across multiple files. In that case the
    convention is so:
        urs2nc('northbeach', 'outfile')
    Where an array of 3 input files produce 4 output files.
"""

from numpy._pytesttester import PytestTester
test = PytestTester(__name__)
del PytestTester


