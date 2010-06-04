""" Modules for performing conversions between file types, or for
    resampling a given file. 
    
    In general, the naming convention follows this rule:
        <file_ext_in>2<file_ext_out>('filename_in', 'filename_out')
        
    for example:
        sww2dem('northbeach.sww', 'outfile.dem')
"""

# Add path of package to PYTHONPATH to allow C-extensions to be loaded
import sys
sys.path += __path__

