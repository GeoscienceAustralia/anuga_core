"""Track IP of data files included in this distribution. 


"""

from os import remove, walk, sep
from os.path import join

def identify_data_files(distro_dir):
    """ Identify potential data files that might violate IP
    """

    print '---------------------------------------------'
    print 'Files that need to be assessed for IP issues:'
    print '---------------------------------------------'

    # Print header
    dirwidth = 72
    print '---------------------------------------------'
    print 'Directory'.ljust(dirwidth), 'File'
    print '---------------------------------------------'

    # Ignore source code files
    extensions_to_ignore = ['.py','.c','.h', '.f'] #,'gif']

    # Ignore certain other files
    files_to_ignore = ['README.txt']

    for dirpath, dirnames, filenames in walk(distro_dir):

        #print 'Searching dir', dirpath
   

        for filename in filenames:


            # Ignore extensions that need no IP check
            ignore = False
            for ext in extensions_to_ignore:
                if filename.endswith(ext):
                    ignore = True

            if filename in files_to_ignore:
                ignore = True

            if ignore is False:
                subdirs = dirpath.split(sep)
                
                print join(subdirs[3:],sep).ljust(dirwidth), filename
        


# FIXME (Ole): Here we could put in a check testing if
# all the files above have a .license file associated with them
# explaining their origins.

