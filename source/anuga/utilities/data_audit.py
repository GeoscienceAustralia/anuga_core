"""Track IP of data files included in this distribution. 
"""

from os import remove, walk, sep
from os.path import join, splitext


def IP_verified(directory):
    """Find and audit potential data files that might violate IP

    This is the public function to be used to ascertain that
    all data in the specified directory tree has been audited according
    to the GA data IP tracking process.

    if IP_verified is False:
        # Stop and take remedial action
        ...
    else:
        # Proceed boldly with confidence
        

    """

    print '---------------------------------------------'
    print 'Files that need to be assessed for IP issues:'
    print '---------------------------------------------'

    # Print header
    dirwidth = 72
    print '---------------------------------------------'
    print 'File'.ljust(dirwidth), 'Status'
    print '---------------------------------------------'

    # Identify data files
    all_files_accounted_for = True
    for dirpath, datafile in identify_datafiles(directory):
        print join(dirpath, datafile) + ': ',

        basename, ext = splitext(datafile)

        # Look for a XML license file with the .lic
        try:
            fid = open(join(dirpath, basename + '.lic'))
        except IOError:
            print 'NO LICENSE FILE'
            all_files_accounted_for = False
        else:
            if license_file_is_valid(fid):
                print 'OK'
            else:
                print 'LICENSE FILE NOT VALID'
                all_files_accounted_for = False
            fid.close()

    # Return result        
    return all_files_accounted_for


def identify_datafiles(root):
    """ Identify files that might contain data
    """

    # Ignore source code files
    extensions_to_ignore = ['.py','.c','.h', '.f'] #, '.gif', '.jpg', '.png']

    # Ignore generated stuff 
    extensions_to_ignore += ['.pyc', '.o', '.so', '~']
    extensions_to_ignore += ['.aux', '.log', '.idx', 'ilg', '.ind',
                             '.bbl', '.blg']

    # Ignore license files themselves
    extensions_to_ignore += ['.lic']    
    

    # Ignore certain other files
    files_to_ignore = ['README.txt']

    # Ignore directories
    directories_to_ignore = ['anuga_work', 'pymetis', 'obsolete_code',
                             'anuga_parallel', 'anuga_viewer',
                             'planning', 'coding_standards',
                             'experimentation',
                             '.svn', 'misc', '.metadata']

    for dirpath, dirnames, filenames in walk(root):

        for ignore in directories_to_ignore:
            if ignore in dirnames:
                dirnames.remove(ignore)  # don't visit ignored directories

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
                yield dirpath, filename


def license_file_is_valid(fid):
    """Check that XML license file is valid
    """

    # TODO

    print fid.read()
    import sys
    from xml.dom import minidom
    doc =  minidom.parse(fid)
