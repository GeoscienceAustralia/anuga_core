"""Implementation of tools to do with system administration made as platform independent as possible.


"""

import sys
import os

def log_to_file(filename, s, verbose=False):
    """Log string to file name
    """

    fid = open(filename, 'a')
    if verbose: print s
    fid.write(s + '\n')
    fid.close()


def get_user_name():
    """Get user name provide by operating system
    """

    if sys.platform == 'win32':
        #user = os.getenv('USERPROFILE')
        user = os.getenv('USERNAME')
    else:
        user = os.getenv('LOGNAME')


    return user    

def get_host_name():
    """Get host name provide by operating system
    """

    if sys.platform == 'win32':
        host = os.getenv('COMPUTERNAME')
    else:
        host = os.uname()[1]


    return host    

def get_revision_number():
    """Get the version number of the SVN
    NOTE: This requires that the command svn is on the system PATH
    (simply aliasing svn to the binary will not work)
    """

    # Create dummy info 
    #info = 'Revision: Version info could not be obtained.'
    #info += 'A command line version of svn must be availbable '
    #info += 'on the system PATH, access to the subversion '
    #info += 'repository is necessary and the output must '
    #info += 'contain a line starting with "Revision:"'
    

    #FIXME (Ole): Change this so that svn info is attempted first.
    # If that fails, try to read a stored file with that same info (this would be created by e.g. the release script). Failing that, throw an exception.

    #FIXME (Ole): Move this and store_version_info to utilities


    try:
        from anuga.stored_version_info import version_info
    except:
	msg = 'No version info stored and command "svn" is not '
	msg += 'recognised on the system PATH.\n\n'
	msg += 'If ANUGA has been installed from a distribution e.g. as '
	msg += 'obtained from SourceForge,\n'
	msg += 'the version info should be '
	msg += 'available in the automatically generated file '
	msg += 'stored_version_info.py\n'
	msg += 'in the anuga root directory.\n'
	msg += 'If run from a Subversion sandpit, '
	msg += 'ANUGA will try to obtain the version info '
	msg += 'by using the command: "svn info".\n'
	msg += 'In this case, make sure svn is accessible on the system path. '
	msg += 'Simply aliasing svn to the binary will not work. '
	msg += 'Good luck!'

        # No file available - try using Subversion
        try:
            # The null stuff is so this section fails quitly.
            # This could cause the svn info command to fail due to
            # the redirection being bad on some platforms.
            # If that occurs then change this code.
            if sys.platform[0:3] == 'win':
                fid = os.popen('svn info 2> null')
            else:
                fid = os.popen('svn info 2>/dev/null')
	
        except:
            raise Exception(msg)
        else:
            #print 'Got version from svn'            
            version_info = fid.read()
	    
	    if version_info == '':
	        raise Exception(msg)    
    else:
        pass
        #print 'Got version from file'

            
    for line in version_info.split('\n'):
        if line.startswith('Revision:'):
            break

    fields = line.split(':')
    msg = 'Keyword "Revision" was not found anywhere in text: %s' %version_info
    assert fields[0].startswith('Revision'), msg            

    try:
        revision_number = int(fields[1])
    except:
        msg = 'Revision number must be an integer. I got %s' %fields[1]
        msg += 'Check that the command svn is on the system path' 
        raise Exception(msg)                
        
    return revision_number


def store_version_info(destination_path='.', verbose=False):
    """Obtain current version from Subversion and store it.
    
    Title: store_version_info()

    Author: Ole Nielsen (Ole.Nielsen@ga.gov.au)

    CreationDate: January 2006

    Description:
        This function obtains current version from Subversion and stores it
        is a Python file named 'stored_version_info.py' for use with
        get_version_info()

        If svn is not available on the system PATH, an Exception is thrown
    """

    # Note (Ole): This function should not be unit tested as it will only
    # work when running out of the sandpit. End users downloading the
    # ANUGA distribution would see a failure.
    #
    # FIXME: This function should really only be used by developers (
    # (e.g. for creating new ANUGA releases), so maybe it should move
    # to somewhere else.
    
    import config

    try:
        fid = os.popen('svn info')
    except:
        msg = 'Command "svn" is not recognised on the system PATH'
        raise Exception(msg)
    else:    
        txt = fid.read()
        fid.close()


        # Determine absolute filename
        if destination_path[-1] != os.sep:
            destination_path += os.sep
            
        filename = destination_path + config.version_filename

        fid = open(filename, 'w')

        docstring = 'Stored version info.\n\n'
        docstring += 'This file provides the version for distributions '
        docstring += 'that are not accessing Subversion directly.\n'
        docstring += 'The file is automatically generated and should not '
        docstring += 'be modified manually.\n'
        fid.write('"""%s"""\n\n' %docstring)
        
        fid.write('version_info = """\n%s"""' %txt)
        fid.close()


        if verbose is True:
            print 'Version info stored to %s' %filename

def safe_crc(string):
    """64 bit safe crc computation.

       See Guido's 64 bit fix at http://bugs.python.org/issue1202            
    """

    from zlib import crc32
    import os

    x = crc32(string)
        
    if os.name == 'posix' and os.uname()[4] == 'x86_64':
        crcval = x - ((x & 0x80000000) << 1)
    else:
        crcval = x
        
    return crcval


def compute_checksum(filename, max_length=2**20):
    """Compute the CRC32 checksum for specified file

    Optional parameter max_length sets the maximum number
    of bytes used to limit time used with large files.
    Default = 2**20 (1MB)
    """

    fid = open(filename, 'rb') # Use binary for portability
    crcval = safe_crc(fid.read(max_length))
    fid.close()

    return crcval

def get_pathname_from_package(package):
    """Get pathname of given package (provided as string)

    This is useful for reading files residing in the same directory as
    a particular module. Typically, this is required in unit tests depending
    on external files.

    The given module must start from a directory on the pythonpath
    and be importable using the import statement.

    Example
    path = get_pathname_from_package('anuga.utilities')

    """

    exec('import %s as x' %package)

    path = x.__path__[0]
    
    return path

    # Alternative approach that has been used at times
    #try:
    #    # When unit test is run from current dir
    #    p1 = read_polygon('mainland_only.csv')
    #except: 
    #    # When unit test is run from ANUGA root dir
    #    from os.path import join, split
    #    dir, tail = split(__file__)
    #    path = join(dir, 'mainland_only.csv')
    #    p1 = read_polygon(path)
        
            


