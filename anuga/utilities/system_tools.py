"""Implementation of tools to do with system administration made as platform independent as possible.


"""

import sys
import os
import string
import urllib
import urllib2
import getpass
import tarfile
import warnings
import pdb

try:
    import hashlib
except ImportError:
    import md5 as hashlib

import anuga.utilities.log as log


def log_to_file(filename, s, verbose=False, mode='a'):
    """Log string to file name
    """

    fid = open(filename, mode)
    if verbose: print s
    fid.write(s + '\n')
    fid.close()


def get_user_name():
    """Get user name provide by operating system
    """

    import getpass
    user = getpass.getuser()

    #if sys.platform == 'win32':
    #    #user = os.getenv('USERPROFILE')
    #    user = os.getenv('USERNAME')
    #else:
    #    user = os.getenv('LOGNAME')


    return user    

def get_host_name():
    """Get host name provide by operating system
    """

    if sys.platform == 'win32':
        host = os.getenv('COMPUTERNAME')
    else:
        host = os.uname()[1]


    return host    

    
    
    
    
def __get_revision_from_svn_entries__():
    """Get a subversion revision number from the .svn/entries file."""

    
    msg = '''
No version info stored and command 'svn' is not recognised on the system PATH.

If ANUGA has been installed from a distribution e.g. as obtained from SourceForge,
the version info should be available in the automatically generated file
'stored_version_info.py' in the anuga root directory.

If run from a Subversion sandpit, ANUGA will try to obtain the version info by
using the command 'svn info'.  In this case, make sure the command line client
'svn' is accessible on the system path.  Simply aliasing 'svn' to the binary will
not work.

If you are using Windows, you have to install the file svn.exe which can be
obtained from http://www.collab.net/downloads/subversion.

Good luck!
'''

    try:
        fd = open(os.path.join('.svn', 'entries'))
    except:
        #raise Exception, msg


        #FIXME SR: Need to fix this up
        # svn 1.7 no longer has a .svn folder in all folders
        # so will need a better way to get revision number
        
        from anuga.utilities.stored_version_info import version_info
        return process_version_info(version_info)

    line = fd.readlines()[3]
    fd.close()
    try:
        revision_number = int(line)
    except:
        msg = ".svn/entries, line 4 was '%s'?" % line.strip()
        raise Exception, msg

    return revision_number

def __get_revision_from_svn_client__():
    """Get a subversion revision number from an svn client."""

    import subprocess

    if sys.platform[0:3] == 'win':
        #print 'On Win'
        try:
            #FIXME SR: This works for python 2.6
            cmd = r'"C:\Program Files\TortoiseSVN\bin\SubWCRev.exe" .'
            version_info = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]

            #print 'Version_Info', version_info
            #fid = os.popen(r'C:\Program Files\TortoiseSVN\bin\SubWCRev.exe')
        except:
            return __get_revision_from_svn_entries__()
        else:
            #version_info = fid.read()
            if version_info == '':
                return __get_revision_from_svn_entries__()


        # split revision number from data
        for line in version_info.split('\n'):
            if line.startswith('Updated to revision '):
                break
            if line.startswith('Last committed at revision'):
                break

        #print line
        fields = line.split(' ')
        msg = 'Keyword "Revision" was not found anywhere in text: %s' % version_info
        assert fields[0].startswith('Updated')  or fields[0].startswith('Last'), msg


        try:
            if fields[0].startswith('Updated'):
                revision_number = int(fields[3])
            if fields[0].startswith('Last'):
                revision_number = int(fields[4])
        except:
            msg = ('Revision number must be an integer. I got "%s" from '
                   '"SubWCRev.exe".' % line)
            raise Exception, msg
    else:                   # assume Linux
        try:
            fid = os.popen('svn info . 2>/dev/null')
        except:
            return __get_revision_from_svn_entries__()
        else:
            version_info = fid.read()
            if version_info == '':
                return __get_revision_from_svn_entries__()

        # Split revision number from data
        for line in version_info.split('\n'):
            if line.startswith('Revision:'):
                break
        fields = line.split(':')
        msg = 'Keyword "Revision" was not found anywhere in text: %s' % version_info
        assert fields[0].startswith('Revision'), msg
        
        
        #if ':' in version_info:
        #    revision_number, _ = version_info.split(':', 1)
        #    msg = ('Some modules have not been checked in. '
        #           'Using last version from repository: %s' % revision_number)
        #    warnings.warn(msg)
        #else:
        #    revision_number = version_info

        try:
            revision_number = int(fields[1])
        except:
            msg = ("Revision number must be an integer. I got '%s' from "
                   "'svn'." % fields[1])
            raise Exception, msg

    return revision_number

    
    
    
    
def get_revision_number():
    """Get the version number of this repository copy.

    Try getting data from stored_version_info.py first, otherwise
    try using SubWCRev.exe (Windows) or svnversion (linux), otherwise
    try reading file .svn/entries for version information, otherwise
    throw an exception.

    NOTE: This requires that the command svn is on the system PATH
    (simply aliasing svn to the binary will not work)
    """

    # try to get revision information from stored_version_info.py
    try:
        return __get_revision_from_svn_client__()
    except:
        try:
            from anuga.stored_version_info import version_info
            return process_version_info(version_info)
        except:
            from anuga.version import version
            return version

def process_version_info(version_info):

    # split revision number from data
    for line in version_info.split('\n'):
        if line.startswith('Revision:'):
            break

    fields = line.split(':')
    msg = 'Keyword "Revision" was not found anywhere in text: %s' % version_info
    assert fields[0].startswith('Revision'), msg

    try:
        revision_number = int(fields[1])
    except:
        msg = ("Revision number must be an integer. I got '%s'.\n"
               'Check that the command svn is on the system path.'
               % fields[1])
        raise Exception, msg

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
    
    import anuga.config as config
    import subprocess

    #txt = subprocess.Popen('svn info', shell=True, stdout=subprocess.PIPE).communicate()[0]
    try:
        #fid = os.popen('svn info')
        #FIXME SR: This works for python 2.6
        txt = subprocess.Popen('svn info', shell=True, stdout=subprocess.PIPE).communicate()[0]
    except:
        msg = 'Command "svn" is not recognised on the system PATH'
        raise Exception(msg)
    else:    
        #txt = fid.read()
        #fid.close()

        if verbose: print 'response ',txt


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
            log.critical('Version info stored to %s' % filename)


def safe_crc(string):
    """64 bit safe crc computation.

    See http://docs.python.org/library/zlib.html#zlib.crc32:

        To generate the same numeric value across all Python versions 
        and platforms use crc32(data) & 0xffffffff.
    """

    from zlib import crc32

    return crc32(string) & 0xffffffff


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
        
    
def clean_line(str, delimiter):
    """Split a string into 'clean' fields.

    str        the string to process
    delimiter  the delimiter string to split 'line' with

    Returns a list of 'cleaned' field strings.

    Any fields that were initially zero length will be removed.
    If a field contains '\n' it isn't zero length.
    """

    return [x.strip() for x in str.strip().split(delimiter) if x != '']


################################################################################
# The following two functions are used to get around a problem with numpy and
# NetCDF files.  Previously, using Numeric, we could take a list of strings and
# convert to a Numeric array resulting in this:
#     Numeric.array(['abc', 'xy']) -> [['a', 'b', 'c'],
#                                      ['x', 'y', ' ']]
#
# However, under numpy we get:
#     numpy.array(['abc', 'xy']) -> ['abc',
#                                    'xy']
#
# And writing *strings* to a NetCDF file is problematic.
#
# The solution is to use these two routines to convert a 1-D list of strings
# to the 2-D list of chars form and back.  The 2-D form can be written to a 
# NetCDF file as before.
#
# The other option, of inverting a list of tag strings into a dictionary with
# keys being the unique tag strings and the key value a list of indices of where
# the tag string was in the original list was rejected because:
#    1. It's a lot of work
#    2. We'd have to rewite the I/O code a bit (extra variables instead of one)
#    3. The code below is fast enough in an I/O scenario
################################################################################

def string_to_char(l):
    '''Convert 1-D list of strings to 2-D list of chars.'''

    if not l:
        return []

    if l == ['']:
        l = [' ']

    maxlen = reduce(max, map(len, l))
    ll = [x.ljust(maxlen) for x in l]
    result = []
    for s in ll:
        result.append([x for x in s])
    return result


def char_to_string(ll):
    '''Convert 2-D list of chars to 1-D list of strings.'''

    return map(string.rstrip, [''.join(x) for x in ll])

################################################################################

def get_vars_in_expression(source):
    '''Get list of variable names in a python expression.'''

    import compiler
    from compiler.ast import Node

    def get_vars_body(node, var_list=[]):
        if isinstance(node, Node):
            if node.__class__.__name__ == 'Name':
                for child in node.getChildren():
                    if child not in var_list:
                        var_list.append(child)
            for child in node.getChildren():
                if isinstance(child, Node):
                    for child in node.getChildren():
                        var_list = get_vars_body(child, var_list)
                    break

        return var_list

    return get_vars_body(compiler.parse(source))


def get_web_file(file_url, file_name, auth=None, blocksize=1024*1024):
    '''Get a file from the web (HTTP).

    file_url:  The URL of the file to get
    file_name: Local path to save loaded file in
    auth:      A tuple (httpproxy, proxyuser, proxypass)
    blocksize: Block size of file reads
    
    Will try simple load through urllib first.  Drop down to urllib2
    if there is a proxy and it requires authentication.

    Environment variable HTTP_PROXY can be used to supply proxy information.
    PROXY_USERNAME is used to supply the authentication username.
    PROXY_PASSWORD supplies the password, if you dare!
    '''

    # Simple fetch, if fails, check for proxy error
    try:
        urllib.urlretrieve(file_url, file_name)
        return (True, auth)     # no proxy, no auth required
    except IOError, e:
        if e[1] == 407:     # proxy error
            pass
        elif e[1][0] == 113:  # no route to host
            print 'No route to host for %s' % file_url
            return (False, auth)    # return False
        else:
            print 'Unknown connection error to %s' % file_url
            return (False, auth)

    # We get here if there was a proxy error, get file through the proxy
    # unpack auth info
    try:
        (httpproxy, proxyuser, proxypass) = auth
    except:
        (httpproxy, proxyuser, proxypass) = (None, None, None)

    # fill in any gaps from the environment
    if httpproxy is None:
        httpproxy = os.getenv('HTTP_PROXY')
    if proxyuser is None:
        proxyuser = os.getenv('PROXY_USERNAME')
    if proxypass is None:
        proxypass = os.getenv('PROXY_PASSWORD')

    # Get auth info from user if still not supplied
    if httpproxy is None or proxyuser is None or proxypass is None:
        print '-'*72
        print ('You need to supply proxy authentication information.')
        if httpproxy is None:
            httpproxy = raw_input('                    proxy server: ')
        else:
            print '         HTTP proxy was supplied: %s' % httpproxy
        if proxyuser is None:
            proxyuser = raw_input('                  proxy username: ') 
        else:
            print 'HTTP proxy username was supplied: %s' % proxyuser
        if proxypass is None:
            proxypass = getpass.getpass('                  proxy password: ')
        else:
            print 'HTTP proxy password was supplied: %s' % '*'*len(proxyuser)
        print '-'*72

    # the proxy URL cannot start with 'http://', we add that later
    httpproxy = httpproxy.lower()
    if httpproxy.startswith('http://'):
        httpproxy = httpproxy.replace('http://', '', 1)

    # open remote file
    proxy = urllib2.ProxyHandler({'http': 'http://' + proxyuser
                                              + ':' + proxypass
                                              + '@' + httpproxy})
    authinfo = urllib2.HTTPBasicAuthHandler()
    opener = urllib2.build_opener(proxy, authinfo, urllib2.HTTPHandler)
    urllib2.install_opener(opener)
    try:
        webget = urllib2.urlopen(file_url)
    except urllib2.HTTPError, e:
        print 'Error received from proxy:\n%s' % str(e)
        print 'Possibly the user/password is wrong.'
        return (False, (httpproxy, proxyuser, proxypass))

    # transfer file to local filesystem
    fd = open(file_name, 'wb')
    while True:
        data = webget.read(blocksize)
        if len(data) == 0:
            break
        fd.write(data)
    fd.close
    webget.close()

    # return successful auth info
    return (True, (httpproxy, proxyuser, proxypass))


def tar_file(files, tarname):
    '''Compress a file or directory into a tar file.'''

    if isinstance(files, basestring):
        files = [files]

    o = tarfile.open(tarname, 'w:gz')
    for file in files:
        o.add(file)
    o.close()


def untar_file(tarname, target_dir='.'):
    '''Uncompress a tar file.'''

    o = tarfile.open(tarname, 'r:gz')
    members = o.getmembers()
    for member in members:
        o.extract(member, target_dir)
    o.close()


def get_file_hexdigest(filename, blocksize=1024*1024*10):
    '''Get a hex digest of a file.'''

    if hashlib.__name__ == 'hashlib':
        m = hashlib.md5()       # new - 'hashlib' module
    else:
        m = hashlib.new()       # old - 'md5' module - remove once py2.4 gone
    fd = open(filename, 'r')
            
    while True:
        data = fd.read(blocksize)
        if len(data) == 0:
            break
        m.update(data)
                                                                
    fd.close()
    return m.hexdigest()


def make_digest_file(data_file, digest_file):
    '''Create a file containing the hex digest string of a data file.'''
    
    hexdigest = get_file_hexdigest(data_file)
    fd = open(digest_file, 'w')
    fd.write(hexdigest)
    fd.close()


def file_length(in_file):
    '''Function to return the length of a file.'''

    fid = open(in_file)
    data = fid.readlines()
    fid.close()
    return len(data)


#### Memory functions
_proc_status = '/proc/%d/status' % os.getpid()
_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}
_total_memory = 0.0
_last_memory = 0.0

def _VmB(VmKey):
    '''private method'''
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to MB
    return (float(v[1]) * _scale[v[2]]) / (1024.0*1024.0)


def MemoryUpdate(print_msg=None,str_return=False):
    '''print memory usage stats in MB.
    '''
    global _total_memory, _last_memory

    _last_memory = _total_memory
    _total_memory = _VmB('VmSize:')

    #if print_msg is not None:
    mem_diff = _total_memory - _last_memory
    return (mem_diff,_total_memory)
    

