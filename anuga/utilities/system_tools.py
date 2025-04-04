"""Implementation of tools to do with system administration made as platform independent as possible.
"""

import sys
import os
from six.moves import urllib
#import urllib.request, urllib.parse, urllib.error
import getpass
import tarfile
import warnings
import platform
import pdb
from functools import reduce

# Record Python version
major_version = int(platform.python_version_tuple()[0])
version = platform.python_version()

try:
    import hashlib
except ImportError:
    import md5 as hashlib

   
def log_to_file(filename, s, verbose=False, mode='a'):
    """Log string to file name
    """

    fid = open(filename, mode)
    if verbose: print(s)
    fid.write(s + '\n')
    fid.close()


def get_user_name():
    """Get user name provide by operating system
    """

    import getpass
    user = getpass.getuser()
    
    return user    

def get_host_name():
    """Get host name provide by operating system
    """

    if sys.platform == 'win32':
        host = os.getenv('COMPUTERNAME')
    else:
        host = os.uname()[1]


    return host    

    
def get_version():
    """Get anuga version number as stored in anuga.__version__
    """

    import anuga
    return anuga.__version__

    
def get_revision_number():
    """Get the (git) sha of this repository copy.
    """
    from anuga import __git_sha__ as revision
    return revision
    

def get_revision_date():
    """Get the (git) revision date of this repository copy.
    """

    from anuga import __git_committed_datetime__ as revision_date
    return revision_date 
  
   
# FIXME (Ole): We should rewrite this to use GIT revision information
# And then update get_revision_number and date to use this if it can't 
# be obtained directly from GIT.
def store_svn_revision_info(destination_path='.', verbose=False):
    """Obtain current revision from Subversion and store it.
    
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
        # FIXME (Ole): Use git module here via get_revision_number/date
        txt = subprocess.Popen('svn info', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
    except:
        txt = 'Revision: 0'
    else:    
        if verbose: print('response ',txt)

        # Determine absolute filename
        if destination_path[-1] != os.sep:
            destination_path += os.sep
            
        filename = destination_path + config.revision_filename

        fid = open(filename, 'w')

        docstring = 'Stored revision info.\n\n'
        docstring += 'This file provides the version for distributions '
        docstring += 'that are not accessing Subversion directly.\n'
        docstring += 'The file is automatically generated and should not '
        docstring += 'be modified manually.\n'
        fid.write('"""%s"""\n\n' %docstring)
        
        fid.write('revision_info = """\n%s"""' %txt)
        fid.close()


        if verbose is True:
            print('Revision info stored to %s' % filename)



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


def get_anuga_pathname():
    """Get pathname of anuga install location 

    Typically, this is required in unit tests depending
    on external files.

    """
    
    import anuga
    import os
    
    return os.path.dirname(anuga.__file__)

    
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

    # Execute import command
    # See https://stackoverflow.com/questions/1463306/how-does-exec-work-with-locals
    exec('import %s as x' % package, globals())

    # # Get and return path
    # return x.__path__[0]
    import os
    return os.path.dirname(x.__file__)



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
    """Convert 1-D list of strings to 2-D list of chars."""

    if not l:
        return []

    if l == ['']:
        l = [' ']


    maxlen = reduce(max, list(map(len, l)))
    ll = [x.ljust(maxlen) for x in l]
    result = []
    for s in ll:
        result.append([x for x in s])
    return result



def char_to_string(ll):
    """Convert 2-D list of chars to 1-D list of strings."""

    # https://stackoverflow.com/questions/23618218/numpy-bytes-to-plain-string
    # bytes_string.decode('UTF-8')

    # We might be able to do this a bit more shorthand as we did in Python2.x
    # i.e return [''.join(x).strip() for x in ll]

    # But this works for now.

    result = []
    for i in range(len(ll)):
        x = ll[i]
        string = ''
        for j in range(len(x)):
            c = x[j]
            if type(c) == str:
                string += c
            else:
                string += c.decode()            

        result.append(string.strip())
        
    return result


################################################################################

def get_vars_in_expression(source):
    """Get list of variable names in a python expression."""

    # https://stackoverflow.com/questions/37993137/how-do-i-detect-variables-in-a-python-eval-expression
    
    import ast
        
    variables = {}
    syntax_tree = ast.parse(source)
    for node in ast.walk(syntax_tree):
        if type(node) is ast.Name:
            variables[node.id] = 0  # Keep first one, but not duplicates
                
    # Only return keys
    result = list(variables.keys()) # Only return keys i.e. the variable names
    result.sort() # Sort for uniqueness
    return result


def get_web_file(file_url, file_name, auth=None, blocksize=1024*1024):
    """Get a file from the web (HTTP).

    file_url:  The URL of the file to get
    file_name: Local path to save loaded file in
    auth:      A tuple (httpproxy, proxyuser, proxypass)
    blocksize: Block size of file reads
    
    Will try simple load through urllib first.  Drop down to urllib2
    if there is a proxy and it requires authentication.

    Environment variable HTTP_PROXY can be used to supply proxy information.
    PROXY_USERNAME is used to supply the authentication username.
    PROXY_PASSWORD supplies the password, if you dare!
    """

    # Simple fetch, if fails, check for proxy error
    try:
        urllib.request.urlretrieve(file_url, file_name)
        return (True, auth)     # no proxy, no auth required
    except IOError as e:
        if e[1] == 407:     # proxy error
            pass
        elif e[1][0] == 113:  # no route to host
            print('No route to host for %s' % file_url)
            return (False, auth)    # return False
        else:
            print('Unknown connection error to %s' % file_url)
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
        print('-'*72)
        print ('You need to supply proxy authentication information.')
        if httpproxy is None:
            httpproxy = input('                    proxy server: ')
        else:
            print('         HTTP proxy was supplied: %s' % httpproxy)
        if proxyuser is None:
            proxyuser = input('                  proxy username: ') 
        else:
            print('HTTP proxy username was supplied: %s' % proxyuser)
        if proxypass is None:
            proxypass = getpass.getpass('                  proxy password: ')
        else:
            print('HTTP proxy password was supplied: %s' % '*'*len(proxyuser))
        print('-'*72)

    # the proxy URL cannot start with 'http://', we add that later
    httpproxy = httpproxy.lower()
    if httpproxy.startswith('http://'):
        httpproxy = httpproxy.replace('http://', '', 1)

    # open remote file
    proxy = urllib.request.ProxyHandler({'http': 'http://' + proxyuser
                                              + ':' + proxypass
                                              + '@' + httpproxy})
    authinfo = urllib.request.HTTPBasicAuthHandler()
    opener = urllib.request.build_opener(proxy, authinfo, urllib.request.HTTPHandler)
    urllib.request.install_opener(opener)
    try:
        webget = urllib.request.urlopen(file_url)
    except urllib.error.HTTPError as e:
        print('Error received from proxy:\n%s' % str(e))
        print('Possibly the user/password is wrong.')
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
    """Compress a file or directory into a tar file."""

    if isinstance(files, str):
        files = [files]

    o = tarfile.open(tarname, 'w:gz')
    for file in files:
        o.add(file)
    o.close()


def untar_file(tarname, target_dir='.'):
    """Uncompress a tar file."""
    

    o = tarfile.open(tarname, 'r:gz')
    members = o.getmembers()
    for member in members:
        o.extract(member, target_dir, filter='data')
    o.close()


def get_file_hexdigest(filename, blocksize=1024*1024*10):
    """Get a hex digest of a file."""

    if hashlib.__name__ == 'hashlib':
        m = hashlib.md5()       # new - 'hashlib' module
    else:
        m = hashlib.new()       # old - 'md5' module - remove once py2.4 gone
    fd = open(filename, 'r')
            
    while True:
        data = fd.read(blocksize)
        if len(data) == 0:
            break
        m.update(data.encode())
                                                                
    fd.close()
    return m.hexdigest()


def make_digest_file(data_file, digest_file):
    """Create a file containing the hex digest string of a data file."""
    
    hexdigest = get_file_hexdigest(data_file)
    fd = open(digest_file, 'w')
    fd.write(hexdigest)
    fd.close()


def file_length(in_file):
    """Function to return the length of a file."""

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
    """private method"""
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
    return (float(v[1]) * _scale[v[2]]) // (1024.0 * 1024.0)


def MemoryUpdate(print_msg=None,str_return=False):
    """print memory usage stats in MB.
    """
    global _total_memory, _last_memory

    _last_memory = _total_memory
    _total_memory = _VmB('VmSize:')

    #if print_msg is not None:
    mem_diff = _total_memory - _last_memory
    return (mem_diff,_total_memory)
    

