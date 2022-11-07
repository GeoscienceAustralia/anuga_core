""" Generic file utilities for creating, parsing deleting
    and naming files in a manner consistent across ANUGA.
"""

import os
import sys
import csv
import numpy as num
import shutil
from . import log


def make_filename(s):
    """Transform argument string into a standard filename

        Convert a possible filename into a standard form.
        s Filename to process.
        The new filename string.
    """

    s = s.strip()
    s = s.replace(' ', '_')
    s = s.replace('(', '')
    s = s.replace(')', '')
    s = s.replace('__', '_')

    return s


def check_dir(path, verbose=None):
    """Check that specified path exists.
    If path does not exist it will be created if possible

    USAGE:
       checkdir(path, verbose):

    ARGUMENTS:
        path -- Directory
        verbose -- Flag verbose output (default: None)

    RETURN VALUE:
        Verified path including trailing separator
    """

    import os.path

    if sys.platform in ['nt', 'dos', 'win32', 'what else?']:
        unix = 0
    else:
        unix = 1

    # add terminal separator, if it's not already there
    if path[-1] != os.sep:
        path = path + os.sep

    # expand ~ or ~username in path
    path = os.path.expanduser(path)

    # create directory if required
    if not (os.access(path, os.R_OK and os.W_OK) or path == ''):
        try:
            exitcode = os.mkdir(path)

            # Change access rights if possible
            if unix:
                exitcode = os.system('chmod 775 ' + path)
            else:
                pass  # FIXME: What about access rights under Windows?

            if verbose:
                log.critical('MESSAGE: Directory %s created.' % path)
        except:
            log.critical('WARNING: Directory %s could not be created.' % path)
            if unix:
                try:
                    path = os.environ['TMPDIR']
                except KeyError:
                    path = '/tmp/'
            else:
                path = 'C:' + os.sep

            log.critical("Using directory '%s' instead" % path)

    return path


def del_dir(path):
    """Recursively delete directory path and all its contents
    """

    if os.path.isdir(path):
        for file in os.listdir(path):
            X = os.path.join(path, file)

            if os.path.isdir(X) and not os.path.islink(X):
                del_dir(X)
            else:
                try:
                    os.remove(X)
                except:
                    log.critical("Could not remove file %s" % X)

        os.rmdir(path)


def rmgeneric(path, func, verbose=False):
    ERROR_STR = """Error removing %(path)s, %(error)s """

    try:
        func(path)
        if verbose:
            log.critical('Removed %s' % path)
    except OSError as xxx_todo_changeme:
        (errno, strerror) = xxx_todo_changeme.args
        log.critical(ERROR_STR % {'path': path, 'error': strerror})


def removeall(path, verbose=False):
    if not os.path.isdir(path):
        return

    for x in os.listdir(path):
        fullpath = os.path.join(path, x)
        if os.path.isfile(fullpath):
            f = os.remove
            rmgeneric(fullpath, f)
        elif os.path.isdir(fullpath):
            removeall(fullpath)
            f = os.rmdir
            rmgeneric(fullpath, f, verbose)


def create_filename(datadir, filename, format, size=None, time=None):
    """Create a standard filename.

    datadir   directory where file is to be created
    filename  filename 'stem'
    format    format of the file, becomes filename extension
    size      size of file, becomes part of filename
    time      time (float), becomes part of filename

    Returns the complete filename path, including directory.

    The containing directory is created, if necessary.
    """

    FN = check_dir(datadir) + filename

    if size is not None:
        FN += '_size%d' % size

    if time is not None:
        FN += '_time%.2f' % time

    FN += '.' + format

    return FN


def get_files(datadir, filename, format, size):
    """Get all file (names) with given name, size and format
    """

    import glob

    dir = check_dir(datadir)
    pattern = dir + os.sep + filename + '_size=%d*.%s' % (size, format)

    return glob.glob(pattern)


def get_all_directories_with_name(look_in_dir='', base_name='', verbose=False):
    '''
    Finds all the directories in a "look_in_dir" which contains a "base_name".

    Returns: a list of strings

    Usage:     iterate_over = get_all_directories_with_name(dir, name)
    then:      for swwfile in iterate_over:
                   do stuff

    Check "export_grids" and "get_maximum_inundation_data" for examples
    '''

    if look_in_dir == "":
        look_in_dir = "."                                  # Unix compatibility

    dir_ls = os.listdir(look_in_dir)
    iterate_over = [x for x in dir_ls if base_name in x]

    if len(iterate_over) == 0:
        msg = 'No files of the base name %s' % base_name
        raise IOError(msg)

    if verbose:
        log.critical('iterate over %s' % iterate_over)

    return iterate_over


def get_all_swwfiles(look_in_dir='', base_name='', verbose=False):
    '''
    Finds all the sww files in a "look_in_dir" which contains a "base_name".
    will accept base_name with or without the extension ".sww"

    Returns: a list of strings

    Usage:     iterate_over = get_all_swwfiles(dir, name)
    then
               for swwfile in iterate_over:
                   do stuff

    Check "export_grids" and "get_maximum_inundation_data" for examples
    '''

    # plus tests the extension
    name, extension = os.path.splitext(base_name)

    if extension != '' and extension != '.sww':
        msg = 'file %s%s must be a NetCDF sww file!' % (base_name, extension)
        raise IOError(msg)

    if look_in_dir == "":
        look_in_dir = "."                                   # Unix compatibility

    dir_ls = os.listdir(look_in_dir)
    iterate_over = [x[:-4] for x in dir_ls if name in x and x[-4:] == '.sww']
    if len(iterate_over) == 0:
        msg = 'No files of the base name %s' % name
        raise IOError(msg)

    if verbose:
        log.critical('iterate over %s' % iterate_over)

    return iterate_over


def get_all_files_with_extension(look_in_dir='',
                                 base_name='',
                                 extension='.sww',
                                 verbose=False):
    '''Find all files in a directory with given stem name.
    Finds all the sww files in a "look_in_dir" which contains a "base_name".

    Returns: a list of strings

    Usage:     iterate_over = get_all_swwfiles(dir, name)
    then
               for swwfile in iterate_over:
                   do stuff

    Check "export_grids" and "get_maximum_inundation_data" for examples
    '''

    # plus tests the extension
    name, ext = os.path.splitext(base_name)

    if ext != '' and ext != extension:
        msg = 'base_name %s must be a file with %s extension!' \
              % (base_name, extension)
        raise IOError(msg)

    if look_in_dir == "":
        look_in_dir = "."                               # Unix compatibility

    dir_ls = os.listdir(look_in_dir)
    iterate_over = [x[:-4]
                    for x in dir_ls if name in x and x[-4:] == extension]

    if len(iterate_over) == 0:
        msg = 'No files of the base name %s in %s' % (name, look_in_dir)
        raise IOError(msg)

    if verbose:
        log.critical('iterate over %s' % iterate_over)

    return iterate_over


def copy_code_files(dir_name, filename1, filename2=None, verbose=False):
    """Copies "filename1" and "filename2" to "dir_name".

    Each 'filename' may be a string or list of filename strings.

    Filenames must be absolute pathnames
    """

    def copy_file_or_sequence(dest, file):

        if isinstance(file, str):
            shutil.copy(file, dir_name)
            if verbose:
                log.critical('File %s copied' % file)
        elif isinstance(file, (tuple, list)):
            for f in file:
                copy_file_or_sequence(dest, f)
        else:
            raise Exception('Unknow argument for file: %s', file)
        


    # check we have a destination directory, create if necessary
    if not os.path.isdir(dir_name):
        if verbose:
            log.critical('Make directory %s' % dir_name)
        os.mkdir(dir_name, 0o777)

    if verbose:
        log.critical('Output directory: %s' % dir_name)

    copy_file_or_sequence(dir_name, filename1)

    if not filename2 is None:
        copy_file_or_sequence(dir_name, filename2)
