""" Generic file utilities for creating, parsing deleting
    and naming files in a manner consistent across ANUGA.
"""


import os, sys
import csv
import numpy as num

##
# @brief Convert a possible filename into a standard form.
# @param s Filename to process.
# @return The new filename string.
def make_filename(s):
    """Transform argument string into a Sexsuitable filename
    """

    s = s.strip()
    s = s.replace(' ', '_')
    s = s.replace('(', '')
    s = s.replace(')', '')
    s = s.replace('__', '_')

    return s


##
# @brief Check that a specified filesystem directory path exists.
# @param path The dirstory path to check.
# @param verbose True if this function is to be verbose.
# @note If directory path doesn't exist, it will be created.
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

            if verbose: log.critical('MESSAGE: Directory %s created.' % path)
        except:
            log.critical('WARNING: Directory %s could not be created.' % path)
            if unix:
                path = '/tmp/'
            else:
                path = 'C:' + os.sep

            log.critical("Using directory '%s' instead" % path)

    return path


##
# @brief Delete directory and all sub-directories.
# @param path Path to the directory to delete.
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


##
# @brief ??
# @param path 
# @param __func__ 
# @param verbose True if this function is to be verbose.
# @note ANOTHER OPTION, IF NEED IN THE FUTURE, Nick B 7/2007
def rmgeneric(path, func, verbose=False):
    ERROR_STR= """Error removing %(path)s, %(error)s """

    try:
        func(path)
        if verbose: log.critical('Removed %s' % path)
    except OSError, (errno, strerror):
        log.critical(ERROR_STR % {'path' : path, 'error': strerror })


##
# @brief Remove directory and all sub-directories.
# @param path Filesystem path to directory to remove.
# @param verbose True if this function is to be verbose.
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


##
# @brief Create a standard filename.
# @param datadir Directory where file is to be created.
# @param filename Filename 'stem'.
# @param format Format of the file, becomes filename extension.
# @param size Size of file, becomes part of filename.
# @param time Time (float), becomes part of filename.
# @return The complete filename path, including directory.
# @note The containing directory is created, if necessary.
def create_filename(datadir, filename, format, size=None, time=None):
    FN = check_dir(datadir) + filename

    if size is not None:
        FN += '_size%d' % size

    if time is not None:
        FN += '_time%.2f' % time

    FN += '.' + format

    return FN


##
# @brief Get all files with a standard name and a given set of attributes.
# @param datadir Directory files must be in.
# @param filename Filename stem.
# @param format Filename extension.
# @param size Filename size.
# @return A list of fielnames (including directory) that match the attributes.
def get_files(datadir, filename, format, size):
    """Get all file (names) with given name, size and format
    """

    import glob

    dir = check_dir(datadir)
    pattern = dir + os.sep + filename + '_size=%d*.%s' % (size, format)

    return glob.glob(pattern)


##
# @brief Find all files in a directory that contain a given string.
# @param look_in_dir Path to the directory to look in.
# @param base_name String that files must contain.
# @param verbose True if this function is to be verbose.
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
        raise IOError, msg

    if verbose: log.critical('iterate over %s' % iterate_over)

    return iterate_over



##
# @brief Find all SWW files in a directory with given stem name.
# @param look_in_dir The directory to look in.
# @param base_name The file stem name.
# @param verbose True if this function is to be verbose.
# @return A list of found filename strings.
# @note Will accept 'base_name' with or without '.sww' extension.
# @note If no files found, raises IOError exception.
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
        raise IOError, msg

    if look_in_dir == "":
        look_in_dir = "."                                   # Unix compatibility

    dir_ls = os.listdir(look_in_dir)
    iterate_over = [x[:-4] for x in dir_ls if name in x and x[-4:] == '.sww']
    if len(iterate_over) == 0:
        msg = 'No files of the base name %s' % name
        raise IOError, msg

    if verbose: log.critical('iterate over %s' % iterate_over)

    return iterate_over


##
# @brief Find all files in a directory that contain a string and have extension.
# @param look_in_dir Path to the directory to look in.
# @param base_name Stem filename of the file(s) of interest.
# @param extension Extension of the files to look for.
# @param verbose True if this function is to be verbose.
# @return A list of found filename strings.
# @note If no files found, raises IOError exception.
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
        raise IOError, msg

    if look_in_dir == "":
        look_in_dir = "."                               # Unix compatibility

    dir_ls = os.listdir(look_in_dir)
    iterate_over = [x[:-4] for x in dir_ls if name in x and x[-4:] == extension]

    if len(iterate_over) == 0:
        msg = 'No files of the base name %s in %s' % (name, look_in_dir)
        raise IOError, msg

    if verbose: log.critical('iterate over %s' % iterate_over)

    return iterate_over


##
# @brief Read a CSV file and convert to a dictionary of {key: column}.
# @param file_name The path to the file to read.
# @param title_check_list List of titles that *must* be columns in the file.
# @return Two dicts: ({key:column}, {title:index}).
# @note WARNING: Values are returned as strings.
def load_csv_as_dict(file_name, title_check_list=None):
    """
    Load in the csv as a dictionary, title as key and column info as value.
    Also, create a dictionary, title as key and column index as value,
    to keep track of the column order.

    Two dictionaries are returned.

    WARNING: Values are returned as strings.
             Do this to change a list of strings to a list of floats
                 time = [float(x) for x in time]
    """

    # FIXME(Ole): Consider dealing with files without headers
    # FIXME(Ole): Consider a wrapper automatically converting text fields
    #             to the right type by trying for: int, float, string
    
    attribute_dic = {}
    title_index_dic = {}
    titles_stripped = [] # List of titles

    reader = csv.reader(file(file_name))

    # Read in and manipulate the title info
    titles = reader.next()
    for i, title in enumerate(titles):
        header = title.strip()
        titles_stripped.append(header)
        title_index_dic[header] = i
    title_count = len(titles_stripped)

    # Check required columns
    if title_check_list is not None:
        for title_check in title_check_list:
            if not title_index_dic.has_key(title_check):
                msg = 'Reading error. This row is not present %s' % title_check
                raise IOError, msg

    # Create a dictionary of column values, indexed by column title
    for line in reader:
        n = len(line) # Number of entries
        if n != title_count:
            msg = 'Entry in file %s had %d columns ' % (file_name, n)
            msg += 'although there were %d headers' % title_count
            raise IOError, msg
        for i, value in enumerate(line):
            attribute_dic.setdefault(titles_stripped[i], []).append(value)

    return attribute_dic, title_index_dic


          
##
# @brief Convert CSV file to a dictionary of arrays.
# @param file_name The path to the file to read.
def load_csv_as_array(file_name):
    """
    Convert CSV files of the form:

    time, discharge, velocity
    0.0,  1.2,       0.0
    0.1,  3.2,       1.1
    ...

    to a dictionary of numeric arrays.


    See underlying function load_csv_as_dict for more details.
    """

    X, _ = load_csv_as_dict(file_name)

    Y = {}
    for key in X.keys():
        Y[key] = num.array([float(x) for x in X[key]])

    return Y



def copy_code_files(dir_name, filename1, filename2, verbose=False):
    """Copies "filename1" and "filename2" to "dir_name".

    Each 'filename' may be a string or list of filename strings.

    Filenames must be absolute pathnames
    """

    ##
    # @brief copies a file or sequence to destination directory.
    # @param dest The destination directory to copy to.
    # @param file A filename string or sequence of filename strings.
    def copy_file_or_sequence(dest, file):
        if hasattr(file, '__iter__'):
            for f in file:
                shutil.copy(f, dir_name)
                if verbose:
                    log.critical('File %s copied' % f)
        else:
            shutil.copy(file, dir_name)
            if verbose:
                log.critical('File %s copied' % file)

    # check we have a destination directory, create if necessary
    if not os.path.isdir(dir_name):
        if verbose:
            log.critical('Make directory %s' % dir_name)
        mkdir(dir_name, 0777)

    if verbose:
        log.critical('Output directory: %s' % dir_name)        

    copy_file_or_sequence(dir_name, filename1)

    if not filename2 is None:
        copy_file_or_sequence(dir_name, filename2)
