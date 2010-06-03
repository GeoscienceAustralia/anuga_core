"""
    A set of functions which extend the capabilities of the Python csv
    module.
    
    These have been left as functions to aviod confusion with the standard
    csv module.
"""

import csv
import numpy as num


def load_csv_as_dict(file_name, title_check_list=None, delimiter=','):
    """
    Load in the csv as a dictionary, title as key and column info as value.
    Also, create a dictionary, title as key and column index as value,
    to keep track of the column order.

    file_name The path to the file to read.
    
    title_check_list List of titles that *must* be columns in the file.

    delimiter is the delimiter used to separate the fields

    return 2 dictionaries: ({key:column}, {title:index}).

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

    reader = csv.reader(file(file_name), delimiter=delimiter)

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
        if n < title_count:
            msg = 'Entry in file %s had %d columns ' % (file_name, n)
            msg += 'although there were %d headers' % title_count
            raise IOError, msg
        for i, value in enumerate(line[:title_count]):  # skip trailing data
            attribute_dic.setdefault(titles_stripped[i], []).append(value)

    return attribute_dic, title_index_dic


          
def load_csv_as_array(file_name, delimiter = ','):
    """
    Convert CSV files of the form:

    time, discharge, velocity
    0.0,  1.2,       0.0
    0.1,  3.2,       1.1
    ...

    to a dictionary of numeric arrays.

    file_name The path to the file to read.
    delimiter is the delimiter used to separate the fields    

    See underlying function load_csv_as_dict for more details.
    """

    X, _ = load_csv_as_dict(file_name, delimiter=delimiter)


    # Return result as a dict of arrays
    ret = {}
    for key in X.keys():
        ret[key] = num.array([float(x) for x in X[key]])
            
    return ret


def load_csv_as_matrix(file_name, delimiter = ','):
    """
    Convert CSV files of the form:

    time, discharge, velocity
    0.0,  1.2,       0.0
    0.1,  3.2,       1.1
    ...

    to a numeric matrix.

    file_name The path to the file to read.
    delimiter is the delimiter used to separate the fields    

    See underlying function load_csv_as_dict for more details.
    """

    X, title_indices = load_csv_as_dict(file_name, delimiter=delimiter)

    col_titles = title_indices.keys()

    # Return result as a 2D array
    ret = num.zeros((len(X[col_titles[0]]), len(title_indices)), float)

    header = []
    for col_title in col_titles:
        index = title_indices[col_title]
        header.append(col_title)
        for i, x in enumerate(X[col_title]):
            ret[i, index] = float(x)

    return header, ret



##
# @brief Store keyword params into a CSV file.
# @param verbose True if this function is to be verbose.
# @param kwargs Dictionary of keyword args to store.
# @note If kwargs dict contains 'file_name' key, that has the output filename.
#       If not, make up a filename in the output directory.
def store_parameters(verbose=False, **kwargs):
    """
    Store "kwargs" into a temp csv file, if "completed" is in kwargs,
    csv file is kwargs[file_name] else it is kwargs[output_dir]+details_temp.csv

    Must have a file_name keyword arg, this is what is writing to.
    might be a better way to do this using CSV module Writer and writeDict.

    writes file to "output_dir" unless "completed" is in kwargs, then
    it writes to "file_name" kwargs
    """

    import types

    # Check that kwargs is a dictionary
    if type(kwargs) != types.DictType:
        raise TypeError

    # is 'completed' in kwargs?
    completed = kwargs.has_key('completed')

    # get file name and removes from dict and assert that a file_name exists
    if completed:
        try:
            file = str(kwargs['file_name'])
        except:
            raise 'kwargs must have file_name'
    else:
        # write temp file in output directory
        try:
            file = str(kwargs['output_dir']) + 'detail_temp.csv'
        except:
            raise 'kwargs must have output_dir'

    # extracts the header info and the new line info
    line = ''
    header = ''
    count = 0
    keys = kwargs.keys()
    keys.sort()

    # used the sorted keys to create the header and line data
    for k in keys:
        header += str(k)
        line += str(kwargs[k])
        count += 1
        if count < len(kwargs):
            header += ','
            line += ','
    header += '\n'
    line += '\n'

    # checks the header info, if the same, then write, if not create a new file
    # try to open!
    try:
        fid = open(file, 'r')
        file_header = fid.readline()
        fid.close()
        if verbose: log.critical('read file header %s' % file_header)
    except:
        msg = 'try to create new file: %s' % file
        if verbose: log.critical(msg)
        #tries to open file, maybe directory is bad
        try:
            fid = open(file, 'w')
            fid.write(header)
            fid.close()
            file_header=header
        except:
            msg = 'cannot create new file: %s' % file
            raise Exception, msg

    # if header is same or this is a new file
    if file_header == str(header):
        fid = open(file, 'a')
        fid.write(line)
        fid.close()
    else:
        # backup plan,
        # if header is different and has completed will append info to
        # end of details_temp.cvs file in output directory
        file = str(kwargs['output_dir']) + 'detail_temp.csv'
        fid = open(file, 'a')
        fid.write(header)
        fid.write(line)
        fid.close()

        if verbose:
            log.critical('file %s', file_header.strip('\n'))
            log.critical('head %s', header.strip('\n'))
        if file_header.strip('\n') == str(header):
            log.critical('they equal')

        msg = 'WARNING: File header does not match input info, ' \
              'the input variables have changed, suggest you change file name'
        log.critical(msg)

