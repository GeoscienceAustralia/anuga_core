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


          
def load_csv_as_array(file_name, delimiter=','):
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

    Y = {}
    for key in X.keys():
        Y[key] = num.array([float(x) for x in X[key]])

    return Y

