#!/usr/bin/env python

'''Routine to add a header to a CSV file.
'''

import os
import tempfile


##
# @brief Add a header line to a text CSV file.
# @param file Path to the file to add header to.
# @param header_list A list of strings - the header to insert.
# @param be_green If True, go easy on memory, but slower.
# @note Checks that the first line of the original file has the same
#       number of fields as the new header line.
# @note Raises exception if above is not true.
# @note The 'be_green' option is not yet implemented.
def add_csv_header(file, header_list, be_green=False):
    '''Add a CSV header line to a text file.'''

    # open the file to insert header into
    try:
        fd = open(file, 'r')
    except IOError, e:
        msg = "Can't open file '%s': %s" % (file, str(e))
        raise Exception, msg
    except:
        msg = "Can't open file '%s'" % file
        raise Exception, msg

    # get a temporary file.
    # must create in same directory as input file, as we rename it.
    input_dir = os.path.dirname(file)
    (tmp_f, tmp_filename) = tempfile.mkstemp(dir=input_dir)
    tmp_fd = os.fdopen(tmp_f, 'w')

    # check the header, create header _string.
    if header_list is None or header_list is []:
        msg = "add_csv_header: header_list can't be None or []"
        raise Exception, msg

    header_string = ','.join(header_list) + '\n'

    # copy header to output file, then input file
    tmp_fd.write(header_string)
    data = fd.readlines()
    columns = data[0].strip().split(',')

    if len(columns) != len(header_list):
        msg = ("add_csv_header: File %s has %d columns but header "
               "has %d columns" % (file, len(columns), len(header_list)))
        raise Exception, msg

    data = ''.join(data)
    tmp_fd.write(data)

    # close and rename all files
    tmp_fd.close()
    fd.close()
    os.rename(tmp_filename, file)


if __name__ == '__main__':
    import sys

    file_data = '1,2,3\n4,5,6\n7,8,9'
    header = ['alpha', 'bravo', 'charlie']
    filename = '/tmp/add_csv_header.csv'

    # create test file
    fd = open(filename, 'w')
    fd.write(file_data)
    fd.close()

    add_csv_header(filename, header)

    # read test file
    fd = open(filename, 'r')
    data = fd.readlines()
    fd.close

    # check if data as expected
    data = ''.join(data)
    header_string = ','.join(header) + '\n'
    expected = header_string + file_data

    msg = 'Expected data:\n%s\ngot:\n%s' % (expected, data)
    assert expected == data, msg

