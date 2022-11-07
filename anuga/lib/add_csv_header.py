#!/usr/bin/env python

'''Routine to add a header to a CSV file.
'''

import os
import tempfile


def add_csv_header(file, header_list, be_green=False):
    """Add a CSV header line to a text file.

    file         Path to the file to add header to
    header_list  A list of strings - the header to insert
    be_green     If True, go easy on memory, but slower

    Checks that the first line of the original file has the same number of
    fields as the new header line.  Raises exception if thi is not true.

    The 'be_green' option is not yet implemented.
    """

    # open the file to insert header into
    try:
        fd = open(file, 'r')
    except IOError as e:
        msg = "Can't open file '%s': %s" % (file, str(e))
        raise Exception(msg)
    except:
        msg = "Can't open file '%s'" % file
        raise Exception(msg)

    # get a temporary file.
    # must create in same directory as input file, as we rename it.
    input_dir = os.path.dirname(file)
    (tmp_f, tmp_filename) = tempfile.mkstemp(dir=input_dir)
    tmp_fd = os.fdopen(tmp_f, 'w')

    # check the header, create header _string.
    if header_list is None or header_list is []:
        msg = "add_csv_header: header_list can't be None or []"
        raise Exception(msg)

    header_string = ','.join(header_list) + '\n'

    # copy header to output file, then input file
    tmp_fd.write(header_string)

    if be_green:
        first_line = True
        for line in fd:
            if first_line:
                first_line = False
                columns = line.strip().split(',')
                if len(columns) != len(header_list):
                    msg = ("add_csv_header: File %s has %d columns but header "
                           "has %d columns" % (file, len(columns), len(header_list)))
                    raise Exception(msg)
            tmp_fd.write(line)
    else:
        data = fd.readlines()
        columns = data[0].strip().split(',')

        if len(columns) != len(header_list):
            msg = ("add_csv_header: File %s has %d columns but header "
                   "has %d columns" % (file, len(columns), len(header_list)))
            raise Exception(msg)

        data = ''.join(data)
        tmp_fd.write(data)

    # close and rename all files
    tmp_fd.close()
    fd.close()
    os.rename(tmp_filename, file)


if __name__ == '__main__':
    import sys
    import os

    file_data = '1,2,3\n4,5,6\n7,8,9'
    header = ['alpha', 'bravo', 'charlie']
    filename = '/tmp/add_csv_header.csv'
    filename2 = '/tmp/add_csv_header2.csv'

######
# Create file and test function.
######

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

######
# Test the 'be_green' option.
######

    # create test file
    fd = open(filename2, 'w')
    fd.write(file_data)
    fd.close()

    add_csv_header(filename2, header, be_green=True)

    # read test file
    fd = open(filename2, 'r')
    data = fd.readlines()
    fd.close

    # check if data as expected
    data = ''.join(data)
    header_string = ','.join(header) + '\n'
    expected = header_string + file_data

    msg = 'Expected data:\n%s\ngot:\n%s' % (expected, data)
    assert expected == data, msg

    os.remove(filename)
    os.remove(filename2)
