#!/usr/bin/env python

"""
The MaxAsc() function.
Takes 1 or more ASC files and generates an output ASC file with
an element-wise maximum.
"""


import sys
import re

# if user does 'from ... import *' only give her/him MaxAsc
__all__ = ['MaxAsc']

HEADER_SIZE = 6

# pattern string used to split multimax data
SpacesPatternString = ' +'

# generate 're' pattern for 'any number of spaces'
SpacesPattern = re.compile(SpacesPatternString)


def MaxAsc(out_file, in_files):
    """
    MaxAsc('output_filename', ['list', 'of', 'filenames'])
    
    The output file is an ASC file with each element being the maximum of
    the corresponding element in all the input ASC files.  The output file
    has the same shape as the input file(s).
    """
    
    # get all file data into memory
    file_data = []
    for f in in_files:
        fd = open(f, 'r')
        data = fd.readlines()
        file_data.append(data)
        fd.close()

    # check we have same number of lines in each file
    num_lines = len(file_data[0])
    for (i, d) in enumerate(file_data):
        if len(d) != num_lines:
            raise RuntimeError("File %s has the wrong number of lines "
                   "(%d, expected %d)." % (in_files[i], len(d), num_lines))
    
    # open the output file
    out_fd = open(out_file, 'w')

    # read header lines, check same, write out
    for i in range(HEADER_SIZE):
        line = file_data[0][i]
        for (j, f) in enumerate(file_data):
            d = f[i]
            if d != line:
                out_fd.close()
                raise RuntimeError("File %s has the wrong header at line %d." % \
                      (in_files[j], i))
        out_fd.write(line)

    # read data lines
    for line_num in range(HEADER_SIZE, num_lines):
        lines = []
        col_file = ''
        columns = None
        for (i, f) in enumerate(file_data):
            data = f[line_num]
            if len(data) > 1:
                data = data.strip()
                data = SpacesPattern.split(data)
                if columns:
                    if columns != len(data):
                        out_fd.close()
                        raise RuntimeError("File %s doesn't have same number of columns "
                               "(%d) as file %s (%d), line %d!?" %
                               (fd.name, len(data), col_file, columns, line_num))
                else:
                    col_file = in_files[i]
                    columns = len(data)
                fdata = [float(value) for value in data]
                lines.append(fdata)
        outline = ''
        for i in range(columns):
            maximum = lines[0][i]
            for d in lines:
                if maximum < d[i]:
                    maximum = d[i]
            outline += ' %10.4e ' % maximum
        out_fd.write('%s\n' % outline)
        
    # close the output file
    out_fd.close()
