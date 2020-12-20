#!/usr/bin/env python

"""CSV file utility routines.
"""

import csv

def read_csv_file(filename, key_col, data_col):
    """Read data from a CSV file, get 'key_col' and 'data_col' columns.

    Returns ((key[0], data[0]), ...).
    """

    # Start reading the CSV file
    data = []
    fd = open(filename, 'r')
    csv_reader = csv.reader(fd)

    # Get header row, calculate required column indices
    h = next(csv_reader)
    header = [x.strip() for x in h]
    if key_col not in header:
        msg = ("Column '%s' not in file %s"
               % (key_col, filename))

        fd.close()
        raise Exception(msg)
    if data_col not in header:
        msg = ("Column '%s' not in file %s"
               % (data_col, filename))

        fd.close()
        raise Exception(msg)

    key_index = header.index(key_col)
    data_index = header.index(data_col)

    # read data, extract columns, save
    result = []
    for line in csv_reader:
        try:
            key_data = line[key_index].strip()
            data_data = line[data_index].strip()
            result.append((key_data, data_data))
        except:
            pass

    fd.close()

    return result


def merge_csv_key_values(file_title_list, output_file,
                         key_col='hours', data_col='stage'):
    """Select key and value columns from 'N' CSV files, write one CSV file.

    file_title_list: a list of (filename, new_data_column_title) values, one
                     for each input file
    output_file:     the output CSV file path
    key_col:         column header string of key column
    data_col:        column header string of data column

    The output file will look like:
        <key_col>,   <new_data_column_title1>, <new_data_column_title2>, ...
        <key_value>, <data1>,                  <data2>,                  ...
        <key_value>, <data1>,                  <data2>,                  ...
        <key_value>, <data1>,                  <data2>,                  ...
        <key_value>, <data1>,                  <data2>,                  ...

    There is an assumption that the <key_value> values are the same across
    all files for the same row.  This is tested in the code below.
    """


    # Get number of input files, check we have 1 or more

    num_files = len(file_title_list)
    if num_files == 0:
        msg = "List 'file_title_list' is empty!?"
        raise Exception(msg)

    # Read data from all files
    file_data = []
    for (filename, title) in file_title_list:
        data = read_csv_file(filename, key_col, data_col)
        file_data.append((filename, title, data))

    # Now, file_data -> [(filename, title, [(k,v), (k,v), ...], ...]
    # Sanity check, check num rows same in all files
    num_rows = None
    for (fn, t, d) in file_data:
        if num_rows is None:
            num_rows = len(d)
        else:
            if num_rows != len(d):
                msg = ('File %s has different number of rows from %s, '
                       'expected %d columns, got %d'
                       % (fn, file_data[0][0], num_rows, len(d)))
                raise Exception(msg)

    # Sanity check, check key values same in same rows
    first_key_values = [v[0] for v in file_data[0][2]]
    for (fn, t, d) in file_data:
        key_values = [v[0] for v in d]
        if key_values != first_key_values:
            msg = ('Key values differ between files %s and %s!?'
                   % (fn, file_data[0][0]))
            raise Exception(msg)

    # Open output file
    out_fd = open(output_file, 'w',newline='')
    out_csv = csv.writer(out_fd)

    # Write column rows to output file
    header = [key_col]
    for (fn, col, d) in file_data:
        header.append(col)
    out_csv.writerow(header)

    # Write data rows to output file
    file_kv_list = [x[2] for x in file_data]
    for i in range(num_rows):
        data_row = [file_kv_list[0][i][0]]
        for file_data in file_kv_list:
            data_row.append(file_data[i][1])
        out_csv.writerow(data_row)

    out_fd.close()
