#!/usr/bin/env python
'''Routine to automatically order boundary points'''

################################################################################
# This routine is used to automatically order boundary points in a CSV file.
#
# usage: order_boundary(infile, outfile)
#
# where infile  is the string path of the input CSV file,
#   and outfile is the string path of the output CSV file.
#
# The input file is expected to have the format:
#        longitude,latitude,index
#        150.0833,-37.5,2758
#        150.1,-37.4667,2765
#        150.1167,-37.3833,2769
#        150.1333,-36.7167,2771
#        150.1333,-36.7667,2774
#        150.15,-36.5833,2777
#        150.15,-36.6333,2780
#        150.15,-36.6833,2783
#        150.15,-36.8167,2787
#
# The first point in the ordered output file is assumed to be the first point in
# the input file (ie, line 2).  The header line is preserved, as are all fields
#following the first two fields in each data line.
################################################################################

import csv


def order_boundary(infile, outfile):
    """Order a CSV file of boundary points.

    infile   path to input filep
    outfile  path to output filep

    Input file will have a header line that must be preserved.
    File format is: (longitude, latitude, <other_fields>)
    Fields after long+lat must be preserved.
    """

    def sort_points(unordered, ordered, id):
        """Sort a list of point tuples.

        unordered  unordered list of points.
        ordered    return list of ordered points.
        id         is index into 'unordered' of point to put into 'ordered'

        This code:
            . moves 'id' point from 'unordered' into 'ordered'
            . finds 'id' of next point to be moved
            . if more points, recurse
        """

        # move 'id' point from unordered to ordered, get x0, y0
        ordered.append(unordered[id])
        x0 = unordered[id][0]
        y0 = unordered[id][1]
        del unordered[id]

        # find nearest point to 'id' point in 'unordered'
        nearest = None
        for i, p in enumerate(unordered):
            x1 = p[0]
            y1 = p[1]
            # get square of distance between the 2 points
            d2 = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)
            if nearest is None or d2 < nearest:
                next = i
                nearest = d2

        # if there are more points in 'unordered', recurse
        if nearest is not None:
            sort_points(unordered, ordered, next)

        return ordered

    # read file into list of tuples, convert first 2 fields to floats.
    fd = open(infile, 'r')
    data = []
    for row in csv.reader(fd):
        try:
            row[0] = float(row[0])
        except:
            pass
        try:
            row[1] = float(row[1])
        except:
            pass
        data.append(tuple(row))
    fd.close()

    # preserve header row
    header = data[0]
    del data[0]

    # order the data
    ordered_data = sort_points(data, [], 0)

    # write ordered data to output file
    fd = open(outfile, 'w',newline="")
    w = csv.writer(fd)
    w.writerow(header)

    for d in ordered_data:
        w.writerow(d)
    fd.close()
