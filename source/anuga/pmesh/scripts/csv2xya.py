#!/usr/bin/python
"""
Convert a .cvs file from gis to a xya file, for pmesh.
(removes the 1st column in a cvs file and renames it from .cvs to .xya)
"""
import os, sys

usage = "usage: %s input_file_name.csv" %         os.path.basename(sys.argv[0])

if (len(sys.argv)< 2) or not (sys.argv[1][-4:]== ".csv"):
    print usage
else:
    read_file_name = sys.argv[1]
    write_file_name = read_file_name[:-4] + ".xya"
    file_obj = open(read_file_name)
    lines = file_obj.readlines()
    file_obj.close()

    
    fd = open(write_file_name,'w')
    for line in lines:
        elements = line.split(',')
        elements.pop(0) # remove the 1st column
        if len(elements) > 0: #if there is a 2nd column, write it
            fd.write(elements.pop(0))
        for element in elements:
            fd.write(',')
            fd.write(element)
