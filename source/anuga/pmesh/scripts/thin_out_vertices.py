#!/usr/bin/python
"""
Script to thin out the number of vertices in an xya file.
Put's a square mesh over the domain and then only keeps
one vertice in each square.
"""
import os, sys
sys.path.append('..')
from mesh import *

usage = "usage: %s file_name square_length" %         os.path.basename(sys.argv[0])

if len(sys.argv) < 3:
    print usage
else:
    file_name_path= sys.argv[1]
    square_length= float(sys.argv[2])
    wallis_name = file_name_path[:-4] + "_" + str(square_length) + ".xya"
    m = loadxyafile(file_name_path," ")
    m.thinoutVertices(square_length)
    print "Saving file " + wallis_name + "."
    m.exportxyafile(wallis_name)
