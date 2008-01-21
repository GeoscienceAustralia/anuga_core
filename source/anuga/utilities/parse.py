import sys
from xml.dom import minidom


def print_tree(n, indent=0):
    while n:
        print " "*indent, n
        print_tree(n.firstChild, indent+4)
        n = n.nextSibling



#doc = minidom.parse(sys.argv[1])

fid = open(sys.argv[1])

#print fid.read()

#doc = minidom.parse(fid)
#print doc
#print_tree(doc)

from data_audit import license_file_is_valid
license_file_is_valid(fid, '.')


