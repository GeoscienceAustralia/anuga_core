from .data_audit import license_file_is_valid
import sys

def print_tree(n, indent=0):
    while n:
        print(" "*indent, n)
        print_tree(n.firstChild, indent+4)
        n = n.nextSibling

fid = open(sys.argv[1])
license_file_is_valid(fid, '.')
fid.close()


