#!/usr/bin/env python3
"""
    Merge a list of .sww files together into a single file.
"""



import argparse
from anuga.anuga_exceptions import ANUGAError
from anuga.utilities.sww_merge import sww_merge_parallel  # Import the required function


parser = argparse.ArgumentParser(description='Merge sww files created from parallel run')
parser.add_argument('-np', type=int, default = 4,
                help='number of processors used to produce sww files')
parser.add_argument('-f', type=str, default="domain",
                help='domain global name')
parser.add_argument('-v', nargs='?', type=bool, const=True, default=False,
                help='verbosity')
parser.add_argument('-delete_old', nargs='?', type=bool, const=True, default=False,
                help='Flag to delete the input files')
args = parser.parse_args()

np = args.np
domain_global_name = args.f
verbose = args.v
delete_old = args.delete_old


try:
    sww_merge_parallel(domain_global_name, np, verbose, delete_old)
except:
    msg = 'ERROR: When merging sww files %s '% domain_global_name
    print(msg)
    raise


