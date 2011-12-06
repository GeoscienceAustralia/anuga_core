#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$06/12/2011 7:02:15 PM$"


from anuga.utilities.sww_merge import sww_merge

if __name__ == "__main__":

    import argparse
    from anuga.anuga_exceptions import ANUGAError


    parser = argparse.ArgumentParser(description='Merge sww files created from parallel run')
    parser.add_argument('-np', type=int, default = 4,
                   help='number of processors used to produce sww files')
    parser.add_argument('-f', type=str, default="domain",
                   help='base sww file name')
    parser.add_argument('-v', nargs='?', type=bool, const=True, default=False,
                   help='verbosity')

    args = parser.parse_args()

    np = args.np
    filebase = args.f
    verbose = args.v

    #print np
    #print filebase
    print verbose
    


    output = filebase+".sww"
    swwfiles = [ filebase+"_P"+str(v)+"_"+str(np)+".sww" for v in range(np)]

    print swwfiles
    
    try:
        sww_merge(swwfiles, output, verbose)
    except:
        msg = 'ERROR: When merging sww files '+" ".join(swwfiles)
        print msg
        raise
