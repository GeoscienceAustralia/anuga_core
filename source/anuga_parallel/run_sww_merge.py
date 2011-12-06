#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$06/12/2011 7:02:15 PM$"


from anuga.utilities.sww_merge import sww_merge

if __name__ == "__main__":

    filename = '100y'
    np = 4;
    verbose = True


    output = filename+".sww"
    swwfiles = [filename+"_P"+str(v)+"_"+str(np)+".sww" for v in range(np)]

    print swwfiles
    
    sww_merge(swwfiles, output, verbose)
