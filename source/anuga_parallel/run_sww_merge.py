#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$06/12/2011 7:02:15 PM$"


from anuga.utilities.sww_merge import sww_merge

if __name__ == "__main__":

    filename = 'merimbula'
    np = 2;
    verbose = True


    output = filename+".sww"
    swwfiles = " ".join(filename+"P_"+str(np)+"_"+str(v)+".sww" for v in range(np))

    sww_merge(swwfiles, output, verbose)
