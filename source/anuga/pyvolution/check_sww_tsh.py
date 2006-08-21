
######################
# Module imports 
import sys
from os import sep, path

import data_manager
from load_mesh.loadASCII import import_mesh_file
from shallow_water import Domain
from pmesh2domain import pmesh_to_domain_instance

def check_sww_tsh(sww_file, tsh_file, verbose = False):
    [xmin, xmax, ymin, ymax, stagemin, stagemax] =  \
           data_manager.extent_sww(sww_file)
    if verbose == True:print "Extent of ", sww_file
    if verbose == True:print "xmin", xmin 
    if verbose == True:print "xmax", xmax 
    if verbose == True:print "ymin", ymin 
    if verbose == True:print "ymax", ymax 
    if verbose == True:print "stagemin", stagemin 
    if verbose == True:print "stagemax", stagemax
    
    domain = pmesh_to_domain_instance(tsh_file, Domain)
    [tsh_xmin, tsh_xmax, tsh_ymin, tsh_ymax] =  domain.get_extent()
    if verbose == True:print "Extent of ", tsh_file
    if verbose == True:print "tsh_xmin", tsh_xmin 
    if verbose == True:print "tsh_xmax", tsh_xmax 
    if verbose == True:print "tsh_ymin", tsh_ymin 
    if verbose == True:print "tsh_ymax", tsh_ymax
    is_subset = xmin < tsh_xmin and xmax > tsh_xmax and  \
                ymin < tsh_ymin and ymax > tsh_ymax 
    if verbose == True:
        if is_subset:
            print "tsh within sww"
        else:
            print "WARNING: tsh NOT within sww"
        
  
#-------------------------------------------------------------
if __name__ == "__main__":
    """
    Load in a mesh and data points with attributes.
    Fit the attributes to the mesh.
    Save a new mesh file.
    """
    import os, sys
    usage = "usage: %s *.sww *.tsh" %os.path.basename(sys.argv[0])

    if len(sys.argv) < 3:
        print usage
    else:
        sww_file = sys.argv[1]
        tsh_file = sys.argv[2]
        check_sww_tsh(sww_file, tsh_file, verbose = True)
