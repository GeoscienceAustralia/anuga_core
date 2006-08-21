"""Given a .tsh, print an sww
"""

######################
# Module imports 
#

import sys
from os import sep, path
sys.path.append('..'+sep+'pyvolution')

from shallow_water import Domain
from pmesh2domain import pmesh_to_domain_instance
import time, os 
from data_manager import get_dataobject   
from utilities.numerical_tools import mean

def tsh2sww(infilename, sww_file_name = None, verbose = False):
    """
    This converts a mesh file (.tsh/.msh) to an .sww file.
    This is usefull to visualise the mesh.
    
    Note: This currently just writes the output file in the input file dir.
    """
    if verbose == True:print 'Creating domain from', infilename
    domain = pmesh_to_domain_instance(infilename, Domain)
    if verbose == True:print "Number of triangles = ", len(domain)

    domain.smooth = True
    domain.format = 'sww'   #Native netcdf visualisation format
        
    file_path, filename = path.split(infilename)
    filename, ext = path.splitext(filename)
    
    if not (sww_file_name == None):
        file_path, filename = path.split(sww_file_name)
        filename, ext = path.splitext(filename)
    domain.filename = filename
        
    domain.reduction = mean
    if verbose == True:print "file_path",file_path
    if file_path == "":file_path = "."
    domain.set_datadir(file_path)

    if verbose == True:
        print "Output written to " + domain.get_datadir() + sep + \
              domain.filename + "." + domain.format  
    sww = get_dataobject(domain)
    sww.store_connectivity()
    sww.store_timestep('stage')

if __name__ == "__main__":
    usage = "usage:python %s pmesh_file_name [verbose|non_verbose]" %         path.basename(sys.argv[0])
    if len(sys.argv) < 2:
        print usage
    else:
        filename = sys.argv[1]
        verbose = False
    if len(sys.argv) > 2:
        if sys.argv[2][0] == "v" or sys.argv[2][0] == "V":
            verbose = True
        tsh2sww(filename, verbose = verbose)
    
