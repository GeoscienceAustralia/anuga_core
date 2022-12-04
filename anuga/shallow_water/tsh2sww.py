"""Given a .tsh, print an sww
"""


######################
# Module imports 
#

import sys
from os import sep, path
sys.path.append('..'+sep+'pyvolution')

from anuga import Domain
from anuga.abstract_2d_finite_volumes.pmesh2domain import pmesh_to_domain_instance
import time, os 
from anuga.file.sww import SWW_file
from anuga.utilities.numerical_tools import mean
import anuga.utilities.log as log


def tsh2sww(infilename, sww_file_name = None, verbose = False):
    """
    This converts a mesh file (.tsh/.msh) to an .sww file.
    This is usefull to visualise the mesh.
    
    Note: This currently just writes the output file in the input file dir.
    """
    if verbose == True: log.critical('Creating domain from %s' % infilename)
    domain = pmesh_to_domain_instance(infilename, Domain)
    if verbose == True: log.critical("Number of triangles = %d" % len(domain))

    domain.smooth = True
    domain.format = 'sww'   #Native netcdf visualisation format
        
    file_path, filename = path.split(infilename)
    filename, ext = path.splitext(filename)
    
    if not (sww_file_name is None):
        file_path, filename = path.split(sww_file_name)
        filename, ext = path.splitext(filename)
    domain.set_name(filename)
        
    domain.reduction = mean
    if verbose == True: log.critical("file_path %s" % file_path)
    if file_path == "":file_path = "."
    domain.set_datadir(file_path)

    if verbose == True:
        log.critical("Output written to %s%s%s.%s" % (domain.get_datadir(), sep, domain.get_name(), domain.format))
    
    sww = SWW_file(domain)
    sww.store_connectivity()
    sww.store_timestep('stage')

if __name__ == "__main__":
    usage = "usage:python %s pmesh_file_name [verbose|non_verbose]" %         path.basename(sys.argv[0])
    if len(sys.argv) < 2:
        print(usage)
    else:
        filename = sys.argv[1]
        verbose = False
    if len(sys.argv) > 2:
        if sys.argv[2][0] == "v" or sys.argv[2][0] == "V":
            verbose = True
        tsh2sww(filename, verbose = verbose)
    
