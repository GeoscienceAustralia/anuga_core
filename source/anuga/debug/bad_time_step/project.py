"""Common filenames and locations for topographic data, meshes and outputs.
Also includes origin for slump scenario.
"""

from os import sep, environ, getenv, getcwd
import sys
from time import localtime, strftime
from os import mkdir, access, F_OK
from anuga.abstract_2d_finite_volumes.util import add_directories, \
     copy_code_files

from anuga.utilities.system_tools import get_user_name

boundary_file = 'boundary.tsm'

class Project:    
    def __init__(self,
                 trunk,
                 outputdir_name = None,
                 home = None):
        
        self.user = get_user_name()
        if home is None:
            home = getenv('INUNDATIONHOME') #Sandpit's parent dir   
        self.home = home
        # Create the structure of where the output directories will go
        scenariodir=add_directories(home, trunk)

        #Derive subdirectories and filenames
        if outputdir_name is None:
            #gets time for dir
            outputdir_name = strftime('%Y%m%d_%H%M%S',localtime()) 
        general_outputdir = scenariodir+sep+'output'+sep
        self.outputdir = general_outputdir+outputdir_name+sep
        self.meshdir = scenariodir+sep+'meshes'+sep

        # creates copy of output dir structure, if it doesn't exist
        if not access(self.meshdir,F_OK):
            mkdir (self.meshdir)       
        if not access(general_outputdir,F_OK):
            mkdir (general_outputdir)
        if not access(self.outputdir,F_OK):
            mkdir (self.outputdir)

        self.codedir = getcwd()+sep


#-------------------------------------------------------------
if __name__ == "__main__":
    p = Project(['eagle'], 'lego','.')
    print p.outputdir
    print p.meshdir
