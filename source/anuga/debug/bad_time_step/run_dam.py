"""

Script for running a breaking wave simulation of Jon Hinwoods wave tank.
Note: this is based on the frinction_ua_flume_2006 structure.


Duncan Gray, GA - 2007 



"""


#----------------------------------------------------------------------------
# Import necessary modules
#----------------------------------------------------------------------------

# Standard modules
import time
from time import localtime, strftime
import sys
from shutil import copy
from os import path, sep
from os.path import dirname  #, basename


# Related major packages
from anuga.shallow_water import Domain, Reflective_boundary, \
                            Dirichlet_boundary,  Time_boundary, \
                            File_boundary, \
                            Transmissive_Momentum_Set_Stage_boundary
from anuga.fit_interpolate.interpolate import interpolate_sww2csv
from anuga.abstract_2d_finite_volumes.util import start_screen_catcher, \
     copy_code_files, file_function
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import File_boundary_time

# Scenario specific imports
import project                 # Definition of file names and polygons
import create_mesh

def elevation_function(x,y):
    from Numeric import zeros, size, Float
    slope = 4 ## Bit of a magic Number
    
    z = zeros(size(x), Float)
    for i in range(len(x)):
        if x[i] < slope:
            z[i] = 0.0
        else:
            z[i] = (x[i]-slope)*0.1
    return z 

def main(friction=0.01, outputdir_name=None, is_trial_run=False):
    basename = 'zz' + str(friction)
    if is_trial_run is True:
        outputdir_name += '_test'
        yieldstep = 0.1
        finaltime = 15.
        maximum_triangle_area=0.01
    else:
        yieldstep = 0.02
        finaltime = 15.1
        maximum_triangle_area=0.0001
        
        maximum_triangle_area=0.001
        
    pro_instance = project.Project([],
                                   outputdir_name=outputdir_name,
                                   home='.')
    print "The output dir is", pro_instance.outputdir
    copy_code_files(pro_instance.outputdir,__file__,
                    dirname(project.__file__) \
                    + sep + project.__name__+'.py')
    copy (pro_instance.codedir + 'run_dam.py',
          pro_instance.outputdir + 'run_dam.py')
    copy (pro_instance.codedir + 'create_mesh.py',
          pro_instance.outputdir + 'create_mesh.py')
    
    mesh_filename = pro_instance.meshdir + basename + '.msh'

    #--------------------------------------------------------------------------
    # Copy scripts to output directory and capture screen
    # output to file
    #--------------------------------------------------------------------------

    # creates copy of code in output dir
    if is_trial_run is False:
        start_screen_catcher(pro_instance.outputdir, rank, pypar.size())

    print 'USER:    ', pro_instance.user
    #-------------------------------------------------------------------------
    # Create the triangular mesh
    #-------------------------------------------------------------------------

    # this creates the mesh
    #gate_position = 12.0
    create_mesh.generate(mesh_filename,
                         maximum_triangle_area=maximum_triangle_area)

    head,tail = path.split(mesh_filename)
    copy (mesh_filename,
          pro_instance.outputdir + tail )
    #-------------------------------------------------------------------------
    # Setup computational domain
    #-------------------------------------------------------------------------
    domain = Domain(mesh_filename, use_cache = False, verbose = True)
   

    print 'Number of triangles = ', len(domain)
    print 'The extent is ', domain.get_extent()
    print domain.statistics()

    
    domain.set_name(basename)
    domain.set_datadir(pro_instance.outputdir)
    domain.set_quantities_to_be_stored(['stage', 'xmomentum', 'ymomentum'])
    domain.set_minimum_storable_height(0.001)
    #domain.set_store_vertices_uniquely(True)  # for writting to sww

    #-------------------------------------------------------------------------
    # Setup initial conditions
    #-------------------------------------------------------------------------

    domain.set_quantity('stage', 0.06)
    domain.set_quantity('friction', friction) 
    domain.set_quantity('elevation', elevation_function)

    
    print 'Available boundary tags', domain.get_boundary_tags()

    # Create boundary function from timeseries provided in file
    function = file_function(project.boundary_file, domain, verbose=True)
    Bts = Transmissive_Momentum_Set_Stage_boundary(domain, function)

    Br = Reflective_boundary(domain)
    #Bd = Dirichlet_boundary([0.3,0,0]) 
    #Bts = Time_boundary(domain, function)
    domain.set_boundary( {'wall': Br, 'wave': Bts} )

    #-------------------------------------------------------------------------
    # Evolve system through time
    #-------------------------------------------------------------------------
    t0 = time.time()

    for t in domain.evolve(yieldstep, finaltime):
    
        domain.write_time()
        print 'That took %.2f seconds' %(time.time()-t0)
        print 'finished'

    points = [[2.8,0.225],  #-1.8m from SWL
              [5.1,0.225],  #0.5m from SWL
              [6.6,0.225],  #2m from SWL
              [6.95,0.255], #2.35m from SWL
              [7.6,0.255],  #3m from SWL
              [8.2,0.255],  #3.5m from SWL
              [9.2,0.255]  #4.5m from SWL
              ]

#     points = [[gate_position - 0.65,0.2],
#               [gate_position - 0.55,0.2],
#               [gate_position - 0.45,0.2],
#               [gate_position - 0.35,0.2],
#               [gate_position - 0.25,0.2]
#               ]

    #-------------------------------------------------------------------------
    # Calculate gauge info
    #-------------------------------------------------------------------------

    if False:
        interpolate_sww2csv(pro_instance.outputdir + basename +".sww",
                            points,
                            pro_instance.outputdir + "depth_manning_"+str(friction)+".csv",
                            pro_instance.outputdir + "velocity_x.csv",
                            pro_instance.outputdir + "velocity_y.csv")
  
    return pro_instance

#-------------------------------------------------------------
if __name__ == "__main__":
    main( is_trial_run = True,
         outputdir_name='Hinwood_draft')
