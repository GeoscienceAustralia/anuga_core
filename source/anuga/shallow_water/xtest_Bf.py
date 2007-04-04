"""

this test compares the gauge data from a boundary sww file which is precomputed and 
stored in this directory, and the same gauge location from a simple model that uses 
the same boundary sww file for its boundary conditions

Basically tests that file_boundary and evolve is working correctly
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

# Standard modules
from os import sep,getcwd, access, F_OK, mkdir, getenv
from os.path import dirname, basename,abspath 
from shutil import copy
import time, sys, os, tempfile
#import sys
#import os
#import tempfile

# Related major packages
from anuga.shallow_water import Domain,Dirichlet_boundary,File_boundary,Transmissive_boundary, Field_boundary
#from anuga.shallow_water.shallow_water_domain import Field_boundary
from Numeric import allclose, array
from anuga.pmesh.mesh_interface import create_mesh_from_regions
from anuga.abstract_2d_finite_volumes.util import start_screen_catcher, copy_code_files, sww2timeseries, get_data_from_file
#from anuga.caching import myhash
# Application specific imports

#------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#------------------------------------------------------------------------------
out_dir = dirname(abspath(__file__))+sep
#print'temp_name'
fileName="temp"
#import sys
#sys.exit()
meshes_dir_name = 'small.tsh' # this will be local

tide = 2.4


#EVOLVED MODEL

#--------------------------------------------------------------------------
# Create the triangular mesh based on overall clipping polygon with a
# tagged
# boundary and interior regions defined in project.py along with
# resolutions (maximal area of per triangle) for each polygon
#--------------------------------------------------------------------------

#all = [[469000,7760000],[470000,7758000],[468000,7758000]]
#all = [[470000,7760000],[469000,7758000],[468000,7760000]]

#all_large = [[474000,7764000],[474000,7756000],[466000,7756000],[466000,7764000]]
#all = [[470000,7764000],[470000,7756000],[466000,7756000],[466000,7764000]]
all = [[465184,7764500],[470397,7764510],[470407,7758988],[465195,7758979]]
# N, W, S, E

create_mesh_from_regions(all,
#                             boundary_tags={'ocean': [0,1,2]},
#                             boundary_tags={'ocean': [ 2],'side': [0, 1]},
                             boundary_tags={'ocean': [ 0],'side': [1, 3], 'back': [2]},
                             maximum_triangle_area=50000,
                             filename=meshes_dir_name,
                             use_cache=False,
#                             verbose=True
                             )

#-------------------------------------------------------------------------
# Setup computational domain
#-------------------------------------------------------------------------
#print 'Setup computational domain'
domain = Domain( meshes_dir_name, verbose=True)

from anuga.shallow_water.data_manager import urs2sww
boundaries_dir_name = 'o_test'

# convert MUX urs files to an SWW file output
#print 'boundary file is: ',boundaries_dir_name
from caching import cache
cache(urs2sww,
          (boundaries_dir_name,
           boundaries_dir_name), 
          {'verbose': False,
           'mint': 9200, 'maxt': 11200,
#           'origin': domain.geo_reference.get_origin(),
           'fail_on_NaN': False},
           verbose = False,
           )
    
#-------------------------------------------------------------------------
# Setup initial conditions
#-------------------------------------------------------------------------

#print 'Start Set quantity'

domain.set_quantity('elevation', -42.3)

#print 'Finished Set quantity'

#------------------------------------------------------
# Set domain parameters
#------------------------------------------------------ 

domain.set_quantity('stage', tide)
domain.set_quantity('friction', 0.01)
domain.set_name(fileName)
domain.set_datadir(out_dir)

#-------------------------------------------------------------------------
# Setup boundary conditions 
#-------------------------------------------------------------------------
#Bf = File_boundary(out_dir + 'o_test_8500_12000.sww',
#                  domain, time_thinning=24, use_cache=True, verbose=True)
Bf = Field_boundary(out_dir + 'o_test.sww',
                  domain, mean_stage=tide, use_cache=True, verbose=False)

#print 'finished reading boundary file'

#Br = Reflective_boundary(domain)
Bt = Transmissive_boundary(domain)
Bd = Dirichlet_boundary([tide,0,0])

#print'set_boundary'

domain.set_boundary({'back': Bt,
                     'side': Bd,
                    'ocean': Bf
                    })
#print'finish set boundary'

#----------------------------------------------------------------------------
# Evolve system through time
#----------------------------------------------------------------------------

#remove mesh file!!! and o_test.sww
t0 = time.time()

for t in domain.evolve(yieldstep = 60, finaltime = 1920):
    domain.write_time()
#    domain.write_boundary_statistics()



#Gets timeseries from boundary sww and evolved sww
home = getenv('INUNDATIONHOME') #Sandpit's parent dir   
#user = get_user_name()
data = 'data'
state = 'western_australia'
scenario_name = 'dampier.sww'

scenario = 'dampier_tsunami_scenario_2006'
#scenario = 'test_dampier'
an = 'anuga'
bo = 'boundaries'

run_time = 'blank'
#run_time = project.run_time
production_dirs = {run_time: 'URS evolved data'#,
                   #'boundaries': 'URS boundary condition'
                   }

topo = 'topographies'
out = 'outputs'
urs = 'urs'
gridded = '1_10000'


#gauge_boundary_filename = 'boundary_gauge_near_top.csv'
gauge_boundary_filename = 'gauges_time_series_b_near_top.csv'
gauge_evolved_filename = 'gauges_time_series_near_top.csv'

boundary_dir_filename = os.path.join(out_dir,gauge_boundary_filename)
#print'boundary_dir_filename',boundary_dir_filename

evolved_dir_filename= os.path.join(out_dir,gauge_evolved_filename)

#print'boundary_dir_filename',boundary_dir_filename
#print'evolved_dir_filename',evolved_dir_filename
                                     
swwfiles = {}
swwfile = out_dir + fileName + '.sww'
swwfiles[swwfile] = run_time
#print"swwfiles",swwfiles,"shallow_water"
        
texname, elev_output = sww2timeseries(swwfiles,
                                      out_dir+sep+"gauges.csv",
                                      production_dirs,
                                      report = False,
                                      plot_quantity = ['stage', 'xmomentum', 'ymomentum'],
                                      surface = False,
                                      time_min = None,
                                      time_max = None,
                                      title_on = False,
                                      use_cache = True,
                                      verbose = False)
                                      
swwfiles = {}
swwfile = out_dir + boundaries_dir_name + '.sww'
swwfiles[swwfile] = run_time
#print"swwfiles",swwfiles,"shallow_water"
        
texname, elev_output = sww2timeseries(swwfiles,
                                      out_dir+sep+"boundary_gauges.csv",
                                      production_dirs,
                                      report = False,
                                      plot_quantity = ['stage', 'xmomentum', 'ymomentum'],
                                      surface = False,
                                      time_min = None,
                                      time_max = None,
                                      title_on = False,
                                      use_cache = True,
                                      verbose = False)

#makes the csv files from the evolved model
e_header, e_data = get_data_from_file(evolved_dir_filename)

e_time = e_data[:,0]
e_stage = e_data[:,1]
e_momentum = e_data[:,2]
e_speed = e_data[:,3]
e_elevation = e_data[:,4]

#print'boundary_dir_filename',boundary_dir_filename
b_header, b_data = get_data_from_file(boundary_dir_filename)

b_time = b_data[:,0]
b_stage = b_data[:,1]
b_momentum = b_data[:,2]
b_speed = b_data[:,3]
b_elevation = b_data[:,4]


# compares the 2 models
j=0
k=0
b_sample = []
e_sample = []
for i in range(len(b_time)):
#    if j<(len(e_time)) and b_time[i] == e_time[j]:
    if j<(len(e_time)) and k<(len(e_time)) and b_time[i] == e_time[j]:
        b_sample.append(float(tide+b_stage[i]))
        e_sample.append(float(e_stage[k]))
#        if k <len(e_time): print 'time e equal b:', b_time[i],i, j, b_stage[i], b_sample[i], e_stage[k],(len(e_time)-1)
        j = j +1
        k = k +1


#print len(b_sample), len(e_stage),e_stage,type(e_stage)
e_stage = e_stage.tolist()
e_stage.pop()
#e_stage.pop()
 
#print len(b_sample), len(e_stage),e_stage


#assert allclose (b_sample, e_stage, 0.2, 0.2) 
print "test successful" 
os.remove(fileName+'.sww')
#print 'evolved_dir_filename',evolved_dir_filename
os.remove('gauges_time_series_near_top.csv')
os.remove('gauges_time_series_b_near_top.csv')
os.remove('gauges_t0.csv')
os.remove('gauges_maxmins.csv')
os.remove(meshes_dir_name)


assert allclose (b_sample, e_sample, 0.5, 0.5)







                             
                             
                                          





