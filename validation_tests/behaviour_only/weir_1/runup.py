"""Runup example from the manual, slightly modified
"""
#---------
#Import Modules
#--------
import anuga

import numpy
from anuga import myid, finalize, distribute, barrier

from math import exp

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

scale_me=1.0

## Set up mesh
boundaryPolygon=[ [0., 0.], [0., 100.], [100.0, 100.0], [100.0, 0.0]]
midResPolygon=[ [30., 30.], [30., 70.], [70., 70.], [70., 30.]]
higherResPolygon=[ [40., 40.], [40., 60.], [60., 60.], [60., 40.]]
# Riverwall = list of lists, each with a set of x,y,z (and optional QFactor) values
#riverWall={ 'centralWall':
#                [ [50., 0.0, -0.0],
#                  [50., 45., -0.0],
#                  [50., 46., -0.2],
#                  [50., 54., -0.2],
#                  [50., 55., -0.0], 
#                  [50., 100.0, -0.0]] 
#          }
riverWall={ 'leftWall':
                [ [50., 0.0, -0.0],
                  [50., 45.5, -0.0]],
             'centralWall': 
                [[50., 45.5, -0.2],
                 [50., 54.5, -0.2]],
             'rightWall':
                [ [50., 54.5, -0.0], 
                  [50., 100.0, -0.0]] 
          }

#riverWall_Par={'centralWall':{'Qfactor':1.0}}
# Try to avoid any shallow-water type solution -- becomes unstable
#riverWall_Par={'centralWall':{'Qfactor':1.0, 's1': 0.999, 's2':0.9999, 'h1':100, 'h2':150}}

# The boundary polygon + riverwall breaks the mesh into multiple regions
# Must define the resolution in these areas with an xy point + maximum area
# Otherwise triangle.c gets confused
regionPtAreas=[ [99., 99., 10.0*10.0*0.5],
                [1., 1., 10.0*10.0*0.5],
                [45., 45., 1.0*1.0*0.5],
                [55., 55., 1.0*1.0*0.5],
                [65., 65., 3.0*3.0*0.5],
                [35., 35., 3.0*3.0*0.5] ]

if myid == 0:
    #==================================================================
    # Create Sequential Domain
    #==================================================================
    anuga.create_mesh_from_regions(boundaryPolygon, 
                             boundary_tags={'left': [0],
                                            'top': [1],
                                            'right': [2],
                                            'bottom': [3]},
                               maximum_triangle_area = 1.0e+20,
                               minimum_triangle_angle = 28.0,
                               filename = 'runup.msh',
                               interior_regions = [ [higherResPolygon, 1.*1.*0.5],
                                                    [midResPolygon, 3.0*3.0*0.5]],
                               breaklines=riverWall.values(),
                               use_cache=False,
                               verbose=True,
                               regionPtArea=regionPtAreas)
    
    domain=anuga.create_domain_from_file('runup.msh')
    
    
    domain.set_flow_algorithm(alg)
    
    
    domain.set_name('runup_riverwall')                         
    domain.set_datadir('.')                         
    domain.set_store_vertices_uniquely()
    
    #------------------
    # Define topography
    #------------------
    
    def topography(x,y):
        return -x/150.*scale_me 
    
    def stagefun(x,y):
        stg=-0.5*scale_me
        return stg 
    
    domain.set_quantity('elevation',topography)     # Use function for elevation
    
    domain.set_quantity('friction',0.03)             # Constant friction
    
    domain.set_quantity('stage', stagefun)              # Constant negative initial stage
else:
    domain = None

#======================================================================
# create Parallel Domain
#======================================================================    
domain = distribute(domain)

domain.riverwallData.create_riverwalls(riverWall)

#--------------------------
# Setup boundary conditions
#--------------------------

Br=anuga.Reflective_boundary(domain)            # Solid reflective wall

def boundaryFun(t):
    output=-0.4*exp(-t/100.)-0.1
    output=min(output,-0.11)
    #output=min(output,-0.3)
    return output  

Bin_tmss = anuga.Transmissive_momentum_set_stage_boundary(domain=domain, function = boundaryFun) 

#----------------------------------------------
# Associate boundary tags with boundary objects
#----------------------------------------------
domain.set_boundary({'left': Br, 'right': Bin_tmss, 'top': Br, 'bottom':Br})

#------------------------------------------------------------------------------
# Produce a documentation of parameters
#------------------------------------------------------------------------------
if myid == 0:
    parameter_file=open('parameters.tex', 'w')
    parameter_file.write('\\begin{verbatim}\n')
    from pprint import pprint
    pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
    parameter_file.write('\\end{verbatim}\n')
    parameter_file.close()

#------------------------------
#Evolve the system through time
#------------------------------

for t in domain.evolve(yieldstep=10.0,finaltime=4000.0):
    if myid == 0 and verbose: print(domain.timestepping_statistics())
    # Print velocity as we go
#     uh=domain.quantities['xmomentum'].centroid_values
#     vh=domain.quantities['ymomentum'].centroid_values
#     depth=domain.quantities['height'].centroid_values
#     depth=depth*(depth>1.0e-06) + 1.0e-06
#     vel=((uh/depth)**2 + (vh/depth)**2)**0.5
#     print('peak speed is', vel.max())


domain.sww_merge(delete_old=True)

finalize()


