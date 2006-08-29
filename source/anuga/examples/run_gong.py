"""Example of shallow water wave equation.

Specific methods pertaining to the 2D shallow water equation
are imported from shallow_water
for use with the generic finite volume framework

A example of running this program is;
python run_tsh.py visualise hill.tsh 0.05 1
"""

######################
# Module imports 
#

import sys
from os import sep, path
sys.path.append('..'+sep+'pyvolution')

from shallow_water import Domain, Reflective_boundary, Dirichlet_boundary,\
     Transmissive_boundary, Time_boundary, Add_value_to_region, File_boundary
from mesh_factory import rectangular
from anuga.pyvolution.pmesh2domain import pmesh_to_domain, pmesh_to_domain_instance

from Numeric import array

import time

######################
# Domain

import sys

usage = "usage: %s ['visual'|'non-visual'] pmesh_file_name yieldstep finaltime" %         path.basename(sys.argv[0])

if len(sys.argv) < 4:
    print usage
else:
    if sys.argv[1][0] == "n" or sys.argv[1][0] == "N":
        visualise = False
    else:    
        visualise = True
    filename = sys.argv[2]
    yieldstep = float(sys.argv[3])
    finaltime = float(sys.argv[4])
    
    print 'Creating domain from', filename
    domain = pmesh_to_domain_instance(filename, Domain)
    print "Number of triangles = ", len(domain)

    domain.checkpoint = False #True
    domain.default_order = 1
    domain.visualise = visualise
    domain.smooth = True
    domain.set_datadir('.')
    domain.starttime = 20000 

    if (domain.visualise):
        domain.store = False  #True    #Store for visualisation purposes
    else:
        domain.store = True  #True    #Store for visualisation purposes
        domain.format = 'sww'   #Native netcdf visualisation format
    
        file_path, filename = path.split(filename)
        filename, ext = path.splitext(filename)
        domain.filename = filename + '_' + '_ys'+ str(yieldstep) + \
                          '_ft' + str(finaltime)
        print "Output being written to " + domain.get_datadir() + sep + \
              domain.filename + "." + domain.format


    #Set friction
    manning = 0.07
    inflow_stage = 10.0
    domain.set_quantity('friction', manning)

    ######################
    # Boundary conditions
    #
    print 'Boundaries'
    reflective = Reflective_boundary(domain)
    Bt = Transmissive_boundary(domain)

    #Constant inflow
    Bd = Dirichlet_boundary(array([inflow_stage, 0.0, 0.0]))

    #Time dependent inflow
    from math import sin, pi
    Bw = Time_boundary(domain=domain,
                       f=lambda x: array([(1 + sin(x*pi/4))*\
                        (inflow_stage*(sin(2.5*x*pi)+0.7)),0,0]))

    boundary_file = 'kermadec2.sww'

    Fb = File_boundary(boundary_file, domain, verbose=True)

    print 'Available boundary tags are', domain.get_boundary_tags()
    print 'The extent is ', domain.get_extent()

    #Set the ocean...
    domain.set_quantity('stage', 0.0)

    #Set boundary conditions
    
    tags = {}
   
    tsunami = Fb
    tags['w1'] = tsunami 
    tags['w2'] = tsunami
    tags['w3'] = tsunami
    tags['w4'] = tsunami
    tags['w5'] = tsunami
    tags['w6'] = tsunami
    tags['w7'] = tsunami
    tags['wave'] = tsunami
    tags['external'] = tsunami
    tags['exterior'] = tsunami

    domain.set_boundary(tags)

    # region tags
    
    domain.check_integrity()

    ######################
    #Evolution
    t0 = time.time()
    for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
        domain.write_time()
    
    print 'That took %.2f seconds' %(time.time()-t0)

    
