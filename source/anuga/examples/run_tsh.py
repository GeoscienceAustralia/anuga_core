"""Example of shallow water wave equation.

Specific methods pertaining to the 2D shallow water equation
are imported from shallow_water
for use with the generic finite volume framework

A example of running this program is;
python run_tsh.py n hill.tsh 0.05 1
"""

######################
# Module imports 
#

from Numeric import array
import time
import sys
from os import sep, path

from anuga.shallow_water import Domain, Reflective_boundary, \
     Dirichlet_boundary, Transmissive_boundary, Time_boundary
from anuga.abstract_2d_finite_volumes.region import Add_value_to_region, \
     Set_region
from anuga.visualiser import RealtimeVisualiser



#from anuga.config import default_datadir 

######################
# Domain

import sys


    ######NEW
def add_x_y(x, y):
    return x+y

    ######NEW

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
    domain = Domain(filename)

    # check if the visualiser will work
    try:
        xx = RealtimeVisualiser(domain)
    except:
        print "Warning: Error in RealtimeVisualiser. Could not visualise."
        visualise = False
    
    print "Number of triangles = ", len(domain)
    print "domain.geo_reference",domain.geo_reference 
    domain.checkpoint = False #True
    domain.default_order = 1
    domain.smooth = True
    domain.set_datadir('.')

    if (visualise):
        domain.store = False  #True    #Store for visualisation purposes
    else:
        domain.store = True  #True    #Store for visualisation purposes
        domain.format = 'sww'   #Native netcdf visualisation format
    
        file_path, filename = path.split(filename)
        filename, ext = path.splitext(filename)
        if domain.smooth is True:
            s = 'smooth'
        else:
            s = 'nonsmooth'        
        domain.set_name(filename + '_' + s + '_ys'+ str(yieldstep) + \
                        '_ft' + str(finaltime))
        print "Output being written to " + domain.get_datadir() + sep + \
              domain.get_name() + "." + domain.format


    #Set friction
    manning = 0.07
    inflow_stage = 10.0

    
    domain.set_quantity('friction', manning)

    #domain.set_quantity('stage', add_x_y)
    #domain.set_quantity('elevation',
    #                        domain.quantities['stage'].vertex_values+ \
    #                        domain.quantities['elevation'].vertex_values)
    #domain.set_quantity('stage', 0.0)
    

    ######################
    # Boundary conditions
    #
    print 'Boundaries'
    reflective = Reflective_boundary(domain)
    Bt = Transmissive_boundary(domain)

    #Constant inflow
    Bd = Dirichlet_boundary(array([3, 0.0, 0.0]))

    #Time dependent inflow
    from math import sin, pi
    Bw = Time_boundary(domain=domain,
                       f=lambda x: array([(1 + sin(x*pi/4))*\
                        (inflow_stage*(sin(2.5*x*pi)+0.7)),0,0]))


    print 'Available boundary tags are', domain.get_boundary_tags()

    #Set boundary conditions
    
    tags = {}
    tags['left'] = Bw
    tags['1'] = Bd
    
    tags['wave'] = Bd
    tags['wave'] = Time_boundary(domain=domain,
                       f=lambda x: array([(1 + sin(x*pi/4))*\
                        (0.15*(sin(2.5*x*pi)+0.7)),0,0]))
    tags['internal'] = None
    tags['levee'] = None
    tags['0'] = reflective
    tags['wall'] = reflective 
    tags['external'] = reflective  
    tags['exterior'] = reflective
    tags['open'] = Bd  
    tags['opening'] = None 

    domain.set_boundary(tags)

    # region tags
    
    domain.set_region(Set_region('slow', 'friction', 20, location='unique vertices'))
    domain.set_region(Set_region('silo', 'elevation', 20, location='unique vertices'))
    domain.set_region(Set_region('wet', 'elevation', 0, location='unique vertices'))
    domain.set_region(Set_region('dry', 'elevation', 2, location='unique vertices'))
    domain.set_region(Add_value_to_region('wet', 'stage', 1.5, location='unique vertices', initial_quantity='elevation'))
    domain.set_region(Add_value_to_region('dry', 'stage', 0, location='unique vertices', initial_quantity='elevation'))
    
    #print domain.quantities['elevation'].vertex_values
    #print domain.quantities['stage'].vertex_values
         
    domain.check_integrity()
    
    # prepare the visualiser
    if visualise is True:
        v = RealtimeVisualiser(domain)
        v.render_quantity_height('elevation', dynamic=False)
        v.render_quantity_height('stage', dynamic=True)
        v.colour_height_quantity('stage', (0.0, 0.0, 0.8))
        v.start()
    ######################
    #Evolution
    t0 = time.time()
    for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
        domain.write_time()
        if visualise is True:
            v.update()
    if visualise is True:
        v.evolveFinished()

    
    print 'That took %.2f seconds' %(time.time()-t0)

    
