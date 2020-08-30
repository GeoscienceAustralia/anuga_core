import numpy
import anuga
from anuga.structures import internal_boundary_functions
from anuga.structures.internal_boundary_functions import pumping_station_function


end_point0 = [49.0,50.0]
end_point1 = [51.0,50.0]

end_points = [end_point0, end_point1]

inlet1_poly = [[[end_point0[0]-10, 45.0],[end_point0[0]-10,55],
                [end_point0[0],55],[end_point0[0],45],[end_point0[0]-10, 45.0]]]
              
inlet2_poly = [[[end_point1[0], 45.0],[end_point1[0],55],
                [end_point1[0]+10,55],[end_point1[0]+10,45],[end_point1[0], 45.0]]]
              

def tobreaklines(riverWall):
    rw_values = list(riverWall.values())
    return [numpy.array(rw_values[0])[:,0:2].tolist()]


boundaryPolygon = [ [0., 0.], [0., 100.], [100.0, 100.0], [100.0, 0.0]]
wallLoc = 50.
# The boundary polygon + riverwall breaks the mesh into multiple regions
# Must define the resolution in these areas with an xy point + maximum area
# Otherwise triangle.c gets confused
length = 2.0
res = length*length*0.5
regionPtAreas = [ [99., 99., res],
                  [1., 1., res],
                  [45, 50, res],
                  [55, 50, res]]

wallHeight=10.
InitialOceanStage=2.
InitialLandStage=6.

riverWall = { 'centralWall':
                           [ [wallLoc, 0.0, wallHeight],
                             [wallLoc, 100.0, wallHeight]] 
                        }

riverWall_Par = {'centralWall':{'Qfactor':1.0}}

domain = anuga.create_domain_from_regions(boundaryPolygon, 
                         boundary_tags={'left': [0],
                                        'top': [1],
                                        'right': [2],
                                        'bottom': [3]},
                           maximum_triangle_area = 10.0,
                           minimum_triangle_angle = 28.0,
                           interior_regions =[ ], #[ [higherResPolygon, 1.*1.*0.5],
                                                  #  [midResPolygon, 3.0*3.0*0.5]],
                           breaklines=tobreaklines(riverWall)+inlet1_poly+inlet2_poly,
                           regionPtArea=regionPtAreas,
                           use_cache=False,
                           verbose=False)


domain.set_name('run_pump')
domain.set_store_vertices_uniquely(True)


#=======================================
# Setup Initial conditions
#=======================================
def topography(x,y):
    return -x/150. 

def stagefun(x,y):
    stg = InitialOceanStage*(x>=wallLoc) + InitialLandStage*(x<wallLoc)
    return stg 


# NOTE: Setting quantities at centroids is important for exactness of tests
domain.set_quantity('elevation',topography,location='centroids')     
domain.set_quantity('stage', stagefun,location='centroids')            


#========================================
# Setup wall down the middle of the domain
#========================================
domain.riverwallData.create_riverwalls(riverWall,riverWall_Par,verbose=False) 


#========================================
# Boundary conditions
# Simple reflective BC all around
#========================================
Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom':Br})


#========================================
# Setup Pump
# (1) First setup the pump characteristics
# (2) Then locate the pump using the operator
#========================================
pump_function = anuga.pumping_station_function(
            domain=domain,
            pump_capacity=100.0,
            hw_to_start_pumping=0.0,
            hw_to_stop_pumping=-1.0,
            initial_pump_rate=100.0, 
            pump_rate_of_increase = 1.0, 
            pump_rate_of_decrease = 1.0, 
            verbose=True)



pump = anuga.Internal_boundary_operator(domain, pump_function,
                                        width = 10.0,
                                        height = 1.0,
                                        apron = 10.0,
                                        end_points=end_points,
                                        verbose=True)


#============================================
# Evolve.
# Monitor the amount of water on each side
# of the wall. The sum should remain constant,
# and the change should be match the pump
# capacity
#============================================
region1 = anuga.Region(domain, polygon=[[0.0,0.0], [50.0,0.0], [50.0, 100.0], [0.0,100.0]])
region2 = anuga.Region(domain, polygon=[[50.0,0.0], [100.0,0.0], [100.0, 100.0], [50.0,100.0]])

for t in domain.evolve(yieldstep=1.0, duration=60):
    domain.print_timestepping_statistics()
    stage = domain.get_quantity('stage')
    elev  = domain.get_quantity('elevation')
    height = stage - elev

    print (anuga.indent + 'Integral1 = ', height.get_integral(region=region1))
    print (anuga.indent + 'Integral2 = ', height.get_integral(region=region2))
    print (anuga.indent + 'Total Integral = ', height.get_integral())
    #pump.print_timestepping_statistics()
    
