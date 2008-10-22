#!/usr/bin/env python


"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment

This is a very simple test of the parallel algorithm using the simplified parallel API
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from Numeric import allclose

from anuga.pmesh.mesh_interface import create_mesh_from_regions
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.polygon import is_inside_polygon

from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.shallow_water import Dirichlet_boundary
from anuga.shallow_water import Time_boundary
from anuga.shallow_water import Transmissive_boundary

from parallel_api import distribute, myid, numprocs


#--------------------------------------------------------------------------
# Setup computational domain on processor 0
#--------------------------------------------------------------------------
if myid == 0 :
    points, vertices, boundary = rectangular_cross(10, 10) # Basic mesh

    #print 'points', points
    #print 'vertices', vertices
    
    domain = Domain(points, vertices, boundary) # Create domain
else:
    domain = None

#--------------------------------------------------------------------------
# Setup initial conditions
#--------------------------------------------------------------------------

def topography(x,y): 
    return -x/2                              # linear bed slope


if myid == 0:
    domain.set_quantity('elevation', topography) # Use function for elevation
    domain.set_quantity('friction', 0.0)         # Constant friction 
    domain.set_quantity('stage', expression='elevation') # Dry initial stage


#--------------------------------------------------------------------------
# Create the parallel domain
#--------------------------------------------------------------------------

domain = distribute(domain, verbose=True)

domain.set_name('runup')                    # Set sww filename
domain.set_datadir('.')                     # Set output dir
domain.set_maximum_allowed_speed(100)       # 
domain.set_quantities_to_be_stored(None)
domain.tight_slope_limiters = 0 # Backwards compatibility (14/4/7)
domain.H0 = 0 # Backwards compatibility (6/2/7)
domain.use_centroid_velocities = 0 # Backwards compatibility (7/5/8)

#------------------------------------------------------------------------------
# Setup boundary conditions
# This must currently happen *after* domain has been distributed
#------------------------------------------------------------------------------

Br = Reflective_boundary(domain)      # Solid reflective wall
Bd = Dirichlet_boundary([-0.2,0.,0.]) # Constant boundary values

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})



#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

interpolation_points = [[0.4,0.5], [0.6,0.5], [0.8,0.5], [0.9,0.5]]
gauge_values = []
local_interpolation_points = []
for i, point in enumerate(interpolation_points):
    gauge_values.append([]) # Empty list for timeseries

    if is_inside_polygon(point, domain.get_boundary_polygon()):

        # FIXME: One point appears on multiple processes
        # Need to get true boundary somehow
        
        #print 'P%d: point=[%f,%f]' %(myid, point[0], point[1])
        local_interpolation_points.append(i)

# Hack before we excluded ghosts.
if numprocs == 2:
    if myid == 0:
        del local_interpolation_points[0]                
        #local_interpolation_points = [1,2,3]
if numprocs == 3:
    if myid == 1:
        del local_interpolation_points[0]
if numprocs == 4:
    if myid == 0:
        del local_interpolation_points[1] #2
        del local_interpolation_points[1] #3                
    if myid == 3:
        del local_interpolation_points[1]



print 'P%d has points = %s' %(myid, local_interpolation_points)


time = []

for t in domain.evolve(yieldstep = 0.1, finaltime = 5.0):
    domain.write_time()

    # Record time series at known points
    time.append(domain.get_time())
    
    stage = domain.get_quantity('stage')
    w = stage.get_values(interpolation_points=interpolation_points)

    print 'P%d:w(%f) = '%(myid,domain.get_time()),w
    
    for i, _ in enumerate(interpolation_points):
        gauge_values[i].append(w[i])

## if myid == 0:
##     for i, (x,y) in enumerate(interpolation_points):

##         try:
##             from pylab import *
##         except:
##             pass
##         else:
##             ion()
##             hold(False)
##             plot(time, gauge_values[i], 'r.-')
##             #time, predicted_gauge_values[i], 'k-')

##             title('Gauge %d (%f,%f)' %(i,x,y))
##             xlabel('time(s)')
##             ylabel('stage (m)')    
##             #legend(('Observed', 'Modelled'), shadow=True, loc='upper left')
##             #savefig('Gauge_%d.png' %i, dpi = 300)

##         raw_input('Next')
        


# Reference from sequential version (also available as a
# unit test in test_shallow_water_domain)
# Added Friday 13 October 2006 by Ole

        G0 = ensure_numeric([-0.20000000000000001, -0.20000000000000001, -0.19920600846161715, -0.19153647344085376, -0.19127622768281194, -0.1770671909675095, -0.16739412133181927, -0.16196038919122191, -0.15621633053131384, -0.15130021599977705, -0.13930978857215484, -0.19349274358263582, -0.19975307598803765, -0.19999897143103357, -0.1999999995532111, -0.19999999999949952, -0.19999999999949952, -0.19999999999949952, -0.19997270012494556, -0.19925805948554556, -0.19934513778450533, -0.19966484196394893, -0.1997352860102084, -0.19968260481750394, -0.19980280797303882, -0.19998804881822749, -0.19999999778075916, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167, -0.19999999999966167])

        G1 = ensure_numeric([-0.29999999999999993, -0.29999588068034899, -0.29250047332330331, -0.28335081844518584, -0.26142206997410805, -0.22656028856329835, -0.21224087216745585, -0.19934324109114465, -0.1889857939783175, -0.18146311603911383, -0.17401078727434263, -0.15419361061257214, -0.16225060576782063, -0.19010941396999181, -0.20901161407004412, -0.21670683975774699, -0.21771386270738891, -0.21481284465869752, -0.21063120869004387, -0.20669243364582401, -0.20320707386714859, -0.19984087691926442, -0.19725417448019505, -0.19633783049069981, -0.19650494599999785, -0.19708316838336942, -0.19779309449413818, -0.19853070294429562, -0.19917342167307153, -0.19964814677795845, -0.19991627610824922, -0.20013162970144974, -0.20029864969405509, -0.20036259676501131, -0.20030682824965193, -0.20016105135750167, -0.19997664501985918, -0.19980185871568762, -0.19966836175417696, -0.19958856744312226, -0.19955954696194517, -0.19956950051110917, -0.19960377086336181, -0.19964885299433241, -0.19969427478531132, -0.19973301547655564, -0.19976121574277764, -0.19977765285688653, -0.19978315117522441, -0.19977994634841739, -0.19977101394878494])
        
        G2 = ensure_numeric([-0.40000000000000002, -0.39077401254732241, -0.33350466136630474, -0.29771023004255281, -0.27605439066140897, -0.25986156218997497, -0.24502185018573647, -0.231792624329521, -0.21981564668803993, -0.20870707082936543, -0.19877739883776599, -0.18980922837977957, -0.17308011674005838, -0.16306400164013773, -0.17798470933304333, -0.1929554075869116, -0.20236705191987037, -0.20695767560655007, -0.20841025876092567, -0.20792102174869989, -0.20655350005579293, -0.20492002526259828, -0.20310627026780645, -0.20105983335287836, -0.19937394565794653, -0.19853917506699659, -0.19836389977624452, -0.19850305023602796, -0.19877764028836831, -0.19910928131034669, -0.19943705712418805, -0.19970344172958865, -0.19991076989870474, -0.20010020127747646, -0.20025937787100062, -0.20035087292905965, -0.20035829921463297, -0.20029606557316171, -0.20019606915365515, -0.20009096093399206, -0.20000371608204368, -0.19994495432920584, -0.19991535665176338, -0.19990981826533513, -0.19992106419898723, -0.19994189853516578, -0.19996624091229293, -0.19998946016985167, -0.20000842303470234, -0.20002144460718174, -0.20002815561337187])
        
        G3 = ensure_numeric([-0.45000000000000001, -0.37631169657400332, -0.33000044342859486, -0.30586045469008522, -0.28843572253009941, -0.27215308978603808, -0.25712951540331219, -0.2431608296216613, -0.23032023651386374, -0.2184546873456619, -0.20735123704254332, -0.19740397194806389, -0.1859829564064375, -0.16675980728362105, -0.16951575032846536, -0.1832860872609344, -0.19485758939241243, -0.20231368291811427, -0.20625610376074754, -0.20758116241495619, -0.20721445402086161, -0.20603406830353785, -0.20450262808396991, -0.2026769581185151, -0.2007401212066364, -0.19931160535777592, -0.19863606301128725, -0.19848511940572691, -0.19860091042948352, -0.19885490669377764, -0.19916542732701112, -0.19946678238611959, -0.19971209594104697, -0.19991912886512292, -0.2001058430788881, -0.20024959409472989, -0.20032160254609382, -0.20031583165752354, -0.20025051539293123, -0.2001556115816068, -0.20005952955420872, -0.1999814429561611, -0.19992977821558131, -0.19990457708664208, -0.19990104785490476, -0.19991257153954825, -0.19993258231880562, -0.19995548502882532, -0.19997700760919687, -0.19999429663503748, -0.20000588800248761])

        


# Only compare those that belong to this process id
G = [G0, G1, G2, G3]

for i in local_interpolation_points:
    msg = 'P%d, point #%d: Computed time series and reference time series are different: %s'\
          %(myid, i, gauge_values[i]-G[i])
    assert allclose(gauge_values[i], G[i]), msg

print 'P%d completed succesfully using points = %s' %(myid, local_interpolation_points)

