#!/usr/bin/env python


"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment

This is a very simple test of the parallel algorithm using the simplified parallel API
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import numpy as num

from anuga.pmesh.mesh_interface import create_mesh_from_regions

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.polygon import is_inside_polygon

from anuga.interface import Domain
from anuga.interface import Reflective_boundary
from anuga.interface import Dirichlet_boundary
from anuga.interface import Time_boundary
from anuga.interface import Transmissive_boundary

from anuga.interface import rectangular_cross

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
    #So full domain is only constructed for the 0th processor (will save on memory)
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
domain = distribute(domain, verbose=False)

#--------------------------------------------------------------------------
# Setup domain parameters
#--------------------------------------------------------------------------
domain.set_name('runup')                    # Set sww filename
domain.set_datadir('.')                     # Set output dir

domain.set_default_order(1)        
domain.set_quantities_to_be_stored(None)
domain.set_maximum_allowed_speed(100) #FIXME (Ole): try to remove this

# FIXME (Ole): Need tests where this is commented out
domain.tight_slope_limiters = 0 # Backwards compatibility (14/4/7)
domain.H0 = 0 # Backwards compatibility (6/2/7)
domain.use_centroid_velocities = 0 # Backwards compatibility (7/5/8)


#------------------------------------------------------------------------------
# Setup boundary conditions
# This must currently happen *AFTER* domain has been distributed
#------------------------------------------------------------------------------

Br = Reflective_boundary(domain)      # Solid reflective wall
Bd = Dirichlet_boundary([-0.2,0.,0.]) # Constant boundary values

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})



#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

interpolation_points = [[0.4,0.51], [0.6,0.51], [0.8,0.51], [0.9,0.51]]
gauge_values = []
tri_ids = []
for i, point in enumerate(interpolation_points):
    gauge_values.append([]) # Empty list for timeseries

    #if is_inside_polygon(point, domain.get_boundary_polygon()):
    try:
        k = domain.get_triangle_containing_point(point)
        #print 'KKK',myid, k, domain.tri_full_flag[k]
        if domain.tri_full_flag[k] == 1:
            tri_ids.append(k)
        else:
            tri_ids.append(-1)            
    except:
        tri_ids.append(-2)

#print myid, tri_ids
#print myid, domain.tri_full_flag[tri_ids[0]]


#print myid, domain.tri_full_flag



print 'P%d has points = %s' %(myid, tri_ids)


time = []

for t in domain.evolve(yieldstep = 0.1, finaltime = 5.0):
    if myid == 0: domain.write_time()

    # Record time series at known points
    time.append(domain.get_time())
    
    stage = domain.get_quantity('stage')

    for i in range(4):
        if tri_ids[i] > -1:
            gauge_values[i].append(stage.centroid_values[tri_ids[i]])




G0 = [-0.19166666666666665, -0.19166666666666665, -0.1908726751282839, -0.18320314010752042, -0.18294289434947858, -0.17186692881708804, -0.16703039399924297, -0.16196038919122194, -0.15621633053131387, -0.15130021599977714, -0.13930978857215476, -0.18515941024930252, -0.19141974265470424, -0.19166563809769938, -0.19166666621987774, -0.19166666666616614, -0.19166666666616614, -0.19166666666616614, -0.19163936679820112, -0.19092472615522754, -0.19101180444286325, -0.19133150862916881, -0.1914019526842341, -0.19134927149082989, -0.19146947464147604, -0.19165471548476315, -0.19166666444742567, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831, -0.19166666666632831]

G1 = [-0.29166666666666669, -0.29160604748111429, -0.28349212083663766, -0.26637122760395054, -0.25078182919408626, -0.22796566989512818, -0.21386588163016981, -0.20291676395026542, -0.19288397535979634, -0.18529863721918491, -0.17833464440180594, -0.16041289672597714, -0.15844675190030885, -0.18498595153509523, -0.20472977209518117, -0.21430915468432699, -0.21656349614594508, -0.2143710033505706, -0.21047592509949745, -0.20690277124897924, -0.20383621964071927, -0.2009312504622697, -0.19832906347760781, -0.19715218778669685, -0.19713511857310209, -0.19760258124272423, -0.19821523236790017, -0.19886735995941604, -0.19945089045017395, -0.19990124934465239, -0.20016164831872438, -0.20035891530554095, -0.20053280091253389, -0.20062879519403856, -0.20061453716771055, -0.20050516254700712, -0.20034457406881809, -0.2001804893235041, -0.20004787348381434, -0.19996313076460828, -0.19992645686226715, -0.19992808844809359, -0.19995462843450704, -0.19999336711919727, -0.20003430609120498, -0.20007059456815027, -0.20009823569348062, -0.20011560727571823, -0.20012297891084577, -0.20012200982886325, -0.2001151929510443]

        
G2 = [-0.39166666666666666, -0.38244067921398905, -0.33350466136630463, -0.29771023004255276, -0.27605439066140891, -0.25986156218997503, -0.24502185018573649, -0.23179262432952102, -0.21981564668803996, -0.20870707082936543, -0.19877739883776596, -0.18980922837977954, -0.1730801167400583, -0.16306400164063853, -0.17798470933316624, -0.19295540736943456, -0.20236705173335867, -0.20695767548229582, -0.20841025868691554, -0.20792102171307641, -0.20655350005113712, -0.2049200252815718, -0.20310627030300929, -0.20105983339290284, -0.19937394568595421, -0.1985391750807548, -0.19836389978207072, -0.19850305023816406, -0.19877764028800521, -0.19910928130800573, -0.19943705712074122, -0.19970344172627957, -0.19991076989513878, -0.20010020127344161, -0.20025937786728282, -0.20035087292654716, -0.20035829921368681, -0.20029606557358018, -0.2001960691549278, -0.20009096093557327, -0.20000371608352077, -0.19994495433034862, -0.1999153566524963, -0.19990981826568235, -0.19992106419902292, -0.19994189853498176, -0.19996624091198109, -0.19998946016949451, -0.20000842303436642, -0.20002144460691332, -0.20002815561319481]
        
G3 = [-0.44166666666666665, -0.37631169657400343, -0.33000044342859486, -0.30586045469008505, -0.28843572253009925, -0.27215308978603797, -0.25712951540331214, -0.2431608296216613, -0.23032023651386366, -0.21845468734566184, -0.20735123704254327, -0.19740397194806383, -0.18598295640643708, -0.16675980728412546, -0.1695157503295662, -0.18328608714846414, -0.19485758921965809, -0.2023136827680253, -0.20625610365424957, -0.20758116234856749, -0.20721445398906116, -0.20603406830101498, -0.20450262810357958, -0.2026769581530177, -0.20074012124164839, -0.19931160538111795, -0.1986360630235105, -0.19848511941161767, -0.1986009104317408, -0.19885490669335967, -0.19916542732465842, -0.19946678238299068, -0.19971209593781808, -0.19991912886144106, -0.20010584307496243, -0.20024959409137078, -0.20032160254398582, -0.2003158316568579, -0.20025051539345248, -0.20015561158283834, -0.20005952955567172, -0.19998144295750273, -0.19992977821660762, -0.19990457708729134, -0.19990104785520091, -0.19991257153955966, -0.19993258231861782, -0.19995548502852453, -0.19997700760886039, -0.19999429663472679, -0.20000588800224417]


# Only compare those that belong to this process id
G = [G0, G1, G2, G3]


success = True

for i in range(4):
    if tri_ids[i] > -1:
        #print 'myid = %g, allclose(gauge_values[%g], G%g) = %g' % (myid, i,i, allclose(gauge_values[i], G[i]))
        success = success and num.allclose(gauge_values[i], G[i])


if success:
    print 'Successful completion on processor ',myid
else:
    print 'Failure on processor ',myid
