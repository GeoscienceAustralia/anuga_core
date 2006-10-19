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
# Setup computational domain
#--------------------------------------------------------------------------
points, vertices, boundary = rectangular_cross(10, 10) # Basic mesh
domain = Domain(points, vertices, boundary) # Create domain

#--------------------------------------------------------------------------
# Setup initial conditions
#--------------------------------------------------------------------------

def topography(x,y): 
    return -x/2                              # linear bed slope

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

# Hack
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
    #domain.write_time()

    
    # Record time series at known points
    time.append(domain.get_time())
    
    stage = domain.get_quantity('stage')
    w = stage.get_values(interpolation_points=interpolation_points)
    
    for i, _ in enumerate(interpolation_points):
        gauge_values[i].append(w[i])


for i, (x,y) in enumerate(interpolation_points):
        
    try:
        from pylab import *
    except:
        pass
    else:
        ion()
        hold(False)
        plot(time, gauge_values[i], 'r.-')
        #time, predicted_gauge_values[i], 'k-')
        
        title('Gauge %d (%f,%f)' %(i,x,y))
        xlabel('time(s)')
        ylabel('stage (m)')    
        #legend(('Observed', 'Modelled'), shadow=True, loc='upper left')
        #savefig('Gauge_%d.png' %i, dpi = 300)
    
        raw_input('Next')
        


# Reference from sequential version (also available as a
# unit test in test_shallow_water_domain)
# Added Friday 13 October 2006 by Ole

G0 = ensure_numeric([-0.20000000000000001, -0.19999681443389281, -0.1986192343695303, -0.19147413648863046, -0.19132688908678019, -0.17642317476621105, -0.167376262630034, -0.16192452887426961, -0.15609171725778803, -0.15127107084302249, -0.14048864340360018, -0.19296484125327093, -0.19997006390580363, -0.19999999999937063, -0.19999999999937063, -0.19999999999938772, -0.19999999999938772, -0.19999999999938772, -0.19999999999938772, -0.19974288463035494, -0.19951636867991712, -0.19966301435195755, -0.19981082259800226, -0.19978575003960128, -0.19992942471933109, -0.19999999931029933, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989, -0.19999999999906989])

G1 = ensure_numeric([-0.29999999999999993, -0.29988962537199199, -0.29293904425532025, -0.28329367722887888, -0.25999146407696289, -0.22613875068011896, -0.21190705052094994, -0.19900707995208217, -0.18876305176191882, -0.18132447501091936, -0.17395459512711151, -0.15562414200985644, -0.16212999953643359, -0.18964422820514618, -0.20871181844346975, -0.21672207791083464, -0.21774940291862779, -0.21482868050219833, -0.21057786776704043, -0.20649663432591045, -0.20294932949211578, -0.19974459897911329, -0.19733648772704043, -0.19641404599824669, -0.19654095699184146, -0.19709942852191994, -0.19780873983410741, -0.19853259125123518, -0.19916495938961168, -0.19965391267799168, -0.19993539587158982, -0.2001383705551133, -0.20029344332295113, -0.20035349748150011, -0.20029886541561631, -0.20015541958920294, -0.19997273066429103, -0.19979879448668514, -0.19966016997024041, -0.19957558009501869, -0.19955725674938532, -0.19958083002853366, -0.19961752462568647, -0.19965296611330258, -0.19968998132634594, -0.19972532942208607, -0.19975372922008239, -0.19977196116929855, -0.19977951443660594, -0.19977792107284789, -0.19976991595502003])

G2 = ensure_numeric([-0.40000000000000002, -0.39011996186687281, -0.33359026016903887, -0.29757449757405952, -0.27594124995715791, -0.25970211955309436, -0.24482929492054245, -0.23156757139219822, -0.21956485769139392, -0.20844522129026694, -0.19856327660654355, -0.18962303467030903, -0.17371085465024955, -0.16429840256208336, -0.17793711732368575, -0.19287799702389993, -0.20236271260796762, -0.20700727993623128, -0.20847704371373174, -0.20796895600687262, -0.20653398626186478, -0.20480656169870676, -0.20295863990994492, -0.20100199602968896, -0.19940642689498472, -0.19858371478015749, -0.19838672154605322, -0.19851093923669558, -0.19878191998909323, -0.19910827645394291, -0.19943514333832094, -0.19971231361970535, -0.19992429278849655, -0.20010744405928019, -0.20025927002359642, -0.20034751667523681, -0.20035504591467249, -0.20029401385620157, -0.20019492358237226, -0.20008934249434918, -0.19999808924091636, -0.19993869218976712, -0.19991589568150098, -0.19991815777945968, -0.19993012995477188, -0.19994576118144997, -0.19996497193815974, -0.19998586151236197, -0.20000487253824847, -0.20001903000364174, -0.20002698661385457])

G3 = ensure_numeric([-0.45000000000000001, -0.37713945714588398, -0.33029565026933816, -0.30598209033945367, -0.28847101155177313, -0.27211191064563195, -0.25701544058818926, -0.24298945948410997, -0.23010402733784807, -0.21820351802867713, -0.20709938367218383, -0.19719881806182216, -0.18568281604361933, -0.16828653906676322, -0.16977310768235579, -0.1832707289594605, -0.19483524345250974, -0.20233480051649216, -0.20630757214159207, -0.20763927857964531, -0.20724458160595791, -0.20599191745446047, -0.20438329669495012, -0.20256105512496606, -0.20071269486729407, -0.19934403619901719, -0.19866860191898347, -0.19849975056296071, -0.19860870923007437, -0.19885838217851401, -0.19916422433758982, -0.19946861981642039, -0.19972267778871666, -0.19993013816258154, -0.20011063428833351, -0.20024891930311628, -0.20031882555219671, -0.20031326268593497, -0.20024881068472311, -0.20015443214902759, -0.20005669097631221, -0.19997542564643309, -0.19992564006223304, -0.19990746148869892, -0.19990923999172872, -0.19991956416813192, -0.19993484556273733, -0.1999538628054662, -0.19997381636620407, -0.19999130900268777, -0.20000388227457688])


# Only compare those that belong to this process id
G = [G0, G1, G2, G3]

for i in local_interpolation_points:
    msg = 'P%d, point #%d: Computed time series and reference time series are different: %s'\
          %(myid, i, gauge_values[i]-G[i])
    assert allclose(gauge_values[i], G[i]), msg

print 'P%d completed succesfully using points = %s' %(myid, local_interpolation_points)

