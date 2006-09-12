"""Example of shallow water wave equation.

Flat bed with rotational wind stress

"""

######################
# Module imports 
#
from mesh_factory import rectangular
from shallow_water import Domain, Dirichlet_boundary, Wind_stress

#Create basic mesh
N = 20
length = 200
points, vertices, boundary = rectangular(N, N, length, length)

#Create shallow water domain
domain = Domain(points, vertices, boundary)
domain.smooth = True
domain.visualise = False
domain.store = True
domain.default_order=2
domain.set_name('wind_rotation')

#Set initial conditions
h = 1.0
domain.set_quantity('elevation', 0.0)
domain.set_quantity('stage', h)
domain.set_quantity('friction', 0.01)


#Variable windfield implemented using functions
def speed(t,x,y):
    """Large speeds halfway between center and edges
    Low speeds at center and edges
    """

    from math import pi
    from Numeric import sqrt, exp, cos

    c = (length/2, length/2)
    r = sqrt((x - c[0])**2 + (y - c[1])**2)/length
    factor = exp( -(r-0.15)**2 )

    #return 9000 * factor
    return 4000 * factor * (cos(t*2*pi/150) + 2)


def phi(t,x,y):
    """Rotating field
    """
    
    from math import pi
    from Numeric import sqrt, exp, cos, arctan2, choose, less

    c = (length/2, length/2)
    xx = (x - c[0])/length
    yy = (y - c[1])/length
    angle = arctan2(yy,xx)

    #Take normal direction (but reverse direction every 50 seconds)
    #if sin(t*2*pi/100) < 0:
    #    sgn = -1
    #else:
    #    sgn = 1
    #angle += sgn*pi/2
    angle -= pi/2    

    #Convert to degrees and return
    return angle/pi*180

    
domain.forcing_terms.append( Wind_stress(speed, phi) )


#Add lateral wind gusts bearing 25 deg
def gust(t,x,y): 
    from math import sin, pi
    from Numeric import zeros, ones, Float

    N = len(x)

    tt = sin(2*pi*t/200)

    if tt > 0.9:
        return 6000*tt*ones(N, Float)
    else:
        return zeros(N, Float)
    
domain.forcing_terms.append(Wind_stress(gust, 25))


######################
# Boundary conditions
#Br = Reflective_boundary(domain)
Bd = Dirichlet_boundary([h, 0, 0])
domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})

######################
#Evolution
for t in domain.evolve(yieldstep = 0.5, finaltime = 1000):
    domain.write_time()

print 'Done'    

