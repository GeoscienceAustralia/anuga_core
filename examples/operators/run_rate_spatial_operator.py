"""Simple water flow example using ANUGA

Water flowing down a channel with a topography that varies with time
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import numpy
from anuga import rectangular_cross
from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary

from anuga import indent


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
length = 24.
width = 5.
dx = dy = 0.2 #.1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
                                               len1=length, len2=width)
domain = Domain(points, vertices, boundary)
domain.set_name() # Output name based on script
print (domain.statistics())


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    """Complex topography defined by a function of vectors x and y."""

    z = -x/100

    N = len(x)
    for i in range(N):
        # Step
        if 2 < x[i] < 4:
            z[i] += 0.4 - 0.05*y[i]

        # Permanent pole
        #if (x[i] - 8)**2 + (y[i] - 2)**2 < 0.4**2:
        #    z[i] += 1


#        # Dam
#        if 12 < x[i] < 13 and y[i] > 3.0:
#            z[i] += 0.4
#
#        if 12 < x[i] < 13 and y[i] < 2.0:
#            z[i] += 0.4


#        # Dam
#        if 12 < x[i] < 13:
#            z[i] += 0.4

            
    return z


def pole_increment(x,y,t):
    """This provides a small increment to a pole located mid stream
    For use with variable elevation data
    """

    z = 0.0*x
    
    if t<10.0:
        return z
    
    N = len(x)
    for i in range(N):
        # Pole 1
        if (x[i] - 12)**2 + (y[i] - 3)**2 < 0.4**2:
            z[i] += 0.1

    for i in range(N):
        # Pole 2
        if (x[i] - 14)**2 + (y[i] - 2)**2 < 0.4**2:
            z[i] += 0.05

    return z


def pole(t):

    if t<10:
        return 0.0
    elif t>15:
        return 0.0
    else:
        return 0.05


domain.set_quantity('elevation', topography)           # elevation is a function
domain.set_quantity('friction', 0.01)                  # Constant friction
domain.set_quantity('stage', expression='elevation')   # Dry initial condition

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = Dirichlet_boundary([0.4, 0, 0])          # Inflow
Br = Reflective_boundary(domain)              # Solid reflective wall
Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow

domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
polygon1 = [ [10.0, 0.0], [11.0, 0.0], [11.0, 5.0], [10.0, 5.0] ]
polygon2 = [ [12.0, 2.0], [13.0, 2.0], [13.0, 3.0], [12.0, 3.0] ]


from anuga.operators.rate_operators import Rate_operator

op1 = Rate_operator(domain, rate=10.0, polygon=polygon2)

area1 = numpy.sum(domain.areas[op1.indices])
Q1 = 10.0*area1
print ('op1 Q ',Q1)

op2 = Rate_operator(domain, rate=10.0, radius=0.5, center=(10.0, 3.0))

area2 = numpy.sum(domain.areas[op2.indices])
Q2 = 10.0*area2
print ('op2 Q ',Q2)


def rain(x,y,t):
    """Function to calculate "rain"
    input x,y should be considered to be numpy arrays
    abd t a scalar
    """
    if t<=4.0:
        return (x+y)*1.0
    else:
        return 0.0*x

factor = 1e-3
op3 = Rate_operator(domain, rate = rain, factor=factor)
Q3 = numpy.sum(op3.get_spatial_rate()*domain.areas)*factor


#op3()
#domain.fractional_step_operators.remove(op3)


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
accum = 0.0
yieldstep = 0.1
finaltime = 5.0
for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
    domain.print_timestepping_statistics()
    domain.print_operator_timestepping_statistics()

    stage = domain.get_quantity('stage')
    elev  = domain.get_quantity('elevation')
    height = stage - elev

    print (indent + 'Integral = ', height.get_integral())
    print (indent + 'Exact accumultion = ', accum)
    
    dd = max(min(yieldstep,4.0-t),0.0)
    accum += (Q1+Q2)*yieldstep + dd*Q3






