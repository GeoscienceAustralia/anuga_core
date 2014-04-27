import os
import random
from math import sqrt, pow, pi 
from channel_domain import *
from numpy import allclose, array, zeros, ones, take, sqrt
from anuga_1d.config import g, epsilon

from anuga_1d.base.generic_mesh import uniform_mesh


print "varying width to mimic rotationally symmetric dam break"

# Define functions for initial quantities



def initialize_plotting(domain, style = '-k',
                        stage_lim = [-1.0, 40.0],
                        velocity_lim = [-5.0, 5.0]):

    import pylab
    pylab.ion()


    x = domain.get_vertices().flatten()

    z = domain.quantities['elevation'].vertex_values.flatten()
    w = domain.quantities['stage'].vertex_values.flatten()
    h = domain.quantities['height'].vertex_values.flatten()
    v = domain.quantities['velocity'].vertex_values.flatten()
    b = domain.quantities['width'].vertex_values.flatten()

    print x.shape
    print z.shape

    #-------------------------------
    # Top plot
    #-------------------------------
    domain.plot1 = pylab.subplot(211)

    domain.zplot, = pylab.plot(x, z)
    domain.wplot, = pylab.plot(x, w, style)

    domain.plot1.set_ylim(stage_lim)
    #pylab.xlabel('Position')
    pylab.ylabel('Stage')


    #-------------------------------
    # Bottom Plot
    #-------------------------------
    domain.plot3 = pylab.subplot(212)

    domain.vplot, = pylab.plot(x, v, style)

    domain.plot3.set_ylim(velocity_lim)

    pylab.xlabel('Position')
    pylab.ylabel('Velocity')


def update_plotting(domain):

    import pylab

    #x = domain.get_vertices().flatten()
    z = domain.quantities['elevation'].vertex_values.flatten()
    w = domain.quantities['stage'].vertex_values.flatten()
    h = domain.quantities['height'].vertex_values.flatten()
    v = domain.quantities['velocity'].vertex_values.flatten()
    b = domain.quantities['width'].vertex_values.flatten()


    domain.zplot.set_ydata(z)
    domain.wplot.set_ydata(w)
    #domain.bplot.set_ydata(b)
    domain.vplot.set_ydata(v)

    pylab.draw()


def hold_plotting(domain,save=None):

    update_plotting(domain)
    import pylab

    pylab.ioff()

    if save != None:
        file = save+".pdf"
        pylab.savefig(file)

    pylab.show()



def finalize_plotting(domain):

    pass





def bed(x):
    y = zeros(len(x),'f')
    return y


def width(x):
    
    return x*2.0*pi



def initial_area(x):

    a_width = width(x)

    y = numpy.where (x <= 50.0,  15.0*a_width, 2.0*a_width )

    return y


import time


# Define cells for finite volume and their size
N = 100

domain = Domain(*uniform_mesh(N, x_0 = 0.0, x_1 = 100.0))


# Set initial values of quantities - default to zero
#domain.set_quantity('stage',6.0)
domain.set_quantity('elevation',bed)
domain.set_quantity('width',width)
domain.set_quantity('area', initial_area)


#domain.setstageflag = True
# Set boundry type, order, timestepping method and limiter
Br = Reflective_boundary(domain)
domain.set_boundary({'left':Br, 'right':Br})
domain.order = 2
domain.set_timestepping_method('rk2')
domain.set_CFL(1.0)
domain.set_limiter("vanleer")


#AreaC = domain.quantities['area'].centroid_values
#BedC = domain.quantities['elevation'].centroid_values
#WidthC = domain.quantities['width'].centroid_values
##
#AreaC[:] = (8.0 - BedC)* WidthC


# Start timer
t0 = time.time()
i=0

#print 'elevation vertex values'
#print domain.quantities['elevation'].vertex_values
#print 'stage vertex values'
#print domain.quantities['stage'].vertex_values
#print 'area vertex values'
#print domain.quantities['area'].vertex_values
#print 'width vertex values'
#print domain.quantities['width'].vertex_values


domain.distribute_to_vertices_and_edges()



# Set final time and yield time for simulation
finaltime = 2.0
yieldstep = 0.2

initialize_plotting(domain, style = '.k', stage_lim = [0.0, 16.0],
                            velocity_lim = [-10.0, 10.0])

for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    domain.write_time()

    update_plotting(domain)

#-------------------------------------------------------------
# Fine grid solution
#-------------------------------------------------------------

domain1 = Domain(*uniform_mesh(1000, x_0 = 0.0, x_1 = 100.0))


# Set initial values of quantities - default to zero
#domain.set_quantity('stage',6.0)
domain1.set_quantity('elevation',bed)
domain1.set_quantity('width',width)
domain1.set_quantity('area', initial_area)


#domain.setstageflag = True
# Set boundry type, order, timestepping method and limiter
Br = Reflective_boundary(domain1)
domain1.set_boundary({'left':Br, 'right':Br})
domain1.order = 2
domain1.set_timestepping_method('rk2')
domain1.set_CFL(1.0)
domain1.set_limiter("vanleer")

finaltime = 2.0
yieldstep = finaltime

initialize_plotting(domain1, stage_lim = [0.0, 16.0],
                            velocity_lim = [-10.0, 10.0])

for t in domain1.evolve(yieldstep = yieldstep, finaltime = finaltime):
    domain1.write_time()

    update_plotting(domain1)


hold_plotting(domain1, save="circular_dam_break_well_balanced")




