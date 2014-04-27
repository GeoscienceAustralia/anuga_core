#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$15/06/2010 5:02:36 PM$"


import os
from math import sqrt, pow, pi

import numpy

from anuga_1d.channel.channel_domain import *
from anuga_1d.config import g, epsilon
from anuga_1d.base.generic_mesh import uniform_mesh


def run_evolve():


    print "Channel Flow 1 Padarn Test"

    # Define functions for initial quantities
    def initial_area(x):
        return 1.4691*width(x)

    def width(x):
        x1=(x/1000)*(x/1000)
        x2=x1*(x/1000)
        x3=x2*(x/1000)
        return 10-64*(x1-2*x2+x3)

    def bed(x):
        y = numpy.ones(len(x),numpy.float)

        return numpy.where( (x<525) & (x>475),y,0.0)


    def initial_discharge(x):
        return 20

    import time

    # Set final time and yield time for simulation
    finaltime = 50.0
    yieldstep = 10.0

    # Length of channel (m)
    L = 1000.0
    # Define the number of cells
    N = 200

    # Create domain with centroid points as defined above
    domain = Domain(*uniform_mesh(N))


    # Set initial values of quantities - default to zero
    domain.set_quantity('area', initial_area)
    domain.set_quantity('width',width)
    domain.set_quantity('elevation',bed)
    domain.set_quantity('discharge',initial_discharge)

    # Set boundry type, order, timestepping method and limiter
    Bd = Dirichlet_boundary([14,20,0,1.4,20/14,9,1.4])
    domain.set_boundary({'left': Bd , 'right' : Bd })


    domain.order = 2
    domain.set_timestepping_method('rk2')
    domain.set_CFL(1.0)
    domain.set_limiter("vanleer")
    #domain.set_limiter("minmod")
    #domain.h0=0.0001

    # Start timer
    t0 = time.time()

    for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
        domain.write_time()



import cProfile
cProfile.run('run_evolve()', 'evolve_prof')


import pstats
p = pstats.Stats('evolve_prof')

#p.strip_dirs().sort_stats(-1).print_stats(20)

p.sort_stats('cumulative').print_stats(30)





