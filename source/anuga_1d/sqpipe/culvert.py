import os
from math import sqrt, pi, sin, cos
import numpy
import time

from anuga_1d.config import g, epsilon
from anuga_1d.base.generic_mesh import uniform_mesh
import anuga_1d.sqpipe.sqpipe_domain as dom

#===============================================================================
# setup problem
#===============================================================================

L_x = 50.0         ## length of channel

def elevation(x):
    z = -x*0.1
    return z

def stage(x):

    w = elevation(x) + 1.0

    return w


def height(x):

    h = stage(x) - elevation(x)

    return h

def width(x):
    return numpy.ones_like(x)

def top(x):

    t = numpy.ones_like(x)*4

    t = numpy.where(x<-40, 20, t)
    t = numpy.where(x>40, 20, t)
    return t

def area(x):
    return height(x)*width(x)

def friction(x):
    return numpy.ones_like(x)*0.01


#===============================================================================

def get_domain():
    N = 20
    print "Evaluating domain with %d cells" %N

    points, boundary = uniform_mesh(N, x_0 = -L_x, x_1 = L_x)

    domain = dom.Domain(points, boundary, bulk_modulus = 100.0)

    domain.set_spatial_order(2)
    domain.set_timestepping_method('rk2')
    domain.set_CFL(0.5)
    domain.set_limiter("vanleer")

    domain.set_beta(1.0)

    domain.set_quantity('area', area)
    domain.set_quantity('stage', stage)
    domain.set_quantity('elevation',elevation)
    domain.set_quantity('width',width)
    domain.set_quantity('top',top)
    domain.set_quantity('friction',friction)

    Br = dom.Reflective_boundary(domain)
    Bt = dom.Transmissive_boundary(domain)

    domain.set_boundary({'left': Br, 'right' : Br})

    return domain

def animate_domain(domain, yieldstep, finaltime):
    import pylab

    pylab.ion()

    x, z, w, h, v, t, s, m, M = get_quantities(domain)

    zplot, wplot, ztplot, hplot, tplot, vplot, splot, Mplot, mplot = make_plots(x, z, w, h, v, t, s, m, M)

    for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):

        domain.write_time()
        x, z, w, h, v, t, s, m, M = get_quantities(domain)
        
        zplot.set_ydata(z)
        ztplot.set_ydata(z+t)
        wplot.set_ydata(w)
        hplot.set_ydata(h)
        tplot.set_ydata(t)
        vplot.set_ydata(v)
        splot.set_ydata(s)
        mplot.set_ydata(m)
        Mplot.set_ydata(M)


        pylab.ion()
        pylab.draw()




def plot_domain(domain):
    import pylab

    x, z, w, h, v, t, s, m, M = get_quantities(domain)

    pylab.ioff()
    pylab.hold(False)

    make_plots(x, z, w, h, v, t, s, m, M)
    
    pylab.show()

def write_domain(domain, outfile):
    x = domain.get_centroids()
    z = domain.get_quantity('elevation', 'centroids')
    w = domain.get_quantity('stage', 'centroids')

    f = open(outfile, 'w')
    for i in range(len(x)):
        f.write("%s %s %s\n" % (x[i], z[i], w[i]))
    f.close()

def get_quantities(domain):
    x = domain.get_centroids()
    z = domain.get_quantity('elevation', 'centroids')
    w = domain.get_quantity('stage', 'centroids')
    h = domain.get_quantity('height', 'centroids')
    v = domain.get_quantity('velocity', 'centroids')
    t = domain.get_quantity('top', 'centroids')
    s = domain.state
    m = domain.get_mass()
    M = m.sum() * numpy.ones_like(x)

    return x, z, w, h, v, t, s, m, M

def make_plots(x, z, w, h, v, t, s, m, M):
    import pylab

    #fig = pylab.gcf()
    #fig.set_size_inches(12,12, forward=True)

    


    plot1 = pylab.subplot(321)
    zplot, = pylab.plot(x, z)
    wplot, = pylab.plot(x, w)
    ztplot, = pylab.plot(x, z+t)
    plot1.set_xlim([-60,60])
    plot1.set_ylim([-10,10])
    pylab.xlabel('Position')
    pylab.ylabel('Stage')

    plot2 = pylab.subplot(322)
    hplot, = pylab.plot(x, h)
    tplot, = pylab.plot(x, t)
    plot2.set_xlim([-60,60])
    plot2.set_ylim([-1,5])
    pylab.xlabel('Position')
    pylab.ylabel('Height')

    plot3 = pylab.subplot(323)
    vplot, = pylab.plot(x, v)
    plot3.set_xlim([-60,60])
    plot3.set_ylim([-6,6])
    pylab.xlabel('Position')
    pylab.ylabel('Velocity')

    plot4 = pylab.subplot(324)
    splot, = pylab.plot(x, s)
    plot4.set_xlim([-60,60])
    plot4.set_ylim([-1,2])
    pylab.xlabel('Position')
    pylab.ylabel('State')

    plot5 = pylab.subplot(325)
    mplot, = pylab.plot(x, m)
    plot5.set_xlim([-60,60])
    plot5.set_ylim([-1,10])
    pylab.xlabel('Position')
    pylab.ylabel('Mass')

    plot6 = pylab.subplot(326)
    Mplot, = pylab.plot(x, M)
    plot6.set_xlim([-60,60])
    plot6.set_ylim([-1,450])
    pylab.xlabel('Position')
    pylab.ylabel('Total Mass')




    return zplot, wplot, ztplot, hplot, tplot, vplot, splot, Mplot, mplot
