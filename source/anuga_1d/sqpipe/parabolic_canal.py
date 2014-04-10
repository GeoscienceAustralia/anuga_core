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


z_infty = 10.0       ## max equilibrium water depth at lowest point.
L_x = 2500.0         ## width of channel
A0 = 0.5*L_x                  ## determines amplitudes of oscillations
omega = sqrt(2*g*z_infty)/L_x ## angular frequency of osccilation

def analytic_canal(x,t):
    u = numpy.zeros_like(x)    ## water velocity
    h = numpy.zeros_like(x)    ## water depth

    ## Define Basin Bathymetry
    z = numpy.zeros_like(x) ## elevation of basin
    w = numpy.zeros_like(x)   ## elevation of water surface

    z[:] = z_infty*(x**2/L_x**2)
    u[:] = -A0*omega*sin(omega*t)
    w[:] = numpy.maximum(z_infty+2*A0*z_infty/L_x*cos(omega*t)*(x/L_x-0.5*A0/(L_x)*cos(omega*t)),z)
    h[:] = numpy.maximum(w-z, 0.0)

    T = 2.0*pi/omega

    return u,h,w,z, T


def stage(x):
    t=0.0
    u,h,w,z,T = analytic_canal(x,t)

    return w

def elevation(x):
    t=0.0
    u,h,w,z,T = analytic_canal(x,t)

    return z


def height(x):
    t=0.0
    u,h,w,z,T = analytic_canal(x,t)

    return h

def width(x):
    return numpy.ones_like(x)

def top(x):
    return 4.0 * numpy.ones_like(x)

def area(x):
    return height(x)*width(x)


#===============================================================================

def get_domain():
    N = 100
    print "Evaluating domain with %d cells" %N
    
    domain = dom.Domain(*uniform_mesh(N, x_0 = -2.0*L_x, x_1 = 2.0*L_x), bulk_modulus = 100.0)

    domain.set_spatial_order(2)
    domain.set_timestepping_method('rk2')
    domain.set_CFL(1.0)
    domain.set_limiter("vanleer")

    domain.set_beta(1.0)

    domain.set_quantity('area', area)
    domain.set_quantity('stage', stage)
    domain.set_quantity('elevation',elevation)
    domain.set_quantity('width',width)
    domain.set_quantity('top',top)

    Br = dom.Reflective_boundary(domain)

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
    
    plot1 = pylab.subplot(321)
    zplot, = pylab.plot(x, z)
    wplot, = pylab.plot(x, w)
    ztplot, = pylab.plot(x, z+t)    
    plot1.set_ylim([-1,50])
    pylab.xlabel('Position')
    pylab.ylabel('Stage')

    plot2 = pylab.subplot(322)
    hplot, = pylab.plot(x, h)
    tplot, = pylab.plot(x, t)
    plot2.set_ylim([-1,20])
    pylab.xlabel('Position')
    pylab.ylabel('Height')

    plot3 = pylab.subplot(323)
    vplot, = pylab.plot(x, v)
    plot3.set_ylim([-30,30])
    pylab.xlabel('Position')
    pylab.ylabel('Velocity')

    plot4 = pylab.subplot(324)
    splot, = pylab.plot(x, s)
    plot4.set_ylim([-1,2])
    pylab.xlabel('Position')
    pylab.ylabel('State')

    plot5 = pylab.subplot(325)
    mplot, = pylab.plot(x, m)
    plot5.set_ylim([-1,1000])
    pylab.xlabel('Position')
    pylab.ylabel('Mass')

    plot6 = pylab.subplot(326)
    Mplot, = pylab.plot(x, M)
    plot6.set_ylim([-1,45000])
    pylab.xlabel('Position')
    pylab.ylabel('Total Mass')

    return zplot, wplot, ztplot, hplot, tplot, vplot, splot, Mplot, mplot
