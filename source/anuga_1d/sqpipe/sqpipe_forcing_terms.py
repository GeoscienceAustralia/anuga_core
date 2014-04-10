#! /usr/bin/python


"""
Forcing terms for the shallow water equations, gravity, friction etc
"""

__author__="Stephen Roberts"
__date__ ="$05/06/2010 5:49:35 PM$"


def gravity_for_loops(domain):
    """Apply gravitational pull in the presence of bed slope
    """

    from anuga_1d.base.util import gradient

    xmom  = domain.quantities['xmomentum'].explicit_update

    Stage = domain.quantities['stage']
    Elevation = domain.quantities['elevation']

    h = Stage.vertex_values - Elevation.vertex_values
    b = Elevation.vertex_values
    w = Stage.vertex_values

    x = domain.get_vertex_coordinates()
    g = domain.g

    for k in range(domain.number_of_elements):
        avg_h = sum( h[k,:] )/2

        #Compute bed slope
        #x0, y0, x1, y1, x2, y2 = x[k,:]
        x0, x1 = x[k,:]
        #z0, z1, z2 = v[k,:]
        b0, b1 = b[k,:]


        #zx, zy = gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2)
        bx = gradient(x0, x1, b0, b1)

        #Update momentum (explicit update is reset to source values)
        xmom[k] += -g*bx*avg_h
        #xmom[k] = -g*bx*avg_h
        #stage[k] = 0.0


def gravity(domain):
    """Apply gravitational pull in the presence of bed slope
    """

    d  = domain.quantities['discharge'].explicit_update

    Elevation = domain.quantities['elevation']
    Height    = domain.quantities['height']

    hc = Height.centroid_values
    z  = Elevation.vertex_values

    x = domain.vertices
    g = domain.g

    x0 = x[:,0]
    x1 = x[:,1]

    z0 = z[:,0]
    z1 = z[:,1]

    zx = (z1-z0)/(x1-x0)

    d[:] = d - g*zx*hc

def gravity_press(domain):
    """Apply gravitational pull in the presence of bed slope
    """

    d  = domain.quantities['discharge'].explicit_update

    Elevation = domain.quantities['elevation']
    Height    = domain.quantities['height']

    hc = Height.centroid_values
    z  = Elevation.vertex_values

    x = domain.vertices
    g = domain.g

    x0 = x[:,0]
    x1 = x[:,1]

    z0 = z[:,0]
    z1 = z[:,1]

    zx = (z1-z0)/(x1-x0)

    d[:] = d - g*zx*hc

def manning_friction(domain):
    """Apply (Manning) friction to water momentum
    """

    from math import sqrt, fabs

    eta = domain.quantities['friction'].centroid_values

    hw = domain.quantities['area'].centroid_values
    w  = domain.quantities['width'].centroid_values
    uhw = domain.quantities['discharge'].centroid_values
    
    h = hw/w
    
    uhw_update = domain.quantities['discharge'].semi_implicit_update
    #ymom_update = domain.quantities['ymomentum'].semi_implicit_update

    N = domain.number_of_elements
    eps = 1.0e-10
    g = domain.g

    #print 'mannings'
    for k in range(N):
        if eta[k] >= eps:
            if h[k] >= eps:
            	#S = -g * eta[k]**2 * sqrt((uh[k]**2 + vh[k]**2))
                S = -g * eta[k]**2 * fabs(uhw[k])
            	S /= h[k]**(7.0/3)
                S /= w[k]

            	#Update momentum
            	uhw_update[k] += S*uhw[k]
            	#ymom_update[k] += S*vh[k]

def linear_friction(domain):
    """Apply linear friction to water momentum

    Assumes quantity: 'linear_friction' to be present
    """

    from math import sqrt

    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].centroid_values
    h = w-z

    uh = domain.quantities['xmomentum'].centroid_values
#    vh = domain.quantities['ymomentum'].centroid_values
    tau = domain.quantities['linear_friction'].centroid_values

    xmom_update = domain.quantities['xmomentum'].semi_implicit_update
#    ymom_update = domain.quantities['ymomentum'].semi_implicit_update

    N = domain.number_of_elements
    eps = domain.minimum_allowed_height
    g = domain.g #Not necessary? Why was this added?

    for k in range(N):
        if tau[k] >= eps:
            if h[k] >= eps:
            	S = -tau[k]/h[k]

            	#Update momentum
            	xmom_update[k] += S*uh[k]
 #           	ymom_update[k] += S*vh[k]



def check_forcefield(f):
    """Check that f is either
    1: a callable object f(t,x,y), where x and y are vectors
       and that it returns an array or a list of same length
       as x and y
    2: a scalar
    """


    if callable(f):
        #N = 3
        N = 2
        #x = ones(3, numpy.float)
        #y = ones(3, numpy.float)
        x = ones(2, numpy.float)
        #y = ones(2, numpy.float)

        try:
            #q = f(1.0, x=x, y=y)
            q = f(1.0, x=x)
        except Exception, e:
            msg = 'Function %s could not be executed:\n%s' %(f, e)
	    #FIXME: Reconsider this semantics
            raise msg

        try:
            q = numpy.array(q, numpy.float)
        except:
            msg = 'Return value from vector function %s could ' %f
            msg += 'not be converted into a numpy array of numpy.floats.\n'
            msg += 'Specified function should return either list or array.'
            raise msg

        #Is this really what we want?
        msg = 'Return vector from function %s ' %f
        msg += 'must have same lenght as input vectors'
        assert len(q) == N, msg

    else:
        try:
            f = float(f)
        except:
            msg = 'Force field %s must be either a scalar' %f
            msg += ' or a vector function'
            raise msg
    return f


