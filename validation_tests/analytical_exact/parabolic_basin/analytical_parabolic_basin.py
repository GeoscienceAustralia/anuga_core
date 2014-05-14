"""
Analytical solution of Thacker for oscillation on
a parabolic basin.
"""
from numpy import zeros,sqrt,sin,cos
from anuga import g

def analytic_cannal(x,t, D0=4., L=10., A=2., g=g):
    omega = sqrt(2*D0*g)/L
    N = len(x)
    u = zeros(N)              ## water velocity
    h = zeros(N)              ## water depth
    ## Define Basin Bathymetry
    z = zeros(N)            ## elevation of basin
    w = zeros(N)              ## elevation of water surface
    for i in range(N):
        z[i] = D0*(x[i]**2/L**2)
        u[i] = -A*omega*sin(omega*t)
        w[i] = D0 + 2*A*D0/L**2 * cos(omega*t)*( x[i] - A/2.*cos(omega*t))
        if w[i] <= z[i] :
            u[i] = 0.0
            w[i] = z[i]
    h = w - z
    return u,h,w,z
