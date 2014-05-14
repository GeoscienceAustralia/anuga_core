"""
Thacker's analytical solution for a paraboloid oscillation
of water on a paraboloid basin

Sudi Mungkasi, ANU 2012
"""
from numpy import zeros, sin, cos, sqrt
from anuga import g


# Set bed elevation
def bed_elevation(x,y,D0=1000.,L=2500.,R0=2000.,g=g):
    n = x.shape[0]
    z = 0*x
    for i in range(n):
        r = sqrt(x[i]*x[i] + y[i]*y[i])
        z[i] = -D0*(1.-r*r/L/L)
    return z

def analytic_sol(x,y,t,D0=1000.,L=2500.,R0=2000.,g=g):
    A = (L**4 - R0**4)/(L**4 + R0**4)
    omega = 2./L*sqrt(2.*g*D0)
    z = bed_elevation(x,y)
    n = x.shape[0]
    w = 0*x
    u = 0*x
    h = 0*x    
    for i in range(n):
        r = sqrt(x[i]*x[i] + y[i]*y[i])
        w[i] = D0*((sqrt(1-A*A))/(1.-A*cos(omega*t))
                -1.-r*r/L/L*((1.-A*A)/((1.-A*cos(omega*t))**2)-1.))
        u[i] = 0.5*omega*r*A*sin(omega*t) / (1.0-A*cos(omega*t))
        if x[i] < 0.0:
            u[i] = -u[i]
        h[i] = w[i] - z[i]
        if w[i] < z[i]:
            w[i] = z[i]
            u[i] = 0.0
            h[i] = 0.0
    return w,u,h
