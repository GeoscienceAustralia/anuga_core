"""
Subcritical flow over a bump.
Ref1: Houghton & Kasahara, Nonlinear shallow fluid flow over an isolated ridge.
Comm. Pure and Applied Math. DOI:10.1002/cpa.3160210103

Ref2: Delestre et al, 2012, SWASHES: a compilation of shallow water
analytic solutions..., Int J Numer Meth Fluids, DOI:10.1002/fld.3741

Sudi Mungkasi, ANU 2012
"""
from numpy import zeros, linspace
from scipy.optimize import fsolve
from pylab import plot, ylim, show
from anuga import g

qA = 4.42  # This is the imposed momentum
hx = 2.0   # This is the water height downstream

def analytic_sol(x):    
    def elevation(x):
        z_b = zeros(len(x))
        for i in range(len(x)):
            if (8.0 <= x[i] <= 12.0):
                z_b[i] = 0.2 - 0.05*(x[i]-10.0)**2.0
            else:
                z_b[i] = 0.0
        return z_b
    z = elevation(x)

    def find_hL(h): #to find the water height at every spatial point
        return h**3 + (zb-qA**2/(2*g*hx**2)-hx)*h**2 + qA**2/(2*g)
    h = zeros(len(x))
    for i in range(len(x)):
        zb = z[i]
        h[i] = fsolve(find_hL, 2.0)
    return h,z

##N = 401
##L = 25.
##x = linspace(0.0,L,N)
##h,z = analytic_sol(x)
##plot(x,h+z, x,z)
##ylim([-0.1, 2.1])
##show()
