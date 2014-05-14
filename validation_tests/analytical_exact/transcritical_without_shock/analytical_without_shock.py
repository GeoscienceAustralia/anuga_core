"""
Transcritical flow over a bump without a shock.
Ref1: Houghton & Kasahara, Nonlinear shallow fluid flow over an isolated ridge.
Comm. Pure and Applied Math. DOI:10.1002/cpa.3160210103

Ref2: Delestre et al, 2012, SWASHES: a compilation of shallow water
analytic solutions..., Int J Numer Meth Fluids, DOI:10.1002/fld.3741

Sudi Mungkasi, ANU 2012
"""
from numpy import zeros, linspace
from scipy.optimize import fsolve
from pylab import plot, show
from anuga import g


q0  = 1.53  # This is the imposed momentum
h_d = 0.66  # This is the water height downstream

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
    zM= max(z)

    def find_hM(hM): #to find the water height at the maxima of the bump
        return h_d**3 + (-q0**2/(2*g*hM**2)-hM-zM)*h_d**2 + q0**2/(2*g)
    hM = fsolve(find_hM, 0.5)

    def find_h(h): #to find the water height at every spatial point after hM is found
        return h**3 + (zb-q0**2/(2*g*hM**2)-hM-zM)*h**2 + q0**2/(2*g)
    h = zeros(len(x))
    for i in range(len(x)):
        zb = z[i]
        #h[i] = fsolve(find_h, 1.0)
        if x[i] < 10:
            h[i] = fsolve(find_h, 1.0)
        else:
            h[i] = fsolve(find_h, 0.4)
    return h, z

##N = 401
##L = 25.
##x = linspace(0.0,L,N)
##h,z=analytic_sol(x)
##plot(x,h+z, x,z)
##plot(x, 1.53/h)
##show()
