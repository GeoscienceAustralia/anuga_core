"""
Supercritical flow over a bump.


Steve Roberts, ANU 2014
"""
from numpy import zeros
from scipy.optimize import fsolve
from anuga import g


qA = 1.0  # This is the imposed momentum
hx = 1.0   # This is the water height downstream

def analytic_sol(x):    

    def elevation(x):
        z_b = zeros(len(x))
        e_w=2.0 # Width of expansion region (where depth changes)
        for i in range(len(x)):
            if (x[i] <= 10.0-e_w/2.0):
                z_b[i] = 0.2
            elif (10.0-e_w/2.0 < x[i] < 10.0+e_w/2.0):
                z_b[i] = 0.2 - 0.2*(x[i]-(10.-e_w/2.0))/e_w
            else:
                z_b[i] = 0.0
        return z_b
    z = elevation(x)

 
    def find_hL(h): #to find the water height at every spatial point
        return h**3 + (zb-qA**2/(2*g*hx**2)-hx)*h**2 + qA**2/(2*g)
    h = zeros(len(x))
    for i in range(len(x)):
        zb = z[i]
        h[i] = fsolve(find_hL, hx)


    return h,z


