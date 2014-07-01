"""
Supercritical flow over a bump.


Steve Roberts, ANU 2014
"""
from numpy import zeros
from scipy.optimize import fsolve
from anuga import g


qA = 10.0  # This is the imposed momentum
hx = 0.5   # This is the water height downstream

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
        h[i] = fsolve(find_hL, 0.5)


    return h,z


