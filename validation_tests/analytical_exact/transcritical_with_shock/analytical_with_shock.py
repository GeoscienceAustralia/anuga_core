"""
Transcritical flow over a bump with a shock.
Ref1: Houghton & Kasahara, Nonlinear shallow fluid flow over an isolated ridge.
Comm. Pure and Applied Math. DOI:10.1002/cpa.3160210103

Ref2: Delestre et al, 2012, SWASHES: a compilation of shallow water
analytic solutions..., Int J Numer Meth Fluids, DOI:10.1002/fld.3741

Sudi Mungkasi, ANU 2012
"""
from numpy import zeros, linspace, array, sqrt
from scipy.optimize import fsolve
from pylab import plot,show,ylim
from anuga import g


qA  = 0.18  # This is the imposed momentum
hx  = 0.33  # This is the water height downstream

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
        return hx**3 + (-qA**2/(2*g*hM**2)-hM-zM)*hx**2 + qA**2/(2*g)
    hM = fsolve(find_hM, 0.5)

    def find_h(h): #to find the water height at every spatial point after hM is found
        return h**3 + (zb-qA**2/(2*g*hM**2)-hM-zM)*h**2 + qA**2/(2*g)

    def find_hL(h): #to find the water height at every spatial point after hM is found
        return h**3 + (zb-qA**2/(2*g*hx**2)-hx)*h**2 + qA**2/(2*g)

    def find_shock(Q):
        h1 = Q[0]
        h2 = Q[1]
        zs = Q[2]
        F = zeros(3)
        F[0] = h1**3 + (zs-qA**2/(2*g*hM**2)-hM-zM)*h1**2 + qA**2/(2*g)
        F[1] = h2**3 + (zs-qA**2/(2*g*hx**2)-hx)*h2**2 + qA**2/(2*g)
        F[2] = (qA**2)*(1.0/h1 - 1.0/h2) + 0.5*g*(h1**2 - h2**2)
        return F

    z_shock = fsolve(find_shock,array([0.2,0.3,0.005]))[2]

    def find_zshock(x):
        return 0.2 - 0.05*(x-10.0)**2 - z_shock
    x_shock = fsolve(find_zshock,11.5)[0]

    h=zeros(len(x))
    for i in range(len(x)):
        if x[i] <= 10.0:
            zb   = z[i]
            h[i] = fsolve(find_h,0.4)
        elif x[i] <= x_shock:
            zb   = z[i]
            h[i] = fsolve(find_h,0.1)
        elif x[i] <= 12.0:
            zb   = z[i]
            h[i] = fsolve(find_hL,0.3)
        else:
            h[i] = hx
    return h,z
        
##N = 10001
##L = 25.
##x = linspace(0.0,L,N)
##h,z = analytic_sol(x)
##plot(x,h+z,'b.-', x,z,'k-')
##ylim([0.0, 0.5])
##show()
