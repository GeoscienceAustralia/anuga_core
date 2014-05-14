"""
MacDonald's steady flow on a short channel.

Ref1: MacDonald I, Baines MJ, Nichols NK, Samuels PG. Analytic benchmark solutions for open-channel
flows. Journal of Hydraulic Engineering November 1997; 123(11):1041-1045.
DOI: 10.1061/(ASCE)0733-9429(1997)123:11(1041).

Ref2: Delestre et al, 2012, SWASHES: a compilation of shallow water
analytic solutions..., Int J Numer Meth Fluids.
DOI:10.1002/fld.3741

Sudi Mungkasi, ANU 2012
"""
from numpy import zeros, linspace, array, sqrt
from scipy.integrate import quad
from pylab import plot,show,ylim
from anuga import g


L     = 100.  # Channel length
q = q_s = 2.    # Steady discharge
q_up  = q_s   # Upstream discharge
#h_down= h_L   # Downstream water height

#Parameters for water height
a1 = 0.674202
a2 = 21.7112
a3 = 14.492
a4 = 1.4305
x_s=200./3   # Specified position of a shock
n  = 0.0328  # Manning's friction coefficient


def water_height(x):
    N = len(x)
    h = zeros(N)
    for i in range(N):
        if x[i] <= x_s:
            h[i] = (4./g)**(1./3) * (4./3 - x[i]/100.) - 9.*x[i]/1000*(x[i]/100 - 2./3)
        else:
            r = x[i]/100. - 2./3.
            h[i] = (4./g)**(1./3) * (a1*r**4 + a1*r**3 - a2*r**2 + a3*r + a4)
    return h

##def diff_water_height(x):
##    N = len(x)
##    hx = zeros(N)
##    for i in range(N):
##        if x[i] <= x_s:
##            hx[i] = (4./g)**(1./3)*(-0.01) - (9./1000)*(x[i]/100 - 2./3) - (9.*x[i]/1000)*0.01
##        else:
##            r  = x[i]/100.-2./3
##            rx = 0.01
##            hx[i] = (4./g)**(1./3)*(a1*4.*r**3*rx + a1*3.*r**2*rx - a2*2.*r*rx + a3*rx)
##    return hx

def bed(x):
    N = len(x)
    Z = zeros((N,2))
    def func1(x):
        h = (4./g)**(1./3) * (4./3 - x/100.) - 9.*x/1000*(x/100 - 2./3)
        hx = (4./g)**(1./3)*(-0.01) - (9./1000)*(x/100 - 2./3) - (9.*x/1000)*0.01                
        return (q**2/(g*h**3) - 1.)*hx - n**2*q*abs(q)/h**(10./3)
    def func2(x):
        r = x/100. - 2./3.
        rx = 0.01
        h = (4./g)**(1./3) * (a1*r**4 + a1*r**3 - a2*r**2 + a3*r + a4)
        hx = (4./g)**(1./3)*(a1*4.*r**3*rx + a1*3.*r**2*rx - a2*2.*r*rx + a3*rx)                
        return (q**2/(g*h**3) - 1.)*hx - n**2*q*abs(q)/h**(10./3)    
    for i in range(N):
        if x[i] <= x_s:
            Z[i] = quad(func1, x[i], x_s)
            Z[i] += quad(func2, x_s, L)
        else:           
            Z[i] = quad(func2, x[i], L)
    #print 'Integral = ', -Z[:,0], ' with error = ', Z[:,1]
    return -1.*Z[:,0]

def analytic_sol(x):
    return water_height(x), bed(x)

##N = 1001
##X = linspace(0.,L,N)
##H = water_height(X)
##Hx= diff_water_height(X)
##Z = bed(X)
##
##plot(X,H+Z, X,Z)
##show()
