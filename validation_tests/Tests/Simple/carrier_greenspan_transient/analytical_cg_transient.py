"""
Transient water flows using ANUGA,
where water driven up a linear sloping beach and time varying boundary.
Ref: Carrier and Greenspan, Journal of Fluid Mechanics, 1958
"""
from scipy import linspace, zeros, array, dot
from math import log, ceil, sqrt
from scipy.optimize import fsolve

eps = 0.2
a_const = 1.5*(1.0+0.9*eps)**0.5

def shoreline(t):
    def f(u):
        numer = 5.0*(2.0/a_const*(u+t))**3.0 - (2.0/a_const*(u+t))**5.0
        denom = (1.0 + 4.0/a_const**2.0 * (u+t)**2.0)**4.0
        return u - 8.0*eps/a_const*numer/denom    
    u = fsolve(f,0.0)  # t must be involved here.
    lam = 2.0/a_const*(u+t)
    x = -0.5*u**2.0 + eps - eps*(1.0 + 3*lam**2.0 - 2.0*lam**4.0)/(1.0+lam**2.0)**3.0
    return x, u

def bed(x):
    return x

def analytical_sol(X,t)    :
    p, v = shoreline(t)
    N_X = len(X)
    Vel = zeros(N_X)
    Stage = zeros(N_X)
    q = zeros(2)
    #q[0] = p
    #q[1] = v
    z = bed(X)

    def fun(q): #Here q=(eta, u)
        fun = zeros(2)
        fun[0] = q[0] + 0.5*q[1]**2.0 - eps*(1.0 - 2.0*(1.25-(2.0/a_const*(q[1]+t))*1.0j) / ((1.0-(2.0/a_const*(q[1]+t))*1.0j)**2.0 + (16.0*(q[0]-x)/a_const**2.0))**1.5 + 1.5*(1.0-(2.0/a_const*(q[1]+t))*1.0j)**2.0/((1.0-(2.0/a_const*(q[1]+t))*1.0j)**2.0+(16.0*(q[0]-x)/a_const**2.0))**2.5).real
        fun[1] = q[1] - 8.0*eps/a_const*(1.0/((1.0-(2.0/a_const*(q[1]+t))*1.0j)**2.0 + (16.0*(q[0]-x)/a_const**2.0))**1.5 - 0.75*(1.0-(2.0/a_const*(q[1]+t))*1.0j)/((1.0-(2.0/a_const*(q[1]+t))*1.0j)**2.0+(16.0*(q[0]-x)/a_const**2.0))**2.5).imag
        return fun

    for m in range(N_X):
        k = N_X - m-1
        x = X[k]
        if X[k] > p:
            Vel[k] = 0.0
            Stage[k] = z[k]
        else:
            q = fsolve(fun,q)
            Stage[k] = q[0]
            Vel[k] = q[1]
    return Stage, z, Vel

"""
T = linspace(0.0, 40.0, 2)  # Given time
X = linspace(-0.5, 0.3, 201)
for i in range(len(T)):
    t=T[i]
    w, z, u = analytical_sol(X,t)
    from pylab import clf,plot,title,xlabel,ylabel,legend,savefig,show,hold,subplot
    hold(False)
    clf()
    plot1 = subplot(211)
    plot(X,w, X,z)
    xlabel('Position')
    ylabel('Stage')
    #plot1.set_ylim([0.0,0.25])
    legend(('Stage', ' Bed'),
           'lower right', shadow=False)        
    plot2 = subplot(212)
    plot(X,u)
    xlabel('Position')
    ylabel('Velocity')
    legend(('Velocity', ' '),
           'upper right', shadow=False)
    filename = "cg_"
    filename += str(i)
    filename += ".png"
    savefig(filename)
    show()
"""
