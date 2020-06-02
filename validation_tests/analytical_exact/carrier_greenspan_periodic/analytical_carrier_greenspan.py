from scipy import sin, cos, sqrt, linspace, pi, dot
from numpy import zeros, array
from scipy.special import jn
from scipy.optimize import fsolve
from gaussPivot import *


#Parameters
g = 9.81

def j0(x):
    return jn(0.0, x)

def j1(x):
    return jn(1.0, x)

def analytic_cg(points, t=0.0, h0=5e2, L=5e4, a=1.0, Tp=900.0):
    #
    #Exact solution of the Carrier-Greenspan periodic solution.
    #Ref1: Carrier and Greenspan, J. Fluid Mech., 1958
    #Ref2: Mungkasi and Roberts, Int. J. Numer. Meth. Fluids, 2012
    #Here, points are the spatial evaluating points, t is time variable, h0 is the water depth at Origin
    #L is the horizontal length considered, a is the amplitude of perturbation at Origin
    #Tp is the period of oscillation.
    #
    #Analytical computations#################################################################
    def shore(t):
        def g(u):
            return u + 2.0*A*pi/T*sin(2.0*pi/T*(t+u))    
        u = fsolve(g,0.0)
        xi = -0.5*u*u + A*cos(2.0*pi/T*(t+u))
        position = 1.0 + xi
        return position, u

    def func(q): #Here q=(w, u)
        f = zeros(2)
        f[0] = q[0] + 0.5*q[1]**2.0 - A*j0(4.0*pi/T*(1.0+q[0]-x)**0.5)*cos(2.0*pi/T*(tim+q[1]))
        f[1] = q[1] + A*j1(4.0*pi/T*(1.0+q[0]-x)**0.5)*sin(2.0*pi/T*(tim+q[1]))/(1+q[0]-x)**0.5
        return f  

    def newtonRaphson2(f,q,tol=1.0e-15): ##1.0e-9 may be too large
        for i in range(50):                                                                                                                                    
            h = 1.0e-4      ##1.0e-4 may be too large.
            n = len(q)
            jac = zeros((n,n))
            if 1.0+q[0]-x<0.0:
                #print "Problem occurs i=",i
                #print "PROBLEM OCCURS.......... 1.0+q[0]-x=",1.0+q[0]-x
                q[0] = x-1.0 + 0.0001
                #print "NOW problem is fixed as  1.0+q[0]-x=",1.0+q[0]-x
            f0 = f(q)
            for i in range(n):
                temp = q[i]
                q[i] = temp + h
                f1 = f(q)
                q[i] = temp
                jac[:,i] = (f1 - f0)/h
            if sqrt(dot(f0,f0)/len(q)) < tol: return q
            dq = gaussPivot(jac,-f0)
            q = q + dq
            if sqrt(dot(dq,dq)) < tol*max(max(abs(q)),1.0): return q
        print('Too many iterations')    
    ##################################################################################
    N = len(points)
    
    eps = a/h0           # Dimensionless amplitude
    T = Tp*sqrt(g*h0)/L  # Dimensionless period
    tim = t*sqrt(g*h0)/L # Dimensionless time
    A = eps/j0(4.0*pi/T) # Dimensionless max horizontal displacement of shore
    
    W = zeros(N) # Dimensional stage
    P = zeros(N) # Dimensional momentum
    Z = zeros(N) # Dimensional bed
    H = zeros(N) # Dimensional height
    U = zeros(N) # Dimensional velocity
    
    q = zeros(2)
    pos_shore, vel_shore = shore(tim) # Dimensionless shore
    #q[0] = pos_shore
    #q[1] = vel_shore

    for m in range(N):
        i = N - m-1
        X = points[i]            # Dimensional
        if X >= pos_shore*L:     # Dimensional
            W[i] = (h0/L)*X - h0 # Dimensional
            U[i] = 0.0           # Dimensional           
            Z[i] = (h0/L)*X - h0 # Dimensional
            H[i] = 0.0           # Dimensional
            P[i] = 0.0           # Dimensional
        else:
            x = X/L              # Dimensionless
            #q = fsolve(func,q)  # Dimensionless (fsolve is too sensitive!)
            q = newtonRaphson2(func,q)  # Dimensionless (better than fsolve)           
            W[i] = q[0]*h0      # Dimensional
            U[i] = q[1]         # It needs dimensionalisation
            Z[i] = (h0/L)*X - h0# Dimensional
            H[i] = W[i] - Z[i]  # Dimensional
            P[i] = H[i]*U[i]    # It needs dimensionalisation
    U = U*sqrt(g*h0) #This is the dimensionalisation
    P = P*sqrt(g*h0) #This is the dimensionalisation
    return W, P, Z, H, U


if __name__ == "__main__":
    points = linspace(0.0, 55000.0, 1000)
    W, P, Z, H, U = analytic_cg(points,t=300.0, h0=5e2, L=5e4, a=1.0, Tp=900.0)

    from pylab import clf,plot,title,xlabel,ylabel,legend,savefig,show,hold,subplot,ion
    hold(False)
    clf()
    plot1 = subplot(311)
    plot(points/1e+4,W,  points/1e+4,Z)
    ylabel('Stage')

    plot2 = subplot(312)
    plot(points/1e+4,P,'b-')
    ylabel('Momentum')      

    plot3 = subplot(313)
    plot(points/1e+4,U,'b-')
    xlabel('Position / 10,000')
    ylabel('Velocity')
    legend(('Analytical Solution','Numerical Solution'),
           'upper center', shadow=False)

    show()
