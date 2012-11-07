"""
This module contains various auxiliary function used by pyvolution.
"""

def mean(x):
    from numpy import sum
    return sum(x)/len(x)


def gradient(x0, x1, q0, q1):

    if q1-q0 != 0:
        a = (q1-q0)/(x1-x0)
    else:
        a = 0
        
    return a

def  minmod(beta_p,beta_m):
    if (abs(beta_p) < abs(beta_m)) & (beta_p*beta_m > 0.0):
        phi = beta_p
    elif (abs(beta_m) < abs(beta_p)) & (beta_p*beta_m > 0.0):
        phi = beta_m
    else:
        phi = 0.0
    return phi

def  minmod_kurganov(a,b,c):
    from numpy import sign
    if sign(a)==sign(b)==sign(c):
        return sign(a)*min(abs(a),abs(b),abs(c))
    else:
        return 0.0

def  maxmod(a,b):
    if (abs(a) > abs(b)) & (a*b > 0.0):
        phi = a
    elif (abs(b) > abs(a)) & (a*b > 0.0):
        phi = b
    else:
        phi = 0.0
    return phi

def vanleer(a,b):
    if abs(a)+abs(b) > 1e-12:
        return (a*abs(b)+abs(a)*b)/(abs(a)+abs(b))
    else:
        return 0.0

def vanalbada(a,b):
    if a*a+b*b > 1e-12:
        return (a*a*b+a*b*b)/(a*a+b*b)
    else:
        return 0.0

def calculate_wetted_area(x1,x2,z1,z2,w1,w2):
    if (w1 > z1) & (w2 < z2) & (z1 <= z2):
        x = ((w2-z1)*(x2-x1)+x1*(z2-z1)-x2*(w2-w1))/(z2-z1+w1-w2)
        A = 0.5*(w1-z1)*(x-x1)
        L = x-x1
    elif (w1 < z1) & (w2 > z2) & (z1 < z2):
        x = ((w2-z1)*(x2-x1)+x1*(z2-z1)-x2*(w2-w1))/(z2-z1+w1-w2)
        A = 0.5*(w2-z2)*(x2-x)
        L = x2-x
    elif (w1 < z1) & (w2 > z2) & (z1 >= z2):
        x = ((w1-z2)*(x2-x1)+x2*(z2-z1)-x1*(w2-w1))/(z2-z1+w1-w2)
        A = 0.5*(w2-z2)*(x2-x)
        L = x2-x
    elif (w1 > z1) & (w2 < z2) & (z1 > z2):
        x = ((w1-z2)*(x2-x1)+x2*(z2-z1)-x1*(w2-w1))/(z2-z1+w1-w2)
        A = 0.5*(w1-z1)*(x-x1)
        L = x-x1
    elif (w1 <= z1) & (w2 <= z2):
        A = 0.0
    elif (w1 == z1) & (w2 > z2) & (z2 < z1):
        A = 0.5*(x2-x1)*(w2-z2)
    elif (w2 == z2) & (w1 > z1) & (z1 < z2):
        A = 0.5*(x2-x1)*(w1-z1)
    return A


def calculate_new_wet_area(x1,x2,z1,z2,A):
    from numpy import sqrt
    min_centroid_height = 1.0e-3
    # Assumes reconstructed stage flat in a wetted cell
    M = (z2-z1)/(x2-x1)
    L = (x2-x1)
    min_area = min_centroid_height*L
    max_area = 0.5*(x2-x1)*abs(z2-z1)
    if A < max_area:
        if (z1 < z2):
            x = sqrt(2*A/M)+x1
            wet_len = x-x1
            wc = z1 + sqrt(M*2*A)
        elif (z2 < z1):
            x = -sqrt(-2*A/M)+x2
            wet_len = x2-x 
            wc = z2+sqrt(-M*2*A)
        else:
            wc = A/L+0.5*(z1+z2)
            wet_len = x2-x1
    else:
        wc = 0.5*(z1+z2)+A/L
        wet_len = x2-x1
            
    return wc,wet_len

def calculate_new_wet_area_analytic(x1,x2,z1,z2,A,t):
    min_centroid_height = 1.0e-3
    # Assumes reconstructed stage flat in a wetted cell
    M = (z2-z1)/(x2-x1)
    L = (x2-x1)
    min_area = min_centroid_height*L
    max_area = 0.5*(x2-x1)*abs(z2-z1)
    w1,uh1 = analytic_cannal(x1,t)
    w2,uh2 = analytic_cannal(x2,t)
    if (w1 > z1) & (w2 < z2) & (z1 <= z2):
        print "test1"
        x = ((w2-z1)*(x2-x1)+x1*(z2-z1)-x2*(w2-w1))/(z2-z1+w1-w2)
        wet_len = x-x1
    elif (w1 < z1) & (w2 > z2) & (z1 < z2):
        print "test2"
        x = ((w2-z1)*(x2-x1)+x1*(z2-z1)-x2*(w2-w1))/(z2-z1+w1-w2)
        wet_len = x2-x
    elif (w1 < z1) & (w2 > z2) & (z1 >= z2):
        print "test3"
        x = ((w1-z2)*(x2-x1)+x2*(z2-z1)-x1*(w2-w1))/(z2-z1+w1-w2)
        wet_len = x2-x
    elif (w1 > z1) & (w2 < z2) & (z1 > z2):
        print "test4"
        x = ((w1-z2)*(x2-x1)+x2*(z2-z1)-x1*(w2-w1))/(z2-z1+w1-w2)
        wet_len = x-x1
    elif (w1 >= z1) & (w2 >= z2):
        print "test5"
        wet_len = x2-x1 
    else: #(w1 <= z1) & (w2 <= z2)
        print "test5"
        if (w1 > z1) | (w2 > z2):
            print "ERROR"
        wet_len = x2-x1        
    return w1,w2,wet_len,uh1,uh2

def analytic_cannal(C,t):
    from numpy import zeros,sqrt,sin,cos

    
    #N = len(C)
    #u = zeros(N)    ## water velocity
    #h = zeros(N)    ## water depth
    x = C
    g = 9.81


    ## Define Basin Bathymetry
    #z_b = zeros(N) ## elevation of basin
    #z = zeros(N)   ## elevation of water surface
    z_infty = 10.0       ## max equilibrium water depth at lowest point.
    L_x = 2500.0         ## width of channel

    A0 = 0.5*L_x                  ## determines amplitudes of oscillations
    omega = sqrt(2*g*z_infty)/L_x ## angular frequency of osccilation

    x1 = A0*cos(omega*t)-L_x # left shoreline
    x2 = A0*cos(omega*t)+L_x # right shoreline
    if (x >=x1) & (x <= x2):
        z_b = z_infty*(x**2/L_x**2) ## or A0*cos(omega*t)\pmL_x
        u = -A0*omega*sin(omega*t)
        z = z_infty+2*A0*z_infty/L_x*cos(omega*t)*(x/L_x-0.5*A0/(L_x)*cos(omega*t))
    else:
       z_b = z_infty*(x**2/L_x**2)
       u=0.0
       z = z_b
    h = z-z_b
    return z,u*h
