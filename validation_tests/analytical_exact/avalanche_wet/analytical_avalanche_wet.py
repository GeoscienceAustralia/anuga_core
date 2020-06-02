"""
Analytical solution to debris avalanche involving a shock wave
Ref: Mungkasi and Roberts, Pure and Applied Geophysics 2012
"""

from numpy import array, zeros, linspace
from numpy import arctan, sqrt, sin, cos, tan

#Parameters
g = 9.81           # gravity
h_0 = 20.0         # depth upstream
h_1 = 10.0         # depth downstream
L = 200.0          # length of stream/domain


#delt = 0.0
friction_slope = 0.05             #tan(delta) # NOTE THAT friction_slope must less than bed_slope
#thet = pi/36.0
bed_slope = 0.1                   #tan(theta) #0.1      # bottom slope, positive if it is increasing bottom.
thet = arctan(bed_slope)
F2 = g*cos(thet)*cos(thet)*friction_slope
c0 = sqrt(g*h_0)                  # sound speed
m = -1.0*g*bed_slope + F2         # auxiliary variable


# Function for analytical solution
def analytical_sol(X,t):
    #The following is for computing the shock velocity
    c0 = sqrt(g*h_0) #right celerity
    c1 = sqrt(g*h_1) #left celerity
    Smin=-21.0
    Smax=20.0
    for i in range(100):
        S=(Smin+Smax)/2.0       
        u2=S-c1*c1/4.0/S*(1.0+sqrt(1.0+8.0*S*S/c1/c1))
        c2=c1*sqrt(0.5*(sqrt(1.0+8.0*S*S/c1/c1)-1.0))
        func = -2.0*c0 -u2 +2.0*c2
        if (func > 0.0):
            Smin=S
        else:
            Smax=S
    if( abs(S) > 99.0):
        print('no convergence')
    #h2= h_1/(1.0-u2/S)
    h2 = c2**2/9.81
    #print "S=",S
    #print "c2=",c2
    #print "u2=",u2
    
    N = len(X)
    u = zeros(N)
    h = zeros(N)
    w = zeros(N)
    z = zeros(N)
    mom = zeros(N)
    if abs(t)>1e-9:
        for i in range(N):
            # Calculate Analytical Solution at time t > 0
            if X[i]-0.5*m*t**2 <= S*t :
                u[i] = m*t 
                h[i] = h_1
            elif X[i]-0.5*m*t**2 <= (u2+c2)*t :
                u[i] = u2+m*t
                h[i] = h2
            elif X[i]-0.5*m*t**2 < c0*t :
                u[i] = 2.0/3.0*((X[i]-0.5*m*t**2)/t - sqrt(g*h_0))  +m*t
                h[i] = 4.0/(9.0*g)*((X[i]-0.5*m*t**2)/(2.0*t) + sqrt(g*h_0))**2
            else:
                u[i] = m*t
                h[i] = h_0 
            z[i] = bed_slope*X[i]
            w[i] = h[i] + z[i]
            mom[i] = u[i]*h[i]
    else:
        for i in range(N):
            # Calculate Analytical Solution at time t > 0
            if X[i] <= 0.0:
                u[i] = 0.0
                h[i] = h_1
            else:
                u[i] = 0.0
                h[i] = h_0
            z[i] = bed_slope*X[i]
            w[i] = h[i] + z[i]
            mom[i] = u[i]*h[i]    
    return u,h,w,z,mom

"""
# Run a simulation
n = 800            # number of cells
cell_len = L/n     # length of each cell
points = zeros(n+1)
for i in range (n+1):
    points[i] = i*cell_len - 0.5*L
X = points
t = 4.0
#The following is for plotting the result.
Uv,Hv,Wv,Zv,Mv = analytical_sol(X,t)
from pylab import clf,plot,title,xlabel,ylabel,legend,savefig,show,hold,subplot
hold(False)
clf()

plot1 = subplot(211)
plot(X,Wv,'b-', X,Zv,'k-')
ylabel('Stage')
#plot1.set_ylim([-50.0,75.0])#([-1.0,21.0])
legend(('analytical solution', 'bed elevation'), 'upper left', shadow=False)

plot2 = subplot(212)
plot(X,Mv,'b-')
xlabel('Position')
ylabel('Momentum')
#plot2.set_ylim([-200.0,10.0])#([-90.0,10.0])   
show()
"""
