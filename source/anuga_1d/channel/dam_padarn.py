import os
from math import sqrt, pi
from channel_domain_Ab import *
from Numeric import allclose, array, zeros, ones, Float, take, sqrt
from config import g, epsilon


h1 = 5.0
h0 = 0.0

def analytical_sol(C,t):
    if t==0:
        t=0.0001
    #t  = 0.0     # time (s)
    # gravity (m/s^2)
    #h1 = 10.0    # depth upstream (m)
    #h0 = 0.0     # depth downstream (m)
    L = 2000.0   # length of stream/domain (m)
    n = len(C)    # number of cells

    u = zeros(n,Float)
    h = zeros(n,Float)
    x = C-3*L/4.0
    

    for i in range(n):
        # Calculate Analytical Solution at time t > 0
        u3 = 2.0/3.0*(sqrt(g*h1)+x[i]/t) 
        h3 = 4.0/(9.0*g)*(sqrt(g*h1)-x[i]/(2.0*t))*(sqrt(g*h1)-x[i]/(2.0*t))
        u3_ = 2.0/3.0*((x[i]+L/2.0)/t-sqrt(g*h1))
        h3_ = 1.0/(9.0*g)*((x[i]+L/2.0)/t+2*sqrt(g*h1))*((x[i]+L/2.0)/t+2*sqrt(g*h1))

        if ( x[i] <= -1*L/2.0+2*(-sqrt(g*h1)*t)):
            u[i] = 0.0
            h[i] = h0
        elif ( x[i] <= -1*L/2.0-(-sqrt(g*h1)*t)):
            u[i] = u3_
            h[i] = h3_

        elif ( x[i] <= -t*sqrt(g*h1) ):
            u[i] = 0.0 
            h[i] = h1 
        elif ( x[i] <= 2.0*t*sqrt(g*h1) ):
            u[i] = u3 
            h[i] = h3 
        else:
            u[i] = 0.0 
            h[i] = h0 

    return h , u*h, u

#def newLinePlot(title='Simple Plot'):
#   import Gnuplot
#    gg = Gnuplot.Gnuplot(persist=0)
#    gg.terminal(postscript)
#    gg.title(title)
#    gg('set data style linespoints') 
#    gg.xlabel('x')
#    gg.ylabel('y')
#    return gg

#def linePlot(gg,x1,y1,x2,y2):
#    import Gnuplot
#    plot1 = Gnuplot.PlotItems.Data(x1.flat,y1.flat,with="linespoints")
#    plot2 = Gnuplot.PlotItems.Data(x2.flat,y2.flat, with="lines 3")
#    g.plot(plot1,plot2)

h2=5.0
k=1

print "TEST 1D-SOLUTION III -- DRY BED"

def stage(x):
    y = zeros(len(x),Float)
    for i in range(len(x)):
        if x[i]<=L/4.0:
            y[i] = h0*width([x[i]])
        elif x[i]<=3*L/4.0:
            y[i] = h2*width([x[i]])
        else:
            y[i] = h0*width([x[i]])
    return y

def width(x):
    return k

 
import time

finaltime = 10.0
yieldstep = finaltime
L = 2000.0     # Length of channel (m)
number_of_cells = [810]
k = 0
widths = [1,2,5]
heights= []
velocities = []

for i in range(len(widths)):
    k=widths[i]
    for i in range(len(number_of_cells)):
        N = int(number_of_cells[i])
        print "Evaluating domain with %d cells" %N
        cell_len = L/N # Origin = 0.0
        points = zeros(N+1,Float)
        for j in range(N+1):
            points[j] = j*cell_len
        
        domain = Domain(points)
    
        domain.set_quantity('area', stage)
        domain.set_quantity('width',width)
        print "width in cell 1",domain.quantities['width'].vertex_values[1]
        domain.set_boundary({'exterior': Reflective_boundary(domain)})
        domain.order = 2
        domain.set_timestepping_method('rk2')
        domain.set_CFL(1.0)
        domain.set_limiter("vanleer")
        #domain.h0=0.0001
 
        t0 = time.time()

        for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):            
            domain.write_time()

        N = float(N)
        HeightC = domain.quantities['height'].centroid_values 
        DischargeC = domain.quantities['discharge'].centroid_values
        C = domain.centroids
        h, uh, u = analytical_sol(C,domain.time)
        h_error = 1.0/(N)*sum(abs(h-HeightC))
        u_error = 1.0/(N)*sum(abs(uh-DischargeC))
        #print "h_error %.10f" %(h_error[k])
        #print "uh_error %.10f"% (uh_error[k])
        k = k+1
        print 'That took %.2f seconds' %(time.time()-t0)
        X = domain.vertices
        heights.append(domain.quantities['height'].vertex_values)
        velocities.append( domain.quantities['velocity'].vertex_values)
        #stage = domain.quantities['stage'].vertex_values
        h, uh, u = analytical_sol(X.flat,domain.time)
        x = X.flat
            
        print "Error in height", h_error
        print "Error in xmom", u_error
        #from pylab import plot,title,xlabel,ylabel,legend,savefig,show,hold,subplot
from pylab import *
import pylab as p
import matplotlib.axes3d as p3
print 'test 2'
#hold(False)
print 'test 3'
plot1 = subplot(211)
print 'test 4'
plot(x,heights[0].flat,x,h)
print 'test 5'
plot1.set_ylim([-1,11])
xlabel('Position')
ylabel('Stage')
legend(('Numerical Solution','Analytic Soltuion'), 'upper right', shadow=True)
plot2 = subplot(212)


plot(x,velocities[0].flat,x,u)
plot2.set_ylim([-35,35])
xlabel('Position')
ylabel('Velocity')

print heights[0].flat-heights[1].flat
   
file = "dry_bed_"
file += str(number_of_cells[i])
file += ".eps"
#savefig(file)
show()
    

