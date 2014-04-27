import os
from math import sqrt, pi
from channel_domain import *
from Numeric import allclose, array, zeros, ones, Float, take, sqrt
from config import g, epsilon


h1 = 5.0
h0 = 0.0

## def analytical_sol(C,t):
    
##     #t  = 0.0     # time (s)
##     # gravity (m/s^2)
##     #h1 = 10.0    # depth upstream (m)
##     #h0 = 0.0     # depth downstream (m)
##     L = 2000.0   # length of stream/domain (m)
##     n = len(C)    # number of cells

##     u = zeros(n,Float)
##     h = zeros(n,Float)
##     x = C-3*L/4.0
    

##     for i in range(n):
##         # Calculate Analytical Solution at time t > 0
##         u3 = 2.0/3.0*(sqrt(g*h1)+x[i]/t) 
##         h3 = 4.0/(9.0*g)*(sqrt(g*h1)-x[i]/(2.0*t))*(sqrt(g*h1)-x[i]/(2.0*t))
##         u3_ = 2.0/3.0*((x[i]+L/2.0)/t-sqrt(g*h1))
##         h3_ = 1.0/(9.0*g)*((x[i]+L/2.0)/t+2*sqrt(g*h1))*((x[i]+L/2.0)/t+2*sqrt(g*h1))

##         if ( x[i] <= -1*L/2.0+2*(-sqrt(g*h1)*t)):
##             u[i] = 0.0
##             h[i] = h0
##         elif ( x[i] <= -1*L/2.0-(-sqrt(g*h1)*t)):
##             u[i] = u3_
##             h[i] = h3_

##         elif ( x[i] <= -t*sqrt(g*h1) ):
##             u[i] = 0.0 
##             h[i] = h1 
##         elif ( x[i] <= 2.0*t*sqrt(g*h1) ):
##             u[i] = u3 
##             h[i] = h3 
##         else:
##             u[i] = 0.0 
##             h[i] = h0 

##     return h , u*h, u

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

h2=10.0

print "TEST 1D-SOLUTION III -- DRY BED"

def stage(x):
    return 5

def width(x):
    return 1

def bot(x):
    return x/500
 
import time

finaltime = 20.0
yieldstep = finaltime
L = 2000.0     # Length of channel (m)
number_of_cells = [810]#,200,500,1000,2000,5000,10000,20000]
h_error = zeros(len(number_of_cells),Float)
uh_error = zeros(len(number_of_cells),Float)
k = 0
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
    domain.set_quantity('elevation',bot)
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
    #h, uh, u = analytical_sol(C,domain.time)
    #h_error[k] = 1.0/(N)*sum(abs(h-StageC))
    #uh_error[k] = 1.0/(N)*sum(abs(uh-XmomC))
    #print "h_error %.10f" %(h_error[k])
    #print "uh_error %.10f"% (uh_error[k])
    k = k+1
    print 'That took %.2f seconds' %(time.time()-t0)
    X = domain.vertices
    HeightQ = domain.quantities['height'].vertex_values
    VelocityQ = domain.quantities['velocity'].vertex_values
    #stage = domain.quantities['stage'].vertex_values
    #h, uh, u = analytical_sol(X.flat,domain.time)
    x = X.flat
    z =bot(x)
    w=HeightQ.flat+z
    
    #from pylab import plot,title,xlabel,ylabel,legend,savefig,show,hold,subplot
    from pylab import *
    import pylab as p
    import matplotlib.axes3d as p3
    print 'test 2'
    hold(False)
    print 'test 3'
    plot1 = subplot(211)
    print 'test 4'
    plot(x,z,x,w)
    print 'test 5'
    plot1.set_ylim([-1,11])
    xlabel('Position')
    ylabel('Stage')
    #legend(( 'Numerical Solution'),
    #       'upper right', shadow=True)
    plot2 = subplot(212)
    plot(x,VelocityQ.flat)
    plot2.set_ylim([-2,2])
    
    xlabel('Position')
    ylabel('Velocity')
   
    file = "dry_bed_"
    file += str(number_of_cells[i])
    file += ".eps"
    #savefig(file)
    show()
    
    
#print "Error in height", h_error
#print "Error in xmom", uh_error
