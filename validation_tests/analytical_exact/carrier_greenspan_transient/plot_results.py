"""
    Quick plot of the Carrier-Greenspan outputs

"""
import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from analytical_cg_transient import *

p_st = util.get_output('carrier_greenspan.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[10]
v2=(p2_st.y==v)

#Dimensional parameters
L   = 5e4         # Length of channel (m)      
h0  = 5e2         # Height at origin when the water is still
g   = 9.81        # Gravity

#Dimensionless solution
W1, Z1, U1 = analytical_sol(p2_st.x[v2]/L, p2_st.time[0]*sqrt(g*h0)/L)
W2, Z2, U2 = analytical_sol(p2_st.x[v2]/L, p2_st.time[1]*sqrt(g*h0)/L)
W3, Z3, U3 = analytical_sol(p2_st.x[v2]/L, p2_st.time[30]*sqrt(g*h0)/L)

#Dimensional solution
W1 = W1*h0                # dimensional
U1 = U1*sqrt(g*h0)        # dimensional
Z1 = Z1*h0                # dimensional
H1 = W1 - Z1              # dimensional
P1 = U1 * H1              # dimensional

W2 = W2*h0                # dimensional
U2 = U2*sqrt(g*h0)        # dimensional
Z2 = Z2*h0                # dimensional
H2 = W2 - Z2              # dimensional
P2 = U2 * H2              # dimensional
      

W3 = W3*h0                # dimensional
U3 = U3*sqrt(g*h0)        # dimensional
Z3 = Z3*h0                # dimensional
H3 = W3 - Z3              # dimensional
P3 = U3 * H3              # dimensional



#Plot the stages##############################################################
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.stage[0,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], p2_st.stage[1,v2], 'b.') 
pyplot.plot(p2_st.x[v2], p2_st.stage[30,v2], 'b.')                    
pyplot.plot(p2_st.x[v2], W1,'r-', label='analytical')
pyplot.plot(p2_st.x[v2], W2,'r-')
pyplot.plot(p2_st.x[v2], W3,'r-')
pyplot.title('Stage at several instants in time')
#pyplot.ylim(-5.0,5.0)
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png')


#Plot the momentums##########################################################
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xmom[0,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], p2_st.xmom[1,v2], 'b.') 
pyplot.plot(p2_st.x[v2], p2_st.xmom[30,v2], 'b.')                    
pyplot.plot(p2_st.x[v2], P1,'r-', label='analytical')
pyplot.plot(p2_st.x[v2], P2,'r-')
pyplot.plot(p2_st.x[v2], P3,'r-')
pyplot.title('Xmomentum at several instants in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('xmom_plot.png')


#Plot the velocities#########################################################
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xvel[0,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], p2_st.xvel[1,v2], 'b.')
pyplot.plot(p2_st.x[v2], p2_st.xvel[30,v2], 'b.')                    
pyplot.plot(p2_st.x[v2], U1,'r-', label='analytical')
pyplot.plot(p2_st.x[v2], U2,'r-')
pyplot.plot(p2_st.x[v2], U3,'r-')
pyplot.title('Xvelocity at several instants in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_plot.png')

