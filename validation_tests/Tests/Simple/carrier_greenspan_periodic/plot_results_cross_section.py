"""
    Quick plot of the Carrier-Greenspan outputs

"""
import anuga.utilities.plot_utils as util
from matplotlib import pyplot as pyplot
from analytical_carrier_greenspan import *

p_st = util.get_output('carrier_greenspan.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[10]
v2=(p2_st.y==v)

W1, P1, Z1, H1, U1 = analytic_cg(p2_st.x[v2], p2_st.time[288], h0=5e2, L=5e4, a=1.0, Tp=900.0)
W2, P2, Z2, H2, U2 = analytic_cg(p2_st.x[v2], p2_st.time[296], h0=5e2, L=5e4, a=1.0, Tp=900.0)
W3, P3, Z3, H3, U3 = analytic_cg(p2_st.x[v2], p2_st.time[304], h0=5e2, L=5e4, a=1.0, Tp=900.0)
W4, P4, Z4, H4, U4 = analytic_cg(p2_st.x[v2], p2_st.time[312], h0=5e2, L=5e4, a=1.0, Tp=900.0)
W5, P5, Z5, H5, U5 = analytic_cg(p2_st.x[v2], p2_st.time[320], h0=5e2, L=5e4, a=1.0, Tp=900.0)
W6, P6, Z6, H6, U6 = analytic_cg(p2_st.x[v2], p2_st.time[328], h0=5e2, L=5e4, a=1.0, Tp=900.0)

#Plot the stages##############################################################
pyplot.clf()
pyplot.ion()
pyplot.plot(p2_st.x[v2], p2_st.stage[288,v2], 'b.', label='numerical') # 0*T/6
pyplot.plot(p2_st.x[v2], p2_st.stage[296,v2], 'b.')                    # 1*T/6
pyplot.plot(p2_st.x[v2], p2_st.stage[304,v2], 'b.')                     # 2*T/6
pyplot.plot(p2_st.x[v2], p2_st.stage[312,v2], 'b.')                     # 3*T/6
pyplot.plot(p2_st.x[v2], p2_st.stage[320,v2], 'b.')                     # 4*T/6
pyplot.plot(p2_st.x[v2], p2_st.stage[328,v2], 'b.')  # 5*T/6

pyplot.plot(p2_st.x[v2], W1,'r-', label='analytical')
pyplot.plot(p2_st.x[v2], W2,'r-')
pyplot.plot(p2_st.x[v2], W3,'r-')
pyplot.plot(p2_st.x[v2], W4,'r-')
pyplot.plot(p2_st.x[v2], W5,'r-')
pyplot.plot(p2_st.x[v2], W6,'r-')


pyplot.title('Stage at several instants in time')
pyplot.ylim(-5.0,5.0)
pyplot.legend(loc=2)
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png')


#Plot the momentums##########################################################
pyplot.clf()
pyplot.ion()
pyplot.plot(p2_st.x[v2], p2_st.xmom[288,v2], 'b.', label='numerical') # 0*T/6
pyplot.plot(p2_st.x[v2], p2_st.xmom[296,v2], 'b.')                    # 1*T/6
pyplot.plot(p2_st.x[v2], p2_st.xmom[304,v2], 'b.')                     # 2*T/6
pyplot.plot(p2_st.x[v2], p2_st.xmom[312,v2], 'b.')                     # 3*T/6
pyplot.plot(p2_st.x[v2], p2_st.xmom[320,v2], 'b.')                     # 4*T/6
pyplot.plot(p2_st.x[v2], p2_st.xmom[328,v2], 'b.')  # 5*T/6

pyplot.plot(p2_st.x[v2], P1,'r-', label='analytical')
pyplot.plot(p2_st.x[v2], P2,'r-')
pyplot.plot(p2_st.x[v2], P3,'r-')
pyplot.plot(p2_st.x[v2], P4,'r-')
pyplot.plot(p2_st.x[v2], P5,'r-')
pyplot.plot(p2_st.x[v2], P6,'r-')

pyplot.title('Xmomentum at several instants in time')
pyplot.legend(loc=4)
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('xmom_plot.png')



#Plot the velocities#########################################################
pyplot.clf()
pyplot.ion()
pyplot.plot(p2_st.x[v2], p2_st.xvel[288,v2], 'b.', label='numerical') # 0*T/6
pyplot.plot(p2_st.x[v2], p2_st.xvel[296,v2], 'b.')                    # 1*T/6
pyplot.plot(p2_st.x[v2], p2_st.xvel[304,v2], 'b.')                     # 2*T/6
pyplot.plot(p2_st.x[v2], p2_st.xvel[312,v2], 'b.')                     # 3*T/6
pyplot.plot(p2_st.x[v2], p2_st.xvel[320,v2], 'b.')                     # 4*T/6
pyplot.plot(p2_st.x[v2], p2_st.xvel[328,v2], 'b.')  # 5*T/6

pyplot.plot(p2_st.x[v2], U1,'r-', label='analytical')
pyplot.plot(p2_st.x[v2], U2,'r-')
pyplot.plot(p2_st.x[v2], U3,'r-')
pyplot.plot(p2_st.x[v2], U4,'r-')
pyplot.plot(p2_st.x[v2], U5,'r-')
pyplot.plot(p2_st.x[v2], U6,'r-')

pyplot.title('Xvelocity at several instants in time')
pyplot.legend(loc=3)
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_plot.png')

