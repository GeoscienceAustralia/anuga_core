"""
    Quick plot of the Carrier-Greenspan outputs

"""
import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from analytical_carrier_greenspan import *

p_st = util.get_output('carrier_greenspan.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[10]
v2=(p2_st.y==v)

x_n = p2_st.x[v2]


n6 = -1
n5 = n6-12
n4 = n5-12
n3 = n4-12
n2 = n3-12
n1 = n2-12

W1, P1, Z1, H1, U1 = analytic_cg(x_n, p2_st.time[n1], h0=5e2, L=5e4, a=1.0, Tp=900.0)
W2, P2, Z2, H2, U2 = analytic_cg(x_n, p2_st.time[n2], h0=5e2, L=5e4, a=1.0, Tp=900.0)
W3, P3, Z3, H3, U3 = analytic_cg(x_n, p2_st.time[n3], h0=5e2, L=5e4, a=1.0, Tp=900.0)
W4, P4, Z4, H4, U4 = analytic_cg(x_n, p2_st.time[n4], h0=5e2, L=5e4, a=1.0, Tp=900.0)
W5, P5, Z5, H5, U5 = analytic_cg(x_n, p2_st.time[n5], h0=5e2, L=5e4, a=1.0, Tp=900.0)
W6, P6, Z6, H6, U6 = analytic_cg(x_n, p2_st.time[n6], h0=5e2, L=5e4, a=1.0, Tp=900.0)

W1_n = p2_st.stage[n1,v2]
W2_n = p2_st.stage[n2,v2]
W3_n = p2_st.stage[n3,v2]
W4_n = p2_st.stage[n4,v2]
W5_n = p2_st.stage[n5,v2]
W6_n = p2_st.stage[n6,v2]

UH1_n = p2_st.xmom[n1,v2]
UH2_n = p2_st.xmom[n2,v2]
UH3_n = p2_st.xmom[n3,v2]
UH4_n = p2_st.xmom[n4,v2]
UH5_n = p2_st.xmom[n5,v2]
UH6_n = p2_st.xmom[n6,v2]

U1_n = p2_st.xvel[n1,v2]
U2_n = p2_st.xvel[n2,v2]
U3_n = p2_st.xvel[n3,v2]
U4_n = p2_st.xvel[n4,v2]
U5_n = p2_st.xvel[n5,v2]
U6_n = p2_st.xvel[n6,v2]

#Plot the stages##############################################################
pyplot.clf()
pyplot.plot(x_n, W1_n, 'b.', label='numerical')  # 0*T/6
pyplot.plot(x_n, W2_n, 'b.')                     # 1*T/6
pyplot.plot(x_n, W3_n, 'b.')                     # 2*T/6
pyplot.plot(x_n, W4_n, 'b.')                     # 3*T/6
pyplot.plot(x_n, W5_n, 'b.')                     # 4*T/6
pyplot.plot(x_n, W6_n, 'b.')                     # 5*T/6

pyplot.plot(x_n, W1,'r-', label='analytical')
pyplot.plot(x_n, W2,'r-')
pyplot.plot(x_n, W3,'r-')
pyplot.plot(x_n, W4,'r-')
pyplot.plot(x_n, W5,'r-')
pyplot.plot(x_n, W6,'r-')


pyplot.title('Stage at several instants in time')
pyplot.ylim(-5.0,5.0)
pyplot.legend(loc=2)
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png')


#Plot the momentums##########################################################
pyplot.clf()
pyplot.plot(x_n, UH1_n, 'b.', label='numerical')  # 0*T/6
pyplot.plot(x_n, UH2_n, 'b.')                     # 1*T/6
pyplot.plot(x_n, UH3_n, 'b.')                     # 2*T/6
pyplot.plot(x_n, UH4_n, 'b.')                     # 3*T/6
pyplot.plot(x_n, UH5_n, 'b.')                     # 4*T/6
pyplot.plot(x_n, UH6_n, 'b.')                     # 5*T/6

pyplot.plot(x_n, P1,'r-', label='analytical')
pyplot.plot(x_n, P2,'r-')
pyplot.plot(x_n, P3,'r-')
pyplot.plot(x_n, P4,'r-')
pyplot.plot(x_n, P5,'r-')
pyplot.plot(x_n, P6,'r-')

pyplot.title('Xmomentum at several instants in time')
pyplot.legend(loc=4)
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('xmom_plot.png')



#Plot the velocities#########################################################
pyplot.clf()
pyplot.plot(x_n, U1_n, 'b.', label='numerical')  # 0*T/6
pyplot.plot(x_n, U2_n, 'b.')                     # 1*T/6
pyplot.plot(x_n, U3_n, 'b.')                     # 2*T/6
pyplot.plot(x_n, U4_n, 'b.')                     # 3*T/6
pyplot.plot(x_n, U5_n, 'b.')                     # 4*T/6
pyplot.plot(x_n, U6_n, 'b.')                     # 5*T/6

pyplot.plot(x_n, U1,'r-', label='analytical')
pyplot.plot(x_n, U2,'r-')
pyplot.plot(x_n, U3,'r-')
pyplot.plot(x_n, U4,'r-')
pyplot.plot(x_n, U5,'r-')
pyplot.plot(x_n, U6,'r-')

pyplot.title('Xvelocity at several instants in time')
pyplot.legend(loc=3)
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_plot.png')



maxx=p2_st.x.max()

minx=p2_st.x.min()

v1 = abs(p2_st.x -minx).argmin()
#v2 = abs(p.x -0.5*(minx+maxx)).argmin()
#v3 = abs(p.x -maxx).argmin()

pyplot.clf()
pyplot.plot(p2_st.time, p2_st.stage[:,v1], label='Left edge of domain')
#pyplot.plot(p.time, p.stage[:,v2], label='Middle of domain')
#pyplot.plot(p.time, p.stage[:,v3], label='Right edge of domain')
pyplot.ylim((-1.0,1.0))
pyplot.legend()
pyplot.xlabel('Time')
pyplot.ylabel('Perturbation')
pyplot.title('Vertical perturbation at the origin over time')
pyplot.savefig('perturbation_at_origin.png')
