"""
    Quick plot of the dam break outputs

"""
import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from analytical_avalanche_dry import *

p_st = util.get_output('avalanche.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[10]
v2=(p2_st.y==v)

x_n = p2_st.x[v2]



u0,h0,w0,z0,p0 = analytical_sol(x_n, p2_st.time[0])
u10,h10,w10,z10,p10 = analytical_sol(x_n, p2_st.time[10])
u30,h30,w30,z30,p30 = analytical_sol(x_n, p2_st.time[30])


w0_n  = p2_st.stage[0,v2]
w10_n = p2_st.stage[10,v2]
w30_n = p2_st.stage[30,v2]

z_n = p2_st.elev[v2]

uh0_n  = p2_st.xmom[0,v2]
uh10_n = p2_st.xmom[10,v2]
uh30_n = p2_st.xmom[30,v2]

u0_n  = p2_st.xvel[0,v2]
u10_n = p2_st.xvel[10,v2]
u30_n = p2_st.xvel[30,v2]



#Plot stages
pyplot.clf()
pyplot.plot(x_n, w0_n,'b.-', label='numerical stage')
pyplot.plot(x_n, w10_n, 'b.-')
pyplot.plot(x_n, w30_n, 'b.-')
pyplot.plot(x_n, w0,'r-', label='analytical stage')
pyplot.plot(x_n, w10,'r-')
pyplot.plot(x_n, w30,'r-')
pyplot.plot(x_n, z_n,'k-', label='bed elevation')
pyplot.title('Stage at several instants in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png', dpi=600)
#pyplot.show()


#Plot xmomentum
pyplot.clf()
pyplot.plot(x_n, uh0_n, 'b.-', label='numerical')
pyplot.plot(x_n, uh10_n, 'b.-')
pyplot.plot(x_n, uh30_n,'b.-')
pyplot.plot(x_n, u0*h0,'r-', label='analytical')
pyplot.plot(x_n, u10*h10,'r-')
pyplot.plot(x_n, u30*h30,'r-')
pyplot.title('Xmomentum at several instants in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('xmom_plot.png')


#Plot velocities
pyplot.clf()
pyplot.plot(x_n, u0_n, 'b.-', label='numerical')
pyplot.plot(x_n, u10_n, 'b.-')
pyplot.plot(x_n, u30_n,'b.-')
pyplot.plot(x_n, u0,'r-', label='analytical')
pyplot.plot(x_n, u10,'r-')
pyplot.plot(x_n, u30,'r-')
pyplot.title('Xvelocity at several instants in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_plot.png')
