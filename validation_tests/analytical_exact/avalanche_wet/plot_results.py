"""
    Quick plot of the avalanche outputs

"""
import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from analytical_avalanche_wet import *

p_st = util.get_output('avalanche.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[20]
v2=(p2_st.y==v)

x_n = p2_st.x[v2]

u0,h0,w0,z0,p0 = analytical_sol(x_n, p2_st.time[0])
u20,h20,w20,z20,p20 = analytical_sol(x_n, p2_st.time[20])
u40,h40,w40,z40,p40 = analytical_sol(x_n, p2_st.time[40])



w0_n  = p2_st.stage[0,v2]
w20_n = p2_st.stage[20,v2]
w40_n = p2_st.stage[40,v2]

z_n = p2_st.elev[v2]

uh0_n  = p2_st.xmom[0,v2]
uh20_n = p2_st.xmom[20,v2]
uh40_n = p2_st.xmom[40,v2]

u0_n  = p2_st.xvel[0,v2]
u20_n = p2_st.xvel[20,v2]
u40_n = p2_st.xvel[40,v2]

#Plot stages
pyplot.clf()
pyplot.plot(x_n, w0_n,'b.-', label='numerical stage')
pyplot.plot(x_n, w20_n, 'b.-')
pyplot.plot(x_n, w40_n, 'b.-')
pyplot.plot(x_n, w0,'r-', label='analytical stage')
pyplot.plot(x_n, w20,'r-')
pyplot.plot(x_n, w40,'r-')
pyplot.plot(x_n, p2_st.elev[v2],'k-', label='bed elevation')
pyplot.title('Stage at several instants in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png', dpi=600)
#pyplot.show()


#Plot xmomentum
pyplot.clf()
pyplot.plot(x_n, uh0_n, 'b.-', label='numerical')
pyplot.plot(x_n, uh20_n, 'b.-')
pyplot.plot(x_n, uh40_n,'b.-')
pyplot.plot(x_n, u0*h0,'r-', label='analytical')
pyplot.plot(x_n, u20*h20,'r-')
pyplot.plot(x_n, u40*h40,'r-')
pyplot.title('Xmomentum at several instants in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('xmom_plot.png')


#Plot velocities
pyplot.clf()
pyplot.plot(x_n, u0_n, 'b.-', label='numerical')
pyplot.plot(x_n, u20_n, 'b.-')
pyplot.plot(x_n, u40_n,'b.-')
pyplot.plot(x_n, u0,'r-', label='analytical')
pyplot.plot(x_n, u20,'r-')
pyplot.plot(x_n, u40,'r-')
pyplot.title('Xvelocity at several instants in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_plot.png')
