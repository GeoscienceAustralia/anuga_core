"""
    Quick plot of the avalanche outputs

"""
import anuga.utilities.plot_utils as util
from matplotlib import pyplot as pyplot
from analytical_avalanche_wet import *

p_st = util.get_output('avalanche.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[20]
v2=(p2_st.y==v)

u0,h0,w0,z0,p0 = analytical_sol(p2_st.x[v2], p2_st.time[0])
u20,h20,w20,z20,p20 = analytical_sol(p2_st.x[v2], p2_st.time[20])
u40,h40,w40,z40,p40 = analytical_sol(p2_st.x[v2], p2_st.time[40])

#Plot stages
pyplot.clf()
#pyplot.ion()
pyplot.plot(p2_st.x[v2], p2_st.stage[0,v2],'b.-', label='numerical stage')
pyplot.plot(p2_st.x[v2], p2_st.stage[20,v2], 'b.-')
pyplot.plot(p2_st.x[v2], p2_st.stage[40,v2], 'b.-')
pyplot.plot(p2_st.x[v2], w0,'r-', label='analytical stage')
pyplot.plot(p2_st.x[v2], w20,'r-')
pyplot.plot(p2_st.x[v2], w40,'r-')
pyplot.plot(p2_st.x[v2], p2_st.elev[v2],'k-', label='bed elevation')
pyplot.title('Stage at several instants in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png', dpi=600)
#pyplot.show()


#Plot xmomentum
pyplot.clf()
#pyplot.ion()
pyplot.plot(p2_st.x[v2], p2_st.xmom[0,v2], 'b.-', label='numerical')
pyplot.plot(p2_st.x[v2], p2_st.xmom[20,v2], 'b.-')
pyplot.plot(p2_st.x[v2], p2_st.xmom[40,v2],'b.-')
pyplot.plot(p2_st.x[v2], u0*h0,'r-', label='analytical')
pyplot.plot(p2_st.x[v2], u20*h20,'r-')
pyplot.plot(p2_st.x[v2], u40*h40,'r-')
pyplot.title('Xmomentum at several instants in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('xmom_plot.png')


#Plot velocities
pyplot.clf()
#pyplot.ion()
pyplot.plot(p2_st.x[v2], p2_st.xvel[0,v2], 'b.-', label='numerical')
pyplot.plot(p2_st.x[v2], p2_st.xvel[20,v2], 'b.-')
pyplot.plot(p2_st.x[v2], p2_st.xvel[40,v2],'b.-')
pyplot.plot(p2_st.x[v2], u0,'r-', label='analytical')
pyplot.plot(p2_st.x[v2], u20,'r-')
pyplot.plot(p2_st.x[v2], u40,'r-')
pyplot.title('Xvelocity at several instants in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_plot.png')
