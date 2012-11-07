"""
    Quick plot of the dam break outputs

"""
import anuga.utilities.plot_utils as util
from matplotlib import pyplot as pyplot
import analytical_dam_break_wet as analytic

p_st = util.get_output('dam_break.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[10]
v2=(p2_st.y==v)


h0 = 1.0
h1 = 10.0

h10,u10 = analytic.vec_dam_break(p2_st.x[v2], p2_st.time[10], h0=h0, h1=h1)
h50,u50 = analytic.vec_dam_break(p2_st.x[v2], p2_st.time[50], h0=h0, h1=h1)
h100,u100 = analytic.vec_dam_break(p2_st.x[v2], p2_st.time[100], h0=h0, h1=h1)

#Plot stages
pyplot.clf()
pyplot.ion()
pyplot.plot(p2_st.x[v2], p2_st.stage[10,v2],'b.', label='numerical')
pyplot.plot(p2_st.x[v2], p2_st.stage[50,v2], 'b.')
pyplot.plot(p2_st.x[v2], p2_st.stage[100,v2], 'b.')
pyplot.plot(p2_st.x[v2], h10,'r-', label='analytical')
pyplot.plot(p2_st.x[v2], h50,'r-')
pyplot.plot(p2_st.x[v2], h100,'r-')
pyplot.title('Stage at several instants in time')
pyplot.legend(loc=3)
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png')


#Plot xmomentums
pyplot.clf()
pyplot.ion()
pyplot.plot(p2_st.x[v2], p2_st.xmom[10,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], p2_st.xmom[50,v2], 'b.')
pyplot.plot(p2_st.x[v2], p2_st.xmom[100,v2],'b.')
pyplot.plot(p2_st.x[v2], u10*h10,'r-', label='analytical')
pyplot.plot(p2_st.x[v2], u50*h50,'r-')
pyplot.plot(p2_st.x[v2], u100*h100,'r-')
pyplot.title('Xmomentum at several instants in time')
pyplot.legend(loc=2)
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('xmom_plot.png')


#Plot velocities
pyplot.clf()
pyplot.ion()
pyplot.plot(p2_st.x[v2], p2_st.xvel[10,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], p2_st.xvel[50,v2], 'b.')
pyplot.plot(p2_st.x[v2], p2_st.xvel[100,v2],'b.')
pyplot.plot(p2_st.x[v2], u10,'r-', label='analytical')
pyplot.plot(p2_st.x[v2], u50,'r-')
pyplot.plot(p2_st.x[v2], u100,'r-')
pyplot.title('Xvelocity at several instants in time')
pyplot.legend(loc=2)
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_plot.png')
