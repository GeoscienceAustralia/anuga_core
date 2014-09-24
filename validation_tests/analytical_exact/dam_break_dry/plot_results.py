"""
    Quick plot of the dam break outputs

"""
import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
import analytical_dam_break_dry as analytic
import numpy

p_st = util.get_output('dam_break.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[10]
v2=numpy.argwhere(p2_st.y==v).flatten()
iv2 = numpy.argsort(p2_st.x[v2])
v2 = v2[iv2]

#p_dev = util.get_output('dam_break.sww', 0.001)
#p2_dev=util.get_centroids(p_dev, velocity_extrapolation=True)

h0 = 1e-9
h1 = 10.0

h10,u10 = analytic.vec_dam_break(p2_st.x[v2], p2_st.time[10], h0=h0, h1=h1)
h50,u50 = analytic.vec_dam_break(p2_st.x[v2], p2_st.time[50], h0=h0, h1=h1)
#h100,u100 = analytic.vec_dam_break(p2_st.x[v2], p2_st.time[100], h0=h0, h1=h1)

#Plot stages
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.stage[10,v2],'b.', label='numerical t=10')
pyplot.plot(p2_st.x[v2], p2_st.stage[50,v2], 'g.', label ='numerical t=50')
#pyplot.plot(p2_st.x[v2], p2_st.stage[100,v2], 'b.')
pyplot.plot(p2_st.x[v2], h10,'r-', label='analytical t=10')
pyplot.plot(p2_st.x[v2], h50,'y-', label='analytical t=50')
#pyplot.plot(p2_st.x[v2], h100,'r.')
pyplot.title('Stage at several instants in time')
pyplot.legend(loc=3)
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png')


#Plot xmomentum
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xmom[10,v2], 'b.', label='numerical t=10')
pyplot.plot(p2_st.x[v2], p2_st.xmom[50,v2], 'g.', label='numerical t=50')
#pyplot.plot(p2_st.x[v2], p2_st.xmom[100,v2],'b.')
pyplot.plot(p2_st.x[v2], u10*h10,'r-', label='analytical t=10')
pyplot.plot(p2_st.x[v2], u50*h50,'y-', label='analytical t=50')
#pyplot.plot(p2_st.x[v2], u100*h100,'r-')
pyplot.title('Xmomentum at several instants in time')
pyplot.legend(loc=2)
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('xmom_plot.png')


#Plot velocities
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xvel[10,v2], 'b.', label='numerical t=10')
pyplot.plot(p2_st.x[v2], p2_st.xvel[50,v2], 'g.', label='numerical t=50')
#pyplot.plot(p2_st.x[v2], p2_st.xvel[100,v2],'b.')
pyplot.plot(p2_st.x[v2], u10,'r-', label='analytical t=10')
pyplot.plot(p2_st.x[v2], u50,'y-', label='analytical t=50')
#pyplot.plot(p2_st.x[v2], u100,'r-')
pyplot.title('Xvelocity at several instants in time')
pyplot.legend(loc=2)
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_plot.png')


