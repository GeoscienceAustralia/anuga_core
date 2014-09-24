"""
    Quick plot of the dam break outputs

"""
import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from numpy import ones, zeros

p_st = util.get_output('immersed_bump.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[10]
v2=(p2_st.y==v)

#p_dev = util.get_output('dam_break.sww', 0.001)
#p2_dev=util.get_centroids(p_dev, velocity_extrapolation=True)

#Plot stages
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.stage[10,v2],'b.', label='numerical stage')
pyplot.plot(p2_st.x[v2], p2_st.elev[v2],'k-', label='discretised bed')
pyplot.plot(p2_st.x[v2], 0.5*ones(len(p2_st.x[v2])),'r-', label='analytical stage')
pyplot.title('Stage at an instant in time')
pyplot.legend(loc='best')
pyplot.ylim([0.0,0.6])
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png')


#Plot xmomentum
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xmom[10,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], zeros(len(p2_st.x[v2])),'r-', label='analytical')
pyplot.title('Xmomentum at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('xmom_plot.png')


#Plot velocities
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xvel[10,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], zeros(len(p2_st.x[v2])),'r-', label='analytical')
pyplot.title('Xvelocity at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_plot.png')
