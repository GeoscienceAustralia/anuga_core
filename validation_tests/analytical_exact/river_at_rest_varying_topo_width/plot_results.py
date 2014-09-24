"""
    Quick plot of the dam break outputs

"""
import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from numpy import ones, zeros, ones_like

p_st = util.get_output('varying_width.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[116]
v2=(p2_st.y==v)


#p_dev = util.get_output('dam_break.sww', 0.001)
#p2_dev=util.get_centroids(p_dev, velocity_extrapolation=True)

#Plot stages
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.stage[-1,v2],'b.', label='numerical stage')
pyplot.plot(p2_st.x[v2], 12.*ones_like(p2_st.x[v2]),'r--', label='analytical stage')
pyplot.plot(p2_st.x[v2], p2_st.elev[v2],'k-', label='discretised bed')
pyplot.title('Stage at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.xlim([0,1500])
pyplot.ylim([0,25])
pyplot.savefig('stage_plot.png')


#Plot xmomentum
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xmom[-1,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], zeros(len(p2_st.x[v2])),'r-', label='analytical')
pyplot.title('Xmomentum at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.xlim([0,1500])
pyplot.savefig('xmom_plot.png')


#Plot velocities
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xvel[-1,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], zeros(len(p2_st.x[v2])),'r-', label='analytical')
pyplot.title('Xvelocity at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.xlim([0,1500])
pyplot.savefig('xvel_plot.png')




# Sweep in y direction

v = p2_st.x[116]
v2=(p2_st.x==v)


#p_dev = util.get_output('dam_break.sww', 0.001)
#p2_dev=util.get_centroids(p_dev, velocity_extrapolation=True)

#Plot stages
pyplot.clf()
pyplot.plot(p2_st.y[v2], p2_st.stage[-1,v2],'b.', label='numerical stage')
pyplot.plot(p2_st.y[v2], 12.*ones_like(p2_st.x[v2]),'r--', label='analytical stage')
pyplot.plot(p2_st.y[v2], p2_st.elev[v2],'k-', label='discretised bed')
pyplot.title('Stage at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Yposition')
pyplot.ylabel('Stage')
pyplot.xlim([-30,30])
pyplot.ylim([-1,26])
pyplot.savefig('stage_plot_y.png')


#Plot xmomentum
pyplot.clf()
pyplot.plot(p2_st.y[v2], p2_st.xmom[-1,v2], 'b.', label='numerical')
pyplot.plot(p2_st.y[v2], zeros(len(p2_st.x[v2])),'r-', label='analytical')
pyplot.title('Xmomentum at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Yposition')
pyplot.ylabel('Xmomentum')
pyplot.xlim([-30,30])
pyplot.savefig('xmom_plot_y.png')


#Plot velocities
pyplot.clf()
pyplot.plot(p2_st.y[v2], p2_st.xvel[-1,v2], 'b.', label='numerical')
pyplot.plot(p2_st.y[v2], zeros(len(p2_st.x[v2])),'r-', label='analytical')
pyplot.title('Xvelocity at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Yposition')
pyplot.ylabel('Xvelocity')
pyplot.xlim([-30,30])
pyplot.savefig('xvel_plot_y.png')
