"""
    Quick plot of the dam break outputs

"""
import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from numpy import ones, zeros

p_st = util.get_output('steep_island.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[10]
v2=(p2_st.y==v)

z=zeros(len(p2_st.x[v2]))
xx = p2_st.x[v2]
for i in range(len(p2_st.x[v2])):
    if 0 <= xx[i] <= 200.0:
        z[i] = max(-0.01*(xx[i]-200) + 4.0, 4.5)
    elif 900.0 <= xx[i] <= 1000.0:
        z[i] = 6.0
    elif 1800.0 <= xx[i] <= 2000.0:
        z[i] = max((4.5/40000)*(xx[i]-1800)**2 + 2.0, 4.5)
    else:
        z[i] = 4.5

#p_dev = util.get_output('dam_break.sww', 0.001)
#p2_dev=util.get_centroids(p_dev, velocity_extrapolation=True)

#Plot stages
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.stage[40,v2],'b.', label='numerical stage')
pyplot.plot(p2_st.x[v2], p2_st.elev[v2],'k-', label='discretised bed')
pyplot.plot(p2_st.x[v2], z,'r-', label='analytical stage')
pyplot.title('Stage at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png')


#Plot xmomentum
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xmom[40,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], zeros(len(p2_st.x[v2])),'r-', label='analytical')
pyplot.title('Xmomentum at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('xmom_plot.png')


#Plot velocities
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xvel[40,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], zeros(len(p2_st.x[v2])),'r-', label='analytical')
pyplot.title('Xvelocity at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_plot.png')
