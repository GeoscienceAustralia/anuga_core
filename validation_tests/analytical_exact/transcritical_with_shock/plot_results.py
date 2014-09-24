"""
    Quick plot for outputs of the transcritical flow with a shock

"""
import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from analytical_with_shock import *
from numpy import ones

p_st = util.get_output('transcritical.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[10]
v2=(p2_st.y==v)

tid = -1

h,z = analytic_sol(p2_st.x[v2])

#Plot the stages##############################################################
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.stage[tid,v2], 'b.-', label='numerical stage') # 0*T/6
pyplot.plot(p2_st.x[v2], h+z,'r-', label='analytical stage')
pyplot.plot(p2_st.x[v2], z,'k-', label='bed elevation')
pyplot.title('Stage at time %s secs'% p2_st.time[tid])
##pyplot.ylim(-5.0,5.0)
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png')


#Plot the momentums##########################################################
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xmom[tid,v2], 'b.-', label='numerical') # 0*T/6
pyplot.plot(p2_st.x[v2], 0.18*ones(len(p2_st.x[v2])),'r-', label='analytical')
pyplot.title('Xmomentum at time %s secs'% p2_st.time[tid])
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('xmom_plot.png')



#Plot the velocities#########################################################
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xvel[tid,v2], 'b.-', label='numerical') # 0*T/6
pyplot.plot(p2_st.x[v2], 0.18/h,'r-', label='analytical')
pyplot.title('Xvelocity at time %s secs'% p2_st.time[tid])
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_plot.png')

