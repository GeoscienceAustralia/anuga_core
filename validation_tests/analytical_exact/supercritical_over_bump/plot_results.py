"""
    Quick plot of the subcritical-flow outputs

"""

import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from analytical_supercritical import analytic_sol
from numpy import ones, arange

import anuga
parser = anuga.create_standard_parser()

parser.add_argument('-tid', type=int, default=-1, help ='timestep id')
args = parser.parse_args()
tid = args.tid
verbose = args.verbose


if verbose: print('Read in swwfile')
p_st = util.get_output('supercritical.sww')
p2_st=util.get_centroids(p_st)

#v = p2_st.y[10]
#v2=(p2_st.y==v)
#v2=(p2_st.y>-1.0)
v2 = arange(len(p2_st.y))

if verbose: print('Calculate analytical solution')
h,z = analytic_sol(p2_st.x[v2])
qexact = 10




#Plot the stages##############################################################
if verbose: print('Create Stage plot')
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.stage[tid,v2], 'b.-', label='numerical stage') # 0*T/6
pyplot.plot(p2_st.x[v2], h+z,'r-', label='analytical stage')
pyplot.plot(p2_st.x[v2], z,'k-', label='bed elevation')
pyplot.title('Stage at time = %s secs'% p2_st.time[tid])
##pyplot.ylim(-5.0,5.0)
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png')


#Plot the momentums##########################################################
if verbose: print('Create Momentum plot')
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xmom[tid,v2], 'b.-', label='numerical') # 0*T/6
pyplot.plot(p2_st.x[v2], qexact*ones(len(p2_st.x[v2])),'r-', label='analytical')
pyplot.title('Xmomentum at time = %s secs'% p2_st.time[tid])
pyplot.legend(loc='best')
##pyplot.ylim([1.52,1.54])
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('xmom_plot.png')



#Plot the velocities#########################################################
if verbose: print('Create Velocity plot')
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xvel[tid,v2], 'b.-', label='numerical') # 0*T/6
pyplot.plot(p2_st.x[v2], qexact/h,'r-', label='analytical')
pyplot.title('Xvelocity at time = %s secs'% p2_st.time[tid])
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('xvel_plot.png')

