"""
    Quick plot of the dam break outputs
    Note: v = p2_st.y[?], where ? is the cell position of the origin.

"""
import anuga.utilities.plot_utils as util
from matplotlib import pyplot as pyplot
from math import sqrt, pi,sin,cos
from numpy import ones
from analytical_paraboloid_basin import analytic_sol

# Get the sww file
p_st = util.get_output('paraboloid.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[11248]#[4898]
v2=(p2_st.y==v)
time_level = 50
w2,u2,h2 = analytic_sol(p2_st.x[v2], v2, p2_st.time[time_level])


#Plot stages
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.stage[time_level,v2],'b.', label='numerical stage')
pyplot.plot(p2_st.x[v2], p2_st.elev[v2],'k-', label='bed elevation')
pyplot.plot(p2_st.x[v2], w2,'r-', label='analytical stage')
pyplot.title('Stage at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('cross_section_stage.png')


#Plot xmomentum
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xmom[time_level,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], u2*h2,'r-', label='analytical')
pyplot.title('Xmomentum at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('cross_section_xmom.png')


#Plot velocities
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xvel[time_level,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], u2,'r-', label='analytical')
pyplot.title('Xvelocity at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('cross_section_xvel.png')


