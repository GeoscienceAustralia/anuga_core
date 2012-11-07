"""
    Quick plot of the dam break outputs

"""
from numpy import zeros
import anuga.utilities.plot_utils as util
from matplotlib import pyplot as pyplot
from math import sqrt, pi, cos, sin
from analytical_parabolic_basin import analytic_cannal

# Get the sww file
p_st = util.get_output('parabola.sww')
p2_st=util.get_centroids(p_st) #(p_st, velocity_extrapolation=True)

index_origin = 3978
v = p2_st.y[index_origin]
v2=(p2_st.y==v)


# Compute the analytical solution
time_level = 200
u,h,w,z = analytic_cannal(p2_st.x[v2], p2_st.time[time_level])

#Plot stages
pyplot.clf()
pyplot.ion()
pyplot.plot(p2_st.x[v2], p2_st.stage[time_level,v2],'b.', label='numerical stage')
pyplot.plot(p2_st.x[v2], p2_st.elev[v2],'k-', label='bed elevation')
pyplot.plot(p2_st.x[v2], w,'r-', label='analytical stage')
pyplot.title('Stage at an instant in time')
pyplot.legend(loc='upper center')
pyplot.xlabel('Xposition')
pyplot.ylabel('Stage')
pyplot.savefig('cross_section_stage.png')


#Plot xmomentum
pyplot.clf()
pyplot.ion()
pyplot.plot(p2_st.x[v2], p2_st.xmom[time_level,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], u*h,'r-', label='analytical')
pyplot.title('Xmomentum at an instant in time')
pyplot.legend(loc='lower right')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('cross_section_xmom.png')


#Plot velocities
pyplot.clf()
pyplot.ion()
pyplot.plot(p2_st.x[v2], p2_st.xvel[time_level,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], u,'r-', label='analytical')
pyplot.title('Xvelocity at an instant in time')
pyplot.legend(loc='upper center')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('cross_section_xvel.png')


