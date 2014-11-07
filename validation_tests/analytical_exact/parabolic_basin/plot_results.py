"""
    Quick plot of the dam break outputs

"""
import numpy
from numpy import zeros
import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from math import sqrt, pi, cos, sin
from analytical_parabolic_basin import analytic_cannal

# Get the sww file
p_st = util.get_output('parabola.sww')
p2_st=util.get_centroids(p_st, velocity_extrapolation=True) #(p_st, velocity_extrapolation=True)

index_origin = 3978
v = p2_st.y[index_origin]
v2=(p2_st.y==v)


# Compute the analytical solution
time_level = 200
u,h,w,z = analytic_cannal(p2_st.x[v2], p2_st.time[time_level])

#Plot stages
pyplot.clf()
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
pyplot.plot(p2_st.x[v2], p2_st.xmom[time_level,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], u*h,'r-', label='analytical')
pyplot.title('Xmomentum at an instant in time')
pyplot.legend(loc='lower right')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xmomentum')
pyplot.savefig('cross_section_xmom.png')


#Plot velocities
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xvel[time_level,v2], 'b.', label='numerical')
pyplot.plot(p2_st.x[v2], u,'r-', label='analytical')
pyplot.title('Xvelocity at an instant in time')
pyplot.legend(loc='upper center')
pyplot.xlabel('Xposition')
pyplot.ylabel('Xvelocity')
pyplot.savefig('cross_section_xvel.png')




# Compute the analytical solution
D0=4.0
L=10.0
A=2.0
g=9.81

omega=numpy.sqrt(2*D0*g)/L
T= 2*numpy.pi/omega
index_origin = 3978

# Free surface over time at the origin of spatial domain.
w = D0 + 2*A*D0/(L**2)*numpy.cos(omega*p2_st.time)*( p2_st.x[index_origin] -A/2.*numpy.cos(omega*p2_st.time))
#w = w*(w>p2_st.elev[index_origin])+p2_st.elev[index_origin]*(w<=p2_st.elev[index_origin]) #not needed

# Velocity at the origin of the spatial domain.
u = -A*omega*numpy.sin(omega*p2_st.time)


#Plot stages
pyplot.clf()
pyplot.plot(p2_st.time,w, color='blue', label='analytical')
pyplot.plot(p2_st.time,p2_st.stage[:,index_origin], color='green', label='numerical')
pyplot.legend(loc='best')
pyplot.title('Stage at the centre of the parabolic basin over time')
pyplot.xlabel('Time')
pyplot.ylabel('Stage')
pyplot.savefig('Stage_centre.png')

#Plot xmomentums
pyplot.clf()
pyplot.plot(p2_st.time,u*w, color='blue', label='analytical') # At the origin of space, w=h as always wet.
pyplot.plot(p2_st.time,p2_st.xmom[:,index_origin], color='green', label='numerical')
pyplot.legend(loc='best')
pyplot.title('Xmomentum at the centre of the parabolic basin over time')
pyplot.xlabel('Time')
pyplot.ylabel('Xmomentum')
pyplot.savefig('Xmom_centre.png')

#Plot velocities
pyplot.clf()
pyplot.plot(p2_st.time,u, color='blue', label='analytical')
pyplot.plot(p2_st.time,p2_st.xvel[:,index_origin], color='green', label='numerical')
pyplot.legend(loc='best')
pyplot.title('Xvelocity at the centre of the parabolic basin over time')
pyplot.xlabel('Time')
pyplot.ylabel('Xvelocity')
pyplot.savefig('Xvel_centre.png')

