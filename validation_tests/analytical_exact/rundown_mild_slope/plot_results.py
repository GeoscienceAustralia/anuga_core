"""
View results of numerical_rundown_channel.py

Simple water flow example using ANUGA: Water flowing down a channel.
It was called "steep_slope" in an old validation test.

Two cases: 'stable' (Froude 1.21) is compared with the steady uniform
solution; 'rollwave' (Froude 2.26) is unstable, and the plots show the roll
waves growing down the slope from the uniform state.
"""
#---------------
# Import Modules
#---------------
import anuga
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from anuga.utilities import plot_utils as util

#-------------------------------
# Define variables of case study
#-------------------------------
bedslope=-0.1
fluxin=20./100. #The momentum flux at the upstream boundary ( = discharge / width)

def analytic(mann):
    """Steady uniform flow: Friction_slope=bedslope and depth*velocity = flux_in."""
    uana= ( mann**(-2.)*abs(bedslope)*fluxin**(4./3.) )**(3./10.) # Velocity
    dana= fluxin/uana # Depth
    return uana, dana

def centreline_profile(p2, xbins, it):
    """Depth averaged over the cells within 1.5 m of the centreline, in x bins."""
    yc = 0.5*(p2.y.min() + p2.y.max())
    band = abs(p2.y - yc) < 1.5
    k = numpy.digitize(p2.x[band], xbins) - 1
    depth = p2.stage[it][band] - p2.elev[band]
    return numpy.array([depth[k == i].mean() for i in range(len(xbins) - 1)])

#==============================================================
# Stable case: n = 0.06
#==============================================================
mann=0.06
uana, dana = analytic(mann)

p=util.get_output('channel_stable.sww', 0.001)
p2=util.get_centroids(p,velocity_extrapolation=True)

# Find an y value close to y==50
tmp=(abs(p2.y-50.)).argmin()
vx=(abs(p2.y - p2.y[tmp])<1.5)

pyplot.clf()
pyplot.plot(p2.x[vx],p2.stage[-1,vx]-p2.elev[vx], 'o', label='numerical')
pyplot.plot((0,100),(dana,dana),label='analytical')
pyplot.ylim([0.7*dana,1.3*dana])
pyplot.xlabel('Xposition m')
pyplot.ylabel('Depth m')
pyplot.title('Stable case: depth down the slope (along y=50.)')
pyplot.legend(loc='best')
pyplot.savefig('depth_x.png')

#--------------------------------------------
# Compare velocity with analytical solution
#--------------------------------------------
# Find an x value close to x==50
tmp=(abs(p2.x-50.)).argmin()
v=(abs(p2.x - p2.x[tmp])<1.5)

pyplot.clf()
pyplot.plot(p2.y[v],p2.stage[-1,v]-p2.elev[v], 'o', label='numerical')
pyplot.plot((0,100),(dana,dana),label='analytical')
pyplot.ylim([0.7*dana,1.3*dana])
pyplot.xlabel('Yposition m')
pyplot.ylabel('Depth m')
pyplot.title('Stable case: depth across the slope (x=50.)')
pyplot.legend(loc='best')
pyplot.savefig('depth_y.png')


pyplot.clf()
pyplot.plot(p2.y[v],p2.xvel[-1,v], 'o', label='numerical')
pyplot.plot((0,100),(uana,uana),label='analytical')
pyplot.ylim([0.7*uana,1.3*uana])
pyplot.xlabel('Yposition along the line x=50')
pyplot.ylabel('Xvelocity m/s')
pyplot.title('Stable case: final Xvelocity around the line x=50.')
pyplot.legend(loc='best')
pyplot.savefig('x_velocity.png')

pyplot.clf()
pyplot.plot(p2.y[v],p2.yvel[-1,v],'o', label='numerical')
pyplot.plot((0,100),(0.0, 0.0),label='analytical')
pyplot.xlabel('Yposition along the line x=50')
pyplot.ylabel('Yvelocity')
pyplot.title('Stable case: final Yvelocity around the line x=50.')
pyplot.legend(loc='best')
pyplot.savefig('y_velocity.png')

#==============================================================
# Roll-wave case: n = 0.03
#==============================================================
mann=0.03
uana, dana = analytic(mann)

p=util.get_output('channel_rollwave.sww', 0.001)
p2=util.get_centroids(p,velocity_extrapolation=True)

length = p2.x.max() + p2.x.min()
xbins = numpy.arange(0.0, length + 0.5, 1.0)
xmid = 0.5*(xbins[1:] + xbins[:-1])
t = numpy.asarray(p2.time)
profiles = numpy.array([centreline_profile(p2, xbins, it) for it in range(len(t))])

# Envelope over the last quarter of the run
late = t >= 0.75*t[-1]

pyplot.clf()
pyplot.fill_between(xmid, profiles[late].min(axis=0), profiles[late].max(axis=0),
                    color='C0', alpha=0.25,
                    label='range over t = %g-%g s' % (t[late][0], t[-1]))
pyplot.plot(xmid, profiles[-1], 'C0-', label='numerical, t = %g s' % t[-1])
pyplot.plot((0,length),(dana,dana),'C1', label='steady uniform (unstable)')
pyplot.ylim([0.0, 2.5*dana])
pyplot.xlabel('Xposition m')
pyplot.ylabel('Depth m')
pyplot.title('Roll-wave case: depth down the slope (centreline)')
pyplot.legend(loc='upper left')
pyplot.savefig('rollwave_depth_x.png')

# Space-time view of the last 30 s, which the simulation stores every 0.5 s
fine = t >= t[-1] - 30.0 - 1.0e-6
pyplot.clf()
pyplot.pcolormesh(xmid, t[fine], profiles[fine]/dana, cmap='RdBu_r', vmin=0.0, vmax=2.0,
                  shading='nearest')
pyplot.colorbar(label='depth / steady uniform depth')
pyplot.xlabel('Xposition m')
pyplot.ylabel('Time s')
pyplot.title('Roll-wave case: centreline depth, last 30 s')
pyplot.savefig('rollwave_xt.png')

# Wave amplitude down the slope: half the late-time range, relative to dana
pyplot.clf()
amp = 0.5*(profiles[late].max(axis=0) - profiles[late].min(axis=0))/dana
pyplot.semilogy(xmid, numpy.maximum(amp, 1.0e-6))
pyplot.xlabel('Xposition m')
pyplot.ylabel('half range of depth / steady uniform depth')
pyplot.title('Roll-wave case: growth of the wave amplitude')
pyplot.grid(True, which='both', alpha=0.3)
pyplot.savefig('rollwave_amplitude.png')
