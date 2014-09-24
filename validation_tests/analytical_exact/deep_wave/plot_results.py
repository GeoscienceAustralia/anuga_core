from anuga.utilities import plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot

p2=util.get_output('data_wave.sww',0.001)
p=util.get_centroids(p2,velocity_extrapolation=True)

maxx=p.x.max()

minx=p.x.min()

v1 = abs(p.x -minx).argmin()
v2 = abs(p.x -0.5*(minx+maxx)).argmin()
v3 = abs(p.x -maxx).argmin()

pyplot.clf()
pyplot.plot(p.time, p.stage[:,v1], label='Left edge of domain')
pyplot.plot(p.time, p.stage[:,v2], label='Middle of domain')
pyplot.plot(p.time, p.stage[:,v3], label='Right edge of domain')
pyplot.ylim((-1.5,2.0))
pyplot.legend()
pyplot.xlabel('Time')
pyplot.ylabel('Amplitude')
pyplot.title('Stage over time at 3 points in space')
pyplot.savefig('wave_atten.png')
#pyplot.show()

pyplot.clf()
pyplot.plot(p.time, p.xmom[:,v1], label='Left edge of domain')
pyplot.plot(p.time, p.xmom[:,v2], label='Middle of domain')
pyplot.plot(p.time, p.xmom[:,v3], label='Right edge of domain')
#pyplot.ylim((-1.,1.0))
pyplot.legend()
pyplot.xlabel('Time')
pyplot.ylabel('Xmomentum')
pyplot.title('Xmomentum over time at 3 points in space')
pyplot.savefig('xmom.png')

pyplot.clf()
pyplot.plot(p.time, p.ymom[:,v1], label='Left edge of domain')
pyplot.plot(p.time, p.ymom[:,v2], label='Middle of domain')
pyplot.plot(p.time, p.ymom[:,v3], label='Right edge of domain')
#pyplot.ylim((-1.,1.0))
pyplot.legend()
pyplot.xlabel('Time')
pyplot.ylabel('Ymomentum')
pyplot.title('Ymomentum over time at 3 points in space')
pyplot.savefig('ymom.png')
