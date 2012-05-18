from anuga.utilities import plot_utils as util
from matplotlib import pyplot as pyplot

p2=util.get_output('data_wave.sww',0.001)
p=util.get_centroids(p2,velocity_extrapolation=True)

maxx=p.x.max()

minx=p.x.min()

v1 = abs(p.x -minx).argmin()
v2 = abs(p.x -0.5*(minx+maxx)).argmin()
v3 = abs(p.x -maxx).argmin()

pyplot.plot(p.time, p.stage[:,v1], label='Left edge of domain')
pyplot.plot(p.time, p.stage[:,v2], label='Middle of domain')
pyplot.plot(p.time, p.stage[:,v3], label='Right edge of domain')
pyplot.ylim((-1.5,2.0))
pyplot.legend()
pyplot.savefig('wave_atten.png')
