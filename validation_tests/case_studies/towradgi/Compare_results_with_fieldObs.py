"""

Compare a Towradgi model run with the 1998 field observations, and write the
figures that results.tex includes:

  Spatial_Depth_and_Error.png  peak depth, with the error at each flood mark
  Error_peakstage.png          modelled against observed levels, and errors
  Pioneer_Bridge_Stage.png     modelled and recorded Pioneer Rd gauge stage

Usage (from this directory):

  python Compare_results_with_fieldObs.py [RUN_DIR_OR_SWW]

With no argument it uses the newest anuga_toml_run output,
OUTPUT/RUN_*_Towradgi_historic_flood/, falling back to the run_towradgi.py
output, MODEL_OUTPUTS/Towradgi_historic_flood.sww.

"""
import glob
import os
import sys

import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from anuga.utilities import plot_utils as util

# The Pioneer Rd logger clock runs 9 h (540 min) ahead of the rainfall/model
# clock. The original 2014 version of this script used the same shift, and a
# least-squares fit of the model to the record gives it too.
LOGGER_SHIFT_MIN = 540.0


def find_sww(arg=None):
    if arg is not None:
        if os.path.isdir(arg):
            return glob.glob(os.path.join(arg, '*.sww'))[0]
        return arg
    runs = sorted(glob.glob('OUTPUT/RUN_*_Towradgi_historic_flood/'
                            'Towradgi_historic_flood.sww'))
    if runs:
        return runs[-1]
    return 'MODEL_OUTPUTS/Towradgi_historic_flood.sww'


sww = find_sww(sys.argv[1] if len(sys.argv) > 1 else None)
print('Using', sww)

p = util.get_output(sww)
p2 = util.get_centroids(p, velocity_extrapolation=True)
# Model time on the rainfall clock (the run may start before t=0 as a burn-in)
model_hours = (p2.time + p.starttime) / 3600.

obs = numpy.genfromtxt('Validation/historic_1998_flood_levels_towradgi_ck.csv',
                       delimiter=',', skip_header=1, dtype=None, encoding=None)
obs_x = numpy.array([r[0] for r in obs])
obs_y = numpy.array([r[1] for r in obs])
obs_name = [str(r[2]) for r in obs]
obs_level = numpy.array([r[3] for r in obs])

pioneer = numpy.genfromtxt('Validation/pioneer_timeseries.txt', skip_header=1)

# Peak stage and depth per centroid. anuga_toml_run tracks the maxima every
# timestep in Max_quantities_*.csv; otherwise use the stored SWW frames.
max_csv = glob.glob(os.path.join(os.path.dirname(sww), 'Max_quantities_*.csv'))
if max_csv:
    mq = numpy.loadtxt(max_csv[0], delimiter=',', skiprows=1)
    cx, cy, peak_stage, peak_depth = mq[:, 0], mq[:, 1], mq[:, 2], mq[:, 3]
else:
    cx, cy = p2.x + p.xllcorner, p2.y + p.yllcorner
    peak_stage = p2.stage.max(axis=0)
    peak_depth = p2.height.max(axis=0)

# Modelled peak at the centroid nearest each flood mark
nearest = numpy.array([((cx - x)**2 + (cy - y)**2).argmin()
                       for x, y in zip(obs_x, obs_y)])
modelled_level = peak_stage[nearest]
err = modelled_level - obs_level

print('Flood marks: %d' % len(err))
print('  RMS error            %.3f m' % numpy.sqrt((err**2).mean()))
print('  mean error (bias)    %+.3f m' % err.mean())
print('  mean absolute error  %.3f m' % numpy.abs(err).mean())
print('  largest error        %.2f m (mark %s)'
      % (numpy.abs(err).max(), obs_name[numpy.abs(err).argmax()]))
print('  within 0.3 m / 0.5 m %d / %d'
      % ((numpy.abs(err) <= 0.3).sum(), (numpy.abs(err) <= 0.5).sum()))

# ---- Pioneer Rd gauge -------------------------------------------------------
gauge = ((p2.x - (obs_x[0] - p.xllcorner))**2 +
         (p2.y - (obs_y[0] - p.yllcorner))**2).argmin()
model_stage = p2.stage[:, gauge]
record_hours = (pioneer[:, 0] - LOGGER_SHIFT_MIN) / 60.
record_stage = pioneer[:, 2]

inside = (record_hours >= model_hours[0]) & (record_hours <= model_hours[-1]) \
    & (record_hours >= 0)
gauge_err = numpy.interp(record_hours[inside], model_hours, model_stage) \
    - record_stage[inside]
print('Pioneer Rd gauge: peak %.2f m at %.2f h (recorded %.2f m at %.2f h), '
      'RMS error %.2f m from t = 0'
      % (model_stage.max(), model_hours[model_stage.argmax()],
         record_stage.max(), record_hours[record_stage.argmax()],
         numpy.sqrt((gauge_err**2).mean())))

fig, ax = pyplot.subplots(figsize=(9, 3.8))
ax.plot(record_hours, record_stage, color='0.15', lw=1.6,
        label='Recorded (logger clock - 9 h)')
ax.plot(model_hours, model_stage, color='#2a6fdb', lw=1.8, label='ANUGA')
if model_hours[0] < 0:
    ax.axvspan(model_hours[0], 0, color='0.9', zorder=0)
    ax.text(model_hours[0] / 2, 4.2, 'burn-in', ha='center', fontsize=8,
            color='0.4')
ax.set_xlim(model_hours[0], model_hours[-1])
ax.set_xlabel('Time from start of rain (h)')
ax.set_ylabel('Stage (m AHD)')
ax.set_title('Pioneer Road bridge gauge')
ax.legend(frameon=False, loc='upper right')
ax.grid(alpha=.3)
fig.tight_layout()
fig.savefig('Pioneer_Bridge_Stage.png', dpi=150)

# ---- Modelled against observed, and the error distribution -------------------
fig, (a1, a2) = pyplot.subplots(1, 2, figsize=(10, 4.2))
lo, hi = 0, 45
a1.scatter(obs_level, modelled_level, s=18, color='#2a6fdb')
a1.plot([lo, hi], [lo, hi], 'k-', lw=0.8)
a1.plot([lo, hi], [lo + 0.5, hi + 0.5], 'k:', lw=0.6)
a1.plot([lo, hi], [lo - 0.5, hi - 0.5], 'k:', lw=0.6)
a1.set_xlabel('Observed flood level (m AHD)')
a1.set_ylabel('Modelled peak level (m AHD)')
a1.set_title('Modelled vs observed (dotted: ±0.5 m)')
a2.hist(err, bins=numpy.arange(-0.8, 0.81, 0.1), color='#2a6fdb',
        edgecolor='white')
a2.axvline(0, color='k', lw=0.8)
a2.set_xlabel('Model minus observed (m)')
a2.set_ylabel('Number of marks')
a2.set_title('Error distribution')
fig.tight_layout()
fig.savefig('Error_peakstage.png', dpi=150)

# ---- Peak depth with the error at each mark ------------------------------------
fig, ax = pyplot.subplots(figsize=(10, 6.0))
wet = peak_depth > 0.01
sc = ax.scatter(cx[wet] - 300000, cy[wet] - 6190000,
                c=numpy.minimum(peak_depth[wet], 4), s=0.3, cmap='Blues',
                vmin=0, vmax=4, rasterized=True)
cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.08, orientation='horizontal',
                  shrink=0.45, anchor=(0.0, 1.0))
cb.set_label('Peak depth (m)')
se = ax.scatter(obs_x - 300000, obs_y - 6190000, c=err, cmap='RdBu_r',
                vmin=-0.8, vmax=0.8, s=36, edgecolor='k', lw=0.5, zorder=3)
fig.colorbar(se, ax=ax, fraction=0.025, pad=0.02).set_label(
    'Model minus observed (m)')
ax.set_aspect('equal')
ax.set_xlabel('Easting - 300000 (m)')
ax.set_ylabel('Northing - 6190000 (m)')
ax.set_title('Peak flood depth and error at the %d surveyed flood marks'
             % len(err))
fig.tight_layout()
fig.savefig('Spatial_Depth_and_Error.png', dpi=150)

print('Wrote Pioneer_Bridge_Stage.png, Error_peakstage.png, '
      'Spatial_Depth_and_Error.png')
