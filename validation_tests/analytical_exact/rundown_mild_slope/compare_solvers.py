"""
Run the roll-wave case of numerical_rundown_channel.py with each flow
algorithm, and plot how each one grows the same inflow disturbance.

All the solvers hold the stable case (see the main report figures), so the
comparison is made where they differ: the growth of the roll waves. The
physical growth is the same for all of them, so differences in the onset and
amplitude come from the scheme's numerical damping.

    python compare_solvers.py [-np N] [-nr] [-v]

-nr reuses the channel_rollwave_<alg>.sww files already on disk.
"""
import os

import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot

import anuga
from anuga.utilities import plot_utils as util

ALGS = ['DE0', 'DE1', 'DE2', 'DE_ader2']

bedslope = -0.1
fluxin = 20./100.
mann = 0.03
uana = (mann**(-2.)*abs(bedslope)*fluxin**(4./3.))**(3./10.)
dana = fluxin/uana


def sww_name(alg):
    return 'channel_rollwave_%s.sww' % alg


def run_solvers(args):
    # Set only for these children: produce_results.py runs the main
    # simulation afterwards from the same process.
    saved = {k: os.environ.get(k) for k in ('RUNDOWN_CASES', 'RUNDOWN_TAG')}
    try:
        for alg in ALGS:
            if args.verbose:
                print(anuga.indent + 'Running the roll-wave case with %s' % alg)
            os.environ['RUNDOWN_CASES'] = 'rollwave'
            os.environ['RUNDOWN_TAG'] = '_' + alg
            args.alg = alg
            res = anuga.run_anuga_script('numerical_rundown_channel.py', args=args)
            assert res == 0, 'roll-wave run with %s failed with return code %d' % (alg, res)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def centreline(sww):
    """Times, x bin centres, and the centreline depth in 1 m x bins."""
    p2 = util.get_centroids(util.get_output(sww, 0.001), velocity_extrapolation=True)
    yc = 0.5*(p2.y.min() + p2.y.max())
    band = abs(p2.y - yc) < 1.5
    length = p2.x.max() + p2.x.min()
    xbins = numpy.arange(0.0, length + 0.5, 1.0)
    k = numpy.digitize(p2.x[band], xbins) - 1
    depth = (p2.stage - p2.elev)[:, band]
    profiles = numpy.array([[d[k == i].mean() for i in range(len(xbins) - 1)]
                            for d in depth])
    return numpy.asarray(p2.time), 0.5*(xbins[1:] + xbins[:-1]), profiles


def plot_solvers():
    results = {alg: centreline(sww_name(alg)) for alg in ALGS}

    # Wave amplitude down the slope: half the range over the last quarter
    pyplot.clf()
    for alg in ALGS:
        t, x, profiles = results[alg]
        late = t >= 0.75*t[-1]
        amp = 0.5*(profiles[late].max(axis=0) - profiles[late].min(axis=0))/dana
        pyplot.semilogy(x, numpy.maximum(amp, 1.0e-5), label=alg)
    pyplot.ylim([1.0e-3, 2.0])
    pyplot.xlabel('Xposition m')
    pyplot.ylabel('half range of depth / steady uniform depth')
    pyplot.title('Roll-wave growth by flow algorithm')
    pyplot.grid(True, which='both', alpha=0.3)
    pyplot.legend(loc='lower right')
    pyplot.savefig('rollwave_solvers.png')

    # Centreline depth in space and time over the last 30 s (stored every
    # 0.5 s), one panel per solver
    fig, axes = pyplot.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True,
                                constrained_layout=True)
    for ax, alg in zip(axes.flat, ALGS):
        t, x, profiles = results[alg]
        fine = t >= t[-1] - 30.0 - 1.0e-6
        im = ax.pcolormesh(x, t[fine], profiles[fine]/dana, cmap='RdBu_r',
                           vmin=0.0, vmax=2.0, shading='nearest')
        ax.set_title(alg)
    for ax in axes[-1]:
        ax.set_xlabel('Xposition m')
    for ax in axes[:, 0]:
        ax.set_ylabel('Time s')
    fig.colorbar(im, ax=axes, label='depth / steady uniform depth', shrink=0.8)
    fig.savefig('rollwave_solvers_xt.png')
    pyplot.close(fig)


if __name__ == '__main__':
    args = anuga.get_args()
    if not getattr(args, 'no_run', False):
        run_solvers(args)
    plot_solvers()
