import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import setup_diagrams          # noqa: E402  the shared flume schematics

setup_diagrams.pier_dambreak('setup.png', 'P1')

if os.path.exists('pier_P1.sww'):
    import numpy as _np
    from netCDF4 import Dataset as _nc
    with _nc('pier_P1.sww') as _f:
        _v = _f.variables['volumes'][:]
        _vx = _f.variables['x'][:].astype(float)
        _vy = _f.variables['y'][:].astype(float)
    _cx = _vx[_v].mean(1)
    _cy = _vy[_v].mean(1)
    _side = _np.hypot(_vx[_v][:, 1] - _vx[_v][:, 0], _vy[_v][:, 1] - _vy[_v][:, 0])
    _zones = [('reservoir', _cx < 0.0),
              ('sand reach', (_cx > 1.0) & (_cx < 2.5)),
              ('at the pier', _np.hypot(_cx - 1.6, _cy - 0.12) < 0.06)]
    _counts = [(n, int(m.sum()), int(round(1000 * _side[m].mean()))) for n, m in _zones]
    setup_diagrams.mesh_figure(
        'pier_P1.sww', 'mesh.png', zoom=(1.4, 1.85, 0.0, 0.24),
        title='Case P1: the mesh, %d triangles' % len(_v),
        zoom_title='detail at the pier, which is meshed as a hole',
        counts=_counts)

"""Figures for the pier dam-break case: bed change maps (Kinect and ANUGA),
centreline bed profiles and the water surface along the centreline."""
import glob
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import anuga
from validate_pier import load_bed, load_wse, model_dz

for npz in sorted(glob.glob('pier_P?_bed.npz')):
    case = npz[5:7]
    r = np.load(npz)
    nev = r['beds'].shape[0]
    x, y, meas = load_bed(case)
    for event in range(1, nev + 1):
        f, _ = model_dz(case, event)
        mod = f(x, y)
        vmax = 25.0
        fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True, sharey=True)
        for a, v, ttl in [(ax[0], meas[:, event - 1], 'Kinect'), (ax[1], mod, 'ANUGA')]:
            sc = a.scatter(x, y, c=v, s=3, cmap='seismic', vmin=-vmax, vmax=vmax, marker='s')
            a.set_aspect('equal'); a.set_title('%s, bed change after dam-break %d: %s (mm)' % (case, event, ttl))
            a.set_ylabel('y (m)')
        ax[1].set_xlabel('x from the gate (m)')
        fig.colorbar(sc, ax=ax, shrink=0.8, label='bed change (mm)')
        fig.savefig('bed_change_%s_event%d.png' % (case, event), dpi=130)
        plt.close(fig)
    # centreline profiles after each event
    fig, ax = plt.subplots(figsize=(11, 4))
    band = np.abs(y - 0.12) < 0.0126
    order = np.argsort(x[band])
    for event in range(1, nev + 1):
        f, _ = model_dz(case, event)
        ax.plot(x[band][order], meas[band, event - 1][order], 'x', ms=3, color='C%d' % (event - 1),
                label='Kinect, after dam-break %d' % event)
        ax.plot(x[band][order], f(x[band], y[band])[order], '-', color='C%d' % (event - 1),
                label='ANUGA, after dam-break %d' % event)
    ax.axvspan(1.6 - 0.015, 1.6 + 0.015, color='0.8', label='pier')
    ax.set_xlabel('x from the gate (m)'); ax.set_ylabel('bed change (mm)'); ax.grid(alpha=0.3)
    ax.set_title('%s: bed change along the centreline (y = 0.12 +- 0.0125 m)' % case)
    ax.legend(fontsize='small', ncol=2)
    fig.tight_layout(); fig.savefig('centreline_%s.png' % case, dpi=130); plt.close(fig)
    # water surface along the centreline during the first dam-break
    if os.path.exists('kinect_%s_wse.csv' % case) and os.path.exists('pier_%s.sww' % case):
        from netCDF4 import Dataset
        t, xw, w = load_wse(case)
        with Dataset('pier_%s.sww' % case) as s:
            ts = s.variables['time'][:]
            xs = s.variables['x'][:]; ys = s.variables['y'][:]
            vols = s.variables['volumes'][:]
            xc = xs[vols].mean(axis=1); yc = ys[vols].mean(axis=1)
            sel = (np.abs(yc - 0.12) < 0.01) & (xc > 1.1) & (xc < 2.1)
            stage = s.variables['stage_c'][:, sel]
            zc = s.variables['elevation_c'][:, sel] if s.variables['elevation_c'].ndim == 2 else s.variables['elevation_c'][sel]
        fig, ax = plt.subplots(figsize=(11, 4))
        times = np.unique(t)
        pick = times[np.linspace(0, len(times) - 1, 5).astype(int)]
        for k, tp in enumerate(pick):
            m = t == tp
            ax.plot(xw[m], w[m], 'x', ms=3, color='C%d' % k, label='Kinect t = %.2f s' % tp)
            i = np.argmin(np.abs(ts - tp))
            o = np.argsort(xc[sel])
            ax.plot(xc[sel][o], 1000 * stage[i][o], '-', color='C%d' % k, label='ANUGA t = %.2f s' % ts[i])
        ax.set_xlabel('x from the gate (m)'); ax.set_ylabel('water surface above the initial bed (mm)')
        ax.set_title('%s: water surface along the centreline, first dam-break' % case); ax.grid(alpha=0.3)
        ax.legend(fontsize='x-small', ncol=5); ax.set_ylim(-5, 60)
        fig.tight_layout(); fig.savefig('wse_%s.png' % case, dpi=130); plt.close(fig)
