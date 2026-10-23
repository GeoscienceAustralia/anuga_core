"""Reduce the Zenodo Kinect archives (Segovia-Burillo et al. 2026, DOI
10.5281/zenodo.17777387) to the small CSV files the case compares against.

Run once, with the archives unpacked under $ZARAGOZA_KINECT_DATA (default
~/Projects/Sediment_Transport/Data/Zaragoza_Kinect_dambreak), as
    python extract_kinect.py P1 P2

Writes, per case:
  kinect_<case>_bed.csv  x_m, y_m, dz1_mm, dz2_mm, dz3_mm: the bed change
                         after each of the three dam-breaks relative to the
                         initial flat sand (sand0), block-averaged to 5 mm
                         inside the measurement mask;
  kinect_<case>_wse.csv  t_s, x_m, wse_mm: the water surface above the
                         initial bed along the flume centreline (y = 0.12 m,
                         averaged over +-1 cm) during the first dam-break,
                         at 10 mm spacing and every 4th frame (~0.13 s).
The Kinect grid is (X, Y) in mm, X along the wall from the gate, Y from the
right wall; its window covers x = 1.12-2.07 m, i.e. the pier reach.
"""
import glob
import os
import sys

import h5py
import numpy as np

ROOT = os.environ.get('ZARAGOZA_KINECT_DATA',
                      os.path.expanduser('~/Projects/Sediment_Transport/Data/Zaragoza_Kinect_dambreak'))
BLOCK = 0.005          # m, bed output spacing
DX_WSE = 0.010         # m, centreline output spacing
Y_CENTRE = 0.12
FRAME_STEP = 4


def load(path, key_prefix):
    with h5py.File(path) as h:
        key = [k for k in h.keys() if k.startswith(key_prefix)][0]
        return h['gX'][()] / 1000.0, h['gY'][()] / 1000.0, h[key][()].astype(float)


def block_mean(x, y, v, good, step):
    """Mean of v over step x step blocks; blocks with < 25% good points dropped."""
    ix = np.floor(x / step).astype(int)
    iy = np.floor(y / step).astype(int)
    keys = ix * 100000 + iy
    out = []
    for k in np.unique(keys[good]):
        sel = (keys == k)
        g = sel & good
        if g.sum() < 0.25 * sel.sum():
            continue
        out.append((x[g].mean(), y[g].mean(), v[g].mean()))
    return np.array(out)


def extract(case):
    d = os.path.join(ROOT, case)
    with h5py.File(os.path.join(d, 'offset.h5')) as h:
        mask = h['mask'][()] > 0.5
    X, Y, s0 = load(os.path.join(d, 'sand0.h5'), 'Dsand')
    cols = {}
    for n in (1, 2, 3):
        _, _, sn = load(os.path.join(d, 'sand%d.h5' % n), 'Dsand')
        cols[n] = sn - s0
    good = mask & np.isfinite(s0)
    for n in cols:
        good &= np.isfinite(cols[n])
    x, y = X.ravel(), Y.ravel()
    g = good.ravel()
    rows = None
    for n in (1, 2, 3):
        b = block_mean(x, y, cols[n].ravel(), g, BLOCK)
        rows = b if rows is None else np.column_stack([rows, b[:, 2]])
    with open('kinect_%s_bed.csv' % case, 'w') as f:
        f.write('# Segovia-Burillo et al. (2026) Kinect data, case %s: bed change after each\n'
                '# dam-break relative to the initial flat sand, 5 mm block means inside the\n'
                '# measurement mask. x from the gate, y from the right wall (m).\n'
                '# x_m, y_m, dz1_mm, dz2_mm, dz3_mm\n' % case)
        for r in rows:
            f.write('%.4f, %.4f, %.2f, %.2f, %.2f\n' % tuple(r))
    print(case, 'bed:', len(rows), 'points')
    # water surface, event 1, centreline
    files = sorted(glob.glob(os.path.join(d, 'water1_h5', 'wse_t*.h5')))
    band = np.abs(Y - Y_CENTRE) <= 0.01
    xs = np.arange(np.floor(X.min() / DX_WSE) * DX_WSE, X.max(), DX_WSE)
    with open('kinect_%s_wse.csv' % case, 'w') as f:
        f.write('# Segovia-Burillo et al. (2026) Kinect data, case %s: water surface above\n'
                '# the initial bed along the flume centreline (y = 0.12 m, +-1 cm) during the\n'
                '# first dam-break; t from the gate opening.\n# t_s, x_m, wse_mm\n' % case)
        nfr = 0
        for fn in files[::FRAME_STEP]:
            t = float(os.path.basename(fn)[5:-3].replace('_', '.'))
            with h5py.File(fn) as h:
                w = h['wse[mm]'][()].astype(float)
            ok = band & mask & np.isfinite(w)
            for xa in xs:
                sel = ok & (np.abs(X - xa) <= DX_WSE / 2)
                if sel.sum() >= 5:
                    f.write('%.3f, %.3f, %.2f\n' % (t, xa, w[sel].mean()))
            nfr += 1
    print(case, 'wse:', nfr, 'frames')


if __name__ == '__main__':
    for case in sys.argv[1:] or ('P1', 'P2'):
        extract(case)
