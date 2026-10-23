"""Shared readers for the van Rijn flume cases."""
import numpy as np


def read_centroid_series(filename, name='sand'):
    """time, concentration, bed, stage, xmom at the centroids, (ntime, ncells),
    and the centroid coordinates (ncells, 2)."""
    from netCDF4 import Dataset
    with Dataset(filename) as f:
        t = f.variables['time'][:].astype(float)
        c = f.variables[name + '_c'][:].astype(float)
        z = f.variables['elevation_c'][:].astype(float)
        w = f.variables['stage_c'][:].astype(float)
        uh = f.variables['xmomentum_c'][:].astype(float)
        x = f.variables['x'][:].astype(float)
        y = f.variables['y'][:].astype(float)
        vols = f.variables['volumes'][:]
    if z.ndim == 1:
        z = np.broadcast_to(z, c.shape).copy()
    xc = x[vols].mean(axis=1)
    yc = y[vols].mean(axis=1)
    return t, c, z, w, uh, np.column_stack([xc, yc])


def read_measured(filename):
    """Two-column CSV with '#' comments -> (col0, col1) arrays."""
    d = np.loadtxt(filename, delimiter=',', comments='#')
    return d[:, 0], d[:, 1]


def along_flume(x, values, bins):
    """Average `values` over cells binned by x; returns bin centres, means."""
    idx = np.digitize(x, bins) - 1
    xm, vm = [], []
    for k in range(len(bins) - 1):
        sel = idx == k
        if sel.any():
            xm.append(x[sel].mean())
            vm.append(values[sel].mean())
    return np.array(xm), np.array(vm)
