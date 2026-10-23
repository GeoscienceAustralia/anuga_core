"""A settling basin with flow: the reference solution.

Sediment-laden water enters a straight channel of uniform steady flow and
its suspended load settles out along the way. The flow is uniform (flat
frictionless bed, depth h, discharge per unit width q = u h fixed by the
inflow boundary), the sediment is a single fraction of settling velocity
v_s with entrainment switched off (critical Shields stress zero), and
deposition is the well-mixed rate

    D = d* v_s c .                                          [D-1]

At steady state the depth-averaged concentration obeys the advection
balance per unit width

    d(q c)/dx = -D = -d* v_s c ,                            [G-3]

so with q constant along the channel

    c(x) = c0 exp(-x / L_s),      L_s = q / (d* v_s) ,

the settling length. The bed rises where the sediment lands, by the
deposited volume over the packing fraction [G-4]:

    dz/dt = D / (1 - lambda) = d* v_s c(x) / (1 - lambda) ,

steady in time once the concentration is, so the bed-rise RATE over a late
interval is the check on the bed update. As the bed rises the depth falls
and the velocity rises to keep q, so L_s does not change; the reference
holds for the whole run.

Sediment mass is conserved: what crossed the inflow boundary, less what
left through the outflow, is either still in suspension or in the bed:

    inflow - outflow = water column + (1 - lambda) * bed volume .
"""
import numpy as np


def settling_length(q, v_s, d_star=1.0):
    """L_s = q / (d* v_s), in metres."""
    return q / (d_star * v_s)


def concentration(x, c0, q, v_s, d_star=1.0):
    """Steady c(x) = c0 exp(-x / L_s), x measured from the inflow."""
    x = np.asarray(x, dtype=float)
    return c0 * np.exp(-x / settling_length(q, v_s, d_star))


def bed_rise_rate(x, c0, q, v_s, d_star=1.0, porosity=0.30):
    """Steady dz/dt = d* v_s c(x) / (1 - lambda), in m/s."""
    return d_star * v_s * concentration(x, c0, q, v_s, d_star) / (1.0 - porosity)


def read_centroid_series(filename, name='silt'):
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
