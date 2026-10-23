"""A migrating bed hump under Grass bedload: the reference solution.

The classic Exner test (Hudson & Sweby 2003; Castro Diaz, Fernandez-Nieto
& Ferreiro 2008): a channel of still-shaped, slowly moving bed carrying a
uniform subcritical flow of depth about 10 m and unit velocity, with a
low sinusoidal hump on the bed. Bedload follows Grass's law

    q_b = A_g u^m,   m = 3                                  [K-6]

and the bed moves by its divergence (Exner),

    (1 - lambda) dz/dt + dq_b/dx = 0 .                      [K-3]

The hump is low against the depth and the Froude number small, so over
the hump the free surface w and the discharge per unit width q barely
change: with w and q held constant, u = q / (w - z) and Exner becomes a
scalar conservation law for z alone,

    dz/dt + c(z) dz/dx = 0,    c(z) = m A_g q^m / ((1 - lambda) (w - z)^(m+1)),

whose solution is carried along characteristics: a point of the bed at
elevation z moves downstream at the speed c(z), faster the higher it is,

    z(x0 + c(z0(x0)) t, t) = z0(x0) .

So the hump migrates and its front steepens, until the crest overtakes the
toe and a shock forms at

    t_shock = 1 / max(-c'(z0) z0'(x0)) ;

the reference here is the smooth solution before that time, evaluated by
mapping a fine grid of x0 forward and interpolating back to the cells.

The frozen free surface is the only approximation: the surface dips by
about (u^2/g) dz/h over the hump, a fraction 1e-3 of the depth here, and
enters c through (w - z)^4, so the reference is good to about 0.5%. The
total bed volume is exact: the flux in at the inflow equals the flux out
at the outflow, both over the same flat bed.
"""
import numpy as np


def initial_bed(x, z_flat=0.1, x0=300.0, width=200.0, height=1.0):
    """z0(x): a sin^2 hump of the given height on a flat bed."""
    x = np.asarray(x, dtype=float)
    z = np.full(x.shape, z_flat)
    inside = (x >= x0) & (x <= x0 + width)
    z[inside] = z_flat + height * np.sin(np.pi * (x[inside] - x0) / width) ** 2
    return z


def wave_speed(z, w, q, A_g, m=3.0, porosity=0.4):
    """c(z) = m A_g q^m / ((1 - lambda) (w - z)^(m + 1))."""
    return m * A_g * q ** m / ((1.0 - porosity) * (w - np.asarray(z, dtype=float)) ** (m + 1.0))


def shock_time(w, q, A_g, m=3.0, porosity=0.4, **hump):
    """The time the smooth solution breaks on the hump's front."""
    x0 = np.linspace(0.0, 1000.0, 200001)
    z0 = initial_bed(x0, **hump)
    c = wave_speed(z0, w, q, A_g, m, porosity)
    dcdx = np.gradient(c, x0)
    return 1.0 / max(-dcdx.min(), 1e-300)


def bed(x, t, w, q, A_g, m=3.0, porosity=0.4, **hump):
    """The reference bed z(x, t) on the points x, by characteristics."""
    x0 = np.linspace(-200.0, 1400.0, 320001)
    z0 = initial_bed(x0, **hump)
    xt = x0 + wave_speed(z0, w, q, A_g, m, porosity) * t
    if np.any(np.diff(xt) <= 0.0):
        raise ValueError('the characteristics have crossed: t is past the shock')
    return np.interp(np.asarray(x, dtype=float), xt, z0)


def read_centroid_series(filename):
    """time, bed, stage, xmom at the centroids (ntime, ncells), and the
    centroid coordinates (ncells, 2)."""
    from netCDF4 import Dataset
    with Dataset(filename) as f:
        t = f.variables['time'][:].astype(float)
        z = f.variables['elevation_c'][:].astype(float)
        w = f.variables['stage_c'][:].astype(float)
        uh = f.variables['xmomentum_c'][:].astype(float)
        x = f.variables['x'][:].astype(float)
        y = f.variables['y'][:].astype(float)
        vols = f.variables['volumes'][:]
    if z.ndim == 1:
        z = np.broadcast_to(z, w.shape).copy()
    xc = x[vols].mean(axis=1)
    yc = y[vols].mean(axis=1)
    return t, z, w, uh, np.column_stack([xc, yc])
