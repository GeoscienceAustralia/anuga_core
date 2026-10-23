"""Settling of suspended sediment in still water: the reference solution.

A tank of still water of depth h0 carries a uniform volumetric concentration
c0 of a single sediment fraction. There is no flow, so the bed shear stress
is zero and nothing is entrained; the only exchange with the bed is
deposition, [D-1] with a constant near-bed profile factor d*:

    D = d* v_s c                                            [D-1]

The conserved variable is the mass per unit area m = h c, so [G-3] with
E = 0 and no advection is

    dm/dt = -D = -d* v_s c ,

and the bed rises by the deposited volume divided by the packing fraction,
[G-4]:

    dz/dt = D / (1 - lambda) .

In still water the free surface does not move, so h = w - z falls as the bed
rises: h(t) = h0 - (z(t) - z0). The system is therefore

    dc/dt = -d* v_s c / h  +  c (dz/dt) / h          (from m = h c)
    dz/dt =  d* v_s c / (1 - lambda)

which has no closed form once the bed feedback is kept. It is integrated
here to machine precision with scipy; that is the reference the numerical
model is compared against. Dropping the feedback (h = h0) gives the familiar
exponential, also provided, which the reference approaches as c0 -> 0:

    c(t) = c0 exp(-d* v_s t / h0)
    z(t) = z0 + h0 c0 (1 - exp(-d* v_s t / h0)) / (1 - lambda)

Sediment mass is conserved between the water column and the bed:

    h c  +  (1 - lambda) (z - z0)  =  h0 c0        for all t.
"""
import numpy as np
from scipy.integrate import solve_ivp


def settling(t, h0, c0, v_s, d_star=1.0, porosity=0.30, z0=0.0):
    """Reference c(t), z(t), h(t) at the times t, with bed feedback.

    Returns three arrays. The concentration is the depth-averaged volumetric
    concentration; z is the bed elevation.
    """
    t = np.asarray(t, dtype=float)
    k = d_star * v_s

    def rhs(_, y):
        m, z = y                        # m = h c is the conserved variable
        h = h0 - (z - z0)
        c = m / h
        D = k * c
        return [-D, D / (1.0 - porosity)]

    sol = solve_ivp(rhs, (0.0, float(t.max())), [h0 * c0, z0], t_eval=t,
                    rtol=1e-12, atol=1e-15, method='DOP853')
    m, z = sol.y
    h = h0 - (z - z0)
    return m / h, z, h


def settling_no_feedback(t, h0, c0, v_s, d_star=1.0, porosity=0.30, z0=0.0):
    """The exponential solution with the depth held at h0."""
    t = np.asarray(t, dtype=float)
    decay = np.exp(-d_star * v_s * t / h0)
    c = c0 * decay
    z = z0 + h0 * c0 * (1.0 - decay) / (1.0 - porosity)
    return c, z


def timescale(h0, v_s, d_star=1.0):
    """The e-folding time h0 / (d* v_s) of the concentration."""
    return h0 / (d_star * v_s)


def read_centroid_series(filename):
    """time, concentration, bed and stage at the centroids, all (ntime, ncells).

    Read straight from the sww rather than through plot_utils.get_centroids,
    which collapses the (time-varying) bed to one slice.
    """
    from netCDF4 import Dataset
    with Dataset(filename) as f:
        t = f.variables['time'][:].astype(float)
        c = f.variables['silt_c'][:].astype(float)
        z = f.variables['elevation_c'][:].astype(float)
        w = f.variables['stage_c'][:].astype(float)
    return t, c, z, w
