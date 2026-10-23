"""Entrainment of bed sediment into still water: the reference solution.

A tank of still water stands over a plane bed of gentle slope S. There is
no flow, so the quadratic-drag closure would give zero shear; the run
instead selects the depth-slope closure [T-7], with the slope FROZEN at the
setup bed, under which the bed shear stress of a cell is a function of its
depth alone:

    tau_b / rho = g h S .                                   [T-7]

The Shields stress

    tau* = (tau_b / rho) / (R g d) = h S / (R d)             [T-3]

exceeds the critical value tau_c* of the sand fraction, and the
non-cohesive entrainment law fires,

    E = v_s E*,   E* = 0.65 gamma0 X / (1 + gamma0 X),   X = tau*/tau_c* - 1
                                                           [E-1], [E-2]

while deposition returns sediment to the bed at the well-mixed rate

    D = d* v_s c .                                          [D-1]

A second fraction of the same grain but a critical stress far above tau*
has X < 0 and is never entrained: that is the threshold control.

With no flow there is no advection, so every cell is a closed system. The
sediment mass per unit area m = h c and the bed elevation z obey

    dm/dt = E(h) - D                                        [G-3]
    dz/dt = -(E(h) - D) / (1 - lambda)                      [G-4]

with the free surface w fixed, h = w - z: as the bed lowers the depth
grows and E with it. That coupled system has no closed form and is
integrated to machine precision here, one ODE per cell; it is the
reference. With the depth held at h0 (a fixed bed) the concentration
relaxes exponentially to c_eq = E(h0)/(d* v_s) with time constant
h0/(d* v_s), also provided.

Sediment mass is conserved between the water column and the bed:

    h c  +  (1 - lambda) (z - z0)  =  0        for all t   (c0 = 0).

Why the slope is frozen. The depth-slope closure normally reads S from the
evolving bed, and once cells erode at different rates the local slopes
steepen, tau_b rises, and erosion accelerates: the closure feeds on the
roughness it creates, and there is no closed-form reference (and, left to
run, the tank scours metres). Freezing S at the setup bed removes that
feedback, keeps S the reach slope the closure was written for, and lets
the bed evolve under it. That is the option a user would choose to run
the anugaSed closure on a moving bed, and it is what the case exercises.
"""
import numpy as np


def entrainment_rate(h, S, v_s, R, d, tau_c_star, gamma0):
    """E in m/s for depth h and bed slope S, [T-7], [T-3], [E-1], [E-2]."""
    h = np.asarray(h, dtype=float)
    S = np.asarray(S, dtype=float)
    tau_star = h * S / (R * d)
    X = tau_star / tau_c_star - 1.0
    gX = gamma0 * np.maximum(X, 0.0)
    return v_s * 0.65 * gX / (1.0 + gX)


def equilibrium_concentration(h, S, v_s, R, d, tau_c_star, gamma0, d_star=1.0):
    """c_eq = E / (d* v_s), the concentration at which deposition balances."""
    return entrainment_rate(h, S, v_s, R, d, tau_c_star, gamma0) / (d_star * v_s)


def timescale(h, v_s, d_star=1.0):
    """The e-folding time h / (d* v_s) of the relaxation."""
    return np.asarray(h, dtype=float) / (d_star * v_s)


def erosion(t, h, S, v_s, R, d, tau_c_star, gamma0, d_star=1.0):
    """Reference c(t) for cells of depth h and bed slope S, from c = 0.

    h and S are arrays with one entry per cell (scalars are broadcast); the
    result has shape (len(t), ncells).
    """
    t = np.asarray(t, dtype=float)
    h = np.atleast_1d(np.asarray(h, dtype=float))
    S = np.atleast_1d(np.asarray(S, dtype=float))
    c_eq = equilibrium_concentration(h, S, v_s, R, d, tau_c_star, gamma0, d_star)
    T = timescale(h, v_s, d_star)
    return c_eq[None, :] * (1.0 - np.exp(-t[:, None] / T[None, :]))


def erosion_with_feedback(t, h0, S, v_s, R, d, tau_c_star, gamma0,
                          d_star=1.0, porosity=0.30, z0=None):
    """Reference c(t), z(t), h(t) with the bed evolving, from c = 0.

    h0 and S are arrays with one entry per cell (scalars are broadcast);
    each result has shape (len(t), ncells). The slope S is constant per
    cell (frozen); the depth h = w - z grows as the bed lowers.
    """
    from scipy.integrate import solve_ivp
    t = np.asarray(t, dtype=float)
    h0 = np.atleast_1d(np.asarray(h0, dtype=float))
    n = len(h0)
    S = np.broadcast_to(np.asarray(S, dtype=float), (n,))
    z0 = np.zeros(n) if z0 is None else np.broadcast_to(np.asarray(z0, dtype=float), (n,))
    w = z0 + h0
    k = d_star * v_s

    def rhs(_, y):
        m = y[:n]
        z = y[n:]
        h = w - z
        net = entrainment_rate(h, S, v_s, R, d, tau_c_star, gamma0) - k * m / h
        return np.concatenate([net, -net / (1.0 - porosity)])

    sol = solve_ivp(rhs, (0.0, float(t.max())), np.concatenate([np.zeros(n), z0]),
                    t_eval=t, rtol=1e-12, atol=1e-15, method='DOP853')
    m = sol.y[:n].T
    z = sol.y[n:].T
    h = w[None, :] - z
    return m / h, z, h


def read_centroid_series(filename, name='sand'):
    """time, concentration, bed, stage and momenta at the centroids.

    All arrays are (ntime, ncells). A fixed bed is stored once in the sww;
    it is broadcast over time here so the caller need not care.
    """
    from netCDF4 import Dataset
    with Dataset(filename) as f:
        t = f.variables['time'][:].astype(float)
        c = f.variables[name + '_c'][:].astype(float)
        z = f.variables['elevation_c'][:].astype(float)
        w = f.variables['stage_c'][:].astype(float)
        uh = f.variables['xmomentum_c'][:].astype(float)
        vh = f.variables['ymomentum_c'][:].astype(float)
    if z.ndim == 1:
        z = np.broadcast_to(z, c.shape).copy()
    return t, c, z, w, uh, vh
