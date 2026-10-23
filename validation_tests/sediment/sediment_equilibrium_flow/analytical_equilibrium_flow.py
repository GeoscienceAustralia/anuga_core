"""Equilibrium suspended load in normal flow: the reference solution.

Water flows down a plane channel of slope S with Manning roughness n at
its normal depth h, so that gravity and friction balance and the flow is
uniform: from Manning's relation

    u = h^(2/3) S^(1/2) / n,       q = u h .                 [F-1]

The quadratic-drag shear closure [T-1] with the Manning friction factor
[T-6] gives

    tau_b / rho = f_c u^2,   f_c = g n^2 h^(-1/3)             [T-1], [T-6]

and in normal flow those two combine to the depth-slope product exactly,

    tau_b / rho = g n^2 h^(-1/3) h^(4/3) S / n^2 = g h S ,

so the Shields stress of a grain of diameter d and submerged specific
gravity R is

    tau* = h S / (R d) .                                     [T-3]

That fixes the non-cohesive entrainment rate

    E = v_s E*,   E* = 0.65 gamma0 X / (1 + gamma0 X),   X = tau*/tau_c* - 1
                                                           [E-1], [E-2]

and, against the well-mixed deposition D = d* v_s c [D-1], the equilibrium
concentration at which the bed exchange vanishes:

    c_eq = E / (d* v_s) = E* / d* .

Clear water enters at the inflow. Along the channel the depth-averaged
concentration obeys the steady advection balance per unit width

    d(q c)/dx = E - d* v_s c ,                               [G-3]

so with q constant it relaxes to c_eq over the settling length,

    c(x) = c_eq (1 - exp(-x / L_s)),      L_s = q / (d* v_s) .

The bed is held fixed (an erodible bed of unlimited supply that is not
lowered): the net exchange E - D is largest at the inflow and would
otherwise scour the reach in a few minutes and take the depth, the stress
and the reference with it. Under a fixed bed the flow and the reference are
exact for the whole run.
"""
import numpy as np


def normal_velocity(h, S, n):
    """Manning: u = h^(2/3) S^(1/2) / n, in m/s."""
    return h ** (2.0 / 3.0) * np.sqrt(S) / n


def shields_stress(h, S, R, d):
    """tau* = h S / (R d) in normal flow under quadratic drag."""
    return h * S / (R * d)


def entrainment_rate(h, S, v_s, R, d, tau_c_star, gamma0):
    """E in m/s, [T-3], [E-1], [E-2]."""
    X = shields_stress(h, S, R, d) / tau_c_star - 1.0
    gX = gamma0 * max(X, 0.0)
    return v_s * 0.65 * gX / (1.0 + gX)


def equilibrium_concentration(h, S, v_s, R, d, tau_c_star, gamma0, d_star=1.0):
    """c_eq = E / (d* v_s), the concentration at which the bed exchange vanishes."""
    return entrainment_rate(h, S, v_s, R, d, tau_c_star, gamma0) / (d_star * v_s)


def settling_length(q, v_s, d_star=1.0):
    """L_s = q / (d* v_s), in metres."""
    return q / (d_star * v_s)


def concentration(x, c_eq, q, v_s, d_star=1.0):
    """Steady c(x) = c_eq (1 - exp(-x / L_s)) from clear water at x = 0."""
    x = np.asarray(x, dtype=float)
    return c_eq * (1.0 - np.exp(-x / settling_length(q, v_s, d_star)))


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
