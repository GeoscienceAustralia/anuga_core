"""Wind stress operator — fractional-step operator for wind-driven momentum forcing."""

import numpy as np
from math import pi
from anuga.config import rho_a, rho_w, eta_w
from anuga.operators.base_operator import Operator
from anuga.utilities.function_utils import evaluate_file_function_all_points


class Wind_stress_operator(Operator):
    """Apply wind stress to water momentum each timestep.

    Wind speed *s* [m/s] and direction *phi* [degrees] drive a surface stress
    that is added to xmomentum and ymomentum.  Both can be scalars or
    functions ``f(t, x, y)`` returning arrays aligned with the centroid
    coordinates.

    The stress formula follows Large & Pond (1981) via::

        S = (eta_w * rho_a / rho_w) * |wind|
        d(xmom)/dt += S * u_wind
        d(ymom)/dt += S * v_wind

    Parameters
    ----------
    domain : anuga.Domain
    speed : float or callable
        Wind speed in m/s.  If callable, signature must be
        ``speed(t, x, y)`` where *x*, *y* are centroid coordinate arrays.
    phi : float or callable
        Wind direction in degrees (standard mathematical convention:
        0° = east, 90° = north).  If callable, same signature as *speed*.
    use_coordinates : bool, optional
        ``True`` (default): callables take ``(t, x, y)``.  ``False``: *speed*
        is a :func:`anuga.file_function` whose quantities are
        ``[speed, angle]``, precomputed at the mesh centroids, and *phi* is
        ignored (see the file example below).
    description, label, logging, verbose : passed to :class:`Operator`.

    Examples
    --------
    Uniform wind::

        W = Wind_stress_operator(domain, speed=10.0, phi=180.0)

    Spatially varying wind functions::

        def my_speed(t, x, y):
            return 20.0 * np.ones_like(x)

        def my_phi(t, x, y):
            return 270.0 * np.ones_like(x)

        W = Wind_stress_operator(domain, my_speed, my_phi)

    Wind field read from a file (``use_coordinates=False``)::

        F = anuga.file_function('wind.sww', domain,
                                quantities=['wind_speed', 'wind_angle'],
                                interpolation_points=domain.get_centroid_coordinates())
        W = Wind_stress_operator(domain, F, use_coordinates=False)

    Here *speed* is the ``file_function`` object itself, its two quantities
    are taken as (speed, angle) in that order, and *phi* is ignored.  The
    interpolation points must be the centroid coordinates.  A time-only
    ``.tms`` file (two attribute columns) works the same way and applies a
    spatially uniform, time-varying wind.
    """

    def __init__(self, domain, speed=0.0, phi=0.0, use_coordinates=True,
                 description=None, label=None, logging=False, verbose=False):

        Operator.__init__(self, domain,
                          description=description, label=label,
                          logging=logging, verbose=verbose)

        self.use_coordinates = bool(use_coordinates)
        if not self.use_coordinates:
            names = getattr(speed, 'quantity_names', None)
            if names is None or len(names) != 2:
                raise ValueError(
                    'Wind_stress_operator(use_coordinates=False) expects a '
                    'file_function with exactly two quantities (speed, angle); '
                    'got %r' % (names,))
            pts = getattr(speed, 'interpolation_points', None)
            if pts is not None and len(pts) != domain.number_of_elements:
                raise ValueError(
                    'Wind_stress_operator(use_coordinates=False): the '
                    'file_function has %d interpolation points but the mesh '
                    'has %d triangles; interpolate at the centroids'
                    % (len(pts), domain.number_of_elements))
        self.speed = speed
        self.phi = phi
        self.const = eta_w * rho_a / rho_w

    def __call__(self):
        domain = self.domain
        t = domain.get_time()
        dt = domain.get_timestep()
        xc = domain.centroid_coordinates
        N = domain.number_of_elements

        if not self.use_coordinates:
            field = evaluate_file_function_all_points(self.speed, t)
            if field.ndim == 1:      # time-only file: uniform over the mesh
                s_vec = np.full(N, field[0])
                phi_vec = np.full(N, field[1])
            else:
                s_vec = field[0]
                phi_vec = field[1]
        elif callable(self.speed):
            s_vec = np.asarray(
                self.speed(t, xc[:, 0], xc[:, 1]), dtype=float).ravel()
        else:
            s_vec = np.full(N, float(self.speed))

        if not self.use_coordinates:
            pass
        elif callable(self.phi):
            phi_vec = np.asarray(
                self.phi(t, xc[:, 0], xc[:, 1]), dtype=float).ravel()
        else:
            phi_vec = np.full(N, float(self.phi))

        phi_rad = phi_vec * (pi / 180.0)
        u = s_vec * np.cos(phi_rad)
        v = s_vec * np.sin(phi_rad)
        S = self.const * np.sqrt(u**2 + v**2)

        self.xmom_c += dt * S * u
        self.ymom_c += dt * S * v

    def parallel_safe(self):
        return True
