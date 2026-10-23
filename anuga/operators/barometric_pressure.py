"""Barometric pressure operator — fractional-step operator for pressure-gradient momentum forcing."""

import numpy as np
from anuga.config import rho_w
from anuga.operators.base_operator import Operator
from anuga.utilities.function_utils import evaluate_file_function_all_points
from anuga.utilities.numerical_tools import gradient


class Barometric_pressure_operator(Operator):
    """Apply barometric pressure gradient forcing to water momentum each timestep.

    For each triangle the horizontal pressure gradient ∇p is computed from
    pressure values at the three vertices and the momentum is updated by::

        d(xmom)/dt += h * dp/dx / rho_w
        d(ymom)/dt += h * dp/dy / rho_w

    where *h* is the local water depth.

    Parameters
    ----------
    domain : anuga.Domain
    pressure : float or callable
        Barometric pressure in Pa.  If callable, signature must be
        ``pressure(t, x, y)`` where *x*, *y* are the global node (vertex)
        coordinate arrays, and the return value must be broadcastable to
        shape ``(number_of_nodes,)``.
    use_coordinates : bool, optional
        ``True`` (default): a callable *pressure* takes ``(t, x, y)``.
        ``False``: *pressure* is a :func:`anuga.file_function` with the single
        quantity ``barometric_pressure`` precomputed at the mesh nodes
        (``interpolation_points=domain.get_nodes()``); see the example below.
    description, label, logging, verbose : passed to :class:`Operator`.

    Examples
    --------
    Uniform pressure::

        P = Barometric_pressure_operator(domain, pressure=101325.0)

    Spatially and temporally varying pressure function::

        def storm_pressure(t, x, y):
            r2 = (x - x0)**2 + (y - y0)**2
            return (p_max - (p_max - p_min) * np.exp(-r2 / R**2)).reshape(1, -1)

        P = Barometric_pressure_operator(domain, storm_pressure, use_coordinates=True)

    Pressure field read from a file::

        F = anuga.file_function('pressure.sww', domain,
                                quantities=['barometric_pressure'],
                                interpolation_points=domain.get_nodes())
        P = Barometric_pressure_operator(domain, F, use_coordinates=False)
    """

    def __init__(self, domain, pressure=101325.0, use_coordinates=True,
                 description=None, label=None, logging=False, verbose=False):

        Operator.__init__(self, domain,
                          description=description, label=label,
                          logging=logging, verbose=verbose)

        self.use_coordinates = bool(use_coordinates)
        if not self.use_coordinates:
            names = getattr(pressure, 'quantity_names', None)
            if names is None or len(names) != 1:
                raise ValueError(
                    'Barometric_pressure_operator(use_coordinates=False) '
                    'expects a file_function with exactly one quantity; '
                    'got %r' % (names,))
            pts = getattr(pressure, 'interpolation_points', None)
            if pts is not None and len(pts) != domain.get_number_of_nodes():
                raise ValueError(
                    'Barometric_pressure_operator(use_coordinates=False): the '
                    'file_function has %d interpolation points but the mesh '
                    'has %d nodes; interpolate at domain.get_nodes()'
                    % (len(pts), domain.get_number_of_nodes()))

        self.pressure = pressure

    def __call__(self):
        domain = self.domain
        t = domain.get_time()
        dt = domain.get_timestep()
        N = domain.number_of_elements

        # Evaluate pressure at each mesh node (unique vertices)
        if not self.use_coordinates:
            field = evaluate_file_function_all_points(self.pressure, t)
            if field.ndim == 1:      # time-only file: uniform over the mesh
                p_nodes = np.full(domain.get_number_of_nodes(), field[0])
            else:
                p_nodes = field[0]
        elif callable(self.pressure):
            node_coords = domain.get_nodes()
            p_nodes = np.asarray(
                self.pressure(t, node_coords[:, 0], node_coords[:, 1]),
                dtype=float).ravel()
        else:
            p_nodes = np.full(domain.get_number_of_nodes(), float(self.pressure))

        # Per-triangle vertex coordinates (3N × 2) and node indices (N × 3)
        xv = domain.get_vertex_coordinates()
        triangles = domain.triangles
        height = self.stage_c - self.elev_c

        for k in range(N):
            p0 = p_nodes[triangles[k, 0]]
            p1 = p_nodes[triangles[k, 1]]
            p2 = p_nodes[triangles[k, 2]]

            k3 = 3 * k
            x0, y0 = xv[k3]
            x1, y1 = xv[k3 + 1]
            x2, y2 = xv[k3 + 2]

            px, py = gradient(x0, y0, x1, y1, x2, y2, p0, p1, p2)

            self.xmom_c[k] += dt * height[k] * px / rho_w
            self.ymom_c[k] += dt * height[k] * py / rho_w

    def parallel_safe(self):
        return True
