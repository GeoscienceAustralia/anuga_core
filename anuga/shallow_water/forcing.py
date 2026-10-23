"""Environmental forcing helpers.

The forcing-function classes that used to live here (``Wind_stress``,
``Rainfall``, ``Inflow``, ``Barometric_pressure`` and their ``_fast``
variants) were deprecated in 4.0 and removed in 4.1.  Use the operators
instead:

* ``anuga.Rate_operator.rainfall()`` / ``anuga.Rate_operator.inflow()``
* ``anuga.Wind_stress_operator``
* ``anuga.Barometric_pressure_operator``

See ``UPGRADING_TO_4.0.md``.  ``Cross_section`` (flow through a polyline)
remains here.
"""

import numpy as num
import anuga.utilities.log as log

from anuga.abstract_2d_finite_volumes.neighbour_mesh import segment_midpoints
from anuga.geospatial_data.geospatial_data import ensure_geospatial


class Cross_section:
    """Class Cross_section - a class to setup a cross section from
    which you can then calculate flow and energy through cross section

    Cross_section(domain, polyline)

    domain:
    polyline: Representation of desired cross section - it may contain
              multiple sections allowing for complex shapes. Assume
              absolute UTM coordinates.
              Format [[x0, y0], [x1, y1], ...]
    verbose:
    """

    def __init__(self,
                 domain,
                 polyline=None,
                 verbose=False):
        """Create an instance of Cross_section.

        domain    domain of interest
        polyline  polyline defining cross section
        verbose   True if this instance is to be verbose
        """

        self.domain = domain
        self.polyline = polyline
        self.verbose = verbose

        import numpy as np

        # Find all intersections and associated triangles.
        self.segments = self.domain.get_intersecting_segments(self.polyline,
                                                              use_cache=True,
                                                              verbose=self.verbose)

        # Get midpoints (kept for get_energy_through_cross_section which needs interpolation)
        self.midpoints = segment_midpoints(self.segments)

        # Make midpoints Geospatial instances
        self.midpoints = ensure_geospatial(self.midpoints, self.domain.geo_reference)

        # Pre-extract per-segment geometry as arrays for vectorised flow computation
        self._tri_ids = np.array([seg.triangle_id for seg in self.segments])
        self._normals = np.array([seg.normal for seg in self.segments])   # (S, 2)
        self._lengths = np.array([seg.length for seg in self.segments])   # (S,)

    def set_verbose(self,verbose=True):
        """Set verbose mode true or flase"""

        self.verbose=verbose

    def get_flow_through_cross_section(self):
        """ Output: Total flow [m^3/s] across cross section.
        """

        # Use centroid values at the known intersected triangle IDs directly —
        # avoids spatial interpolation since midpoints lie inside these triangles.
        uh = self.domain.quantities['xmomentum'].centroid_values[self._tri_ids]
        vh = self.domain.quantities['ymomentum'].centroid_values[self._tri_ids]

        normal_mom = uh * self._normals[:, 0] + vh * self._normals[:, 1]
        return float((normal_mom * self._lengths).sum())


    def get_energy_through_cross_section(self, kind='total'):
        """Obtain average energy head [m] across specified cross section.

        Output:
            E: Average energy [m] across given segments for all stored times.

        The average velocity is computed for each triangle intersected by
        the polyline and averaged weighted by segment lengths.

        The typical usage of this function would be to get average energy of
        flow in a channel, and the polyline would then be a cross section
        perpendicular to the flow.

        #FIXME (Ole) - need name for this energy reflecting that its dimension
        is [m].
        """

        from anuga.config import g, epsilon, velocity_protection as h0

        # Get interpolated values
        stage = self.domain.get_quantity('stage')
        elevation = self.domain.get_quantity('elevation')
        xmomentum = self.domain.get_quantity('xmomentum')
        ymomentum = self.domain.get_quantity('ymomentum')

        w = stage.get_values(interpolation_points=self.midpoints, use_cache=True)
        z = elevation.get_values(interpolation_points=self.midpoints, use_cache=True)
        uh = xmomentum.get_values(interpolation_points=self.midpoints,
                                  use_cache=True)
        vh = ymomentum.get_values(interpolation_points=self.midpoints,
                                  use_cache=True)
        h = w-z                # Depth

        # Compute total length of polyline for use with weighted averages
        total_line_length = 0.0
        for segment in self.segments:
            total_line_length += segment.length

        # Compute and sum flows across each segment
        average_energy = 0.0
        for i in range(len(w)):
            # Average velocity across this segment
            if h[i] > epsilon:
                # Use protection against degenerate velocities
                u = uh[i]/(h[i] + h0/h[i])
                v = vh[i]/(h[i] + h0/h[i])
            else:
                u = v = 0.0

            speed_squared = u*u + v*v
            kinetic_energy = 0.5*speed_squared/g

            if kind == 'specific':
                segment_energy = h[i] + kinetic_energy
            elif kind == 'total':
                segment_energy = w[i] + kinetic_energy
            else:
                msg = 'Energy kind must be either "specific" or "total".'
                msg += ' I got %s' %kind

            # Add to weighted average
            weigth = self.segments[i].length/total_line_length
            average_energy += segment_energy*weigth

        return average_energy
