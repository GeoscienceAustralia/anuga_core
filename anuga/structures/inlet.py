
import anuga.geometry.polygon
from anuga.geometry.polygon import inside_polygon, is_inside_polygon, line_intersect
from anuga.config import velocity_protection, g
from anuga import Region

import math


import numpy as num

def level_stages_to_average(stages, elevations, areas,
                            old_average_depth, new_average_depth):
    """Move an inlet to a new average depth by filling or drawing down to a level.

    Returns the new per-cell stages. The volume change is
    ``(new_average_depth - old_average_depth) * total_area``, applied with the
    "well-mixed reservoir" model the structure operators assume — water finds
    its level:

    * adding volume raises the LOWEST stages to a common level (a high cell is
      never pulled down);
    * removing volume lowers the HIGHEST stages to a common level (a low cell is
      never lifted), clamping each cell at its bed; if the inlet holds less
      water than asked for, it is drained dry and no more.

    This replaces writing a uniform DEPTH across the inlet, which on a sloping
    bed tilts the water surface to follow the bed and disturbs a lake at rest by
    half the bed elevation range across the inlet (issue #229). Leveling is well
    balanced — a level surface with no transfer is untouched — while keeping the
    smoothing property of the old write that the culvert feedback loop relies on
    (a shape-preserving shift was tried first and destabilised discharge into an
    initially dry inlet; see the issue). On a flat bed, filling a dry or level
    inlet reproduces the old uniform-depth result exactly.

    The level is found by bisection on the monotone volume-of-water-above /
    below-a-level function; 100 iterations pins it to the last bit of a double.

    Parameters
    ----------
    stages, elevations, areas : array
        Current per-cell stage, bed elevation and area.
    old_average_depth : float
        The average depth the caller's ``new_average_depth`` was derived from.
        In parallel this is the GLOBAL average over the whole inlet while the
        arrays are one rank's share, so it is passed in rather than recomputed.
    new_average_depth : float
        Target average depth.
    """

    stages = num.asarray(stages, dtype=float)
    elevations = num.asarray(elevations, dtype=float)
    areas = num.asarray(areas, dtype=float)

    total_area = float(areas.sum())
    volume = (new_average_depth - old_average_depth) * total_area

    # Exactly zero transfer must be exactly a no-op: this is the lake-at-rest /
    # well-balancedness case, and it must not be disturbed even at roundoff.
    if volume == 0.0 or total_area <= 0.0:
        return stages.copy()

    if volume > 0.0:
        # Fill: raise the lowest stages to a common level L, where the volume
        # added below L, A(L) = sum a_i * max(L - s_i, 0), equals `volume`.
        # A is continuous, increasing, and A(max_s + volume/total_area) >=
        # total_area * (volume/total_area) = volume, so L lies in [lo, hi].
        lo = float(stages.min())
        hi = float(stages.max()) + volume / total_area
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if float(num.dot(areas, num.maximum(mid - stages, 0.0))) < volume:
                lo = mid
            else:
                hi = mid
        return num.maximum(stages, hi)

    # Drawdown: lower the highest stages to a common level L, clamping at the
    # bed. Removed volume R(L) = sum a_i * min(depth_i, max(s_i - L, 0)) is
    # continuous and decreasing, with R(min_z) = all the water and R(max_s) = 0.
    to_remove = -volume
    depths = num.maximum(stages - elevations, 0.0)
    water = float(num.dot(areas, depths))
    if to_remove >= water:
        # Asked for more than the inlet holds: drain it dry, no further.
        return num.where(depths > 0.0, elevations, stages).astype(float)

    lo = float(elevations.min())
    hi = float(stages.max())
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        removed = float(num.dot(
            areas, num.minimum(depths, num.maximum(stages - mid, 0.0))))
        if removed > to_remove:
            lo = mid
        else:
            hi = mid
    return num.maximum(elevations, num.minimum(stages, lo))


class Inlet:
    """Contains information associated with each inlet
    """

    def __init__(self, domain, poly, verbose=False):

        self.domain = domain
        self.domain_bounding_polygon = self.domain.get_boundary_polygon()
        self.verbose = verbose


        # poly can be either a line, polygon or a region
        if isinstance(poly,Region):
            self.region = poly
        else:
            self.region = Region(domain,poly=poly,expand_polygon=True)

        self.triangle_indices = self.region.indices

        self.compute_area()


    def compute_area(self):

        # Compute inlet area as the sum of areas of triangles identified
        # by line. Must be called after compute_inlet_triangle_indices().
        if len(self.triangle_indices) == 0:
            region = 'Inlet line=%s' % (self.inlet_line)
            msg = 'No triangles have been identified in region '
            raise Exception(msg)

#        self.area = 0.0
#        for j in self.triangle_indices:
#            self.area += self.domain.areas[j]

        self.area = num.sum(self.domain.areas[self.triangle_indices])

        msg = 'Inlet exchange area has area = %f' % self.area
        assert self.area > 0.0



    def get_poly(self):

        return self.poly

    def get_area(self):

        return self.area


    def get_areas(self):

        # Must be called after compute_inlet_triangle_indices().
        return self.domain.areas.take(self.triangle_indices)


    def get_stages(self):

        return self.domain.quantities['stage'].centroid_values.take(self.triangle_indices)
        # self.domain.quantities['stage'].centroid_values would be called gpu_stage_centroid_values in the gpu_interface

    def get_average_stage(self):

        return num.sum(self.get_stages()*self.get_areas())/self.area

    def get_elevations(self):

        return self.domain.quantities['elevation'].centroid_values.take(self.triangle_indices)

    def get_average_elevation(self):

        return num.sum(self.get_elevations()*self.get_areas())/self.area


    def get_xmoms(self):

        return self.domain.quantities['xmomentum'].centroid_values.take(self.triangle_indices)


    def get_average_xmom(self):

        return num.sum(self.get_xmoms()*self.get_areas())/self.area


    def get_ymoms(self):

        return self.domain.quantities['ymomentum'].centroid_values.take(self.triangle_indices)


    def get_average_ymom(self):

        return num.sum(self.get_ymoms()*self.get_areas())/self.area


    def get_depths(self):

        return self.get_stages() - self.get_elevations()


    def get_total_water_volume(self):

       return num.sum(self.get_depths()*self.get_areas())


    def get_average_depth(self):

        return self.get_total_water_volume()/self.area


    def get_velocities(self):

            depths = self.get_depths()
            u = self.get_xmoms()*depths/(depths*depths + velocity_protection)
            v = self.get_ymoms()*depths/(depths*depths + velocity_protection)

            return u, v


    def get_xvelocities(self):

            depths = self.get_depths()
            return self.get_xmoms()*depths/(depths*depths + velocity_protection)

    def get_yvelocities(self):

            depths = self.get_depths()
            return self.get_ymoms()*depths/(depths*depths + velocity_protection)


    def get_average_speed(self):

            u, v = self.get_velocities()

            average_u = num.sum(u*self.get_areas())/self.area
            average_v = num.sum(v*self.get_areas())/self.area

            return math.sqrt(average_u**2 + average_v**2)


    def get_average_velocity_head(self):

        return 0.5*self.get_average_speed()**2/g


    def get_average_total_energy(self):

        return self.get_average_velocity_head() + self.get_average_stage()


    def get_average_specific_energy(self):

        return self.get_average_velocity_head() + self.get_average_depth()



    def set_depths(self,depth):

        self.domain.quantities['stage'].centroid_values.put(self.triangle_indices, self.get_elevations() + depth)


    def set_stages(self,stage):

        self.domain.quantities['stage'].centroid_values.put(self.triangle_indices, stage)

    def set_average_momenta(self, xmom, ymom):
        """Write inlet-average momenta (xmom, ymom), distributed per cell ∝ depth.

        The structure operators compute ONE new momentum value per inlet. The old
        write applied it uniformly, which was consistent with the old uniform-
        depth write (uniform depth + uniform momentum = uniform velocity). With
        the leveling write-back the depths across the inlet differ, and a uniform
        momentum over a nearly-dry cell is an enormous velocity — enough to
        collapse the global timestep. Weighting by depth keeps the VELOCITY field
        uniform instead: cell i gets ``m * depth_i / average_depth``, so the
        area-weighted average momentum is exactly ``m`` and velocity is bounded
        by ``m / average_depth`` everywhere. On uniform depth this reduces to the
        old uniform write.
        """

        if self.area <= 0.0:
            # A rank that takes part in a parallel inlet for its enquiry
            # point alone holds none of the inlet's triangles: nothing to
            # write, and nothing to divide by.
            return
        depths = num.maximum(self.get_stages() - self.get_elevations(), 0.0)
        average_depth = float(num.dot(depths, self.get_areas())) / self.area
        if average_depth <= 0.0:
            weights = 0.0 * depths
        else:
            weights = depths / average_depth
        self.set_xmoms(xmom * weights)
        self.set_ymoms(ymom * weights)

    def set_average_depth(self, new_average_depth, old_average_depth=None):
        """Fill or draw down the inlet to give `new_average_depth`.

        The structure operators' way of writing back a transfer. Unlike
        set_depths(), which flattens the DEPTH and so tilts the water surface on
        a sloping bed, this levels the STAGE — water finds its level — and
        clamps cells at their bed rather than driving them to negative depth.
        See level_stages_to_average().
        """

        if old_average_depth is None:
            old_average_depth = self.get_average_depth()

        new_stages = level_stages_to_average(
            self.get_stages(), self.get_elevations(), self.get_areas(),
            old_average_depth, new_average_depth)
        self.set_stages(new_stages)


    def set_xmoms(self,xmom):

        self.domain.quantities['xmomentum'].centroid_values.put(self.triangle_indices, xmom)


    def set_ymoms(self,ymom):

        self.domain.quantities['ymomentum'].centroid_values.put(self.triangle_indices, ymom)


    def set_elevations(self,elevation):

        self.domain.quantities['elevation'].centroid_values.put(self.triangle_indices, elevation)

    def set_stages_evenly(self,volume):
        """ Distribute volume of water over
        inlet exchange region so that stage is level
        """

        assert volume >= 0.0

        areas = self.get_areas()
        stages = self.get_stages()
        depths = self.get_depths()

        elevations = self.get_elevations()

        #print('elevation')
        #print(elevations)

        stages_order = stages.argsort()

        # accumulate areas of cells ordered by stage
        summed_areas = num.cumsum(areas[stages_order])

        # accumulate the volume need to fill cells
        summed_volume = num.zeros_like(areas)
        summed_volume[1:] = num.cumsum(summed_areas[:-1]*num.diff(stages[stages_order]))

        index = num.nonzero(summed_volume<=volume)[0][-1]

        # calculate stage needed to fill chosen cells with given volume of water
        depth = (volume - summed_volume[index])/summed_areas[index]
        stages[stages_order[0:index+1]] = stages[stages_order[index]]+depth

        #print('stages')
        #print(stages)
        self.set_stages(stages)




    def set_depths_evenly(self,volume):
        """ Distribute volume over all exchange
        cells with equal depth of water
        """

        new_depth = self.get_average_depth() + (volume/self.get_area())
        self.set_depths(new_depth)

