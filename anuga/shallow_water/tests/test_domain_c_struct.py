"""Tests for the Python -> C domain struct boundary (sw_domain_openmp_ext).

Two regressions from the C extension audit:

* array shapes are validated when the struct is (re)filled, so a wrong-sized
  array raises instead of becoming an out-of-bounds pointer in the kernels
  (issue #325);
* the flux substep bookkeeping is per domain, so two domains evolving in
  the same process do not interleave each other's counters (issue #326).
"""
import numpy as num
import pytest

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.shallow_water.boundaries import Reflective_boundary


def _make_domain(method):
    points, vertices, boundary = rectangular_cross(6, 6, len1=10.0, len2=10.0)
    domain = Domain(points, vertices, boundary)
    domain.set_compute_mode('legacy')
    domain.set_timestepping_method(method)
    domain.set_quantity('elevation', lambda x, y: 0.02 * x)
    domain.set_quantity('friction', 0.03)
    # A sloped free surface so there is flow and a non-trivial boundary flux
    domain.set_quantity('stage', lambda x, y: 1.0 + 0.1 * (x > 5.0))
    domain.set_boundary({'left': Reflective_boundary(domain),
                         'right': Reflective_boundary(domain),
                         'top': Reflective_boundary(domain),
                         'bottom': Reflective_boundary(domain)})
    return domain


class TestShapeValidation:

    def test_wrong_length_mesh_array_is_rejected(self):
        domain = _make_domain('euler')
        domain.update_domain_c_struct()   # the correct shapes pass

        domain.areas = domain.areas[:-1].copy()
        with pytest.raises(ValueError, match=r'domain\.areas has shape'):
            domain.update_domain_c_struct()

    def test_wrong_shape_quantity_array_is_rejected(self):
        domain = _make_domain('euler')
        q = domain.quantities['stage']
        q.centroid_values = num.zeros(domain.number_of_elements + 3)
        with pytest.raises(ValueError, match=r'stage\.centroid_values has shape'):
            domain.update_domain_c_struct()

    def test_wrong_boundary_length_is_rejected(self):
        domain = _make_domain('euler')
        q = domain.quantities['xmomentum']
        q.boundary_values = num.zeros(domain.boundary_length + 1)
        with pytest.raises(ValueError, match=r'xmomentum\.boundary_values'):
            domain.update_domain_c_struct()

    def test_none_arrays_are_allowed(self):
        # Quantities with reduced storage have no update arrays; the C
        # struct receives NULL for them and the kernels never touch them.
        domain = _make_domain('euler')
        assert domain.quantities['friction'].explicit_update is None
        domain.update_domain_c_struct()


class TestPerDomainSubstepCounters:

    @staticmethod
    def _flowing_domain(method):
        # A Dirichlet inflow on the left so the boundary mass flux is
        # non-zero and the substep slot it lands in is observable.
        domain = _make_domain(method)
        from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import Dirichlet_boundary
        domain.set_boundary({'left': Dirichlet_boundary([2.0, 1.0, 0.0])})
        domain.update_boundary()
        return domain

    def test_second_substep_survives_a_call_from_another_domain(self):
        """rk2 makes two flux calls per step. If another domain's flux call
        lands between them, the second one must still be treated as the
        second substep: its boundary flux sum goes to slot 1 and, as for any
        non-first substep, the timestep passed in is returned unchanged.

        With the old module-level static counters the intervening call
        (from a domain with a different timestep_fluxcalls) reset the phase,
        so the second call was treated as a new first substep: slot 0 was
        overwritten, slot 1 stayed empty, and the flux timestep was
        recomputed mid-step.
        """
        a_solo = self._flowing_domain('rk2')
        a_solo.compute_fluxes()
        a_solo.compute_fluxes()
        flux = a_solo.boundary_flux_sum[0]
        assert flux != 0.0
        assert a_solo.boundary_flux_sum[1] == flux
        assert a_solo.flux_timestep == a_solo.evolve_max_timestep

        a = self._flowing_domain('rk2')
        b = self._flowing_domain('euler')
        a.compute_fluxes()            # A substep 0
        b.compute_fluxes()            # B, one flux call per step
        a.compute_fluxes()            # A substep 1 -- must not be reset by B

        assert a.boundary_flux_sum[0] == flux
        assert a.boundary_flux_sum[1] == flux
        assert a.flux_timestep == a.evolve_max_timestep

        # B is a one-call-per-step domain; every call is substep 0 for it.
        b.compute_fluxes()
        assert b.boundary_flux_sum[0] != 0.0
        assert b.flux_timestep < b.evolve_max_timestep
