"""Sediment: the legacy (mode 1) and unified (mode 2) paths agree.

Both compute modes call the same sediment code in gpu/core_kernels.c, so a
disagreement is a mapping or binding fault rather than a physics one. The
characteristic fault is a `D->member` load inside an `omp target` region, which
reads a host address on the device and makes the work silently not happen --
and a second is a struct member bound in one Cython extension but not the
other, which leaves an enable flag reading uninitialised memory.

Every mode 1 / mode 2 comparison lives here rather than in the per-topic
modules, for two reasons: they all need the same module-level guard, and each
mode-2 domain costs the NVHPC runtime, so concentrating them keeps the count
in one process low.
"""
import os
import warnings

import numpy as np
import pytest

import anuga
from anuga import Reflective_boundary, rectangular_cross_domain

LEN = 500.0

_gpu_error = None
_gpu_avail = None


def gpu_available():
    global _gpu_error, _gpu_avail
    if _gpu_avail is not None:
        return _gpu_avail
    try:
        from anuga.shallow_water.sw_domain_gpu_ext import init_gpu_domain  # noqa: F401
        _gpu_avail = True
    except Exception as e:
        _gpu_avail = False
        _gpu_error = '%s: %s' % (type(e).__name__, e)
    return _gpu_avail


# Mirrors test_DE_gpu_omp.py: on a GPU-offload build the NVHPC OpenMP-target
# runtime aborts a process that creates many mode-2 domains, so this file skips
# in a normal in-process run and is opted back in by the isolated runner, which
# sets ANUGA_GPU_TESTS_ISOLATED=1.
if (gpu_available() and anuga.gpu_offload_supported()
        and not os.environ.get('ANUGA_GPU_TESTS_ISOLATED')):
    _skip_reason = (
        "GPU-offload build: run this file via "
        "anuga/shallow_water/tests/run_gpu_tests_isolated.sh (one fresh process "
        "per class) - running it in one process aborts the NVHPC OpenMP-target "
        "runtime. Set ANUGA_GPU_TESTS_ISOLATED=1 to force in-process collection.")
    warnings.warn(_skip_reason, stacklevel=1)
    pytest.skip(_skip_reason, allow_module_level=True)

pytestmark = pytest.mark.skipif(not gpu_available(),
                                reason='GPU OpenMP interface not available')


def _require_mode_2(d):
    """Comparing mode 1 against a silent mode-1 fallback is a false green."""
    if getattr(d, 'multiprocessor_mode', None) != 2:
        pytest.fail('mode 2 did not engage: multiprocessor_mode=%r'
                    % getattr(d, 'multiprocessor_mode', None))
    return d


def channel(mode=1, slope=0.01, depth=1.0, n_manning=0.03, dt=1.0,
            nxy=(20, 10), length=(LEN, LEN / 2)):
    d = rectangular_cross_domain(nxy[0], nxy[1], len1=length[0], len2=length[1])
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', lambda x, y: -slope * x)
    d.set_quantity('stage', lambda x, y: -slope * x + depth)
    d.set_quantity('friction', n_manning)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.evolve_max_timestep = dt
    if mode != 1:
        d.set_multiprocessor_mode(mode)
        _require_mode_2(d)
    return d


def _both(configure, finaltime, read):
    """Run the same configuration in both modes and return what `read` gives."""
    out = []
    for mode in (1, 2):
        d = configure(mode)
        d.evolve_to_end(finaltime=finaltime)
        out.append(read(d))
    return out


def test_deposition_agrees():
    def configure(mode):
        d = rectangular_cross_domain(10, 10, len1=LEN, len2=LEN)
        d.set_flow_algorithm('DE0')
        d.set_low_froude(0)
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 1.0)
        d.set_boundary({t: Reflective_boundary(d)
                        for t in d.get_boundary_tags()})
        d.evolve_max_timestep = 1.0
        if mode != 1:
            d.set_multiprocessor_mode(mode)
            _require_mode_2(d)
        d.add_sediment_fraction(name='sand', diameter=1.0e-4, initial_concentration=0.05)
        return d

    a, b = _both(configure, 30.0, lambda d: d.get_tracer('sand').copy())
    assert np.abs(a - b).max() < 1e-8


def test_entrainment_agrees():
    def configure(mode):
        d = channel(mode=mode, nxy=(10, 10), length=(LEN, LEN))
        d.add_sediment_fraction(name='sand', diameter=1.0e-4, tau_c_star=0.04,
                         initial_concentration=0.0)
        return d

    a, b = _both(configure, 40.0, lambda d: d.get_tracer('sand').copy())
    assert np.abs(a - b).max() < 1e-8


def test_the_rouse_near_bed_ratio_agrees():
    def configure(mode):
        d = rectangular_cross_domain(10, 10, len1=LEN, len2=LEN)
        d.set_flow_algorithm('DE0')
        d.set_low_froude(0)
        d.store = False
        d.set_quantity('elevation', lambda x, y: -0.01 * x)
        d.set_quantity('stage', lambda x, y: -0.01 * x + 1.0)
        d.set_quantity('friction', 0.03)
        d.set_boundary({t: Reflective_boundary(d)
                        for t in d.get_boundary_tags()})
        d.evolve_max_timestep = 1.0
        d.sediment_d_star_mode = 1
        if mode != 1:
            d.set_multiprocessor_mode(mode)
            _require_mode_2(d)
        d.add_sediment_fraction(name='s', diameter=1e-4, initial_concentration=0.02)
        return d

    a, b = _both(configure, 30.0, lambda d: d.get_tracer('s').copy())
    assert np.abs(a - b).max() < 1e-8


def test_the_wilson_friction_closure_agrees():
    def configure(mode):
        d = channel(mode=mode, nxy=(20, 10))
        d.set_sediment_friction('wilson', bed='gravel', grain_size=0.05)
        d.sediment_d_star_mode = 1
        d.add_sediment_fraction(name='s', diameter=1e-4, initial_concentration=0.01)
        return d

    a, b = _both(configure, 20.0, lambda d: d.get_tracer('s').copy())
    assert np.abs(a - b).max() < 1e-8


def test_the_depth_slope_shear_closure_agrees():
    def configure(mode):
        d = channel(mode=mode, nxy=(20, 10), length=(200.0, 100.0))
        d.set_shear_closure('depth_slope')
        d.add_sediment_fraction(name='sand', diameter=1e-4, initial_concentration=0.0)
        return d

    a, b = _both(configure, 20.0, lambda d: d.get_tracer('sand').copy())
    assert np.abs(a - b).max() < 1e-8


def test_the_frozen_capped_depth_slope_agrees():
    """The frozen slope lives in a host array uploaded at map time and the
    cap is a scalar; both must reach the device."""
    def configure(mode):
        d = channel(mode=mode, nxy=(20, 10), length=(200.0, 100.0))
        d.set_shear_closure('depth_slope', max_slope=0.005, freeze_slope=True)
        d.add_sediment_fraction(name='sand', diameter=1e-4, initial_concentration=0.0)
        return d

    a, b = _both(configure, 20.0, lambda d: (d.get_tracer('sand').copy(),
                                             d.quantities['elevation'].centroid_values.copy()))
    assert np.abs(a[0] - b[0]).max() < 1e-8
    assert np.abs(a[1] - b[1]).max() < 1e-8


def test_the_cohesive_erosion_route_agrees():
    def configure(mode):
        d = channel(mode=mode, nxy=(20, 10), length=(200.0, 100.0))
        d.set_bed_material('cohesive', tau_crit=0.088)
        d.add_sediment_fraction(name='silt', diameter=6.5e-5, initial_concentration=0.0)
        return d

    a, b = _both(configure, 20.0, lambda d: d.get_tracer('silt').copy())
    assert np.abs(a - b).max() < 1e-8


def test_exner_bed_evolution_agrees():
    def configure(mode):
        d = rectangular_cross_domain(42, 21, len1=2.0, len2=1.0)
        d.set_flow_algorithm('DE0')
        d.set_low_froude(0)
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 1.0)
        d.set_quantity('friction', 0.0)
        d.set_boundary({t: Reflective_boundary(d)
                        for t in d.get_boundary_tags()})
        d.sediment_c_max = 0.80
        d.sediment_c_pack = 0.80
        d.sediment_porosity = 0.28
        if mode != 1:
            d.set_multiprocessor_mode(mode)
            _require_mode_2(d)
        x, y = d.centroid_coordinates[:, 0], d.centroid_coordinates[:, 1]
        c0 = 0.7 * np.exp(-5.0 * (x - 0.9) ** 2 - 50.0 * (y - 0.5) ** 2)
        d.add_sediment_fraction(name='sand', diameter=0.01, rho_s=2400.0, tau_c_star=0.0,
                         initial_concentration=c0)
        return d

    a, b = _both(configure, 20.0,
                 lambda d: d.quantities['elevation'].centroid_values.copy())
    assert np.abs(a - b).max() < 1e-8


def test_bedload_agrees():
    def configure(mode):
        d = channel(mode=mode, slope=0.02, depth=0.5, n_manning=0.025, dt=0.5,
                    nxy=(20, 10), length=(100.0, 50.0))
        d.add_sediment_fraction(name='gravel', diameter=5e-3, tau_c_star=0.0,
                         initial_concentration=0.0)
        d.set_bedload('wong_parker_eq24')
        return d

    a, b = _both(configure, 30.0,
                 lambda d: d.quantities['elevation'].centroid_values.copy())
    assert np.abs(a - b).max() < 1e-8


def test_the_erodible_base_and_bedload_agree():
    """[L-5] adds three struct members and a device mapping, so this is where a
    binding bound in one extension but not the other would show."""
    def configure(mode):
        d = rectangular_cross_domain(30, 8, len1=60.0, len2=16.0)
        d.set_flow_algorithm('DE0')
        d.set_low_froude(0)
        d.store = False
        d.set_quantity('elevation', lambda x, y: -0.01 * x)
        d.set_quantity('stage', lambda x, y: -0.01 * x + 0.6)
        d.set_quantity('xmomentum', 1.2)
        d.set_quantity('friction', 0.03)
        d.set_boundary({t: Reflective_boundary(d)
                        for t in d.get_boundary_tags()})
        d.set_sediment_parameters(porosity=0.3)
        d.set_bedload('wong_parker_eq24')
        for name, diameter in (('fine', 1.0e-4), ('coarse', 6.0e-4)):
            d.add_sediment_fraction(name=name, diameter=diameter)
        d.set_erodible_base(depth=0.02)
        if mode != 1:
            d.set_multiprocessor_mode(mode)
            _require_mode_2(d)
        return d

    def read(d):
        return (d.quantities['elevation'].centroid_values.copy(),
                d.tracer_conserved_values.copy())

    (z1, m1), (z2, m2) = _both(configure, 25.0, read)
    assert np.abs(z1 - z2).max() < 1e-10
    assert np.abs(m1 - m2).max() < 1e-10


def test_angle_of_repose_relaxation_agrees():
    """The relaxation sweeps are Jacobi precisely so that the answer does not
    depend on thread order; if it did, this is where it would show."""
    def configure(mode):
        d = rectangular_cross_domain(40, 12, len1=60.0, len2=16.0)
        d.set_flow_algorithm('DE0')
        d.set_low_froude(0)
        d.store = False
        x, y = d.centroid_coordinates[:, 0], d.centroid_coordinates[:, 1]
        r = np.sqrt((x - 30.0) ** 2 + (y - 8.0) ** 2)
        d.set_quantity('elevation',
                       np.where(r < 8.0, 6.0 * (1.0 - r / 8.0), 0.0),
                       location='centroids')
        d.set_quantity('stage', -1.0)
        d.set_quantity('friction', 0.03)
        d.set_boundary({t: Reflective_boundary(d)
                        for t in d.get_boundary_tags()})
        d.set_sediment_parameters(porosity=0.3)
        d.add_sediment_fraction(name='sand', diameter=2.0e-4)
        d.set_angle_of_repose(30.0, max_sweeps=400)
        if mode != 1:
            d.set_multiprocessor_mode(mode)
            _require_mode_2(d)
        return d

    a, b = _both(configure, 1.0,
                 lambda d: d.quantities['elevation'].centroid_values.copy())
    assert np.abs(a - b).max() < 1e-10


@pytest.mark.parametrize('composition', [False, True])
def test_two_classes_on_an_erodible_base_agree(composition):
    """Two classes over an erodible base, with and without bed composition.

    Guards the explicit device maps on the sediment source loop: with them
    implicit, a longer kernel made nvc map the per-class parameters as
    zero-length sections, so the device read uninitialised memory and mode 2
    drifted 0.25 m from mode 1 here, differently on every run."""
    def configure(mode):
        d = channel(mode=mode, nxy=(10, 10), length=(LEN, LEN))
        d.set_sediment_parameters(porosity=0.3)
        d.add_sediment_fraction(name='fine', diameter=1.0e-4, tau_c_star=0.04,
                                initial_concentration=0.0)
        d.add_sediment_fraction(name='coarse', diameter=6.0e-4, tau_c_star=0.04,
                                initial_concentration=0.0)
        d.set_erodible_base(depth=0.5)
        if composition:
            d.set_bed_composition({'fine': 0.5, 'coarse': 0.5}, active_layer=0.005)
        return d

    def read(d):
        out = [d.quantities['elevation'].centroid_values.copy(),
               d.get_tracer('fine').copy(), d.get_tracer('coarse').copy()]
        if composition:
            out.append(d.sediment_bed_active.copy())
        return out

    a, b = _both(configure, 40.0, read)
    for x, y in zip(a, b):
        assert np.abs(x - y).max() < 1e-8
