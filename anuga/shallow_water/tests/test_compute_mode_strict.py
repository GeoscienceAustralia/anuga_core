"""Strict compute mode: fail instead of falling back (audit issue #331).

Three silent downgrades used to be warnings or prints only: the MPI fallback
from 'unified' to 'legacy', running on the CPU when GPU offload is enabled but
no device was found, and boundary types the device path cannot evaluate. With
``strict=True`` (or ``ANUGA_STRICT_COMPUTE_MODE=1``) each is a RuntimeError,
raised at setup rather than mid-evolve.
"""
import warnings

import numpy as num
import pytest

import anuga
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain, GPU_BOUNDARY_TYPES
from anuga.shallow_water import shallow_water_domain as swd
from anuga.shallow_water.boundaries import (Reflective_boundary,
                                            Transmissive_momentum_set_stage_boundary)


def _domain(**kw):
    points, vertices, boundary = rectangular_cross(4, 4, len1=4.0, len2=4.0)
    domain = Domain(points, vertices, boundary)
    domain.set_compute_mode('legacy')          # start neutral; tests switch
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 1.0)
    return domain


def _all_reflective(domain):
    Br = Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


def _one_host_only_boundary(domain):
    Br = Reflective_boundary(domain)
    Bh = Transmissive_momentum_set_stage_boundary(domain, function=lambda t: 1.0)
    assert Bh.__class__.__name__ not in GPU_BOUNDARY_TYPES
    domain.set_boundary({'left': Bh, 'right': Br, 'top': Br, 'bottom': Br})


class TestStrictFlag:

    def test_default_is_not_strict(self, monkeypatch):
        monkeypatch.delenv('ANUGA_STRICT_COMPUTE_MODE', raising=False)
        d = _domain()
        assert d.compute_mode_strict is False

    def test_environment_variable_sets_the_default(self, monkeypatch):
        monkeypatch.setenv('ANUGA_STRICT_COMPUTE_MODE', '1')
        d = _domain()
        assert d.compute_mode_strict is True

    def test_explicit_argument_sticks_across_calls(self, monkeypatch):
        monkeypatch.delenv('ANUGA_STRICT_COMPUTE_MODE', raising=False)
        d = _domain()
        d.set_compute_mode('legacy', strict=True)
        d.set_compute_mode('legacy')               # strict=None keeps the setting
        assert d.compute_mode_strict is True
        d.set_compute_mode('legacy', strict=False)
        assert d.compute_mode_strict is False


class TestUnsupportedBoundary:

    def test_non_strict_warns_once_at_setup_and_runs(self):
        d = _domain()
        _one_host_only_boundary(d)
        with pytest.warns(UserWarning, match='device path cannot evaluate'):
            d.set_compute_mode('unified')
        assert d._gpu_all_on_gpu is False
        assert d._gpu_cpu_boundary_types == [('left', 'Transmissive_momentum_set_stage_boundary')]
        # The classification is cached: no second warning on evolve.
        with warnings.catch_warnings():
            warnings.simplefilter('error', UserWarning)
            for _ in d.evolve(yieldstep=0.05, finaltime=0.1):
                pass

    def test_strict_refuses_at_setup(self):
        d = _domain()
        _one_host_only_boundary(d)
        with pytest.raises(RuntimeError, match='strict mode'):
            d.set_compute_mode('unified', strict=True)

    def test_strict_is_silent_when_every_boundary_runs_on_the_device(self):
        d = _domain()
        _all_reflective(d)
        with warnings.catch_warnings():
            warnings.simplefilter('error', UserWarning)
            d.set_compute_mode('unified', strict=True)
        assert d._gpu_all_on_gpu is True


class TestFlatherIsDeviceEvaluated:

    def test_flather_boundary_is_classified_for_the_device(self):
        # The class name used to be misspelt in most of the classification
        # copies, so a Flather boundary fell to the host path there.
        d = _domain()
        Br = Reflective_boundary(d)
        Bf = anuga.Flather_external_stage_zero_velocity_boundary(d, function=lambda t: 1.0)
        d.set_boundary({'left': Bf, 'right': Br, 'top': Br, 'bottom': Br})
        d.set_compute_mode('unified', strict=True)   # strict: any host fallback raises
        assert d._gpu_all_on_gpu is True
        assert d._gpu_flather_boundaries == [Bf]


def _device_id_of_a_unified_domain():
    d = _domain()
    _all_reflective(d)
    d.set_compute_mode('unified')
    return d.gpu_interface.gpu_dom.device_id


class TestNoDeviceUnderOffload:

    def test_strict_refuses_the_host_fallback(self, monkeypatch):
        # Pretend offload was enabled process-wide on a machine with no
        # device: strict must refuse rather than run on the host through the
        # OpenMP target fallback. Only meaningful where there is no device.
        if _device_id_of_a_unified_domain() >= 0:
            pytest.skip('a GPU device is present; the no-device path cannot be exercised here')
        d = _domain()
        _all_reflective(d)
        monkeypatch.setattr(swd, 'gpu_offload_enabled', lambda: True)
        with pytest.raises(RuntimeError, match='no GPU device was found'):
            d.set_compute_mode('unified', strict=True)

    def test_non_strict_runs_on_the_host(self, monkeypatch):
        d = _domain()
        _all_reflective(d)
        monkeypatch.setattr(swd, 'gpu_offload_enabled', lambda: True)
        d.set_compute_mode('unified')
        assert d.gpu_interface is not None


class TestMPIFallback:

    def test_strict_refuses_the_legacy_fallback(self, monkeypatch):
        d = _domain()
        _all_reflective(d)
        monkeypatch.setattr(anuga, 'numprocs', 2)
        monkeypatch.setattr(Domain, '_mode2_mpi_available', lambda self: False)
        with pytest.raises(RuntimeError, match='refusing to fall back'):
            d.set_compute_mode('unified', strict=True)

    def test_non_strict_warns_and_falls_back(self, monkeypatch):
        d = _domain()
        _all_reflective(d)
        monkeypatch.setattr(anuga, 'numprocs', 2)
        monkeypatch.setattr(anuga, 'myid', 0)
        monkeypatch.setattr(Domain, '_mode2_mpi_available', lambda self: False)
        with pytest.warns(UserWarning, match='Falling back'):
            d.set_compute_mode('unified')
        assert d.compute_mode == 'legacy'
        assert d.requested_compute_mode == 'unified'


class TestSetGpuOffloadStrict:

    def test_strict_raises_when_offload_is_unsupported(self, monkeypatch):
        monkeypatch.setattr(swd, 'gpu_offload_supported', lambda: False)
        with pytest.raises(RuntimeError, match='strict mode'):
            swd.set_gpu_offload(True, verbose=False, strict=True)

    def test_environment_variable_makes_it_strict(self, monkeypatch):
        monkeypatch.setattr(swd, 'gpu_offload_supported', lambda: False)
        monkeypatch.setenv('ANUGA_STRICT_COMPUTE_MODE', 'yes')
        with pytest.raises(RuntimeError):
            swd.set_gpu_offload(True, verbose=False)

    def test_non_strict_warns_and_returns_false(self, monkeypatch):
        monkeypatch.setattr(swd, 'gpu_offload_supported', lambda: False)
        monkeypatch.delenv('ANUGA_STRICT_COMPUTE_MODE', raising=False)
        with pytest.warns(UserWarning):
            assert swd.set_gpu_offload(True, verbose=False) is False
