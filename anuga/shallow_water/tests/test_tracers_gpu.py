"""Generic passive tracers: the legacy (mode 1) and unified (mode 2) paths agree.

Both compute modes call the same tracer code in `gpu/core_kernels.c`, so a
disagreement between them is a mapping or binding fault rather than a physics
one -- most often a `D->member` load inside an `omp target` region, which reads
a host address on the device and makes the work silently not happen.

The properties tested here are exactly the ones such a fault breaks: the two
modes must agree on the hydrodynamics, on each tracer's concentration and
conserved mass, and on the total; and two tracers must stay distinct on the
device rather than one aliasing the other.

Deliberately NO `sync_from_device()` after the mode-2 evolve. The comparison
reads what an ordinary user's script would see when it inspects centroid values
after `evolve`, which is the thing that has to be right.
"""
import os
import warnings

import numpy as np
import pytest

import anuga
from anuga import Reflective_boundary, rectangular_cross_domain

LEN = 500.0
NXY = 20
FINALTIME = 20.0

_gpu_error = None
_gpu_avail = None


def gpu_available():
    """Check if the GPU OpenMP interface is importable."""
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


# On a GPU-offload build (nvc, -Dgpu_offload=true) the NVHPC OpenMP-target
# runtime aborts the process once many mode-2 domains have been created in it,
# so files that build them are skipped in a normal in-process run and are
# executed one class per fresh process by run_gpu_tests_isolated.sh, which sets
# ANUGA_GPU_TESTS_ISOLATED=1. On a CPU build the omp-target regions run on the
# host and this file is fine in one process. Mirrors test_DE_gpu_omp.py.
if (gpu_available() and anuga.gpu_offload_supported()
        and not os.environ.get('ANUGA_GPU_TESTS_ISOLATED')):
    _skip_reason = (
        "GPU-offload build: run this file via "
        "anuga/shallow_water/tests/run_gpu_tests_isolated.sh (one fresh process "
        "per class) - running it in one process aborts the NVHPC OpenMP-target "
        "runtime. Set ANUGA_GPU_TESTS_ISOLATED=1 to force in-process collection.")
    warnings.warn(_skip_reason, stacklevel=1)
    pytest.skip(_skip_reason, allow_module_level=True)


def _build(mode):
    d = rectangular_cross_domain(NXY, NXY, len1=LEN, len2=LEN)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', lambda x, y: np.where(x < LEN / 2, 2.0, 0.5))
    d.set_quantity('xmomentum', 0.0)
    d.set_quantity('ymomentum', 0.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    # One uniform tracer (tests consistency) and one structured (tests that the
    # device sees per-slot data rather than aliasing slot 0).
    d.add_tracer('uniform', beta=1.0)
    d.add_tracer('wedge', beta=1.0)
    x = d.centroid_coordinates[:, 0]
    d.set_tracer('uniform', 1.0)
    d.set_tracer('wedge', np.where(x < LEN / 2, 1.0, 0.0))
    d.set_multiprocessor_mode(mode)
    return d


@pytest.fixture(scope='module')
def evolved():
    """One mode-1 and one mode-2 run of the same problem."""
    if not gpu_available():
        pytest.skip('GPU OpenMP interface not available: %s' % _gpu_error)

    cpu = _build(1)
    cpu.evolve_to_end(finaltime=FINALTIME)

    gpu = _build(2)
    if getattr(gpu, 'multiprocessor_mode', None) != 2:
        # Comparing mode 1 against a silent mode-1 fallback would be a false
        # green, so refuse rather than pass.
        pytest.fail('mode 2 did not engage: multiprocessor_mode=%r compute_mode=%r'
                    % (getattr(gpu, 'multiprocessor_mode', None),
                       getattr(gpu, 'compute_mode', None)))

    from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device
    sync_to_device(gpu.gpu_interface.gpu_dom)
    gpu.evolve_to_end(finaltime=FINALTIME)
    return cpu, gpu


@pytest.mark.parametrize('name', ['stage', 'xmomentum', 'ymomentum'])
def test_hydrodynamics_agree(evolved, name):
    cpu, gpu = evolved
    a = cpu.quantities[name].centroid_values
    b = gpu.quantities[name].centroid_values
    assert np.abs(a - b).max() < 1e-10


@pytest.mark.parametrize('slot,name', [(0, 'uniform'), (1, 'wedge')])
def test_tracer_fields_agree(evolved, slot, name):
    cpu, gpu = evolved
    ma = cpu.tracer_conserved_values[slot]
    mb = gpu.tracer_conserved_values[slot]
    assert np.abs(ma - mb).max() < 1e-10, 'conserved m disagrees for %r' % name

    ha = (cpu.quantities['stage'].centroid_values
          - cpu.quantities['elevation'].centroid_values)
    hb = (gpu.quantities['stage'].centroid_values
          - gpu.quantities['elevation'].centroid_values)
    wet = (ha > 1e-3) & (hb > 1e-3)
    assert np.abs(ma[wet] / ha[wet] - mb[wet] / hb[wet]).max() < 1e-10

    tot_a = float((ma * cpu.areas).sum())
    tot_b = float((mb * gpu.areas).sum())
    assert abs(tot_a - tot_b) < 1e-10 * max(abs(tot_a), 1.0)


def test_the_two_tracers_stay_distinct_on_the_device(evolved):
    """Guards against slot 1 aliasing slot 0 in the device mapping."""
    _, gpu = evolved
    assert not np.allclose(gpu.tracer_conserved_values[0],
                           gpu.tracer_conserved_values[1])


def test_depth_is_consistent_with_stage_in_both_modes(evolved):
    """A tracer must not corrupt the hydrodynamic state it rides on."""
    for nm, d in zip(('mode 1', 'mode 2'), evolved):
        h = (d.quantities['stage'].centroid_values
             - d.quantities['elevation'].centroid_values)
        assert h.min() > -1e-10, '%s: negative depth %.3e' % (nm, h.min())


# The draining case: a flat bed with the water already moving toward an open
# right wall, so a large fraction of the tracer is flushed out. The point is to
# move enough that a boundary integral which silently read zero could not still
# look "close enough" -- a gently draining domain loses a fraction of a percent
# and would not discriminate.
DRAIN_FINALTIME = 60.0
DRAIN_DEPTH = 2.0
DRAIN_SPEED = 2.0


def _drain_wedge(d):
    """A tracer field that varies in space, so a permutation cannot hide in it."""
    return 0.2 + 0.8 * d.centroid_coordinates[:, 0] / LEN


def _build_draining(mode, field=None):
    """Same problem, but the right wall is open so tracer leaves the domain.

    `field` gives a spatially varying initial concentration; the default is
    uniform, which is what the budget tests want (the total is what matters)
    but useless for anything comparing cell by cell.
    """
    d = rectangular_cross_domain(NXY, NXY, len1=LEN, len2=LEN)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', DRAIN_DEPTH)
    d.set_quantity('xmomentum', DRAIN_DEPTH * DRAIN_SPEED)
    d.set_quantity('ymomentum', 0.0)
    Br = Reflective_boundary(d)
    Bt = anuga.Transmissive_boundary(d)
    d.set_boundary({'left': Br, 'top': Br, 'bottom': Br, 'right': Bt})
    d.add_tracer('c', beta=1.0)
    d.set_tracer('c', 1.0 if field is None else field(d))
    d.set_multiprocessor_mode(mode)
    return d


@pytest.fixture(scope='module')
def drained():
    if not gpu_available():
        pytest.skip('GPU OpenMP interface not available: %s' % _gpu_error)

    cpu = _build_draining(1)
    cpu.evolve_to_end(finaltime=DRAIN_FINALTIME)

    gpu = _build_draining(2)
    if getattr(gpu, 'multiprocessor_mode', None) != 2:
        pytest.fail('mode 2 did not engage: multiprocessor_mode=%r'
                    % getattr(gpu, 'multiprocessor_mode', None))
    from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device
    sync_to_device(gpu.gpu_interface.gpu_dom)
    gpu.evolve_to_end(finaltime=DRAIN_FINALTIME)
    return cpu, gpu


def test_the_tracer_budget_balances_on_the_device(drained):
    """The mode-2 boundary accounting must be real, not a silent zero.

    `tracer_boundary_flux` is written inside the flux kernel's `omp target`
    region, so it is subject to the same trap the module docstring describes:
    load its base pointer off `D` inside the region and the accumulation
    silently does nothing, leaving an integral of exactly 0.0 while the mass
    still changes. Asserting the budget balances catches that.
    """
    for nm, d in zip(('mode 1', 'mode 2'), drained):
        change, flux, disc = d.check_tracer_conservation('c')
        assert change < -0.05 * d.get_tracer_mass('c'), \
            '%s: too little tracer left for this to be a real test' % nm
        assert flux != 0.0, '%s: boundary flux integral is exactly zero' % nm
        assert abs(disc) < 1e-10 * abs(change), \
            '%s: budget does not balance (change %g, flux %g)' % (nm, change, flux)


def test_the_two_modes_agree_on_the_budget(drained):
    cpu, gpu = drained
    c1, f1, _ = cpu.check_tracer_conservation('c')
    c2, f2, _ = gpu.check_tracer_conservation('c')
    assert abs(c1 - c2) < 1e-10 * max(abs(c1), 1.0), 'mass change disagrees'
    assert abs(f1 - f2) < 1e-10 * max(abs(f1), 1.0), 'flux integral disagrees'


def _build_reordered(mode, seed=3):
    """As _build_draining, but with the triangles renumbered first (#277).

    reorder() permutes the tracer blocks in place; the device mapping happens
    afterwards, at evolve time, so mode 2 must see the permuted data.
    """
    # A varying field, and mode set below -- after the reorder.
    d = _build_draining(mode=1, field=_drain_wedge)
    rng = np.random.default_rng(seed)
    d.reorder(rng.permutation(len(d)))
    d.set_multiprocessor_mode(mode)
    return d


def _sorted_by_centroid(d):
    xy = d.centroid_coordinates
    return d.get_tracer('c')[np.lexsort((xy[:, 1], xy[:, 0]))]


def test_a_reordered_domain_gives_the_same_answer_on_the_device():
    if not gpu_available():
        pytest.skip('GPU OpenMP interface not available: %s' % _gpu_error)

    ref = _build_draining(1, field=_drain_wedge)
    ref.evolve_to_end(finaltime=DRAIN_FINALTIME)

    gpu = _build_reordered(2)
    if getattr(gpu, 'multiprocessor_mode', None) != 2:
        pytest.fail('mode 2 did not engage: multiprocessor_mode=%r'
                    % getattr(gpu, 'multiprocessor_mode', None))
    from anuga.shallow_water.sw_domain_gpu_ext import sync_to_device
    sync_to_device(gpu.gpu_interface.gpu_dom)
    gpu.evolve_to_end(finaltime=DRAIN_FINALTIME)

    a = _sorted_by_centroid(ref)
    b = _sorted_by_centroid(gpu)
    assert np.ptp(a) > 1e-6, 'field too flat to discriminate'
    assert np.abs(a - b).max() < 1e-10, \
        'reordering changed the answer on the device'


def _build_for_max(mode, seed=5):
    """A plume that washes out, so the running max differs from the final field."""
    d = rectangular_cross_domain(NXY, NXY, len1=LEN, len2=LEN)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', DRAIN_DEPTH)
    d.set_quantity('xmomentum', DRAIN_DEPTH * DRAIN_SPEED)
    Br = Reflective_boundary(d)
    Bt = anuga.Transmissive_boundary(d)
    d.set_boundary({'left': anuga.Dirichlet_boundary(
                        [DRAIN_DEPTH, DRAIN_DEPTH * DRAIN_SPEED, 0.0]),
                    'top': Br, 'bottom': Br, 'right': Bt})
    d.add_tracer('c', beta=1.0)
    x = d.centroid_coordinates[:, 0]
    d.set_tracer('c', np.where(x < LEN / 4, 1.0, 0.0))
    d.set_tracer_boundary('c', 'left', 0.0)
    from anuga.operators.collect_max_quantities_operator import (
        Collect_max_quantities_operator)
    op = Collect_max_quantities_operator(d)
    d.set_multiprocessor_mode(mode)
    return d, op


def test_the_tracer_maxima_agree_between_the_modes():
    """The max kernel derives c = m/h on the device, as the host path does.

    Reading tracer_centroid_values there would be doubly wrong: stale by a
    step, and -- since D is not mapped inside an omp target region -- a host
    address on the device.
    """
    if not gpu_available():
        pytest.skip('GPU OpenMP interface not available: %s' % _gpu_error)

    cpu, cpu_op = _build_for_max(1)
    cpu.evolve_to_end(finaltime=DRAIN_FINALTIME)

    gpu, gpu_op = _build_for_max(2)
    if getattr(gpu, 'multiprocessor_mode', None) != 2:
        pytest.fail('mode 2 did not engage: multiprocessor_mode=%r'
                    % getattr(gpu, 'multiprocessor_mode', None))
    from anuga.shallow_water.sw_domain_gpu_ext import (sync_to_device,
                                                       get_max_tracers_gpu)
    sync_to_device(gpu.gpu_interface.gpu_dom)
    gpu.evolve_to_end(finaltime=DRAIN_FINALTIME)
    get_max_tracers_gpu(gpu.gpu_interface.gpu_dom, gpu_op.max_tracer)

    a, b = cpu_op.max_tracer[0], gpu_op.max_tracer[0]
    assert np.ptp(a) > 1e-6, 'the maximum field is flat; this proves nothing'
    assert np.abs(a - b).max() < 1e-10, \
        'the modes disagree on the tracer maximum by %g' % np.abs(a - b).max()


def test_the_device_maximum_is_never_below_the_final_field():
    if not gpu_available():
        pytest.skip('GPU OpenMP interface not available: %s' % _gpu_error)

    gpu, op = _build_for_max(2)
    from anuga.shallow_water.sw_domain_gpu_ext import (sync_to_device,
                                                       get_max_tracers_gpu)
    sync_to_device(gpu.gpu_interface.gpu_dom)
    gpu.evolve_to_end(finaltime=DRAIN_FINALTIME)
    get_max_tracers_gpu(gpu.gpu_interface.gpu_dom, op.max_tracer)
    assert np.all(op.max_tracer[0] >= gpu.get_tracer('c') - 1e-12)
