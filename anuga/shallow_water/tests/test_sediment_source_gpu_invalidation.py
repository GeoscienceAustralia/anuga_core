"""An operator may invalidate the GPU interface mid-timestep.

apply_fractional_steps decides `gpu_mode` before running the operators, then
syncs the host-side work of any CPU-only operator back to the device
afterwards. An operator that calls set_tracer_source in between sets
domain.gpu_interface = None -- deliberately, since rebinding the source array
means the C struct and the device mapping must be rebuilt -- so that trailing
sync used to dereference None.

A manufactured-solution source does this on every step, which is how it was
found: the sediment MMS tests failed under ANUGA_DEFAULT_COMPUTE_MODE=unified
with 'NoneType' object has no attribute 'sync_to_device'.
"""

import numpy as num
import pytest

import anuga
from anuga.operators.base_operator import Operator


class _ResetsTheInterface(Operator):
    """A CPU-only fractional-step operator that rebinds the tracer source."""

    def __init__(self, domain):
        Operator.__init__(self, domain)
        self.calls = 0

    def __call__(self):
        self.calls += 1
        n = len(self.domain)
        # Only the FIRST call invalidates the interface (that is when the
        # source array is allocated); calling every step is what a
        # manufactured source does and costs nothing after that.
        self.domain.set_tracer_source('mms', num.full(n, 1.0e-4))

    def parallel_safe(self):
        return True

    def statistics(self):
        return 'resets the interface'

    def timestepping_statistics(self):
        return ''


def _domain(mode):
    d = anuga.rectangular_cross_domain(6, 6, len1=60.0, len2=60.0)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    # A sediment CLASS, not a bare tracer: tracer_external_source is applied
    # by the sediment source kernel, which does nothing with no classes
    # registered -- so a plain tracer would never show the source at all.
    d.add_sediment_fraction(name='mms', diameter=1.0e-4, initial_concentration=0.0)
    d.store = False
    op = _ResetsTheInterface(d)
    d.set_multiprocessor_mode(mode)
    return d, op


@pytest.mark.parametrize('mode', [1, 2])
def test_it_does_not_crash_when_the_interface_is_invalidated(mode):
    """The regression this file is for: it used to raise AttributeError."""
    d, op = _domain(mode)
    for _ in d.evolve(yieldstep=0.5, finaltime=2.0):
        pass
    assert op.calls > 0, 'the operator never ran, so nothing was exercised'


@pytest.mark.parametrize('mode', [1, 2])
def test_the_source_reaches_the_solver(mode):
    """It used to reach it in mode 1 only (#288)."""
    d, _ = _domain(mode)
    for _ in d.evolve(yieldstep=0.5, finaltime=2.0):
        pass
    assert d.get_tracer('mms').max() > 0.0, \
        'the tracer source never took effect'


def test_the_two_modes_agree_on_the_source():
    """The oracle: mode 2 gave exactly 0.0 while mode 1 gave 1.5e-4.

    Three things had to be true at once, and each was a separate defect:
    the source array had to stop being reallocated mid-run (it is now a
    first-class tracer array); it had to be mapped to the device at all;
    and Sediment_transport_operator had to stop running on the device inside
    apply_fractional_steps' host-coherent window, where the trailing
    sync_to_device overwrites whatever the device computed.
    """
    d1, _ = _domain(1)
    for _ in d1.evolve(yieldstep=0.5, finaltime=2.0):
        pass
    d2, _ = _domain(2)
    for _ in d2.evolve(yieldstep=0.5, finaltime=2.0):
        pass

    a = d1.get_tracer('mms')
    b = d2.get_tracer('mms')
    assert a.max() > 0.0, 'the source did nothing even in legacy mode'
    # 1e-8 is the tolerance test_sediment_gpu.py uses for the same
    # mode-1-vs-mode-2 concentration comparison. A one-step staleness would
    # show up as a fraction of the total, i.e. around 1e-5 here, so this is
    # tight enough to catch the failure it exists for.
    assert num.abs(a - b).max() < 1e-8, \
        'the modes disagree by %g' % num.abs(a - b).max()


def test_a_time_varying_source_is_not_stale_on_the_device():
    """set_tracer_source writes the HOST array; the device copy needs pushing.

    A constant source would pass even if the push were missing, since the
    mapped values would happen to be right.
    """
    class _Ramp(_ResetsTheInterface):
        def __call__(self):
            self.calls += 1
            n = len(self.domain)
            self.domain.set_tracer_source('mms',
                                          num.full(n, 1.0e-4 * self.calls))

    out = []
    for mode in (1, 2):
        d = anuga.rectangular_cross_domain(6, 6, len1=60.0, len2=60.0)
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 1.0)
        b = anuga.Reflective_boundary(d)
        d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
        d.add_sediment_fraction(name='mms', diameter=1.0e-4, initial_concentration=0.0)
        d.store = False
        _Ramp(d)
        d.set_multiprocessor_mode(mode)
        for _ in d.evolve(yieldstep=0.5, finaltime=2.0):
            pass
        out.append(d.get_tracer('mms'))

    assert out[0].max() > 0.0
    assert num.abs(out[0] - out[1]).max() < 1e-8, \
        'a ramping source diverges between the modes by %g -- the device copy ' \
        'is stale' % num.abs(out[0] - out[1]).max()
