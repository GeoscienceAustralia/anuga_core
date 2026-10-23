"""Parallel_Inlet_operator must contribute its volume to the mass balance ONCE.

Regression guard for issue #193.

``domain.fractional_step_volume_integral`` is a per-rank LOCAL accumulator that
``Domain.get_fractional_step_volume_integral()`` sums across ranks with an MPI
allreduce.  The parallel inlet's ``volume`` is a GLOBAL quantity — the master computes
it for the whole inlet and broadcasts it to the other participants — so every rank in
``self.procs`` used to add the full value and the allreduce multiplied the inlet's
contribution by the number of participating ranks:

    np=1 -> 40.0 (correct)    np=2 -> 80.0 (2x)    np=4 -> 160.0 (4x)

against a true volume of Q*dt = 40.0.  That silently corrupted the
``Water_volume_statistics`` mass balance in every parallel run with an inlet, in BOTH
compute modes (this was never a GPU-specific bug).

Cannot be caught serially: at np=1 there is only the master, so the bug is invisible.
Hence a real MPI test.
"""

import os
import sys
import unittest

import numpy as num
import pytest

import anuga
from anuga import Reflective_boundary, rectangular_cross_domain
from anuga import Inlet_operator

try:
    from anuga.parallel import distribute, myid, numprocs, finalize
    from mpi4py import MPI
except ImportError:
    pass


Q = 20.0            # m^3/s prescribed inflow
TIMESTEP = 2.0      # s
EXPECTED = Q * TIMESTEP     # the inlet adds exactly this much water, once


def run_inlet_fsvi(compute_mode):
    """Apply the inlet once and return the globally-summed volume integral."""

    if myid == 0:
        domain = rectangular_cross_domain(30, 30, len1=100.0, len2=100.0)
        domain.set_flow_algorithm('DE0')
        domain.set_quantity('elevation', -10.0)
        domain.set_quantity('stage', 1.0)
    else:
        domain = None

    domain = distribute(domain, verbose=False)
    domain.set_multiprocessor_mode(compute_mode)
    domain.set_name('parallel_inlet_fsvi')
    domain.store = False

    Br = Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    domain.distribute_to_vertices_and_edges()

    # Inlet line spanning the middle of the domain, so it deliberately straddles
    # partition boundaries and is shared by several ranks — that is the case the bug
    # needed.
    inlet_line = [[20.0, 50.0], [80.0, 50.0]]
    op = Inlet_operator(domain, inlet_line, Q, label='fsvi_inlet', verbose=False)

    domain.fractional_step_volume_integral = 0.0
    domain.timestep = TIMESTEP
    if op is not None:
        op()

    # This is what Water_volume_statistics reports.
    return domain.get_fractional_step_volume_integral()


def main(verbose=False):
    for compute_mode in (1, 2):
        fsvi = run_inlet_fsvi(compute_mode)

        if myid == 0 and verbose:
            print(f'mode {compute_mode}, np={numprocs}: fsvi={fsvi} expected={EXPECTED}')

        msg = (f'mode {compute_mode}, np={numprocs}: inlet contributed {fsvi} to the '
               f'fractional-step volume integral, expected {EXPECTED}. A multiple of '
               f'the expected value means each participating rank added the global '
               f'volume instead of the master alone (issue #193).')
        assert num.allclose(fsvi, EXPECTED, rtol=0, atol=1.0e-10), msg


@pytest.mark.skipif('mpi4py' not in sys.modules,
                    reason='requires the mpi4py module')
class Test_parallel_inlet_fsvi(unittest.TestCase):
    def test_parallel_inlet_fsvi(self):
        # np=3: enough ranks to share the inlet, so the x(nranks) bug is exposed.
        cmd = anuga.mpicmd(os.path.abspath(__file__), numprocs=3)
        result = os.system(cmd)
        assert result == 0, 'Parallel inlet volume-integral test failed'


if __name__ == '__main__':
    main(verbose=True)
    finalize()
