"""
Unit tests for the three parallel CG primitives in Kinematic_viscosity_operator,
designed to be run inside MPI processes via mpiexec.

Exit code 0 = all assertions passed; non-zero = at least one failure.

Usage:
  python run_parallel_kv_unit_tests.py               # serial (numprocs == 1)
  mpiexec -np 3 python run_parallel_kv_unit_tests.py # parallel
"""

import sys
import numpy as num
import anuga
from anuga import rectangular_cross_domain, Reflective_boundary
from anuga import distribute, myid, numprocs, finalize
from anuga.operators.kinematic_viscosity_operator import Kinematic_viscosity_operator
from anuga.utilities.parallel_abstraction import global_except_hook

sys.excepthook = global_except_hook

failures = []


def check(cond, msg):
    if not cond:
        failures.append(f'[rank {myid}] FAIL: {msg}')


# -----------------------------------------------------------------------
# Build and distribute a small domain
# -----------------------------------------------------------------------
domain = rectangular_cross_domain(6, 6)
domain.set_quantity('elevation', 0.0)
domain.set_quantity('stage', expression='1.0 + 0.5*x')
domain.set_quantity('xmomentum', expression='x - x*x')
domain.set_quantity('ymomentum', expression='y - y*y')

if numprocs > 1:
    domain = distribute(domain)

Br = Reflective_boundary(domain)
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

kv = Kinematic_viscosity_operator(domain, diffusivity='height',
                                   use_triangle_areas=True, verbose=False)
kv.dt = 0.05   # set a nominal timestep so matrices are meaningful


# -----------------------------------------------------------------------
# Test 1: _exchange_ghost_vector
#
# Fill v with the *global* triangle index of each local triangle.
# After exchange, every ghost slot must contain its own global index
# (sent by whichever rank owns that triangle).
# -----------------------------------------------------------------------
if numprocs > 1:
    v = num.zeros(kv.n, dtype=float)
    for i in range(kv.n_full):
        v[i] = float(domain.tri_l2g[i])

    kv._exchange_ghost_vector(v)

    for i in range(kv.n_full, kv.n):
        expected = float(domain.tri_l2g[i])
        check(abs(v[i] - expected) < 0.5,
              f'_exchange_ghost_vector: ghost local={i} global={int(expected)} '
              f'got {v[i]:.1f}')
else:
    # Serial: exchange is a no-op — just confirm the method exists and returns
    v = num.ones(kv.n, dtype=float)
    kv._exchange_ghost_vector(v)
    check(num.all(v == 1.0), '_exchange_ghost_vector changed values in serial mode')


# -----------------------------------------------------------------------
# Test 2: _distributed_dot
#
# A vector of ones on the n_full owned triangles dotted with itself equals
# n_full on each rank; the global Allreduce sum should equal the total
# number of triangles in the full mesh.
# -----------------------------------------------------------------------
v_ones = num.ones(kv.n_full, dtype=float)
got = kv._distributed_dot(v_ones, v_ones)

if numprocs > 1:
    import anuga.utilities.parallel_abstraction as pypar
    from mpi4py.MPI import SUM
    n_full_arr  = num.array([kv.n_full], dtype=num.int64)
    n_full_total = num.zeros(1, dtype=num.int64)
    pypar.comm.Allreduce(n_full_arr, n_full_total, op=SUM)
    expected_dot = float(n_full_total[0])
else:
    expected_dot = float(kv.n_full)

check(abs(got - expected_dot) < 1e-10,
      f'_distributed_dot: expected {expected_dot}, got {got}')


# -----------------------------------------------------------------------
# Test 3: _parabolic_matvec_distributed with dt=0 is the identity
#
# When dt == 0 the operator P = I - 0*A = I, so P*d == d for any d.
# -----------------------------------------------------------------------
kv.dt = 0.0
kv.update_elliptic_matrix(kv.diffusivity)
kv.update_elliptic_boundary_term(
    domain.quantities['xvelocity'])

rng = num.random.default_rng(42 + myid)
d_full = rng.standard_normal(kv.n_full)
Pd = kv._parabolic_matvec_distributed(d_full)

check(num.allclose(Pd, d_full, rtol=1e-12, atol=1e-12),
      f'_parabolic_matvec_distributed with dt=0: max diff = '
      f'{num.max(num.abs(Pd - d_full)):.3e}')


# -----------------------------------------------------------------------
# Test 4: _parabolic_solve_distributed convergence self-check
#
# Choose rhs = P * x_true for a known x_true and verify that the CG
# recovers x_true.  Uses a small dt (0.01) so P is well-conditioned.
# -----------------------------------------------------------------------
kv.dt = 0.01
kv.update_elliptic_matrix(kv.diffusivity)
kv.update_elliptic_boundary_term(domain.quantities['xvelocity'])

rng2 = num.random.default_rng(7)
x_true = rng2.standard_normal(kv.n_full)
rhs = kv._parabolic_matvec_distributed(x_true)

x_sol = kv._parabolic_solve_distributed(rhs, num.zeros(kv.n_full), 5000, 1e-8, 1e-8)

check(num.allclose(x_sol, x_true, rtol=1e-4, atol=1e-6),
      f'_parabolic_solve_distributed self-check: max diff = '
      f'{num.max(num.abs(x_sol - x_true)):.3e}')


# -----------------------------------------------------------------------
# Report and exit
# -----------------------------------------------------------------------
finalize()

if failures:
    for msg in failures:
        print(msg, file=sys.stderr)
    sys.exit(1)
