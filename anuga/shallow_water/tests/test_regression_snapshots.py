"""
Numerical regression snapshots for core ANUGA operations.

These tests capture the exact numerical output of key computations
(evolve, compute_fluxes, extrapolate) so that any unintended change
in numerical results is caught immediately.

Run with --force-regen to regenerate baselines after an intentional change:
    pytest anuga/shallow_water/tests/test_regression_snapshots.py --force-regen
"""

import numpy as np
import pytest

import anuga


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_dam_break_domain(algorithm='DE0', friction=0.0, cells_x=10, cells_y=5):
    """Simple dam-break domain: left half stage=1, right half stage=0."""
    domain = anuga.rectangular_cross_domain(
        cells_x, cells_y, len1=10.0, len2=5.0)
    # These baselines are legacy-recorded numerical references; mode-2 ('unified')
    # differs at the ~1e-6 level from a different reduction/eval order. Pin legacy
    # so the snapshots are deterministic under any ANUGA_DEFAULT_COMPUTE_MODE
    # (mode-2 fidelity is covered by the mode1-vs-mode2 tests in test_DE_gpu_omp.py).
    domain.set_compute_mode('legacy')
    domain.set_name('regression_db')
    domain.set_store(False)
    domain.set_flow_algorithm(algorithm)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('friction', friction)
    domain.set_quantity('stage', lambda x, y: np.where(x < 5.0, 1.0, 0.0))
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    return domain


def _make_thacker_domain(algorithm='DE0'):
    """Planar surface in a parabolic bowl — has an analytical solution."""
    a = 3000.0     # bowl radius (m)
    h0 = 10.0      # central depth (m)
    eta0 = 2.0     # initial surface tilt
    g = 9.8

    def bed(x, y):
        return h0 * ((x / a) ** 2 + (y / a) ** 2 - 1.0)

    def stage(x, y):
        return eta0 * (2 * x / a) - h0 * ((x / a) ** 2 + (y / a) ** 2 - 1.0)

    domain = anuga.rectangular_cross_domain(
        20, 20, len1=2 * a, len2=2 * a,
        origin=(-a, -a))
    # Legacy-recorded baseline — pin legacy (see _make_dam_break_domain).
    domain.set_compute_mode('legacy')
    domain.set_name('regression_thacker')
    domain.set_store(False)
    domain.set_flow_algorithm(algorithm)
    domain.set_quantity('elevation', bed)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', stage)
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    return domain


# ---------------------------------------------------------------------------
# Evolve snapshots
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_dam_break_DE0_stage_snapshot(num_regression):
    """Dam break (DE0, no friction) — stage centroids after 1 s."""
    domain = _make_dam_break_domain(algorithm='DE0')
    for _ in domain.evolve(yieldstep=0.5, duration=1.0):
        pass
    stage = domain.get_quantity('stage').centroid_values
    num_regression.check({'stage': stage}, default_tolerance=dict(atol=1e-6, rtol=1e-5))


@pytest.mark.slow
def test_dam_break_DE1_stage_snapshot(num_regression):
    """Dam break (DE1, no friction) — stage centroids after 1 s."""
    domain = _make_dam_break_domain(algorithm='DE1')
    for _ in domain.evolve(yieldstep=0.5, duration=1.0):
        pass
    stage = domain.get_quantity('stage').centroid_values
    num_regression.check({'stage': stage}, default_tolerance=dict(atol=1e-6, rtol=1e-5))


@pytest.mark.slow
def test_dam_break_with_friction_snapshot(num_regression):
    """Dam break (DE0, Manning n=0.03) — stage and xmomentum after 1 s."""
    domain = _make_dam_break_domain(algorithm='DE0', friction=0.03)
    for _ in domain.evolve(yieldstep=0.5, duration=1.0):
        pass
    num_regression.check({
        'stage': domain.get_quantity('stage').centroid_values,
        'xmomentum': domain.get_quantity('xmomentum').centroid_values,
    }, default_tolerance=dict(atol=1e-6, rtol=1e-5))


@pytest.mark.slow
def test_thacker_planar_surface_snapshot(num_regression):
    """Thacker planar surface in a bowl — stage after 0.5 s (no friction)."""
    domain = _make_thacker_domain(algorithm='DE0')
    for _ in domain.evolve(yieldstep=0.5, duration=0.5):
        pass
    stage = domain.get_quantity('stage').centroid_values
    num_regression.check({'stage': stage}, default_tolerance=dict(atol=1e-4, rtol=1e-4))


# ---------------------------------------------------------------------------
# Extrapolation snapshot
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_extrapolation_edge_values_snapshot(num_regression):
    """Verify edge values produced by second-order extrapolation are stable."""
    domain = _make_dam_break_domain(algorithm='DE0', cells_x=6, cells_y=3)
    # Advance one step to get past initialisation
    for _ in domain.evolve(yieldstep=0.1, duration=0.1):
        pass
    domain.distribute_to_vertices_and_edges()
    num_regression.check({
        'stage_edge':     domain.get_quantity('stage').edge_values.ravel(),
        'xmom_edge':      domain.get_quantity('xmomentum').edge_values.ravel(),
    }, default_tolerance=dict(atol=1e-8, rtol=1e-6))


# ---------------------------------------------------------------------------
# Timestep / CFL snapshot
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_timestep_sequence_snapshot(num_regression):
    """Record the sequence of adaptive timesteps for a dam break (DE0).

    A change here indicates the CFL condition or flux computation has changed.
    """
    domain = _make_dam_break_domain(algorithm='DE0')
    timesteps = []
    for t in domain.evolve(yieldstep=0.25, duration=1.0):
        timesteps.append(float(domain.get_timestep()))
    num_regression.check(
        {'timesteps': np.array(timesteps)},
        default_tolerance=dict(atol=1e-8, rtol=1e-6))
