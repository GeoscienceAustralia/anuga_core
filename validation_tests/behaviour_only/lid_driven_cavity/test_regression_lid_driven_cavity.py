"""Regression test for lid_driven_cavity: checks numerical outputs against stored baseline."""

import pytest
from conftest import ensure_sww, sww_metrics

pytestmark = pytest.mark.slow

_SWW = 'dimensional_lid_driven.sww'


def test_lid_driven_cavity_regression(simdir, num_regression):
    ensure_sww(_SWW, 'numerical_lid_driven_cavity.py')
    metrics = sww_metrics(_SWW)
    num_regression.check(metrics, default_tolerance=dict(atol=1e-3, rtol=1e-3))
