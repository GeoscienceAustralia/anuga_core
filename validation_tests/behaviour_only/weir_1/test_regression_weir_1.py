"""Regression test for weir_1: checks numerical outputs against stored baseline."""

import pytest
from conftest import ensure_sww, sww_metrics

pytestmark = pytest.mark.slow

_SWW = 'runup_riverwall.sww'


def test_weir_1_regression(simdir, num_regression):
    ensure_sww(_SWW, 'runup.py')
    metrics = sww_metrics(_SWW)
    num_regression.check(metrics, default_tolerance=dict(atol=1e-3, rtol=1e-3))
