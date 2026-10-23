"""Regression test for bridge_hecras: sww metrics + correlation with HEC-RAS gauges."""

import pytest
from anuga.utilities import plot_utils as util
from conftest import ensure_sww, sww_metrics, hecras_correlation_metrics

pytestmark = pytest.mark.slow

_SWW    = 'channel_floodplain1.sww'
_GAUGES = 'hecras_bridge_test/gauges.csv'


def test_bridge_hecras_regression(simdir, num_regression):
    ensure_sww(_SWW, 'channel_floodplain1.py')

    metrics = sww_metrics(_SWW)
    p  = util.get_output(_SWW, 0.001)
    pc = util.get_centroids(p, velocity_extrapolation=True)
    metrics.update(hecras_correlation_metrics(pc, _GAUGES, x_channel=15.0))

    # Stage metrics here are mildly mesh-sensitive: final_max_stage shifts ~9e-3
    # between mesh-generation paths, above the 1e-3 default. Use 1e-2 (combined
    # tolerance ~2e-2) so this behaviour check is robust to that drift while
    # still catching gross regressions.
    num_regression.check(metrics, default_tolerance=dict(atol=1e-2, rtol=1e-2))
