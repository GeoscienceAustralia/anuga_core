"""Regression test for lateral_weir_hecras: sww metrics + correlation with HEC-RAS gauges.

NOTE ON EXPECTED LOW CORRELATION
---------------------------------
The HEC-RAS correlation for this case is intentionally low (mean ~0.72, min ~-0.34).
This is a documented, known physics difference — NOT a bug:

  ANUGA transfers both mass AND momentum over the lateral riverwall.
  HEC-RAS transfers mass only, leaving the spilled water with no downstream
  momentum component.

The consequence is that HEC-RAS predicts a much stronger backwater effect
immediately upstream of the weirs (stations y=700-800m), while ANUGA shows
little backwater there because the spilled water carries its downstream velocity
into the side channels.  The authors confirmed this by running ANUGA with
momentum transfer over the weir forced to zero: the two models then agreed well.

A future code change that raises this correlation toward 1.0 would be a red flag
— it would indicate that ANUGA has stopped transferring momentum through riverwalls.

Reference: validation_tests/behaviour_only/lateral_weir_hecras/results.tex
"""

import pytest
from anuga.utilities import plot_utils as util
from conftest import ensure_sww, sww_metrics, hecras_correlation_metrics

pytestmark = pytest.mark.slow

_SWW    = 'channel_floodplain1.sww'
_GAUGES = 'hecras_riverwall_anugaTest/gauges.csv'

_LOW_CORRELATION_NOTE = """
  NOTE: Low HEC-RAS correlation is EXPECTED for this case.
  ANUGA transfers mass+momentum over the riverwall; HEC-RAS transfers mass only.
  This causes a larger backwater in HEC-RAS upstream of the weir (y=700-800m).
  See module docstring and results.tex for full explanation.
"""


def test_lateral_weir_hecras_regression(simdir, num_regression):
    ensure_sww(_SWW, 'channel_floodplain1.py')

    metrics = sww_metrics(_SWW)
    p  = util.get_output(_SWW, 0.001)
    pc = util.get_centroids(p, velocity_extrapolation=True)
    metrics.update(hecras_correlation_metrics(pc, _GAUGES, x_channel=20.0))

    print(_LOW_CORRELATION_NOTE)
    print(f"  mean_hecras_correlation = {metrics['mean_hecras_correlation'][0]:.3f}  (baseline ~0.72)")
    print(f"  min_hecras_correlation  = {metrics['min_hecras_correlation'][0]:.3f}  (baseline ~-0.34)")
    print(f"  mean_hecras_rms_error   = {metrics['mean_hecras_rms_error'][0]:.3f} m")

    num_regression.check(metrics, default_tolerance=dict(atol=1e-3, rtol=1e-3))
