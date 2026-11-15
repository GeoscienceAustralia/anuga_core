"""Tests of the rainfall loss models and of ``losses=`` on the rain operators."""

import numpy as np
import pytest

import anuga
from anuga.rain import (ARR_rate_operator, Hyetograph_rate_operator,
                        Initial_continuing_loss, Initial_proportional_loss,
                        Raster_rate_operator, SCS_curve_number)

MM = 1.0e-3


def flat_domain(mode=None):
    """A 40 m x 20 m flat, dry, walled domain."""
    domain = anuga.rectangular_cross_domain(4, 2, len1=40.0, len2=20.0)
    if mode is not None:
        domain.set_compute_mode(mode)
    domain.set_store(False)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 0.0)
    domain.set_boundary({'left': anuga.Reflective_boundary(domain),
                         'right': anuga.Reflective_boundary(domain),
                         'top': anuga.Reflective_boundary(domain),
                         'bottom': anuga.Reflective_boundary(domain)})
    return domain


class Storm:
    """Uniform rain of 10, 20 and 20 mm in three 10-minute slices."""
    time_step = 600.0
    times = np.array([0.0, 600.0, 1200.0])
    depths = np.array([10.0, 20.0, 20.0]) * MM

    def extract_data_at_locations(self, locations):
        return np.tile(self.depths[:, None], (1, len(locations)))


class Pattern:
    """50 mm as 10, 20 and 20 mm in three 10-minute steps."""
    Tstep = 10
    Ev_dur = 30
    Ev_dep = 50.0
    Tplot = np.array([0.0, 10.0, 20.0, 30.0])
    Rplot = np.array([0.0, 10.0, 20.0, 20.0])


# CN = 80: S = 25400/80 - 254 = 63.5 mm, Ia = 0.2 S = 12.7 mm.
# Runoff after 10, 30 and 50 mm of rain: 0, 17.3^2/80.8 and 37.3^2/100.8 mm.
SCS_RUNOFF = np.array([0.0, 17.3**2/80.8, 37.3**2/100.8]) * MM
SCS_EXCESS = np.diff(np.concatenate([[0.0], SCS_RUNOFF]))

# IL = 15 mm, CL = 2.5 mm/h. Slice 1 (10 mm) all goes to the initial loss.
# Slice 2 (20 mm): 5 mm fills it, and the last 3/4 of the slice loses
# 2.5 * (1/6 h) * 3/4 = 0.3125 mm, leaving 14.6875 mm.
# Slice 3 (20 mm): 20 - 2.5/6 mm.
ILCL_EXCESS = np.array([0.0, 14.6875, 20.0 - 2.5/6.0]) * MM


def excesses(losses, domain):
    losses.bind(domain)
    return np.array([losses.excess(np.full(domain.number_of_elements, d), 600.0)
                     for d in Storm.depths])


def test_scs_curve_number_against_hand_calculation():
    domain = flat_domain()
    losses = SCS_curve_number(curve_number=80.0)
    out = excesses(losses, domain)
    assert np.allclose(out, SCS_EXCESS[:, None], rtol=1e-12, atol=0.0)
    assert np.allclose(losses.cumulative_rain, 50.0*MM)
    assert np.allclose(losses.cumulative_loss, 50.0*MM - SCS_RUNOFF[-1])


def test_scs_curve_number_does_not_depend_on_how_the_storm_is_sliced():
    domain = flat_domain()
    one = SCS_curve_number(curve_number=65.0)
    one.bind(domain)
    whole = one.excess(np.full(domain.number_of_elements, 50.0*MM), 1800.0)
    many = SCS_curve_number(curve_number=65.0)
    assert np.allclose(excesses(many, domain).sum(axis=0), whole, rtol=1e-12)


def test_curve_number_100_loses_nothing():
    domain = flat_domain()
    out = excesses(SCS_curve_number(curve_number=100.0), domain)
    assert np.allclose(out, Storm.depths[:, None], rtol=1e-12)


def test_initial_continuing_loss_against_hand_calculation():
    domain = flat_domain()
    losses = Initial_continuing_loss(initial=15.0, continuing=2.5)
    out = excesses(losses, domain)
    assert np.allclose(out, ILCL_EXCESS[:, None], rtol=1e-12, atol=0.0)
    assert np.all(losses.initial_remaining == 0.0)


def test_zero_initial_and_continuing_loss_loses_nothing():
    domain = flat_domain()
    out = excesses(Initial_continuing_loss(initial=0.0, continuing=0.0), domain)
    assert np.allclose(out, Storm.depths[:, None], rtol=1e-12)


def test_rain_below_the_continuing_loss_gives_no_excess_and_no_rain_is_fine():
    domain = flat_domain()
    losses = Initial_continuing_loss(initial=0.0, continuing=6.0)   # 1 mm per slice
    losses.bind(domain)
    N = domain.number_of_elements
    assert np.all(losses.excess(np.full(N, 0.5*MM), 600.0) == 0.0)
    assert np.all(losses.excess(np.zeros(N), 600.0) == 0.0)


def test_parameters_can_vary_from_cell_to_cell():
    domain = flat_domain()
    x = domain.get_centroid_coordinates(absolute=True)[:, 0]
    by_function = SCS_curve_number(curve_number=lambda x, y: np.where(x < 20.0, 100.0, 80.0))
    out = excesses(by_function, domain).sum(axis=0)
    assert np.allclose(out[x < 20.0], 50.0*MM)
    assert np.allclose(out[x >= 20.0], SCS_RUNOFF[-1])

    by_array = SCS_curve_number(curve_number=np.where(x < 20.0, 100.0, 80.0))
    assert np.allclose(excesses(by_array, domain).sum(axis=0), out)


def test_bad_parameters_are_refused():
    domain = flat_domain()
    for cn in (0.0, 101.0, np.nan, np.ones(3)):
        with pytest.raises(ValueError):
            SCS_curve_number(curve_number=cn).bind(domain)
    with pytest.raises(ValueError):
        Initial_continuing_loss(initial=-1.0, continuing=0.0).bind(domain)


def test_a_loss_model_serves_one_operator():
    losses = SCS_curve_number(curve_number=80.0)
    Raster_rate_operator(flat_domain(), Storm(), losses=losses)
    with pytest.raises(RuntimeError):
        Raster_rate_operator(flat_domain(), Storm(), losses=losses)


def run_storm(make_operator, mode):
    domain = flat_domain(mode)
    op = make_operator(domain)
    volumes = [domain.get_water_volume()
               for _ in domain.evolve(yieldstep=600.0, finaltime=1800.0)]
    return op, np.diff(volumes)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
@pytest.mark.parametrize('rain', ['raster', 'arr', 'hyetograph'])
@pytest.mark.parametrize('model', ['scs', 'ilcl'])
def test_operator_adds_the_excess_and_counts_the_loss(model, rain, mode):
    area = 40.0*20.0

    def make_operator(domain):
        if model == 'scs':
            losses = SCS_curve_number(curve_number=80.0)
        else:
            losses = Initial_continuing_loss(initial=15.0, continuing=2.5)
        if rain == 'raster':
            return Raster_rate_operator(domain, Storm(), losses=losses)
        if rain == 'hyetograph':
            return Hyetograph_rate_operator(
                domain, times=[0.0, 600.0, 1200.0, 1800.0],
                depths=[10.0, 20.0, 20.0], losses=losses)
        depth = np.full(domain.number_of_elements, 50.0)
        return ARR_rate_operator(domain, depth, Pattern(), losses=losses)

    op, added = run_storm(make_operator, mode)
    expected = (SCS_EXCESS if model == 'scs' else ILCL_EXCESS)*area

    assert np.allclose(added, expected, rtol=1e-9, atol=1e-12)
    assert op.cumulative_influx == pytest.approx(expected.sum(), rel=1e-9)
    assert op.cumulative_loss == pytest.approx(50.0*MM*area - expected.sum(), rel=1e-9)


def test_without_losses_nothing_is_counted_as_lost():
    op, added = run_storm(lambda domain: Raster_rate_operator(domain, Storm()), 'legacy')
    assert np.allclose(added, Storm.depths*800.0, rtol=1e-9)
    assert op.cumulative_loss == 0.0


def test_loss_is_counted_on_the_operators_cells_only():
    half = [[0.0, 0.0], [20.0, 0.0], [20.0, 20.0], [0.0, 20.0]]
    op, added = run_storm(
        lambda domain: Raster_rate_operator(
            domain, Storm(), polygon=half,
            losses=SCS_curve_number(curve_number=80.0)), 'legacy')
    assert np.allclose(added, SCS_EXCESS*400.0, rtol=1e-9, atol=1e-12)
    assert op.cumulative_loss == pytest.approx((50.0*MM - SCS_RUNOFF[-1])*400.0, rel=1e-9)


def test_a_slice_the_run_steps_over_still_reaches_the_loss_model():
    domain = flat_domain('legacy')
    losses = SCS_curve_number(curve_number=80.0)
    op = Raster_rate_operator(domain, Storm(), losses=losses)
    op._set_slice(2)
    assert np.allclose(losses.cumulative_rain, 50.0*MM)
    assert np.allclose(op.rate, SCS_EXCESS[2]/600.0)


def test_going_back_to_an_earlier_slice_is_refused():
    domain = flat_domain('legacy')
    op = Raster_rate_operator(domain, Storm(),
                              losses=SCS_curve_number(curve_number=80.0))
    op._set_slice(1)
    with pytest.raises(RuntimeError):
        op._set_slice(0)


def test_initial_proportional_loss_against_hand_calculation():
    # IL = 15 mm, 30% lost after that: 0, then 0.7*15 mm, then 0.7*20 mm.
    domain = flat_domain()
    losses = Initial_proportional_loss(proportional=0.3, initial=15.0)
    out = excesses(losses, domain)
    expected = np.array([0.0, 10.5, 14.0])*MM
    assert np.allclose(out, expected[:, None], rtol=1e-12, atol=0.0)


def test_proportional_loss_alone_is_a_runoff_coefficient():
    domain = flat_domain()
    out = excesses(Initial_proportional_loss(proportional=0.25), domain)
    assert np.allclose(out, 0.75*Storm.depths[:, None], rtol=1e-12)
    with pytest.raises(ValueError):
        Initial_proportional_loss(proportional=1.5).bind(flat_domain())


def test_parameters_from_a_polygon_map():
    """The quantity-setting helper gives a land-use map as f(x, y)."""
    from anuga.utilities.quantity_setting_functions import (
        composite_quantity_setting_function)
    domain = flat_domain()
    paved = [[0.0, 0.0], [20.0, 0.0], [20.0, 20.0], [0.0, 20.0]]
    curve_number = composite_quantity_setting_function(
        [[paved, 100.0], ['All', 80.0]], domain, verbose=False)
    out = excesses(SCS_curve_number(curve_number=curve_number), domain).sum(axis=0)
    x = domain.get_centroid_coordinates()[:, 0]
    assert np.allclose(out[x < 20.0], 50.0*MM)
    assert np.allclose(out[x > 20.0], SCS_RUNOFF[-1])


def test_hyetograph_takes_depths_or_intensities_and_uneven_intervals():
    area = 800.0
    # 10 mm in 10 minutes, then 20 mm in 20 minutes (60 mm/h both)
    times = [0.0, 600.0, 1800.0]

    def volume(**rain):
        domain = flat_domain('legacy')
        op = Hyetograph_rate_operator(domain, times, **rain)
        assert op.time_step is None
        for _ in domain.evolve(yieldstep=600.0, finaltime=2400.0):
            pass
        return domain.get_water_volume()

    assert volume(depths=[10.0, 20.0]) == pytest.approx(30.0*MM*area, rel=1e-9)
    assert volume(intensities=[60.0, 60.0]) == pytest.approx(30.0*MM*area, rel=1e-9)


def test_hyetograph_refuses_bad_input():
    domain = flat_domain()
    for kwargs in (dict(times=[0.0, 600.0]),
                   dict(times=[0.0, 600.0], depths=[1.0], intensities=[1.0]),
                   dict(times=[0.0, 600.0], depths=[1.0, 2.0]),
                   dict(times=[600.0, 0.0], depths=[1.0]),
                   dict(times=[0.0, 600.0], depths=[-1.0])):
        with pytest.raises(ValueError):
            Hyetograph_rate_operator(domain, **kwargs)
