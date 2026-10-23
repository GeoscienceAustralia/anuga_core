"""Structures carry tracers with the water they move (inlet_tracers.py).

Two basins separated by a ridge, joined only by a structure. Before this, a
culvert moved water and no tracer: the dye piled up upstream and clean water
arrived downstream. Now the cells that lose water keep their concentration,
and the tracer they give up arrives with the water, in proportion to each
downstream cell's gain.
"""
import numpy as np
import pytest

import anuga
from anuga import Reflective_boundary, rectangular_cross_domain

L, W = 40.0, 10.0


def two_basins(mode='legacy', c_up=0.02, c_down=0.0, stage_up=2.0, stage_down=0.5):
    d = rectangular_cross_domain(40, 10, len1=L, len2=W)
    d.set_flow_algorithm('DE1')
    d.set_compute_mode(mode)
    d.store = False
    ridge = lambda x, y: np.where(np.abs(x - L / 2) < 2.0, 5.0, 0.0)
    d.set_quantity('elevation', ridge)
    # never below the bed: a negative depth on the ridge would be topped up by
    # the solver and spoil the water budget the tests check against
    d.set_quantity('stage', lambda x, y: np.maximum(
        np.where(x < L / 2, stage_up, stage_down), ridge(x, y)))
    d.set_quantity('friction', 0.01)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.add_tracer('dye')
    x = d.centroid_coordinates[:, 0]
    d.set_tracer('dye', np.where(x < L / 2, c_up, c_down))
    return d


def culvert(d, kind='boyd_box'):
    ends = [[L / 2 - 5.0, W / 2], [L / 2 + 5.0, W / 2]]
    if kind == 'boyd_box':
        return anuga.Boyd_box_operator(d, losses=1.5, width=1.0, height=1.0,
                                       end_points=ends, manning=0.013)
    if kind == 'boyd_pipe':
        return anuga.Boyd_pipe_operator(d, losses=1.5, diameter=1.0,
                                        end_points=ends, manning=0.013)
    if kind == 'weir':
        return anuga.Weir_orifice_trapezoid_operator(
            d, losses=1.5, width=1.0, height=1.0, z1=2.0, z2=2.0,
            end_points=ends, manning=0.013)
    raise ValueError(kind)


def mass(d, side=None):
    m = d.tracer_conserved_values[0] * d.areas
    if side is None:
        return float(m.sum())
    x = d.centroid_coordinates[:, 0]
    return float(m[(x < L / 2) if side == 'up' else (x > L / 2)].sum())


def water(d, side):
    h = (d.quantities['stage'].centroid_values
         - d.quantities['elevation'].centroid_values) * d.areas
    x = d.centroid_coordinates[:, 0]
    return float(h[(x < L / 2) if side == 'up' else (x > L / 2)].sum())


@pytest.mark.parametrize('kind', ['boyd_box', 'boyd_pipe', 'weir'])
@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_the_dye_goes_through_with_the_water(kind, mode):
    d = two_basins(mode)
    op = culvert(d, kind)
    m_up0, w_down0 = mass(d, 'up'), water(d, 'down')
    m0 = mass(d)
    d.evolve_to_end(finaltime=30.0)
    moved_water = water(d, 'down') - w_down0
    assert moved_water > 1.0, 'the structure moved no water'
    # conserved exactly, and what left upstream arrived downstream
    assert mass(d) == pytest.approx(m0, rel=1e-12)
    assert mass(d, 'down') == pytest.approx(m_up0 - mass(d, 'up'), rel=1e-9)
    assert mass(d, 'down') > 0.0
    # it arrives at the upstream concentration: the inflow cells are at c_up
    # and keep it. (Not exactly: if the flow briefly reverses late in the run,
    # water goes back at the lower downstream concentration, which lifts the
    # net tracer per net volume slightly above c_up.)
    assert mass(d, 'down') / moved_water == pytest.approx(0.02, rel=5e-3)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_a_uniform_tracer_stays_uniform(mode):
    """Water moved at concentration c into water at c is still at c: the
    rule neither concentrates nor dilutes."""
    d = two_basins(mode, c_up=0.01, c_down=0.01)
    culvert(d)
    d.evolve_to_end(finaltime=30.0)
    h = d.quantities['stage'].centroid_values - d.quantities['elevation'].centroid_values
    wet = h > 0.01
    assert np.allclose(d.tracer_centroid_values[0][wet], 0.01, rtol=1e-8)


def test_no_transfer_leaves_the_tracer_alone(monkeypatch):
    """A fully blocked culvert moves no water, so the tracer rule must be an
    exact no-op -- even with a concentration front running through the
    inflow region, where mixing the region as a pool would show."""
    def run(rule_on):
        d = two_basins('legacy')
        x = d.centroid_coordinates[:, 0]
        d.set_tracer('dye', np.where(x < 15.0, 0.02, 0.0))   # front at the inlet
        op = anuga.Boyd_box_operator(d, losses=1.5, width=1.0, height=1.0,
                                     blockage=1.0, manning=0.013,
                                     end_points=[[L / 2 - 5.0, W / 2], [L / 2 + 5.0, W / 2]])
        if not rule_on:
            from anuga.structures.inlet_tracers import StructureTracers
            monkeypatch.setattr(StructureTracers, 'apply', lambda self, state, global_sum=None: None)
        d.evolve_to_end(finaltime=10.0)
        monkeypatch.undo()
        return op, d.tracer_conserved_values.copy()
    op, with_rule = run(True)
    assert op.accumulated_flow == 0.0
    _, without = run(False)
    assert np.array_equal(with_rule, without)


def test_the_operator_reports_what_it_moved():
    d = two_basins('legacy')
    op = culvert(d)
    m_down0 = mass(d, 'down')
    d.evolve_to_end(finaltime=30.0)
    assert op._structure_tracers.total_moved['dye'] == pytest.approx(
        mass(d, 'down') - m_down0, rel=1e-9)


def test_the_two_compute_modes_agree():
    """Mode 1 against mode 2. The two compute modes do not run bit-identical
    arithmetic on every build (the GPU, and the Windows compilers, differ in
    the last bits from Linux gcc), and the culvert's feedback loop grows that
    to ~1e-4 of tracer mass per unit area in 30 s, with or without tracers. So
    the check is exact conservation in both modes and a field that agrees to
    the order the flow itself does."""
    out, totals = [], []
    for mode in ('legacy', 'unified'):
        d = two_basins(mode)
        culvert(d)
        m0 = mass(d)
        d.evolve_to_end(finaltime=30.0)
        out.append(d.tracer_conserved_values[0].copy())
        totals.append((m0, mass(d)))
    for m0, m1 in totals:
        assert m1 == pytest.approx(m0, rel=1e-12)
    assert np.abs(out[0] - out[1]).max() < 1e-3
