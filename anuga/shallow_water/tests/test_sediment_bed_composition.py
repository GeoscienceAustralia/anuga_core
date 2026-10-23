"""Bed composition: a per-class active layer over a substrate (Hirano).

Without it every class is entrained from one shared bed, so a fine class keeps
coming out of a coarse bed for as long as the bed lasts (the runaway mud of the
Wax Lake runs). With it, entrainment of a class is scaled by its share of the
active layer, a class cannot leave a cell that does not hold it, and the layers
exchange so that each class is conserved exactly.

The tests: one class reproduces the shared-bed answer; each class is conserved
in a closed box; a class absent from the bed is never entrained; a bed whose
fines are being removed armours; legacy and unified agree; the API refuses the
configurations it cannot honour.
"""
import numpy as np
import pytest

from anuga import Reflective_boundary, rectangular_cross_domain

POROSITY = 0.3
FINE, COARSE = ('fine', 1.0e-4), ('coarse', 5.0e-3)


def build(classes=(('sand', 2.0e-4),), base_depth=0.5, mode='legacy',
          depth=0.6, speed=1.2, n_x=20):
    d = rectangular_cross_domain(n_x, 6, len1=40.0, len2=12.0)
    d.set_compute_mode(mode)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', depth)
    d.set_quantity('xmomentum', speed * depth)
    d.set_quantity('ymomentum', 0.0)
    d.set_quantity('friction', 0.03)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.set_sediment_parameters(porosity=POROSITY)
    for name, diameter in classes:
        d.add_sediment_fraction(name=name, diameter=diameter)
    if base_depth is not None:
        d.set_erodible_base(depth=base_depth)
    return d


def class_totals(d):
    """Per class: suspended + active + substrate, as solid volume [m^3]."""
    a = d.areas
    out = []
    for s in range(d.n_sediment_classes):
        susp = float((d.tracer_conserved_values[s] * a).sum())
        bed = float(((d.sediment_bed_active[s] + d.sediment_bed_substrate[s]) * a).sum())
        out.append(susp + bed)
    return np.array(out)


# ------------------------------------------------------------- equivalence

def test_one_class_reproduces_the_shared_bed():
    """With a single class the share is always 1 and the per-class cap is the
    cell's whole thickness, so the answer must be the shared-bed answer."""
    ref = build()
    ref.evolve_to_end(finaltime=20.0)

    d = build()
    d.set_bed_composition({'sand': 1.0}, active_layer=0.01)
    d.evolve_to_end(finaltime=20.0)

    np.testing.assert_allclose(d.quantities['elevation'].centroid_values,
                               ref.quantities['elevation'].centroid_values,
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(d.tracer_conserved_values[0],
                               ref.tracer_conserved_values[0], rtol=0, atol=1e-12)


def test_the_layers_hold_the_erodible_thickness():
    d = build(classes=(FINE, COARSE), base_depth=0.5)
    d.set_bed_composition({'fine': 0.3, 'coarse': 0.7}, active_layer=0.01)
    held = (d.sediment_bed_active + d.sediment_bed_substrate).sum(axis=0)
    np.testing.assert_allclose(held / (1.0 - POROSITY), 0.5, rtol=1e-12)
    np.testing.assert_allclose(d.sediment_bed_active.sum(axis=0) / (1.0 - POROSITY),
                               0.01, rtol=1e-12)
    F = d.get_bed_composition('active')
    np.testing.assert_allclose(F['fine'], 0.3)


def test_elevation_and_layers_stay_consistent():
    """The Exner update and the layer update are the same exchange, so the
    erodible thickness and the layer content must agree after a run."""
    d = build(classes=(FINE, COARSE), base_depth=0.5)
    d.set_bed_composition({'fine': 0.5, 'coarse': 0.5}, active_layer=0.01)
    d.evolve_to_end(finaltime=20.0)
    held = (d.sediment_bed_active + d.sediment_bed_substrate).sum(axis=0)
    np.testing.assert_allclose(held / (1.0 - POROSITY), d.erodible_thickness(),
                               rtol=0, atol=1e-10)


# ------------------------------------------------------------ conservation

def test_each_class_is_conserved_in_a_closed_box():
    d = build(classes=(FINE, COARSE), base_depth=0.5)
    d.set_bed_composition({'fine': 0.4, 'coarse': 0.6}, active_layer=0.01,
                          substrate={'fine': 0.8, 'coarse': 0.2})
    before = class_totals(d)
    d.evolve_to_end(finaltime=30.0)
    after = class_totals(d)
    np.testing.assert_allclose(after, before, rtol=1e-10)
    assert d.tracer_conserved_values[0].sum() > 0.0      # fine really moved


def test_a_class_absent_from_the_bed_is_not_entrained():
    d = build(classes=(FINE, COARSE), base_depth=0.5)
    d.set_bed_composition({'fine': 0.0, 'coarse': 1.0}, active_layer=0.01)
    d.evolve_to_end(finaltime=20.0)
    assert np.abs(d.tracer_conserved_values[0]).max() == 0.0


def test_without_composition_the_same_bed_does_supply_it():
    """Makes the test above non-vacuous: the flow does entrain fines."""
    d = build(classes=(FINE, COARSE), base_depth=0.5)
    d.evolve_to_end(finaltime=20.0)
    assert d.tracer_conserved_values[0].max() > 0.0


# --------------------------------------------------------------- armouring

def test_a_bed_losing_its_fines_armours():
    """The coarse class (gravel) does not move under this flow. Removing fines from a
    50/50 active layer that is refilled from a 50/50 substrate leaves the
    surface progressively coarser, so far less fine is entrained than from
    the shared bed, which never runs short."""
    gravel = ('coarse', 0.05)      # Shields ~0.02 here, below the threshold
    shared = build(classes=(FINE, gravel), base_depth=0.5)
    shared.evolve_to_end(finaltime=60.0)

    d = build(classes=(FINE, gravel), base_depth=0.5)
    d.set_bed_composition({'fine': 0.5, 'coarse': 0.5}, active_layer=0.002)
    d.evolve_to_end(finaltime=60.0)

    assert np.abs(d.tracer_conserved_values[1]).max() == 0.0     # coarse stays put
    fine_shared = float((shared.tracer_conserved_values[0] * shared.areas).sum())
    fine_comp = float((d.tracer_conserved_values[0] * d.areas).sum())
    assert fine_comp < 0.5 * fine_shared
    assert d.get_bed_composition('active')['fine'].mean() < 0.25


# ------------------------------------------------------------------ modes

def test_legacy_and_unified_agree():
    runs = {}
    for mode in ('legacy', 'unified'):
        d = build(classes=(FINE, COARSE), base_depth=0.5, mode=mode)
        d.set_bed_composition({'fine': 0.5, 'coarse': 0.5}, active_layer=0.005)
        d.evolve_to_end(finaltime=20.0)
        runs[mode] = d
    a, b = runs['legacy'], runs['unified']
    np.testing.assert_allclose(b.sediment_bed_active, a.sediment_bed_active,
                               rtol=1e-9, atol=1e-14)
    np.testing.assert_allclose(b.tracer_conserved_values[:2], a.tracer_conserved_values[:2],
                               rtol=1e-9, atol=1e-14)


# -------------------------------------------------------------------- API

def test_needs_an_erodible_base():
    d = build(base_depth=None)
    with pytest.raises(RuntimeError, match='set_erodible_base'):
        d.set_bed_composition({'sand': 1.0})


def test_fractions_must_sum_to_one():
    d = build(classes=(FINE, COARSE))
    with pytest.raises(ValueError, match='sum to'):
        d.set_bed_composition({'fine': 0.5, 'coarse': 0.3})


def test_unknown_class_is_refused():
    d = build()
    with pytest.raises(ValueError, match='unknown sediment class'):
        d.set_bed_composition({'mud': 1.0})


def test_no_fraction_after_composition():
    d = build()
    d.set_bed_composition({'sand': 1.0})
    with pytest.raises(RuntimeError, match='before set_bed_composition'):
        d.add_sediment_fraction(name='mud', diameter=2.0e-5)


def test_none_switches_it_off():
    d = build()
    d.set_bed_composition({'sand': 1.0})
    d.set_bed_composition(None)
    assert d.sediment_bed_composition == 0 and d.sediment_bed_active is None
