"""Suspended sediment must reach the sww file (issue #274, gap 2).

A grain size IS a tracer -- `add_sediment_fraction` calls `add_tracer`, so
class `s` occupies tracer slot `s` -- which means the tracer writer stores it
with no sediment-specific code at all. `test_tracers_sww.py` says as much in
its docstring, but says it about `add_tracer`; nothing checked the claim
through the entry point a sediment model actually uses.

That matters because the gap #274 reports is precisely an assumption of
coverage: 211 sediment references in shallow_water_domain.py and none in
sww.py, so suspended concentration was invisible in the output while the run
computed it correctly.

The companion file test_sediment_sww_storage.py covers gap 3, the evolving bed.
"""

import numpy as num
import pytest

import anuga

netCDF4 = pytest.importorskip('netCDF4')

D50 = 1.0e-4


def _domain(n=6):
    d = anuga.rectangular_cross_domain(n, n, len1=100.0, len2=100.0)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_quantity('friction', 0.03)
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    return d


def _run(d, tmp_path, name='sediment_sww', finaltime=2.0):
    d.set_name(name)
    d.set_datadir(str(tmp_path))
    for _ in d.evolve(yieldstep=1.0, finaltime=finaltime):
        pass
    return tmp_path / (name + '.sww')


def _vars(path):
    with netCDF4.Dataset(str(path), 'r') as fid:
        return set(fid.variables.keys())


def _read(path, name):
    with netCDF4.Dataset(str(path), 'r') as fid:
        return num.array(fid.variables[name][:])


def test_a_sediment_class_is_written(tmp_path):
    d = _domain()
    d.add_sediment_fraction(name='sand', diameter=D50, initial_concentration=0.02)
    names = _vars(_run(d, tmp_path))
    assert 'sand_c' in names, \
        'suspended concentration is missing from the sww: %s' % sorted(names)


def test_it_is_written_every_timestep(tmp_path):
    """Static storage would record the initial field and miss the transport."""
    d = _domain()
    d.add_sediment_fraction(name='sand', diameter=D50, initial_concentration=0.02)
    sww = _run(d, tmp_path)
    n_times = len(_read(sww, 'time'))
    sand = _read(sww, 'sand_c')
    assert sand.ndim == 2, 'sand_c is static, not per timestep'
    assert sand.shape == (n_times, len(d)), \
        'expected (%d, %d), got %r' % (n_times, len(d), sand.shape)


def test_the_stored_values_are_the_concentrations(tmp_path):
    d = _domain()
    d.add_sediment_fraction(name='sand', diameter=D50, initial_concentration=0.02)
    sww = _run(d, tmp_path)
    assert num.allclose(_read(sww, 'sand_c')[-1], d.get_tracer('sand'),
                        atol=1e-12), \
        'the last slice does not match the final concentration field'


def test_two_classes_do_not_alias(tmp_path):
    """Class s must land in slot s; a shared row would make them identical."""
    d = _domain()
    d.add_sediment_fraction(name='sand', diameter=D50, initial_concentration=0.02)
    d.add_sediment_fraction(name='silt', diameter=2.0e-5, initial_concentration=0.005)
    sww = _run(d, tmp_path)

    sand = _read(sww, 'sand_c')
    silt = _read(sww, 'silt_c')
    assert not num.allclose(sand, silt), 'the two classes wrote the same values'
    assert num.allclose(sand[-1], d.get_tracer('sand'), atol=1e-12)
    assert num.allclose(silt[-1], d.get_tracer('silt'), atol=1e-12)


def test_a_sediment_class_and_a_plain_tracer_coexist(tmp_path):
    """They share the tracer block, so a slot mix-up would show here."""
    d = _domain()
    d.add_sediment_fraction(name='sand', diameter=D50, initial_concentration=0.02)
    d.add_tracer('salinity', initial_value=0.03)
    sww = _run(d, tmp_path)

    names = _vars(sww)
    assert {'sand_c', 'salinity_c'} <= names
    assert num.allclose(_read(sww, 'salinity_c')[-1], d.get_tracer('salinity'),
                        atol=1e-12)
    assert num.allclose(_read(sww, 'sand_c')[-1], d.get_tracer('sand'),
                        atol=1e-12)


def test_settling_is_visible_in_the_file(tmp_path):
    """The end-to-end point: the file shows the concentration changing.

    Storing a field that never moves would satisfy every check above and still
    be useless for validating a sediment model.
    """
    d = _domain()
    d.add_sediment_fraction(name='sand', diameter=D50, initial_concentration=0.02)
    sww = _run(d, tmp_path, finaltime=4.0)
    sand = _read(sww, 'sand_c')
    moved = num.abs(sand[-1] - sand[0]).max()
    assert moved > 1e-8, 'concentration is identical at every timestep (%g)' % moved


def test_store_tracers_false_opts_out(tmp_path):
    d = _domain()
    d.add_sediment_fraction(name='sand', diameter=D50, initial_concentration=0.02)
    d.store_tracers = False
    assert 'sand_c' not in _vars(_run(d, tmp_path))


def test_the_bed_and_the_suspended_load_are_both_recorded(tmp_path):
    """#274 whole: gap 2 and gap 3 have to hold at the same time.

    A sediment run's output is only usable if BOTH the suspended
    concentration and the moving bed are in it -- either one alone gives a
    picture that looks complete and is not.
    """
    d = _domain()
    d.add_sediment_fraction(name='sand', diameter=D50, initial_concentration=0.02)
    assert d.sediment_bed_evolution, 'this test assumes an evolving bed'
    sww = _run(d, tmp_path, finaltime=4.0)

    sand = _read(sww, 'sand_c')
    elev = _read(sww, 'elevation_c')
    assert sand.ndim == 2, 'suspended concentration is not per timestep'
    assert elev.ndim == 2, 'the bed is stored statically; its motion is lost'
    assert num.abs(elev[-1] - elev[0]).max() > 0.0, \
        'the bed never moved, so this proves nothing about storing it'
