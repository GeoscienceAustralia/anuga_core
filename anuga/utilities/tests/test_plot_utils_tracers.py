"""Tracers must be reachable from the sww, not just present in it.

#276 got tracers written; they arrived as `<name>_c`, which is exactly what a
quantity's centroid variable looks like, so `plot_utils` -- which reads a fixed
list of quantities -- could not surface them. The data was in the file and
invisible to ANUGA's own reader.

The file now names its tracers in a `tracer_names` attribute rather than
leaving every reader to guess from a suffix.
"""

import numpy as num
import pytest

import anuga
import anuga.utilities.plot_utils as util

netCDF4 = pytest.importorskip('netCDF4')


def _run(tmp_path, tracers=(('salinity', 0.03), ('dye', 0.5)), name='tr'):
    d = anuga.rectangular_cross_domain(6, 6, len1=60.0, len2=60.0)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    for tname, value in tracers:
        d.add_tracer(tname, initial_value=value)
    d.set_name(name)
    d.set_datadir(str(tmp_path))
    for _ in d.evolve(yieldstep=0.5, finaltime=1.0):
        pass
    return d, str(tmp_path / (name + '.sww'))


# --- the file says what its tracers are ------------------------------------

def test_the_file_records_its_tracer_names(tmp_path):
    _, sww = _run(tmp_path)
    with netCDF4.Dataset(sww, 'r') as fid:
        assert fid.tracer_names == 'salinity dye', \
            'the writer did not record the tracer names in order'


def test_a_file_with_no_tracers_records_an_empty_list(tmp_path):
    """Empty, not absent: it distinguishes "none" from "written before this"."""
    _, sww = _run(tmp_path, tracers=(), name='notr')
    with netCDF4.Dataset(sww, 'r') as fid:
        assert fid.tracer_names == ''
    assert util.get_tracer_names(sww) == []


def test_get_tracer_names_reads_them_in_order(tmp_path):
    _, sww = _run(tmp_path)
    assert util.get_tracer_names(sww) == ['salinity', 'dye']


def test_an_old_file_without_the_attribute_still_works(tmp_path):
    """Falls back to "every <name>_c that is not a known quantity"."""
    _, sww = _run(tmp_path)
    with netCDF4.Dataset(sww, 'a') as fid:
        del fid.tracer_names
    assert sorted(util.get_tracer_names(sww)) == ['dye', 'salinity']


def test_the_fallback_does_not_mistake_a_quantity_for_a_tracer(tmp_path):
    """stage_c and friends are <name>_c too -- the reason for the attribute."""
    _, sww = _run(tmp_path, tracers=(), name='q')
    with netCDF4.Dataset(sww, 'a') as fid:
        del fid.tracer_names
    assert util.get_tracer_names(sww) == [], \
        'a quantity was reported as a tracer'


# --- and the readers surface them ------------------------------------------

def test_get_centroids_exposes_the_tracers(tmp_path):
    """The path that matters: tracers are centroid quantities."""
    d, sww = _run(tmp_path)
    pc = util.get_centroids(sww)
    assert sorted(pc.tracers) == ['dye', 'salinity']
    assert num.allclose(pc.tracers['salinity'][-1], d.get_tracer('salinity'))
    assert num.allclose(pc.tracers['dye'][-1], d.get_tracer('dye'))


def test_get_output_exposes_the_tracers(tmp_path):
    d, sww = _run(tmp_path)
    p = util.get_output(sww)
    assert sorted(p.tracers) == ['dye', 'salinity']
    assert num.allclose(p.tracers['salinity'][-1], d.get_tracer('salinity'))


def test_the_tracers_are_a_dict_not_attributes(tmp_path):
    """So a tracer name can never shadow stage, vel or timeSlices."""
    _, sww = _run(tmp_path)
    pc = util.get_centroids(sww)
    assert isinstance(pc.tracers, dict)
    assert not hasattr(pc, 'salinity'), \
        'tracers were set as attributes, where they could collide'


def test_a_tracer_cannot_be_named_after_a_quantity():
    """Both go to <name>_c, so the tracer would overwrite it in the sww.

    Found by writing these tests: a tracer called 'stage' really did replace
    the stage in the output, with no error.
    """
    d = anuga.rectangular_cross_domain(4, 4)
    with pytest.raises(ValueError, match='quantity'):
        d.add_tracer('stage')
    with pytest.raises(ValueError, match='quantity'):
        d.add_tracer('elevation')


def test_a_tracer_cannot_take_a_reserved_max_name():
    d = anuga.rectangular_cross_domain(4, 4)
    with pytest.raises(ValueError, match='reserved'):
        d.add_tracer('max_depth')


@pytest.mark.parametrize('slices', ['all', 'last', 'max'])
def test_the_time_slicing_options_work(tmp_path, slices):
    _, sww = _run(tmp_path)
    pc = util.get_centroids(sww, timeSlices=slices)
    assert sorted(pc.tracers) == ['dye', 'salinity']
    for v in pc.tracers.values():
        assert v.shape[0] >= 1


def test_no_tracers_means_an_empty_dict(tmp_path):
    _, sww = _run(tmp_path, tracers=(), name='none')
    assert util.get_centroids(sww).tracers == {}
    assert util.get_output(sww).tracers == {}
