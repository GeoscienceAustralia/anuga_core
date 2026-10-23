"""Tracers must reach the sww file.

A tracer is not a Quantity -- it lives in the (n_tracers, N) block so the C
kernel can stride it -- so it takes the same route as a flag-3 quantity: one
dynamic `<name>_c` variable per tracer, centroids only, nothing interpolated
that the solver never computed.

Grain sizes are tracers (Domain.add_sediment_fraction calls add_tracer), so
they are covered by the same path.
"""

import numpy as num
import pytest

import anuga

netCDF4 = pytest.importorskip('netCDF4')


def _domain(n=6):
    d = anuga.rectangular_cross_domain(n, n)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    return d


def _run(d, tmp_path, name='tracer_sww', finaltime=1.0):
    d.set_name(name)
    d.set_datadir(str(tmp_path))
    for _ in d.evolve(yieldstep=0.5, finaltime=finaltime):
        pass
    return tmp_path / (name + '.sww')


def test_no_tracers_means_no_tracer_variables(tmp_path):
    # A domain without tracers must look exactly as it did before. Note
    # domain.store_centroids defaults to True, so stage_c and friends are
    # already written -- what must not appear is a tracer variable.
    d = _domain()
    with netCDF4.Dataset(_run(d, tmp_path)) as ds:
        assert 'stage_c' in ds.variables, 'baseline: centroids are stored'
        assert not [v for v in ds.variables
                    if v.endswith('_c') and v[:-2] not in
                    ('stage', 'elevation', 'friction', 'xmomentum',
                     'ymomentum', 'height', 'xvelocity', 'yvelocity')], \
            'no unexpected _c variables on a tracer-free domain'


def test_a_tracer_is_written_every_timestep(tmp_path):
    d = _domain()
    d.add_tracer('salinity', initial_value=0.02)
    with netCDF4.Dataset(_run(d, tmp_path)) as ds:
        assert 'salinity_c' in ds.variables, \
            'a registered tracer must appear in the sww'
        v = ds.variables['salinity_c']
        assert v.dimensions == ('number_of_timesteps', 'number_of_volumes')
        assert v.shape[0] > 1, 'tracers are dynamic, not stored once'
        assert v.shape[1] == len(d)


def test_the_stored_values_are_the_concentrations(tmp_path):
    d = _domain()
    d.add_tracer('salinity', initial_value=0.02)
    path = _run(d, tmp_path)
    with netCDF4.Dataset(path) as ds:
        first = num.asarray(ds.variables['salinity_c'][0, :])
    # Uniform initial concentration, still water: c stays 0.02 everywhere.
    assert num.allclose(first, 0.02, atol=1e-6), \
        'stored values should be c, not m = h*c'


def test_two_tracers_do_not_alias(tmp_path):
    d = _domain()
    d.add_tracer('a', initial_value=0.10)
    d.add_tracer('b', initial_value=0.90)
    with netCDF4.Dataset(_run(d, tmp_path)) as ds:
        a = num.asarray(ds.variables['a_c'][0, :])
        b = num.asarray(ds.variables['b_c'][0, :])
    assert num.allclose(a, 0.10, atol=1e-6)
    assert num.allclose(b, 0.90, atol=1e-6)


def test_store_tracers_false_opts_out(tmp_path):
    d = _domain()
    d.add_tracer('salinity', initial_value=0.02)
    d.store_tracers = False
    with netCDF4.Dataset(_run(d, tmp_path)) as ds:
        assert 'salinity_c' not in ds.variables


def test_tracer_output_does_not_disturb_the_usual_quantities(tmp_path):
    d = _domain()
    d.add_tracer('salinity', initial_value=0.02)
    with netCDF4.Dataset(_run(d, tmp_path)) as ds:
        for expected in ('stage', 'elevation', 'xmomentum', 'ymomentum'):
            assert expected in ds.variables, \
                '%s went missing once tracers were stored' % expected


def test_a_transported_tracer_actually_changes_in_the_file(tmp_path):
    """The property a modeller depends on: the field moves in the output."""
    d = _domain(10)
    d.set_quantity('elevation', lambda x, y: -0.1 * x)
    d.set_quantity('stage', lambda x, y: num.maximum(0.5 - 0.1 * x, 0.05))
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    # A blob on one side, so transport shows up as a change in the field.
    d.add_tracer('dye', initial_value=0.0)
    idx = d.get_centroid_coordinates()[:, 0] < 0.3
    c = d.get_tracer('dye')
    c[idx] = 1.0
    d.set_tracer('dye', c)

    path = _run(d, tmp_path, name='tracer_move', finaltime=3.0)
    with netCDF4.Dataset(path) as ds:
        arr = num.asarray(ds.variables['dye_c'][:])
    assert arr.shape[0] > 2
    assert not num.allclose(arr[0], arr[-1]), \
        'the tracer field never changed, so nothing was transported or stored'
