"""Mode 2 leaves the host vertex and edge values to their first reader.

At a yield a 'unified' domain brings back only the centroid values; the host
extrapolation to vertices and edges, which costs several GPU steps on a large
mesh, runs when something reads ``vertex_values`` or ``edge_values``. These
tests check that readers see exactly what the old eager extrapolation gave,
that the SWW file is unchanged, and that the extrapolation is skipped at
yields nothing reads at.
"""
import os

import numpy as num
import pytest

import anuga
from anuga import Reflective_boundary, rectangular_cross_domain
from anuga.shallow_water.shallow_water_domain import Domain


def dam_break(mode='unified', name=None, tracer=False):
    """Deep water on the left third, a dry bed with a bump elsewhere: moving
    wet/dry fronts, so the dry-cell parts of the extrapolation matter."""
    d = rectangular_cross_domain(12, 6, len1=60.0, len2=30.0)
    d.set_flow_algorithm('DE1')
    d.set_compute_mode(mode)
    d.set_quantity('elevation', lambda x, y: 0.5 * num.exp(-((x - 40.0) ** 2 + (y - 15.0) ** 2) / 50.0))
    d.set_quantity('stage', lambda x, y: num.where(x < 20.0, 1.0, 0.0))
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    if tracer:
        d.add_tracer('dye')
        d.set_tracer('dye', 0.5)
    if name is None:
        d.store = False
    else:
        d.set_name(name)
    return d


def count_refreshes(monkeypatch):
    calls = []
    orig = Domain._refresh_host_values

    def counted(self):
        calls.append(self.get_time())
        return orig(self)

    monkeypatch.setattr(Domain, '_refresh_host_values', counted)
    return calls


def test_a_yield_leaves_the_vertex_and_edge_values_to_the_first_read():
    d = dam_break()
    for t in d.evolve(yieldstep=0.5, finaltime=2.0):
        if t == 0.0:
            continue
        assert d._host_values_stale
        stage = d.quantities['stage']
        lazy_v = stage.vertex_values.copy()     # this read extrapolates
        assert not d._host_values_stale
        lazy_e = {n: q.edge_values.copy() for n, q in d.quantities.items()
                  if q.edge_values is not None}
        lazy_c = {n: q.centroid_values.copy() for n, q in d.quantities.items()}

        # The eager path: sync, protect and extrapolate everything
        d.distribute_to_vertices_and_edges()
        assert num.array_equal(lazy_v, stage.vertex_values)
        for n, e in lazy_e.items():
            assert num.array_equal(e, d.quantities[n].edge_values), n
        for n, c in lazy_c.items():
            assert num.array_equal(c, d.quantities[n].centroid_values), n


def test_the_centroids_at_a_yield_are_those_of_the_eager_path():
    """Heights, dry momenta and tracer concentrations c = m/h are brought up
    to date at the yield itself, not by the (skipped) extrapolation."""
    d = dam_break(tracer=True)
    for t in d.evolve(yieldstep=0.5, finaltime=2.0):
        if t == 0.0:
            continue
        lazy = {n: q.centroid_values.copy() for n, q in d.quantities.items()}
        lazy_c = d.get_tracer('dye').copy()
        assert d._host_values_stale             # nothing above read the edges
        d.distribute_to_vertices_and_edges()
        for n, c in lazy.items():
            assert num.array_equal(c, d.quantities[n].centroid_values), n
        assert num.array_equal(lazy_c, d.get_tracer('dye'))


def test_nothing_is_extrapolated_at_yields_nothing_reads_at(monkeypatch):
    calls = count_refreshes(monkeypatch)
    d = dam_break()
    for _ in d.evolve(yieldstep=0.25, finaltime=2.0):
        pass
    assert calls == []


def test_a_host_boundary_does_not_set_off_the_extrapolation(monkeypatch):
    """A boundary the device cannot evaluate is evaluated on the host from
    the boundary cells' edges, which come from the device."""
    class Host_reflective_boundary(Reflective_boundary):
        pass                                    # not in GPU_BOUNDARY_TYPES

    calls = count_refreshes(monkeypatch)
    d = dam_break()
    d.set_boundary({t: Host_reflective_boundary(d) for t in d.get_boundary_tags()})
    for _ in d.evolve(yieldstep=0.25, finaltime=1.0):
        pass
    assert not d._gpu_all_on_gpu
    assert calls == []


def test_the_sww_file_is_unchanged(tmp_path, monkeypatch):
    """The default SWW output (centroids, and vertex values averaged from
    them) is bit for bit what the eager extrapolation gave, and writing it
    does not extrapolate."""
    from netCDF4 import Dataset

    monkeypatch.chdir(tmp_path)
    files = {}
    for how in ('lazy', 'eager'):
        d = dam_break(name=how, tracer=True)
        if how == 'eager':
            d._prepare_yield = d.distribute_to_vertices_and_edges
        calls = count_refreshes(monkeypatch)
        for _ in d.evolve(yieldstep=0.25, finaltime=2.0):
            pass
        if how == 'lazy':
            assert calls == []
        files[how] = os.path.join(d.get_datadir(), how + '.sww')

    with Dataset(files['lazy']) as a, Dataset(files['eager']) as b:
        assert set(a.variables) == set(b.variables)
        for v in a.variables:
            assert num.array_equal(a[v][:], b[v][:]), v


def test_unsmoothed_output_extrapolates_on_output_steps_only(tmp_path, monkeypatch):
    """Without centroid averaging the SWW writer reads the vertex values, so
    each stored frame extrapolates once, and yields between output steps do
    not."""
    monkeypatch.chdir(tmp_path)
    calls = count_refreshes(monkeypatch)
    d = dam_break(name='unsmoothed')
    d.set_using_centroid_averaging(False)
    for _ in d.evolve(yieldstep=0.25, outputstep=1.0, finaltime=2.0):
        pass
    # the frames at t = 0, 1 and 2, and none of the 6 yields between them
    assert calls == pytest.approx([0.0, 1.0, 2.0])


def test_legacy_mode_extrapolates_at_every_yield(monkeypatch):
    calls = count_refreshes(monkeypatch)
    d = dam_break(mode='legacy')
    for _ in d.evolve(yieldstep=0.5, finaltime=1.0):
        assert not d._host_values_stale
    assert calls == []


def test_an_old_pickled_quantity_keeps_its_edge_values():
    """Checkpoints from before edge_values was a property."""
    d = rectangular_cross_domain(2, 2)
    q = d.quantities['stage']
    state = dict(q.__dict__)
    state['edge_values'] = state.pop('_edge_values')
    restored = anuga.Quantity.__new__(anuga.Quantity)
    restored.__setstate__(state)
    assert restored.edge_values is q.edge_values
