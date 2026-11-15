#!/usr/bin/env python
"""Meshes sized from a DEM: D8 routing, HAND, size fields, pit removal."""

import heapq
import os
import tempfile

import numpy as num
import pytest

import anuga
from anuga.pmesh.dem_flow_ext import height_above_drainage, route
from anuga.pmesh.dem_mesh import (_finer, elevation_from_dem, read_dem,
                                  remove_mesh_pits, size_field_from_dem)
from anuga.pmesh.size_field import Size_field, area_to_edge

pytest.importorskip('meshpy')


def valley(n=81, cell=5.0, x0=1000.0, y0=2000.0):
    """A V-shaped valley draining south: a channel along x = middle.

    Returns the raster tuple (x, y, Z), Z[i, j] at (x[i], y[j]).
    """
    x = x0 + cell * num.arange(n)
    y = y0 + cell * num.arange(n)
    X, Y = num.meshgrid(x, y, indexing='ij')
    xc = x0 + cell * (n // 2)
    Z = 0.05 * num.abs(X - xc) + 0.01 * (Y - y0)
    return x, y, Z


def reference_route(z, epsilon=1e-3):
    """Plain Python priority-flood + D8, for checking the kernel."""
    nr, nc = z.shape
    dr = [-1, -1, 0, 1, 1, 1, 0, -1]
    dc = [0, 1, 1, 1, 0, -1, -1, -1]
    zf = num.full((nr, nc), num.inf)
    heap = []
    for r in range(nr):
        for c in range(nc):
            if r in (0, nr - 1) or c in (0, nc - 1):
                zf[r, c] = z[r, c]
                heapq.heappush(heap, (z[r, c], r, c))
    while heap:
        v, r, c = heapq.heappop(heap)
        for k in range(8):
            rr, cc = r + dr[k], c + dc[k]
            if 0 <= rr < nr and 0 <= cc < nc and zf[rr, cc] == num.inf:
                zf[rr, cc] = max(z[rr, cc], v + epsilon)
                heapq.heappush(heap, (zf[rr, cc], rr, cc))
    return zf


# ---------------------------------------------------------------------------
# Routing kernel
# ---------------------------------------------------------------------------

def test_route_fills_pit_and_every_cell_drains():
    rng = num.random.default_rng(0)
    n = 40
    r, c = num.mgrid[0:n, 0:n]
    z = 20.0 - 0.3 * r + 2.0 * num.sin(c / 5.0) + rng.normal(0, 0.3, (n, n))
    z[15:21, 15:21] -= 10.0                          # a pit to fill
    z = num.ascontiguousarray(z)

    rec, order, acc = route(z, 1e-3)

    # Every interior cell has a receiver, and following receivers always
    # descends the filled surface, so they all leave the grid
    zf = reference_route(z)
    interior = num.zeros((n, n), bool)
    interior[1:-1, 1:-1] = True
    assert num.all(rec.reshape(n, n)[interior] >= 0)
    flat = zf.ravel()
    has = rec >= 0
    assert num.all(flat[rec[has]] < flat[has])

    # Each cell counts itself once: the cells leaving the grid carry them all
    assert acc.min() == 1.0
    assert acc.ravel()[rec < 0].sum() == n * n
    # a receiver always comes before its cell in the order
    pos = num.empty(n * n, int)
    pos[order] = num.arange(len(order))
    assert num.all(pos[rec[has]] < pos[has])


def test_route_nodata_and_walls():
    x, y, Z = valley(n=31)
    z = num.ascontiguousarray(Z.T)            # [row = y, col = x]
    z[10:13, 12:19] = num.inf                 # a wall across the channel
    z[0, :5] = num.nan                        # nodata

    rec, order, acc = route(z, 1e-3)
    wall = num.isinf(z).ravel()
    nodata = num.isnan(z).ravel()
    assert num.all(acc.ravel()[wall | nodata] == 0)
    assert not num.isin(rec, num.flatnonzero(wall | nodata)).any()
    # Water goes around the wall: everything routed still leaves the grid
    assert acc.ravel()[rec < 0].sum() == (~wall & ~nodata).sum()


def test_height_above_drainage_in_a_valley():
    x, y, Z = valley()
    z = num.ascontiguousarray(Z.T)
    rec, order, acc = route(z, 1e-3)
    hand = height_above_drainage(z, rec, order, acc, 200.0)

    channel = len(x) // 2
    # The channel itself is at 0, and the sides rise at the side slope
    rows = slice(5, len(y) - 30)              # well upstream of the outlet
    assert num.allclose(hand[rows, channel], 0.0)
    off = 10
    assert num.allclose(hand[rows, channel + off], 0.05 * 5.0 * off, atol=0.06)


# ---------------------------------------------------------------------------
# DEM input
# ---------------------------------------------------------------------------

def test_read_dem_asc_and_tuple():
    x, y, Z = valley(n=11)
    with tempfile.TemporaryDirectory() as tmp:
        name = os.path.join(tmp, 'dem.asc')
        with open(name, 'w') as f:
            f.write('ncols 11\nnrows 11\nxllcorner %g\nyllcorner %g\n'
                    'cellsize 5\nNODATA_value -9999\n' % (x[0], y[0]))
            for j in range(len(y) - 1, -1, -1):        # north row first
                f.write(' '.join('%.6f' % v for v in Z[:, j]) + '\n')
        xa, ya, Za = read_dem(name)
    assert num.allclose(xa, x) and num.allclose(ya, y)
    assert num.allclose(Za, Z, atol=1e-5)

    with pytest.raises(ValueError):
        read_dem((x, y, Z[:-1]))


def test_size_field_is_fine_along_the_channel():
    x, y, Z = valley()
    field = size_field_from_dem((x, y, Z), maximum_triangle_area=400.0,
                                streams=[(2000.0, 10.0)], band=1.0)
    channel = x[len(x) // 2]
    ymid = y[len(y) // 2]
    assert field(channel, ymid) == pytest.approx(10.0)
    assert field(x[0], ymid) == pytest.approx(400.0)
    # graded: the edge grows no faster than 0.5 per unit distance
    xs = num.linspace(channel, x[0], 50)
    edge = area_to_edge(field(xs, num.full_like(xs, ymid)))
    assert num.all(num.diff(edge) <= 0.5 * abs(xs[1] - xs[0]) + 5.0)
    # the drainage area comes with it
    xd, yd, A = field.drainage_area
    assert A.shape == Z.shape and A.max() > 2000.0


def test_size_field_floodplain_and_regions():
    x, y, Z = valley()
    channel = x[len(x) // 2]
    ymid = y[len(y) // 2]
    # within 0.6 m of the channel: |x - channel| <= 12 m
    field = size_field_from_dem((x, y, Z), maximum_triangle_area=400.0,
                                floodplain=[(0.6, 2000.0, 20.0)], grade=0.0)
    assert field(channel + 10.0, ymid) == pytest.approx(20.0)
    assert field(channel + 30.0, ymid) == pytest.approx(400.0)

    box = [[x[2], y[2]], [x[8], y[2]], [x[8], y[8]], [x[2], y[8]]]
    field = size_field_from_dem((x, y, Z), maximum_triangle_area=400.0,
                                regions=[(box, 5.0)], grade=0.0)
    assert field(x[5], y[5]) == pytest.approx(5.0)
    assert field(x[40], y[40]) == pytest.approx(400.0)


def test_finer_combines_two_fields():
    x, y, Z = valley(n=21)
    field = size_field_from_dem((x, y, Z), maximum_triangle_area=400.0)
    other = Size_field(num.array([[50.0, 600.0]]), x[0] - 2.5, y[0] - 2.5,
                       50.0)
    both = _finer(field, other)
    assert isinstance(both, Size_field)
    assert both(x[2], y[2]) == pytest.approx(50.0)
    assert both(x[-1], y[2]) == pytest.approx(400.0)
    # a scalar-only callable works too
    both = _finer(field, lambda px, py: 7.0 if px < x[10] else 900.0)
    assert both(x[2], y[2]) == pytest.approx(7.0)
    assert both(x[-1], y[2]) == pytest.approx(400.0)


# ---------------------------------------------------------------------------
# Meshes and domains
# ---------------------------------------------------------------------------

def test_create_domain_from_dem_defaults():
    x, y, Z = valley()
    domain = anuga.create_domain_from_dem(
        (x, y, Z), maximum_triangle_area=200.0, streams=[(2000.0, 10.0)])
    assert set(domain.get_boundary_tags()) == {'bottom', 'right', 'top',
                                               'left'}
    # the mesh covers the DEM's cell centres
    assert domain.get_area() == pytest.approx(
        (x[-1] - x[0]) * (y[-1] - y[0]))
    # finer near the channel
    cen = domain.get_centroid_coordinates(absolute=True)
    near = num.abs(cen[:, 0] - x[len(x) // 2]) < 10.0
    assert domain.areas[near].mean() < 0.25 * domain.areas[~near].mean()
    # elevation is the DEM's: on each sloping side it is linear, and the
    # four-point mean of a plane is its centroid value
    z = domain.quantities['elevation'].centroid_values
    xc = x[len(x) // 2]
    exact = 0.05 * num.abs(cen[:, 0] - xc) + 0.01 * (cen[:, 1] - y[0])
    v = domain.get_vertex_coordinates(absolute=True).reshape(-1, 3, 2)
    one_side = num.all(v[:, :, 0] >= xc, axis=1) | \
        num.all(v[:, :, 0] <= xc, axis=1)
    assert num.allclose(z[one_side], exact[one_side], atol=1e-6)


def test_create_pmesh_from_dem_with_holes_and_polygon():
    x, y, Z = valley()
    outline = [[x[2], y[2]], [x[-3], y[2]], [x[-3], y[-3]], [x[2], y[-3]]]
    hole = [[x[50], y[50]], [x[55], y[50]], [x[55], y[55]], [x[50], y[55]]]
    tags = {'south': [0], 'east': [1], 'north': [2], 'west': [3]}
    pmesh = anuga.create_pmesh_from_dem(
        (x, y, Z), 200.0, bounding_polygon=outline, boundary_tags=tags,
        streams=[(2000.0, 10.0)], interior_holes=[hole])
    domain = anuga.Domain(anuga.pmesh_to_mesh(pmesh))
    assert set(domain.get_boundary_tags()) == {'south', 'east', 'north',
                                               'west', 'interior'}
    area = (x[-3] - x[2]) * (y[-3] - y[2]) - (x[55] - x[50]) ** 2
    assert domain.get_area() == pytest.approx(area)
    assert hasattr(pmesh.size_field, 'drainage_area')


# ---------------------------------------------------------------------------
# Pit removal
# ---------------------------------------------------------------------------

def pit_domain():
    """A plane falling to the west, with a 0.5 m pit and a 3 m pit whose
    ways out do not cross."""
    domain = anuga.rectangular_cross_domain(20, 10, len1=200.0, len2=100.0)
    domain.set_quantity('elevation', lambda x, y: 0.01 * x,
                        location='centroids')
    z = domain.quantities['elevation'].centroid_values.copy()
    cen = domain.centroid_coordinates
    shallow = num.hypot(cen[:, 0] - 55.0, cen[:, 1] - 55.0) < 8.0
    deep = num.hypot(cen[:, 0] - 145.0, cen[:, 1] - 20.0) < 8.0
    z[shallow] -= 0.5
    z[deep] -= 3.0
    domain.set_quantity('elevation', z, location='centroids')
    return domain, z, shallow, deep


def pits(domain):
    """Cells lower than all their neighbours, other than the outlets."""
    z = domain.quantities['elevation'].centroid_values
    nb = domain.neighbours
    lower = num.where(nb >= 0, z[num.maximum(nb, 0)], num.inf).min(axis=1)
    pit = lower > z
    pit[[vol for (vol, _), tag in domain.boundary.items()
         if tag == 'left']] = False
    return pit


def test_fill_removes_all_pits():
    domain, z0, shallow, deep = pit_domain()
    stats = remove_mesh_pits(domain, ['left'], method='fill')
    z = domain.quantities['elevation'].centroid_values
    assert stats['pits_left'] == 0 and not pits(domain).any()
    assert num.all(z >= z0 - 1e-12)
    assert stats['lowered'] == 0 and stats['raised'] > 0
    assert stats['volume_raised'] > 0.0


def test_breach_cuts_instead_and_respects_max_depth():
    domain, z0, shallow, deep = pit_domain()
    stats = remove_mesh_pits(domain, ['left'], method='breach')
    z = domain.quantities['elevation'].centroid_values
    assert not pits(domain).any()
    assert stats['lowered'] > 0
    # the pits keep their beds: nothing in them is raised, and their lowest
    # cells stay within a few epsilon (1 mm) steps of where they were
    for pit in (shallow, deep):
        assert num.all(z[pit] <= z0[pit] + 1e-12)
        assert z[pit].min() == pytest.approx(z0[pit].min(), abs=0.01)
    uncapped = stats['volume_lowered']

    # capped at 1 m: the 0.5 m pit is breached, the 3 m one filled
    domain, z0, shallow, deep = pit_domain()
    stats = remove_mesh_pits(domain, ['left'], method='breach', max_depth=1.0)
    z = domain.quantities['elevation'].centroid_values
    # The uncapped breach cut a 3 m trench all the way to the outlet
    assert stats['volume_lowered'] < 0.2 * uncapped
    assert not pits(domain).any()
    assert num.all(z0 - z <= 1.0 + 1e-9)
    assert num.allclose(z[shallow], z0[shallow])
    assert num.any(z[deep] > z0[deep] + 0.5)


def test_where_limits_the_change():
    domain, z0, shallow, deep = pit_domain()
    where = ~shallow
    stats = remove_mesh_pits(domain, ['left'], method='fill', where=where)
    z = domain.quantities['elevation'].centroid_values
    assert num.array_equal(z[shallow], z0[shallow])
    assert stats['pits_left'] == 0          # counted only where it may act


def test_pit_arguments():
    domain, _, _, _ = pit_domain()
    with pytest.raises(ValueError):
        remove_mesh_pits(domain, ['nowhere'])
    with pytest.raises(ValueError):
        remove_mesh_pits(domain, ['left'], method='dig')
    with pytest.raises(ValueError):
        anuga.create_domain_from_dem(valley(n=11), 200.0, pits='fill')
