#!/usr/bin/env python
"""Mesh resolution from a size field (Triangle's user test)."""

import os
import tempfile
import unittest

import numpy as num
import pytest

import anuga
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.pmesh.mesh_interface import create_pmesh_from_regions
from anuga.pmesh.size_field import (Size_field, size_field_from_regions,
                                    make_refinement_func, area_to_edge,
                                    edge_to_area)

pytest.importorskip('meshpy')


SQUARE = [[0.0, 0.0], [100.0, 0.0], [100.0, 100.0], [0.0, 100.0]]
TAGS = {'bottom': [0], 'right': [1], 'top': [2], 'left': [3]}


def generate(size_field=None, polygon=SQUARE, maximum_triangle_area=500.0,
             **kwargs):
    """Triangulate; return absolute centroids, areas and the Pmesh."""
    m = create_pmesh_from_regions(polygon, TAGS,
                                  maximum_triangle_area=maximum_triangle_area,
                                  size_field=size_field, verbose=False,
                                  **kwargs)
    m.generate_mesh(verbose=False)
    return centroids_areas(m) + (m,)


def centroids_areas(m):
    v = num.asarray(m.tri_mesh.vertices, float)
    t = num.asarray(m.tri_mesh.triangles, int)
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    areas = 0.5 * num.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
                          - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))
    centroids = (a + b + c) / 3.0
    centroids[:, 0] += m.geo_reference.get_xllcorner()
    centroids[:, 1] += m.geo_reference.get_yllcorner()
    return centroids, areas


class Test_size_field(unittest.TestCase):

    def test_callable_field_is_respected(self):
        """No triangle exceeds an analytic field at its centroid."""

        def field(x, y):
            return 2.0 + 0.05 * x * x          # fine at the left, coarse right

        centroids, areas, _ = generate(field)
        target = field(centroids[:, 0], centroids[:, 1])
        assert num.all(areas <= target * (1 + 1e-12))

        # and the mesh really is graded by it
        left = areas[centroids[:, 0] < 20.0]
        right = areas[centroids[:, 0] > 80.0]
        assert left.mean() < 0.1 * right.mean()

        # with many more triangles than without the field
        _, plain, _ = generate(None)
        assert len(areas) > 5 * len(plain)

    def test_field_composes_with_maximum_area(self):
        """The field refines further; it never relaxes the maximum area."""

        centroids, areas, _ = generate(lambda x, y: 1.0e6,
                                       maximum_triangle_area=50.0)
        assert areas.max() <= 50.0 * (1 + 1e-12)

    def test_field_adds_no_segments(self):
        """Only the bounding polygon's segments are in the mesh."""

        field = size_field_from_regions(
            SQUARE, [([[40, 40], [60, 40], [60, 60], [40, 60]], 2.0)], 500.0)
        _, _, m = generate(field)

        v = num.asarray(m.tri_mesh.vertices, float)
        on_segment = v[num.unique(m.tri_mesh.segments)]
        on_boundary = ((on_segment[:, 0] == 0.0) | (on_segment[:, 0] == 100.0)
                       | (on_segment[:, 1] == 0.0) | (on_segment[:, 1] == 100.0))
        assert num.all(on_boundary)

    def test_callable_gets_absolute_coordinates(self):
        """The mesh is relative to its geo reference; the field is not."""

        x0, y0 = 500000.0, 6800000.0
        polygon = [[x0 + x, y0 + y] for x, y in SQUARE]
        seen = []

        def field(x, y):
            seen.append((x, y))
            return 2.0 if x < x0 + 50.0 else 400.0

        centroids, areas, m = generate(field, polygon=polygon)
        assert m.geo_reference.get_xllcorner() == x0

        seen = num.array(seen)
        assert seen[:, 0].min() >= x0 and seen[:, 0].max() <= x0 + 100.0
        assert seen[:, 1].min() >= y0 and seen[:, 1].max() <= y0 + 100.0

        west = centroids[:, 0] < x0 + 50.0
        assert areas[west].max() <= 2.0 * (1 + 1e-12)
        assert areas[~west].mean() > 10.0 * areas[west].mean()

    def test_bad_callable_raises(self):
        """A failing or non-positive field is an error, not a crash or hang."""

        def broken(x, y):
            raise RuntimeError('boom')

        with self.assertRaisesRegex(RuntimeError, 'boom'):
            generate(broken)
        with self.assertRaisesRegex(ValueError, 'must be positive'):
            generate(lambda x, y: 0.0)
        with self.assertRaisesRegex(ValueError, 'must be positive'):
            generate(lambda x, y: float('nan'))

    def test_raster_field_matches_callable(self):
        """A Size_field gives the same mesh as the function it rasterises."""

        x0, y0 = 300000.0, 7000000.0
        polygon = [[x0 + x, y0 + y] for x, y in SQUARE]
        cell = 1.0
        xc = (num.arange(100) + 0.5) * cell
        area = num.tile(2.0 + 0.05 * xc * xc, (100, 1))
        raster = Size_field(area, x0, y0, cell)

        def field(x, y):
            return raster(x, y)

        _, areas_r, m_r = generate(raster, polygon=polygon)
        _, areas_c, m_c = generate(field, polygon=polygon)
        assert len(areas_r) == len(areas_c)
        assert num.allclose(m_r.tri_mesh.vertices, m_c.tri_mesh.vertices)

    def test_size_field_lookup(self):
        area = num.array([[1.0, 2.0, 3.0],
                          [4.0, 5.0, 6.0]])        # row 0 is the south edge
        f = Size_field(area, 10.0, 20.0, 5.0)

        assert f(10.1, 20.1) == 1.0
        assert f(24.9, 20.1) == 3.0
        assert f(10.1, 29.9) == 4.0
        assert f(17.5, 27.5) == 5.0
        # outside: nearest cell
        assert f(-100.0, -100.0) == 1.0
        assert f(1000.0, 1000.0) == 6.0
        assert f(17.5, -100.0) == 2.0
        # arrays
        assert num.allclose(f([10.1, 24.9], [20.1, 29.9]), [1.0, 6.0])

        # the user test Triangle is given agrees, in mesh-relative coordinates
        too_big = f.refinement_func(10.0, 20.0)
        tri = [(7.0, 7.0), (8.0, 7.0), (7.5, 8.5)]       # centroid (7.5, 7.5)
        assert too_big(tri, 5.1)
        assert not too_big(tri, 4.9)
        far = [(-900.0, -900.0), (-899.0, -900.0), (-899.5, -899.0)]
        assert too_big(far, 1.1)
        assert not too_big(far, 0.9)

    def test_size_field_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            Size_field(num.zeros((3, 3)), 0.0, 0.0, 1.0)
        with self.assertRaises(ValueError):
            Size_field(num.array([[1.0, num.nan]]), 0.0, 0.0, 1.0)
        with self.assertRaises(ValueError):
            Size_field(num.ones(4), 0.0, 0.0, 1.0)
        with self.assertRaises(ValueError):
            Size_field(num.ones((2, 2)), 0.0, 0.0, 0.0)
        with self.assertRaises(TypeError):
            make_refinement_func(42.0)

    def test_size_field_from_regions(self):
        inner = [[40.0, 40.0], [60.0, 40.0], [60.0, 60.0], [40.0, 60.0]]
        grade = 0.5
        f = size_field_from_regions(SQUARE, [(inner, 2.0)], 100.0,
                                    grade=grade, cellsize=1.0)

        # the region's area inside it, the maximum far away
        assert num.isclose(f(50.0, 50.0), 2.0)
        assert num.isclose(f(41.5, 58.5), 2.0)
        assert num.isclose(f(1.0, 1.0), 100.0)
        assert f.area.min() >= 2.0 * (1 - 1e-12)
        assert f.area.max() <= 100.0 * (1 + 1e-12)

        # graded: the edge length grows at `grade` per metre from the region
        h = area_to_edge(f.area)
        assert num.abs(num.diff(h, axis=0)).max() <= grade * f.cellsize + 1e-9
        assert num.abs(num.diff(h, axis=1)).max() <= grade * f.cellsize + 1e-9
        h_inner = float(area_to_edge(2.0))
        assert num.isclose(area_to_edge(f(70.5, 50.5)),
                           h_inner + grade * 10.0, atol=grade * f.cellsize)

        # ungraded: a step at the outline
        g = size_field_from_regions(SQUARE, [(inner, 2.0)], 100.0,
                                    grade=0.0, cellsize=1.0)
        assert num.isclose(g(50.0, 50.0), 2.0)
        assert num.isclose(g(65.0, 50.0), 100.0)

        assert num.isclose(edge_to_area(area_to_edge(7.0)), 7.0)

    def test_size_field_from_regions_overlap_and_slivers(self):
        """Regions may overlap (finest wins) and be narrower than a cell."""

        coarse = [[20.0, 20.0], [80.0, 20.0], [80.0, 80.0], [20.0, 80.0]]
        fine = [[40.0, 40.0], [60.0, 40.0], [60.0, 60.0], [40.0, 60.0]]
        for regions in ([(coarse, 50.0), (fine, 2.0)],
                        [(fine, 2.0), (coarse, 50.0)]):
            f = size_field_from_regions(SQUARE, regions, 500.0, cellsize=1.0)
            assert num.isclose(f(50.0, 50.0), 2.0)
            assert f(22.0, 22.0) <= 50.0 * (1 + 1e-12)

        # a strip 0.2 wide on a 5 m raster contains no cell centre
        strip = [[10.0, 51.1], [90.0, 51.1], [90.0, 51.3], [10.0, 51.3]]
        f = size_field_from_regions(SQUARE, [(strip, 2.0)], 500.0,
                                    cellsize=5.0)
        assert num.isclose(f(10.0, 51.2), 2.0)
        assert num.isclose(f(90.0, 51.2), 2.0)

    def test_rasterise_polygon(self):
        """The scanline fill agrees with a point-in-polygon test."""
        from matplotlib.path import Path
        from anuga.pmesh.size_field import _rasterise_polygon

        rng = num.random.default_rng(1)
        xll, yll, cell, nrow, ncol = 3.0, -7.0, 0.7, 150, 170
        xc = xll + (num.arange(ncol) + 0.5) * cell
        yc = yll + (num.arange(nrow) + 0.5) * cell
        X, Y = num.meshgrid(xc, yc)
        centres = num.column_stack([X.ravel(), Y.ravel()])

        # a ragged star, a concave L, and one hanging off the raster
        angle = num.sort(rng.uniform(0.0, 2.0 * num.pi, 60))
        radius = rng.uniform(10.0, 40.0, 60)
        star = num.column_stack([60.0 + radius * num.cos(angle),
                                 45.0 + radius * num.sin(angle)])
        ell = num.array([[10.0, 0.0], [50.0, 0.0], [50.0, 10.0], [20.0, 10.0],
                         [20.0, 60.0], [10.0, 60.0]])
        off = star + [70.0, -50.0]

        for polygon in (star, star[::-1], ell, off):
            expected = Path(polygon).contains_points(centres).reshape(X.shape)
            row0, col0, inside = _rasterise_polygon(polygon, xll, yll, cell,
                                                    nrow, ncol)
            mask = num.zeros((nrow, ncol), bool)
            mask[row0:row0 + inside.shape[0],
                 col0:col0 + inside.shape[1]] = inside
            assert expected.sum() > 100
            assert num.array_equal(mask, expected)

        # entirely outside the raster
        row0, col0, inside = _rasterise_polygon(star + [1.0e4, 1.0e4], xll,
                                                yll, cell, nrow, ncol)
        assert inside.size == 0

    def test_size_field_from_regions_geo_reference(self):
        """Polygons relative to a geo reference give an absolute field."""

        geo = Geo_reference(56, 1000.0, 2000.0)
        inner = [[40.0, 40.0], [60.0, 40.0], [60.0, 60.0], [40.0, 60.0]]
        f = size_field_from_regions(SQUARE, [(inner, 2.0)], 100.0,
                                    cellsize=1.0, poly_geo_reference=geo)
        assert num.isclose(f(1050.0, 2050.0), 2.0)
        assert num.isclose(f(1001.0, 2001.0), 100.0)

    def test_default_cellsize_is_bounded(self):
        f = size_field_from_regions(SQUARE, [], 500.0)
        assert num.isclose(f.cellsize, float(area_to_edge(500.0)))
        f = size_field_from_regions(SQUARE, [], 500.0, max_cells=16)
        assert f.cellsize >= 25.0
        assert num.allclose(f.area, 500.0)

    def test_field_resolves_region_better_than_plain(self):
        """A graded field from a region puts the fine cells in the region."""

        inner = [[40.0, 40.0], [60.0, 40.0], [60.0, 60.0], [40.0, 60.0]]
        field = size_field_from_regions(SQUARE, [(inner, 2.0)], 500.0)
        centroids, areas, _ = generate(field)

        inside = ((num.abs(centroids[:, 0] - 50.0) < 9.0)
                  & (num.abs(centroids[:, 1] - 50.0) < 9.0))
        corner = (centroids[:, 0] < 15.0) & (centroids[:, 1] < 15.0)
        assert areas[inside].max() <= 2.0 * (1 + 1e-12)
        assert areas[corner].mean() > 20.0

    def test_written_mesh_file_uses_field(self):
        """With a filename the field is applied when the file is written."""

        field = size_field_from_regions(
            SQUARE, [([[40, 40], [60, 40], [60, 60], [40, 60]], 2.0)], 500.0)
        fd, name = tempfile.mkstemp(suffix='.msh')
        os.close(fd)
        try:
            m = create_pmesh_from_regions(SQUARE, TAGS,
                                          maximum_triangle_area=500.0,
                                          filename=name, size_field=field,
                                          verbose=False)
            domain = anuga.create_domain_from_file(name)
        finally:
            os.remove(name)
        assert len(domain) == len(m.tri_mesh.triangles)
        assert len(domain) > 200

    def test_create_domain_from_regions(self):
        """size_field and relax_passes reach the triangulation of a domain."""

        field = size_field_from_regions(
            SQUARE, [([[40, 40], [60, 40], [60, 60], [40, 60]], 2.0)], 500.0)

        plain = anuga.create_domain_from_regions(
            SQUARE, TAGS, maximum_triangle_area=500.0)
        fielded = anuga.create_domain_from_regions(
            SQUARE, TAGS, maximum_triangle_area=500.0, size_field=field)
        relaxed = anuga.create_domain_from_regions(
            SQUARE, TAGS, maximum_triangle_area=500.0, size_field=field,
            relax_passes=30)

        assert len(fielded) > 5 * len(plain)
        assert len(relaxed) == len(fielded)
        assert relaxed.radii.min() > fielded.radii.min()
        relaxed.check_integrity()
        assert num.isclose(relaxed.areas.sum(), 100.0 * 100.0)
        assert set(relaxed.get_boundary_tags()) == set(TAGS)

        basic = anuga.create_basic_mesh_from_regions(
            SQUARE, TAGS, maximum_triangle_area=500.0, size_field=field)
        assert basic.number_of_triangles == len(fielded)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_size_field)
    runner = unittest.TextTestRunner()
    runner.run(suite)
