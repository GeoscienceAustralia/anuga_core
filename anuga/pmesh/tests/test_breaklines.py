#!/usr/bin/env python
"""clean_breaklines and the mesh quality check (sliver triangles)."""

import unittest
import warnings

import numpy as num
import pytest

import anuga
from anuga.pmesh.breaklines import clean_breaklines

pytest.importorskip('meshpy')
pytest.importorskip('shapely')

from shapely.geometry import LineString, Point  # noqa: E402


SQUARE = [[0.0, 0.0], [1000.0, 0.0], [1000.0, 1000.0], [0.0, 1000.0]]
TAGS = {'bottom': [0], 'right': [1], 'top': [2], 'left': [3]}

# Three lines meeting badly: two cross at a shallow angle, and a third ends
# 0.3 m short of the first. Meshed as given they make sliver triangles.
CROSSING = [[[100.0, 500.0], [900.0, 500.0]],
            [[150.0, 470.0], [850.0, 530.0]],
            [[500.0, 900.0], [500.0, 500.3]]]


def make_domain(breaklines=None, interior_regions=None, area=2000.0):
    return anuga.create_domain_from_regions(
        SQUARE, TAGS, maximum_triangle_area=area,
        interior_regions=interior_regions, breaklines=breaklines,
        use_cache=False, verbose=False)


def min_separation(lines):
    """Smallest distance between two lines, more than 25 m from a vertex
    they share (the joining segment itself runs up to the other line)."""
    smallest = num.inf
    for i, a in enumerate(lines):
        for b in lines[i + 1:]:
            shared = [p for p in a if p in b]
            la, lb = LineString(a), LineString(b)
            for s in num.linspace(0.0, la.length, 400):
                p = la.interpolate(s)
                if all(p.distance(Point(q)) > 25.0 for q in shared):
                    smallest = min(smallest, lb.distance(p))
    return smallest


class Test_mesh_quality(unittest.TestCase):

    def test_report_on_uniform_mesh(self):
        domain = anuga.rectangular_cross_domain(4, 4, len1=40.0, len2=40.0)
        report = domain.check_mesh_quality()
        assert report['number_of_triangles'] == 64
        assert num.isclose(report['area_min'], 25.0)
        assert num.isclose(report['area_median'], 25.0)
        assert num.isclose(report['radius_ratio'], 1.0)
        assert len(report['failed']) == 0
        assert len(report['smallest']) == 10
        assert 'Mesh quality: 64 triangles' in report['message']

    def test_limits(self):
        domain = anuga.rectangular_cross_domain(4, 4, len1=40.0, len2=40.0)

        # every triangle is 25 m^2: a 30 m^2 limit fails them all
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            report = domain.check_mesh_quality(minimum_triangle_area=30.0)
        assert len(report['failed']) == 64
        assert len(w) == 1 and 'below the limit' in str(w[0].message)

        with self.assertRaises(anuga.MeshQualityError):
            domain.check_mesh_quality(minimum_triangle_area=30.0,
                                      action='raise')
        with self.assertRaises(anuga.MeshQualityError):
            domain.check_mesh_quality(minimum_radius=100.0, action='raise')

        with warnings.catch_warnings():
            warnings.simplefilter('error')
            report = domain.check_mesh_quality(minimum_triangle_area=30.0,
                                               action='none')
            assert len(report['failed']) == 64
            report = domain.check_mesh_quality(minimum_triangle_area=20.0,
                                               minimum_radius=1.0,
                                               action='raise')
            assert len(report['failed']) == 0

        with self.assertRaises(ValueError):
            domain.check_mesh_quality(action='ignore')

    def test_failed_triangles_are_smallest_first_and_located(self):
        domain = make_domain(breaklines=CROSSING)
        report = domain.check_mesh_quality(minimum_triangle_area=5.0,
                                           action='none')
        failed = report['failed']
        assert len(failed) > 0
        radii = domain.get_radii()
        assert num.all(num.diff(radii[failed]) >= 0.0)
        assert failed[0] == num.argmin(radii)
        # the slivers sit where the lines meet, near y = 500
        centroids = domain.get_centroid_coordinates(absolute=True)
        assert num.all(num.abs(centroids[failed, 1] - 500.0) < 20.0)


class Test_clean_breaklines(unittest.TestCase):

    def test_raw_lines_make_slivers_cleaned_lines_do_not(self):
        raw = make_domain(breaklines=CROSSING)
        assert raw.get_areas().min() < 1.0
        with self.assertRaises(anuga.MeshQualityError):
            raw.check_mesh_quality(minimum_triangle_area=5.0, action='raise')

        lines, regions = clean_breaklines(CROSSING, gap=20.0)
        assert regions is None
        domain = make_domain(breaklines=lines)
        report = domain.check_mesh_quality(minimum_triangle_area=5.0,
                                           action='raise')
        assert report['area_min'] > 20.0
        assert report['radius_ratio'] < 20.0

    def test_lines_are_kept_apart_and_joined(self):
        lines, _ = clean_breaklines(CROSSING, gap=20.0)

        # the longest line is kept whole (it gains the join vertices)
        assert lines[0][0] == [100.0, 500.0]
        assert lines[0][-1] == [900.0, 500.0]

        # away from shared vertices nothing is closer than about the gap
        assert min_separation(lines) > 0.9 * 20.0

        # every other piece ends on a vertex of an earlier line
        for i, line in enumerate(lines[1:], start=1):
            earlier = [p for other in lines[:i] for p in other]
            assert line[0] in earlier or line[-1] in earlier

        # the crossing line is cut in two, one piece on each side
        sides = [num.sign(line[len(line) // 2][1] - 500.0)
                 for line in lines[1:]]
        assert -1.0 in sides and 1.0 in sides

    def test_no_join_leaves_a_gap(self):
        lines, _ = clean_breaklines(CROSSING, gap=20.0, join=False)
        smallest = min(LineString(a).distance(LineString(b))
                       for i, a in enumerate(lines) for b in lines[i + 1:])
        assert smallest > 0.9 * 20.0

    def test_short_pieces_dropped_and_vertices_thinned(self):
        wiggly = [[100.0 + i, 200.0 + 0.2 * (i % 2)] for i in range(301)]
        short = [[600.0, 800.0], [620.0, 800.0]]
        lines, _ = clean_breaklines([wiggly, short], gap=10.0)
        assert len(lines) == 1
        assert lines[0][0] == wiggly[0] and lines[0][-1] == wiggly[-1]
        steps = num.hypot(*num.diff(num.array(lines[0]), axis=0).T)
        assert steps.min() >= 6.0 - 1.0e-9
        assert len(lines[0]) < 60

    def test_line_that_only_just_comes_within_the_gap(self):
        # The shorter line runs 30 m from the longer one except at one vertex,
        # where it comes 19.5 m from it: a stretch about a metre long is cut
        # out. The two pieces then end a metre apart, and the second must end
        # on the first's vertex, not be joined to it by a metre-long segment.
        long_line = [[100.0, 500.0], [900.0, 500.0]]
        dipping = [[200.0, 530.0], [480.0, 530.0], [500.0, 519.5],
                   [520.0, 530.0], [800.0, 530.0]]
        lines, _ = clean_breaklines([long_line, dipping], gap=20.0)
        assert len(lines) == 3
        for line in lines:
            steps = num.hypot(*num.diff(num.array(line), axis=0).T)
            assert steps.min() > 0.5 * 20.0, steps.min()
        domain = make_domain(breaklines=lines)
        domain.check_mesh_quality(minimum_triangle_area=5.0, action='raise')

    def test_self_crossing_line(self):
        loop = [[200.0, 200.0], [600.0, 600.0], [600.0, 200.0],
                [200.0, 600.0]]
        lines, _ = clean_breaklines([loop], gap=20.0)
        assert len(lines) >= 2
        domain = make_domain(breaklines=lines)
        domain.check_mesh_quality(minimum_triangle_area=5.0, action='raise')

    def test_line_crossing_an_interior_region(self):
        region = [[300.0, 300.0], [700.0, 300.0], [700.0, 700.0],
                  [300.0, 700.0]]
        # crosses the outline twice, and runs 2 m outside it along the top
        lines = [[[100.0, 520.0], [900.0, 480.0]],
                 [[320.0, 702.0], [680.0, 702.0]]]
        interior = [[region, 500.0]]

        raw = make_domain(breaklines=lines, interior_regions=interior)
        assert raw.get_areas().min() < 5.0

        new_lines, new_regions = clean_breaklines(lines, interior, SQUARE,
                                                  gap=20.0)
        # the line beside the outline is gone, the crossing one is in pieces
        # joined to the outline, which has gained those vertices
        outline = new_regions[0][0]
        assert new_regions[0][1] == 500.0
        assert len(outline) > 4
        for p in region:
            assert p in outline
        assert all(abs(line[0][1] - 702.0) > 1.0 for line in new_lines)
        joined = [p for line in new_lines for p in (line[0], line[-1])
                  if p in outline]
        assert len(joined) >= 2
        assert interior[0][0] == region            # input not modified

        domain = make_domain(breaklines=new_lines, interior_regions=new_regions)
        report = domain.check_mesh_quality(minimum_triangle_area=5.0,
                                           action='raise')
        assert report['area_min'] > 20.0

    def test_bounding_polygon_is_kept_clear_and_unchanged(self):
        lines = [[[-50.0, 400.0], [1050.0, 400.0]]]
        new_lines, regions = clean_breaklines(lines, None, SQUARE, gap=25.0)
        assert regions is None
        assert len(new_lines) == 1
        xs = [p[0] for p in new_lines[0]]
        assert min(xs) >= 24.0 and max(xs) <= 976.0
        assert min(xs) < 26.0 and max(xs) > 974.0

    def test_riverwalls_keep_their_heights(self):
        # two crossing walls with sloping crests, and a breakline that the
        # shorter wall would otherwise have to give way to
        walls = {'long': [[100.0, 500.0, 10.0], [900.0, 500.0, 18.0]],
                 'short': [[500.0, 200.0, 5.0], [500.0, 800.0, 11.0]]}
        par = {'short': {'Qfactor': 0.5}}
        lines = [[[100.0, 300.0], [950.0, 300.0]]]

        new_lines, regions, new_walls, new_par = clean_breaklines(
            lines, None, SQUARE, riverwalls=walls, riverwallPar=par, gap=20.0)
        assert regions is None

        # the long wall is whole; the short one is in two pieces either side
        assert sorted(new_walls) == ['long', 'short', 'short_2']
        assert new_par == {'short': {'Qfactor': 0.5},
                           'short_2': {'Qfactor': 0.5}}
        assert new_walls['long'][0] == [100.0, 500.0, 10.0]
        assert new_walls['long'][-1] == [900.0, 500.0, 18.0]

        # heights follow the original crests: z = 10 + (x - 100) / 100 on
        # the long wall, z = 5 + (y - 200) / 100 on the short one
        for x, y, z in new_walls['long']:
            assert y == 500.0 and num.isclose(z, 10.0 + (x - 100.0) / 100.0)
        for name in ('short', 'short_2'):
            for x, y, z in new_walls[name]:
                assert x == 500.0 and num.isclose(z, 5.0 + (y - 200.0) / 100.0)

        # both pieces of the short wall end on the same vertex of the long
        # one, so the barrier is closed (the crest heights differ there)
        on_long = [p[:2] for p in new_walls['long']]
        ends = [p[:2] for name in ('short', 'short_2')
                for p in (new_walls[name][0], new_walls[name][-1])]
        assert ends.count([500.0, 500.0]) == 2
        assert [500.0, 500.0] in on_long

        # walls come first: the breakline is the one cut, around the short wall
        assert len(new_lines) == 2
        assert all(len(p) == 2 for line in new_lines for p in line)
        short_xy = [p[:2] for p in new_walls['short'] + new_walls['short_2']]
        assert [500.0, 300.0] in short_xy
        assert sum([500.0, 300.0] in line for line in new_lines) == 2

        # and the whole lot meshes cleanly with every wall on mesh edges
        wall_lines = [[p[:2] for p in wall] for wall in new_walls.values()]
        domain = make_domain(breaklines=new_lines + wall_lines)
        domain.check_mesh_quality(minimum_triangle_area=5.0, action='raise')
        domain.riverwallData.create_riverwalls(new_walls, new_par,
                                               verbose=False)
        edge_length = domain.edgelengths.flatten()[
            domain.riverwallData.riverwall_edges]
        # each wall edge is listed once for each of its two triangles
        wall_length = sum(LineString(w).length for w in wall_lines)
        assert num.isclose(edge_length.sum(), 2.0 * wall_length, rtol=1.0e-6)

    def test_riverwall_arguments(self):
        with self.assertRaises(ValueError):
            clean_breaklines([], riverwalls={'w': [[0.0, 0.0], [100.0, 0.0]]})
        # no walls given: two values back; walls given: four
        assert len(clean_breaklines(CROSSING, gap=20.0)) == 2
        out = clean_breaklines([], riverwalls={}, gap=20.0)
        assert out == ([], None, {}, {})

    def test_bad_arguments(self):
        with self.assertRaises(ValueError):
            clean_breaklines(CROSSING, gap=0.0)
        # degenerate lines are ignored
        lines, _ = clean_breaklines([[[1.0, 1.0]], [[2.0, 2.0], [2.0, 2.0]]])
        assert lines == []


if __name__ == '__main__':
    unittest.main()
