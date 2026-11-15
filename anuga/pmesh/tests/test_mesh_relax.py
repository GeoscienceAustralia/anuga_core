#!/usr/bin/env python
"""Worst-cell relaxation of a triangulation."""

import unittest

import numpy as num
import pytest

from anuga.pmesh.mesh_interface import create_pmesh_from_regions
from anuga.pmesh.mesh_relax import relax_worst_cells, segment_sliders
from anuga.pmesh.mesh_relax_ext import cell_quality

pytest.importorskip('meshpy')


SQUARE = [[0.0, 0.0], [100.0, 0.0], [100.0, 100.0], [0.0, 100.0]]
TAGS = {'bottom': [0], 'right': [1], 'top': [2], 'left': [3]}


def make_pmesh(maximum_triangle_area=4.0, **kwargs):
    m = create_pmesh_from_regions(SQUARE, TAGS,
                                  maximum_triangle_area=maximum_triangle_area,
                                  verbose=False, **kwargs)
    return m


def arrays(m):
    v = num.array(m.tri_mesh.vertices, dtype=float, order='C')
    t = num.ascontiguousarray(m.tri_mesh.triangles, dtype=num.int64)
    return v, t


def signed_areas(v, t):
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    return 0.5 * ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
                  - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))


class Test_mesh_relax(unittest.TestCase):

    def test_cell_quality(self):
        # equilateral triangle of edge 2: radii == inradius == 1/sqrt(3)
        v = num.array([[0.0, 0.0], [2.0, 0.0], [1.0, num.sqrt(3.0)]])
        t = num.array([[0, 1, 2]], dtype=num.int64)
        assert num.allclose(cell_quality(v, t, 1.0), 1.0 / num.sqrt(3.0))
        # the same triangle read clockwise is inverted
        assert cell_quality(v, t, -1.0)[0] == -1.0
        assert num.allclose(cell_quality(v, t[:, ::-1].copy(), -1.0),
                            1.0 / num.sqrt(3.0))

        # 3-4-5 triangle: inradius 1; radii is a third of the shortest
        # median... the shortest centroid to edge-midpoint distance
        v = num.array([[0.0, 0.0], [4.0, 0.0], [0.0, 3.0]])
        g = v.mean(axis=0)
        mids = 0.5 * (v + num.roll(v, -1, axis=0))
        radii = num.sqrt(((mids - g) ** 2).sum(axis=1)).min()
        assert num.allclose(cell_quality(v, t, 1.0), min(radii, 1.0))

        # a sliver has good "radii" for its area but a tiny inradius
        v = num.array([[0.0, 0.0], [10.0, 0.0], [5.0, 0.01]])
        assert cell_quality(v, t, 1.0)[0] < 0.01

    def test_relaxation_improves_worst_cell(self):
        m = make_pmesh()
        m.generate_mesh(verbose=False)
        v0, t0 = arrays(m)
        segments = num.array(m.tri_mesh.segments)
        q0 = cell_quality(v0, t0, 1.0)

        before, after = m.relax_worst_cells(passes=30, slide=False)
        v1, t1 = arrays(m)
        q1 = cell_quality(v1, t1, 1.0)

        assert num.isclose(before, q0.min())
        assert num.isclose(after, q1.min())
        assert after > 1.2 * before

        # topology is untouched
        assert num.array_equal(t0, t1)
        assert num.array_equal(segments, m.tri_mesh.segments)
        # no vertex on a segment moved, and some others did
        fixed = num.unique(segments)
        assert num.array_equal(v0[fixed], v1[fixed])
        assert num.any(v0 != v1)
        # no triangle is inverted and the domain is still covered exactly
        assert signed_areas(v1, t1).min() > 0.0
        assert num.isclose(signed_areas(v1, t1).sum(), 100.0 * 100.0)

    def test_generate_mesh_relax_passes(self):
        """relax_passes on generate_mesh is the same as relaxing afterwards."""

        m1 = make_pmesh()
        m1.generate_mesh(verbose=False)
        m1.relax_worst_cells(passes=10)

        m2 = make_pmesh()
        m2.generate_mesh(verbose=False, relax_passes=10)

        m3 = make_pmesh(relax_passes=10)       # remembered on the Pmesh
        m3.generate_mesh(verbose=False)

        assert num.array_equal(m1.tri_mesh.vertices, m2.tri_mesh.vertices)
        assert num.array_equal(m1.tri_mesh.vertices, m3.tri_mesh.vertices)

        m0 = make_pmesh()
        m0.generate_mesh(verbose=False)
        assert not num.array_equal(m0.tri_mesh.vertices, m1.tri_mesh.vertices)

    def test_breakline_vertices_do_not_move(self):
        breakline = [[20.0, 30.0], [80.0, 70.0]]
        for slide in (False, True):
            m = make_pmesh(breaklines=[breakline])
            m.generate_mesh(verbose=False)
            v0, _ = arrays(m)
            m.relax_worst_cells(passes=30, slide=slide)
            v1, t1 = arrays(m)

            # vertices on the line  2x - 3y + 50 = 0  between its ends
            on_line = ((num.abs(2.0 * v0[:, 0] - 3.0 * v0[:, 1] + 50.0)
                        < 1e-9) & (v0[:, 0] >= 20.0) & (v0[:, 0] <= 80.0))
            assert on_line.sum() > 10
            if slide:
                # inserted vertices may slide, but stay on the line and
                # between its ends, which do not move
                w = v1[on_line]
                assert num.abs(2.0 * w[:, 0] - 3.0 * w[:, 1] + 50.0).max() \
                    < 1e-9
                assert w[:, 0].min() >= 20.0 and w[:, 0].max() <= 80.0
                n_in = m.number_of_input_vertices
                assert num.array_equal(v0[:n_in], v1[:n_in])
            else:
                assert num.array_equal(v0[on_line], v1[on_line])
            assert signed_areas(v1, t1).min() > 0.0

    def test_function_leaves_input_alone(self):
        m = make_pmesh()
        m.generate_mesh(verbose=False)
        v0, t0 = arrays(m)
        keep = v0.copy()
        fixed = num.unique(m.tri_mesh.segments)

        v1, before, after = relax_worst_cells(v0, t0, fixed, passes=5)
        assert num.array_equal(v0, keep)
        assert v1 is not v0
        assert after >= before

        # no passes: nothing happens
        v2, _, _ = relax_worst_cells(v0, t0, fixed, passes=0)
        assert num.array_equal(v2, keep)

        # everything fixed: nothing can move
        v3, before, after = relax_worst_cells(v0, t0,
                                              num.arange(len(v0)), passes=5)
        assert num.array_equal(v3, keep)
        assert before == after

    def test_clockwise_triangles(self):
        m = make_pmesh()
        m.generate_mesh(verbose=False)
        v0, t0 = arrays(m)
        fixed = num.unique(m.tri_mesh.segments)

        v1, b1, a1 = relax_worst_cells(v0, t0, fixed, passes=10)
        v2, b2, a2 = relax_worst_cells(v0, t0[:, ::-1].copy(), fixed,
                                       passes=10)
        assert num.isclose(b1, b2)
        assert a2 > b2
        assert signed_areas(v2, t0).min() > 0.0

    def test_segment_sliders(self):
        # A straight segment 0-3 split at 4 and 5, then a corner at 3 and a
        # second segment 3-6 split at 7. Vertices 0-3 and 6 are input.
        v = num.array([[0.0, 0.0], [9.0, 9.0], [9.0, 0.0], [6.0, 0.0],
                       [2.0, 0.0], [4.0, 0.0], [6.0, 4.0], [6.0, 2.0]])
        seg = num.array([[0, 4], [4, 5], [5, 3], [3, 7], [7, 6]])
        def as_set(rows):                 # the order of a and b is free
            return {(int(i), frozenset((int(a), int(b)))) for i, a, b in rows}

        rows = segment_sliders(v, seg, first_inserted=4)
        assert as_set(rows) == {(4, frozenset((0, 5))), (5, frozenset((3, 4))),
                                (7, frozenset((3, 6)))}
        # treat everything as input: nothing slides
        assert len(segment_sliders(v, seg, first_inserted=8)) == 0
        # a bend at an inserted vertex: it may not slide
        v[5] = [4.0, 0.5]
        assert as_set(segment_sliders(v, seg, 4)) == {(7, frozenset((3, 6)))}

    def test_sliding_keeps_segment_shape(self):
        hole = [[40.0, 40.0], [60.0, 40.0], [60.0, 42.0], [40.0, 42.0]]
        breakline = [[20.0, 30.0], [80.0, 70.0]]
        m = make_pmesh(maximum_triangle_area=8.0, interior_holes=[hole],
                       breaklines=[breakline])
        m.generate_mesh(verbose=False)
        v0, t0 = arrays(m)
        segments = num.array(m.tri_mesh.segments)
        n_in = m.number_of_input_vertices
        rows = segment_sliders(v0, segments, n_in)
        assert len(rows) > 10

        before, after = m.relax_worst_cells(passes=30)
        v1, t1 = arrays(m)
        assert after >= before
        assert num.array_equal(t0, t1)
        assert num.array_equal(segments, m.tri_mesh.segments)

        # input vertices (corners, breakline ends) are where they were
        assert num.array_equal(v0[:n_in], v1[:n_in])
        # segment vertices that moved are sliders, still on their segment
        on_seg = num.unique(segments)
        moved = on_seg[num.any(v0[on_seg] != v1[on_seg], axis=1)]
        assert len(moved) > 0
        assert set(moved) <= set(rows[:, 0])
        i, a, b = rows[num.isin(rows[:, 0], moved)].T
        ab, ai = v0[b] - v0[a], v1[i] - v0[a]
        cross = ab[:, 0] * ai[:, 1] - ab[:, 1] * ai[:, 0]
        assert num.abs(cross).max() < 1e-9 * (ab ** 2).sum(axis=1).max()
        # the hole is still a hole of the same area, and nothing inverted
        assert signed_areas(v1, t1).min() > 0.0
        assert num.isclose(signed_areas(v1, t1).sum(),
                           100.0 * 100.0 - 20.0 * 2.0)

    def test_sliding_helps_between_buildings(self):
        """Dense buildings leave few free vertices: there, sliding the points
        the generator put on the buildings' sides is what can still move."""
        rng = num.random.default_rng(0)
        holes, y = [], 10.0
        while y < 70.0:
            x, depth = 10.0, rng.uniform(8.0, 14.0)
            while x < 85.0:
                w = rng.uniform(6.0, 14.0)
                d = depth + rng.uniform(-1.0, 1.0)
                if rng.integers(2):
                    w = rng.uniform(1.9, 3.0)          # a narrow wing
                holes.append([[x, y], [x + w, y], [x + w, y + d], [x, y + d]])
                x += w + rng.uniform(2.0, 5.0)
            y += depth + rng.uniform(3.0, 6.0)

        result = {}
        for slide in (False, True):
            m = make_pmesh(maximum_triangle_area=40.0, interior_holes=holes)
            m.generate_mesh(verbose=False)
            v0, _ = arrays(m)
            _, result[slide] = m.relax_worst_cells(passes=30, slide=slide)
            v1, t1 = arrays(m)
            n_in = m.number_of_input_vertices
            assert num.array_equal(v0[:n_in], v1[:n_in])
            assert signed_areas(v1, t1).min() > 0.0
        assert result[True] >= result[False]

    def test_relax_before_generate_raises(self):
        m = make_pmesh()
        with self.assertRaises(Exception):
            m.relax_worst_cells()


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_mesh_relax)
    runner = unittest.TextTestRunner()
    runner.run(suite)
