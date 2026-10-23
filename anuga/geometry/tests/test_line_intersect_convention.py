"""line_intersect selects by geometry, not by rounding (#231).

A triangle is selected when the segment passes through or along it with
positive length; a point contact (a vertex, or an edge touched from
outside at one point) is not a selection. Along a shared edge both
triangles are selected. The same answer must come out whichever way the
last bit of a coordinate fell, and from whichever compiler built the
extension.
"""
import numpy as np

import anuga
from anuga import rectangular_cross_domain
from anuga.geometry.polygon import line_intersect


def mesh(dx=5.0, dy=2.5):
    """40 x 20 rectangles of dx x dy, four triangles each of area dx*dy/4."""
    d = rectangular_cross_domain(40, 20, len1=40 * dx, len2=20 * dy)
    return d, d.get_vertex_coordinates(absolute=True), d.areas


def select(vc, line):
    return sorted(line_intersect(vc, np.array(line, dtype=float)).tolist())


def test_a_segment_along_a_mesh_line_selects_both_sides_symmetrically():
    d, vc, areas = mesh()
    # x = 60 is a mesh line; y from 20 to 30 spans 4 rectangles of height 2.5
    idx = select(vc, [[60.0, 20.0], [60.0, 30.0]])
    assert len(idx) == 8, idx                      # one triangle per side per rectangle
    x = d.centroid_coordinates[idx, 0]
    assert (x < 60.0).sum() == 4 and (x > 60.0).sum() == 4
    assert np.isclose(areas[idx].sum(), 8 * 3.125)
    # every selected triangle has an edge on x = 60: two vertices at x = 60
    for k in idx:
        assert (np.abs(vc[3 * k:3 * k + 3, 0] - 60.0) < 1e-12).sum() == 2


def test_the_selection_does_not_move_with_the_last_bit():
    """Within rounding of a mesh line is on it, on either side."""
    d, vc, areas = mesh()
    base = select(vc, [[60.0, 20.0], [60.0, 30.0]])
    for eps in (1e-12, -1e-12, 1e-11, -1e-11):
        assert select(vc, [[60.0 + eps, 20.0], [60.0 + eps, 30.0]]) == base, eps


def test_an_off_grid_segment_selects_what_it_crosses():
    d, vc, areas = mesh()
    # x = 61 crosses 4 rectangles; in each it cuts the left, top and bottom
    # triangles (never the right one, which starts at x = 62.5)
    idx = select(vc, [[61.0, 20.0], [61.0, 30.0]])
    assert len(idx) == 12, idx
    assert np.isclose(areas[idx].sum(), 12 * 3.125)


def test_point_contacts_are_not_selections():
    d, vc, areas = mesh()
    # A segment ending exactly at the mesh vertex (60, 20) from the left
    # touches the triangles east of x = 60 only at that point
    idx = select(vc, [[55.0, 20.0], [60.0, 20.0]])
    assert len(idx) == 2, idx                      # the two along y = 20 in that rectangle
    assert np.all(d.centroid_coordinates[idx, 0] < 60.0)
    # A segment across a rectangle's diagonal vertex (its centre) passes
    # through two of its four triangles and touches the other two at a point
    idx = select(vc, [[60.0, 21.25], [65.0, 21.25]])   # centre (62.5, 21.25)
    assert len(idx) == 2, idx


def test_a_segment_inside_one_triangle_selects_it():
    d, vc, areas = mesh()
    # the left triangle of rectangle [60,65]x[20,22.5] has vertices
    # (60,20), (62.5,21.25), (60,22.5); this segment lies strictly inside it
    idx = select(vc, [[60.5, 21.0], [61.0, 21.5]])
    assert len(idx) == 1, idx


def test_expected_count_along_a_long_mesh_line():
    """count = 2 sides x (length / rectangle width) for a line on a mesh line."""
    d, vc, areas = mesh()
    idx = select(vc, [[60.0, 25.0], [140.0, 25.0]])   # y = 25 is a mesh line
    assert len(idx) == 2 * 16, len(idx)
    assert np.isclose(areas[idx].sum(), 32 * 3.125)


def test_a_zero_length_segment_selects_nothing():
    """The side edges of a structure's apron polygon when apron = 0. Treated
    as a point test they pulled in the triangle across the channel bank at
    each end of a bridge's exchange line and changed its discharge."""
    d, vc, areas = mesh()
    assert select(vc, [[60.0, 20.0], [60.0, 20.0]]) == []
    assert select(vc, [[61.0, 21.0], [61.0, 21.0]]) == []      # strictly inside a triangle


def test_a_zero_apron_polygon_selects_what_its_line_does():
    """Region(poly=[p0, p1, p1, p0], expand_polygon=True) must equal
    Region(line=[p0, p1]): a bridge with apron = 0 built exactly that."""
    from anuga import Region, rectangular_cross_domain
    d = rectangular_cross_domain(40, 20, len1=200.0, len2=50.0)
    line = [[60.01, 19.99], [139.99, 19.99]]          # a hair off the mesh lines
    by_line = Region(d, poly=line, expand_polygon=True)
    degenerate = np.array([line[0], line[1], line[1], line[0]])
    by_poly = Region(d, poly=degenerate, expand_polygon=True)
    assert sorted(by_poly.indices.tolist()) == sorted(by_line.indices.tolist())
    assert len(by_line.indices) > 0

