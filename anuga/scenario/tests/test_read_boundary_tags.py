"""
Unit tests for anuga.scenario.read_boundary_tags_line_shapefile

Tests cover:
  - parse_ogr_info_text  (pure Python text parser — no GIS deps)
  - read_boundary_tags_line_shapefile  (CSV path only — no fiona required)
  - read_boundary_tags_line_shapefile  (shapefile / fiona path)
"""
import os
import tempfile
import shutil
import unittest

try:
    import fiona
    HAS_FIONA = True
except ImportError:
    HAS_FIONA = False

try:
    from anuga.scenario.read_boundary_tags_line_shapefile import (
        parse_ogr_info_text,
        read_boundary_tags_line_shapefile,
    )
    HAS_MODULE = True
    SKIP_REASON = ''
except ImportError as _e:
    HAS_MODULE = False
    SKIP_REASON = str(_e)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ogr_lines(tag_attr, features):
    """Build a synthetic ogrinfo -al output as a list of strings.

    Each element of *features* is (tag_value, coords) where coords is a list
    of (x, y) tuples forming one LINESTRING.
    """
    lines = []
    for i, (tag_val, coords) in enumerate(features):
        lines.append('OGRFeature(boundary):' + str(i))
        lines.append('  ' + tag_attr + ' (String) = ' + tag_val)
        coord_str = ','.join('%s %s' % (x, y) for x, y in coords)
        lines.append('  LINESTRING (' + coord_str + ')')
        lines.append('')          # blank line between features
    return lines


# ---------------------------------------------------------------------------
# Tests for parse_ogr_info_text
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestParseOgrInfoText(unittest.TestCase):

    TAG = 'bndryTag'

    def test_single_feature_tag_extracted(self):
        lines = _ogr_lines(self.TAG, [('ocean', [(0, 0), (1, 0)])])
        result = parse_ogr_info_text(lines, self.TAG)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0][1], 'ocean')

    def test_single_feature_coordinates_parsed(self):
        lines = _ogr_lines(self.TAG, [('ocean', [(0.0, 1.0), (2.0, 3.0), (4.0, 5.0)])])
        result = parse_ogr_info_text(lines, self.TAG)
        coords = result[0][0]
        self.assertEqual(len(coords), 3)
        self.assertAlmostEqual(coords[0][0], 0.0)
        self.assertAlmostEqual(coords[0][1], 1.0)
        self.assertAlmostEqual(coords[2][0], 4.0)
        self.assertAlmostEqual(coords[2][1], 5.0)

    def test_two_features_both_returned(self):
        lines = _ogr_lines(self.TAG, [
            ('ocean', [(0, 0), (1, 0)]),
            ('land',  [(1, 0), (1, 1)]),
        ])
        result = parse_ogr_info_text(lines, self.TAG)
        self.assertEqual(len(result), 2)
        tags = [r[1] for r in result]
        self.assertIn('ocean', tags)
        self.assertIn('land', tags)

    def test_non_ogr_lines_are_ignored(self):
        # Extra header lines before the first feature
        lines = [
            'INFO: Open of `test.shp` succeeded.',
            'Layer name: test',
        ] + _ogr_lines(self.TAG, [('river', [(5, 5), (6, 6)])])
        result = parse_ogr_info_text(lines, self.TAG)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0][1], 'river')

    def test_missing_tag_raises(self):
        # ogrinfo lines that contain an OGRFeature but cut off before the tag
        lines = ['OGRFeature(boundary):0']
        with self.assertRaises(Exception):
            parse_ogr_info_text(lines, self.TAG)

    def test_missing_linestring_raises(self):
        # Tag is present but no LINESTRING follows it
        lines = [
            'OGRFeature(boundary):0',
            '  bndryTag (String) = ocean',
            # no LINESTRING line
        ]
        with self.assertRaises(Exception):
            parse_ogr_info_text(lines, self.TAG)

    def test_wrong_tag_attribute_name_raises(self):
        lines = _ogr_lines('wrongAttr', [('ocean', [(0, 0), (1, 0)])])
        # Using 'bndryTag' as tag_attribute but lines have 'wrongAttr'
        with self.assertRaises(Exception):
            parse_ogr_info_text(lines, 'bndryTag')

    def test_three_features_order_preserved(self):
        features = [
            ('a', [(0, 0), (1, 0)]),
            ('b', [(1, 0), (1, 1)]),
            ('c', [(1, 1), (0, 1)]),
        ]
        lines = _ogr_lines(self.TAG, features)
        result = parse_ogr_info_text(lines, self.TAG)
        self.assertEqual([r[1] for r in result], ['a', 'b', 'c'])


# ---------------------------------------------------------------------------
# Tests for read_boundary_tags_line_shapefile (CSV path)
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestReadBoundaryTagsCSV(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _write_csv(self, name, points):
        """Write an x,y CSV polygon file and return its path."""
        path = os.path.join(self.tmpdir, name)
        with open(path, 'w') as f:
            f.write('x,y\n')
            for x, y in points:
                f.write('%s,%s\n' % (x, y))
        return path

    def test_csv_single_tag_returns_polygon_and_tags(self):
        pts = [(0, 0), (1, 0), (1, 1), (0, 1)]
        csv = self._write_csv('poly.csv', pts)
        explicit_tags = [{'tag': 'ocean', 'edges': [0, 1, 2, 3]}]
        poly, tags = read_boundary_tags_line_shapefile(csv, explicit_tags=explicit_tags)
        self.assertEqual(len(poly), 4)
        self.assertIn('ocean', tags)
        self.assertEqual(sorted(tags['ocean']), [0, 1, 2, 3])

    def test_csv_multiple_tags_edges_split(self):
        pts = [(0, 0), (1, 0), (1, 1), (0, 1)]
        csv = self._write_csv('poly.csv', pts)
        explicit_tags = [
            {'tag': 'ocean', 'edges': [0, 1]},
            {'tag': 'land',  'edges': [2, 3]},
        ]
        poly, tags = read_boundary_tags_line_shapefile(csv, explicit_tags=explicit_tags)
        self.assertIn('ocean', tags)
        self.assertIn('land', tags)
        self.assertEqual(sorted(tags['ocean']), [0, 1])
        self.assertEqual(sorted(tags['land']), [2, 3])

    def test_csv_duplicate_tag_edges_merged(self):
        pts = [(0, 0), (1, 0), (1, 1), (0, 1)]
        csv = self._write_csv('poly.csv', pts)
        explicit_tags = [
            {'tag': 'ocean', 'edges': [0, 1]},
            {'tag': 'ocean', 'edges': [2]},
        ]
        _, tags = read_boundary_tags_line_shapefile(csv, explicit_tags=explicit_tags)
        self.assertIn('ocean', tags)
        self.assertEqual(sorted(tags['ocean']), [0, 1, 2])

    def test_csv_without_explicit_tags_raises(self):
        pts = [(0, 0), (1, 0), (1, 1)]
        csv = self._write_csv('poly.csv', pts)
        with self.assertRaises(ValueError):
            read_boundary_tags_line_shapefile(csv)  # no explicit_tags

    def test_nonexistent_file_raises(self):
        with self.assertRaises(ValueError):
            read_boundary_tags_line_shapefile(
                os.path.join(self.tmpdir, 'no_such.csv'),
                explicit_tags=[{'tag': 'x', 'edges': [0]}])

    def test_polygon_coordinates_match_csv(self):
        pts = [(10.0, 20.0), (30.0, 20.0), (30.0, 40.0), (10.0, 40.0)]
        csv = self._write_csv('poly.csv', pts)
        explicit_tags = [{'tag': 't', 'edges': [0]}]
        poly, _ = read_boundary_tags_line_shapefile(csv, explicit_tags=explicit_tags)
        xs = [p[0] for p in poly]
        ys = [p[1] for p in poly]
        self.assertAlmostEqual(min(xs), 10.0)
        self.assertAlmostEqual(max(xs), 30.0)
        self.assertAlmostEqual(min(ys), 20.0)
        self.assertAlmostEqual(max(ys), 40.0)


# ---------------------------------------------------------------------------
# Helpers for shapefile creation
# ---------------------------------------------------------------------------

def _write_shapefile(path, features, tag_attr='bndryTag'):
    """Write a line shapefile using fiona.

    *features* is a list of (tag_value, coords) where coords is a list of
    (x, y) tuples forming one LINESTRING.  Returns *path*.
    """
    schema = {
        'geometry': 'LineString',
        'properties': {tag_attr: 'str'},
    }
    with fiona.open(path, 'w', driver='ESRI Shapefile', schema=schema) as dst:
        for tag_val, coords in features:
            dst.write({
                'geometry': {'type': 'LineString', 'coordinates': coords},
                'properties': {tag_attr: tag_val},
            })
    return path


# ---------------------------------------------------------------------------
# Tests for read_boundary_tags_line_shapefile (shapefile / fiona path)
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE and HAS_FIONA, SKIP_REASON or 'fiona not available')
class TestReadBoundaryTagsShapefile(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _shp(self, name='boundary.shp'):
        return os.path.join(self.tmpdir, name)

    # --- Two-segment polygon ---

    def test_two_segments_polygon_has_correct_point_count(self):
        # Square split into two L-shaped segments that chain together.
        # seg0: (0,0)→(1,0)→(1,1)  and  seg1: (1,1)→(0,1)→(0,0)
        # Assembly produces a 4-point polygon (closure point dropped).
        shp = _write_shapefile(self._shp(), [
            ('east', [(0, 0), (1, 0), (1, 1)]),
            ('west', [(1, 1), (0, 1), (0, 0)]),
        ])
        poly, _ = read_boundary_tags_line_shapefile(shp)
        self.assertEqual(len(poly), 4)

    def test_two_segments_both_tags_present(self):
        shp = _write_shapefile(self._shp(), [
            ('east', [(0, 0), (1, 0), (1, 1)]),
            ('west', [(1, 1), (0, 1), (0, 0)]),
        ])
        _, tags = read_boundary_tags_line_shapefile(shp)
        self.assertIn('east', tags)
        self.assertIn('west', tags)

    def test_two_segments_edge_counts(self):
        # 'east' segment has 2 edges (3 points), 'west' has 2 edges
        shp = _write_shapefile(self._shp(), [
            ('east', [(0, 0), (1, 0), (1, 1)]),
            ('west', [(1, 1), (0, 1), (0, 0)]),
        ])
        _, tags = read_boundary_tags_line_shapefile(shp)
        self.assertEqual(len(list(tags['east'])), 2)
        self.assertEqual(len(list(tags['west'])), 2)

    def test_two_segments_coordinates_match(self):
        shp = _write_shapefile(self._shp(), [
            ('east', [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]),
            ('west', [(1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]),
        ])
        poly, _ = read_boundary_tags_line_shapefile(shp)
        xs = [p[0] for p in poly]
        ys = [p[1] for p in poly]
        self.assertAlmostEqual(min(xs), 0.0)
        self.assertAlmostEqual(max(xs), 1.0)
        self.assertAlmostEqual(min(ys), 0.0)
        self.assertAlmostEqual(max(ys), 1.0)

    # --- Four-segment square ---

    def test_four_segments_polygon_has_four_points(self):
        # A unit square split into four single-edge segments.
        shp = _write_shapefile(self._shp(), [
            ('bottom', [(0, 0), (1, 0)]),
            ('right',  [(1, 0), (1, 1)]),
            ('top',    [(1, 1), (0, 1)]),
            ('left',   [(0, 1), (0, 0)]),
        ])
        poly, _ = read_boundary_tags_line_shapefile(shp)
        self.assertEqual(len(poly), 4)

    def test_four_segments_all_tags_present(self):
        shp = _write_shapefile(self._shp(), [
            ('bottom', [(0, 0), (1, 0)]),
            ('right',  [(1, 0), (1, 1)]),
            ('top',    [(1, 1), (0, 1)]),
            ('left',   [(0, 1), (0, 0)]),
        ])
        _, tags = read_boundary_tags_line_shapefile(shp)
        for name in ('bottom', 'right', 'top', 'left'):
            self.assertIn(name, tags)

    def test_four_segments_each_tag_has_one_edge(self):
        shp = _write_shapefile(self._shp(), [
            ('bottom', [(0, 0), (1, 0)]),
            ('right',  [(1, 0), (1, 1)]),
            ('top',    [(1, 1), (0, 1)]),
            ('left',   [(0, 1), (0, 0)]),
        ])
        _, tags = read_boundary_tags_line_shapefile(shp)
        for name in ('bottom', 'right', 'top', 'left'):
            self.assertEqual(len(list(tags[name])), 1,
                             msg='Expected 1 edge for tag ' + name)

    def test_four_segments_edge_indices_cover_all_edges(self):
        shp = _write_shapefile(self._shp(), [
            ('bottom', [(0, 0), (1, 0)]),
            ('right',  [(1, 0), (1, 1)]),
            ('top',    [(1, 1), (0, 1)]),
            ('left',   [(0, 1), (0, 0)]),
        ])
        _, tags = read_boundary_tags_line_shapefile(shp)
        all_edges = sorted(e for edges in tags.values() for e in edges)
        self.assertEqual(all_edges, [0, 1, 2, 3])

    # --- Single closed LineString ---

    def test_single_closed_segment(self):
        # A single LineString where first == last point is treated as a
        # closed polygon and the duplicate closure point is stripped.
        shp = _write_shapefile(self._shp(), [
            ('ocean', [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]),
        ])
        poly, tags = read_boundary_tags_line_shapefile(shp)
        self.assertEqual(len(poly), 4)
        self.assertIn('ocean', tags)
        self.assertEqual(len(list(tags['ocean'])), 4)

    # --- Custom tag attribute ---

    def test_custom_tag_attribute_name(self):
        shp = _write_shapefile(self._shp(), [
            ('open',   [(0, 0), (1, 0), (1, 1)]),
            ('closed', [(1, 1), (0, 1), (0, 0)]),
        ], tag_attr='myTag')
        poly, tags = read_boundary_tags_line_shapefile(shp, tag_attribute='myTag')
        self.assertEqual(len(poly), 4)
        self.assertIn('open', tags)
        self.assertIn('closed', tags)

    # --- Error cases ---

    def test_wrong_tag_attribute_raises(self):
        shp = _write_shapefile(self._shp(), [
            ('ocean', [(0, 0), (1, 0), (1, 1)]),
            ('land',  [(1, 1), (0, 1), (0, 0)]),
        ], tag_attr='bndryTag')
        # Using the wrong attribute name should raise ValueError
        with self.assertRaises(ValueError):
            read_boundary_tags_line_shapefile(shp, tag_attribute='nonexistent')

    def test_unchained_segments_raises(self):
        # Two segments whose endpoints do not connect — assembly should fail
        shp = _write_shapefile(self._shp(), [
            ('a', [(0, 0), (1, 0)]),
            ('b', [(5, 5), (6, 6)]),  # starts nowhere near (1,0)
        ])
        with self.assertRaises(Exception):
            read_boundary_tags_line_shapefile(shp)


if __name__ == '__main__':
    unittest.main()
