"""GeoPackage in and out (#39, #40)."""
import os
import tempfile
import unittest

import numpy as np
import pytest

import anuga

fiona = pytest.importorskip('fiona')
pytest.importorskip('shapely')

SQUARE = [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]
TRIANGLE = [[20.0, 0.0], [30.0, 0.0], [25.0, 8.0]]

BUILDINGS_CSV = """easting,northing,id,floors
422664.22,870785.46,2,0
422672.48,870780.14,2,0
422668.17,870772.62,2,0
422660.35,870777.17,2,0
422664.22,870785.46,2,0
422661.30,871215.06,3,1
422667.50,871215.70,3,1
422668.30,871204.86,3,1
422662.21,871204.33,3,1
422661.30,871215.06,3,1
"""


class Test_gpkg(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def path(self, name):
        return os.path.join(self.tmp, name)

    def test_polygons_round_trip_with_attributes_and_crs(self):
        gpkg = self.path('regions.gpkg')
        anuga.polygons2gpkg([SQUARE, TRIANGLE], gpkg,
                            attributes=[{'name': 'square', 'k': 1, 'z': 0.5},
                                        {'name': 'triangle', 'k': 2, 'z': 1.5}],
                            crs=32756, layer='regions')
        assert anuga.gpkg_layers(gpkg) == ['regions']
        polys, attrs = anuga.gpkg2polygons(gpkg)
        assert len(polys) == 2
        assert np.allclose(polys[0], SQUARE)          # not closed, ANUGA style
        assert np.allclose(polys[1], TRIANGLE)
        assert attrs[0] == {'name': 'square', 'k': 1, 'z': 0.5}
        assert attrs[1]['name'] == 'triangle' and attrs[1]['k'] == 2
        with fiona.open(gpkg) as src:
            assert src.crs.to_epsg() == 32756
            assert src.schema['geometry'] == 'Polygon'
            assert src.schema['properties']['k'].startswith('int')
        closed, _ = anuga.gpkg2polygons(gpkg, closed=True)
        assert len(closed[0]) == 5 and np.allclose(closed[0][0], closed[0][-1])

    def test_a_closed_input_ring_is_not_doubled(self):
        gpkg = self.path('closed.gpkg')
        anuga.polygons2gpkg([SQUARE + [SQUARE[0]]], gpkg)
        polys, _ = anuga.gpkg2polygons(gpkg)
        assert np.allclose(polys[0], SQUARE)

    def test_multipolygon_gives_one_entry_per_part(self):
        from shapely.geometry import MultiPolygon, Polygon, mapping
        gpkg = self.path('multi.gpkg')
        schema = {'geometry': 'MultiPolygon', 'properties': {'name': 'str'}}
        with fiona.open(gpkg, 'w', driver='GPKG', layer='m', schema=schema) as dst:
            dst.write({'geometry': mapping(MultiPolygon([Polygon(SQUARE), Polygon(TRIANGLE)])),
                       'properties': {'name': 'both'}})
        polys, attrs = anuga.gpkg2polygons(gpkg)
        assert len(polys) == 2
        assert attrs == [{'name': 'both'}, {'name': 'both'}]

    def test_polylines_and_layer_selection(self):
        gpkg = self.path('model.gpkg')
        anuga.polygons2gpkg([SQUARE], gpkg, layer='regions')
        wall = [[0.0, 5.0], [10.0, 5.0], [20.0, 7.0]]
        anuga.polygons2gpkg([wall], gpkg, layer='riverwalls', geometry_type='LineString')
        assert sorted(anuga.gpkg_layers(gpkg)) == ['regions', 'riverwalls']
        lines, _ = anuga.gpkg2polygons(gpkg, layer='riverwalls')
        assert np.allclose(lines[0], wall)
        polys, _ = anuga.gpkg2polygons(gpkg, layer='regions')
        assert np.allclose(polys[0], SQUARE)
        # rewriting a layer replaces it and keeps the other
        anuga.polygons2gpkg([TRIANGLE], gpkg, layer='regions')
        polys, _ = anuga.gpkg2polygons(gpkg, layer='regions')
        assert len(polys) == 1 and np.allclose(polys[0], TRIANGLE)
        assert sorted(anuga.gpkg_layers(gpkg)) == ['regions', 'riverwalls']

    def test_polygon_csv_files_round_trip(self):
        """The Merewether convention: one x,y file per house, no header."""
        d = self.path('houses')
        os.makedirs(d)
        anuga.geometry.polygon.write_polygon(SQUARE, os.path.join(d, 'house000.csv'))
        anuga.geometry.polygon.write_polygon(TRIANGLE, os.path.join(d, 'house001.csv'))
        gpkg = self.path('houses.gpkg')
        files = anuga.polygon_csv_files2gpkg(d, gpkg, crs=32756)
        assert len(files) == 2
        polys, attrs = anuga.gpkg2polygons(gpkg)
        assert [a['name'] for a in attrs] == ['house000', 'house001']
        assert np.allclose(polys[0], SQUARE)
        out = self.path('houses_back')
        paths = anuga.gpkg2polygon_csv_files(gpkg, out)
        assert [os.path.basename(p) for p in paths] == ['house000.csv', 'house001.csv']
        assert np.allclose(anuga.read_polygon(paths[1]), TRIANGLE)

    def test_building_csv_round_trip(self):
        """The easting,northing,id,floors convention of load_csv_as_polygons."""
        csv_in = self.path('buildings.csv')
        with open(csv_in, 'w') as f:
            f.write(BUILDINGS_CSV)
        gpkg = self.path('buildings.gpkg')
        anuga.building_csv2gpkg(csv_in, gpkg, crs=32756)
        polys, attrs = anuga.gpkg2polygons(gpkg)
        assert [a['id'] for a in attrs] == [2, 3]
        assert [a['floors'] for a in attrs] == [0.0, 1.0]
        assert len(polys[0]) == 4                     # closing point dropped
        csv_out = self.path('buildings_back.csv')
        anuga.gpkg2building_csv(gpkg, csv_out)
        polygons, values = anuga.load_csv_as_polygons(csv_out, value_name='floors')
        assert sorted(polygons.keys()) == ['2', '3']
        assert float(values['3']) == 1.0
        heights_in, = [anuga.load_csv_as_building_polygons(csv_in)[1]]
        heights_out = anuga.load_csv_as_building_polygons(csv_out)[1]
        assert heights_out == heights_in
        assert np.allclose(polygons['2'][:4],
                           anuga.load_csv_as_polygons(csv_in, value_name='floors')[0]['2'][:4])

    def test_bad_input_is_rejected(self):
        gpkg = self.path('bad.gpkg')
        with self.assertRaises(ValueError):
            anuga.polygons2gpkg([SQUARE], gpkg, attributes=[{}, {}])
        with self.assertRaises(ValueError):
            anuga.polygons2gpkg([[[0.0, 0.0]]], gpkg)
        with self.assertRaises(ValueError):
            anuga.polygon_csv_files2gpkg(self.path('nothing'), gpkg)
