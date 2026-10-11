"""
Unit tests for anuga.scenario.parse_input_data_toml

Tests cover the three module-level helpers (_load_toml, _glob_files,
_parse_quantity_entries) and all private parsing methods of ProjectDataTOML
via its public constructor.
"""
import math
import os
import shutil
import tempfile
import textwrap
import warnings
import unittest
from unittest.mock import patch, MagicMock

try:
    from anuga.scenario.parse_input_data_toml import (
        _load_toml, _glob_files, _parse_quantity_entries, ProjectDataTOML)
    from anuga.scenario.parse_input_data import ProjectData
    HAS_MODULE = True
    SKIP_REASON = ''
except ImportError as _e:
    HAS_MODULE = False
    SKIP_REASON = str(_e)


# ---------------------------------------------------------------------------
# Helpers shared across test cases
# ---------------------------------------------------------------------------

_MINIMAL_PROJECT = """\
    [project]
    scenario = "test_scenario"
    output_base_directory = "OUTPUT/"
    yieldstep = 60.0
    finaltime = 3600.0
    projection_information = -55
    flow_algorithm = "DE0"

    [mesh]
    bounding_polygon = "extent.shp"
    default_res = 1000000.0
"""


def _write_toml(path, content):
    """Write *content* (a str) to *path* as a UTF-8 binary file (tomllib requirement)."""
    with open(path, 'wb') as fh:
        fh.write(textwrap.dedent(content).encode())


def _touch(path):
    """Create an empty file at *path*."""
    open(path, 'w').close()


def _tp(path):
    """Return a TOML-safe path string.

    TOML basic strings treat backslash as an escape character, so Windows
    paths must use forward slashes (which Windows also accepts).
    """
    return str(path).replace('\\', '/')


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestLoadToml(unittest.TestCase):
    """Tests for the _load_toml helper."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _toml_path(self, name='config.toml'):
        return os.path.join(self.tmpdir, name)

    def test_simple_key_value(self):
        p = self._toml_path()
        _write_toml(p, '[project]\nscenario = "test"\n')
        result = _load_toml(p)
        self.assertEqual(result['project']['scenario'], 'test')

    def test_integer_and_float(self):
        p = self._toml_path()
        _write_toml(p, '[s]\na = 42\nb = 3.14\n')
        r = _load_toml(p)
        self.assertEqual(r['s']['a'], 42)
        self.assertAlmostEqual(r['s']['b'], 3.14)

    def test_boolean_values(self):
        p = self._toml_path()
        _write_toml(p, '[s]\nflag_true = true\nflag_false = false\n')
        r = _load_toml(p)
        self.assertIs(r['s']['flag_true'], True)
        self.assertIs(r['s']['flag_false'], False)

    def test_array_of_tables(self):
        p = self._toml_path()
        _write_toml(p, '[[items]]\nname = "a"\n[[items]]\nname = "b"\n')
        r = _load_toml(p)
        self.assertEqual(len(r['items']), 2)
        self.assertEqual(r['items'][1]['name'], 'b')

    def test_missing_file_raises(self):
        with self.assertRaises((FileNotFoundError, OSError)):
            _load_toml(os.path.join(self.tmpdir, 'nonexistent.toml'))


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestGlobFiles(unittest.TestCase):
    """Tests for the _glob_files helper."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_empty_list_returns_empty(self):
        self.assertEqual(_glob_files([]), [])

    def test_single_exact_match(self):
        f = os.path.join(self.tmpdir, 'line.csv')
        _touch(f)
        result = _glob_files([f])
        self.assertEqual(result, [f])

    def test_glob_pattern_matches_multiple(self):
        for name in ('a.csv', 'b.csv', 'c.csv'):
            _touch(os.path.join(self.tmpdir, name))
        pattern = os.path.join(self.tmpdir, '*.csv')
        result = _glob_files([pattern])
        self.assertEqual(len(result), 3)
        self.assertTrue(all(r.endswith('.csv') for r in result))

    def test_multiple_patterns_concatenated(self):
        for name in ('x.csv', 'y.txt'):
            _touch(os.path.join(self.tmpdir, name))
        result = _glob_files([
            os.path.join(self.tmpdir, '*.csv'),
            os.path.join(self.tmpdir, '*.txt'),
        ])
        self.assertEqual(len(result), 2)

    def test_pattern_no_match_raises(self):
        with self.assertRaises(FileNotFoundError):
            _glob_files([os.path.join(self.tmpdir, '*.csv')])

    def test_results_are_sorted(self):
        for name in ('c.csv', 'a.csv', 'b.csv'):
            _touch(os.path.join(self.tmpdir, name))
        result = _glob_files([os.path.join(self.tmpdir, '*.csv')])
        self.assertEqual(result, sorted(result))


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestParseQuantityEntries(unittest.TestCase):
    """Tests for the _parse_quantity_entries helper."""

    def _parse(self, entries):
        info = []
        return _parse_quantity_entries(entries, info), info

    def test_empty_entries(self):
        (data, clip), _ = self._parse([])
        self.assertEqual(data, [])
        self.assertEqual(clip, [])

    def test_single_all_polygon_numeric_value(self):
        (data, clip), _ = self._parse([{'polygon': 'All', 'value': 5.0}])
        self.assertEqual(data, [['All', 5.0]])
        self.assertEqual(clip[0][0], float('-inf'))
        self.assertEqual(clip[0][1], float('inf'))

    def test_string_value_preserved(self):
        (data, clip), _ = self._parse([{'polygon': 'Extent', 'value': 'dem.asc'}])
        self.assertEqual(data[0][1], 'dem.asc')

    def test_explicit_clip_range(self):
        entry = {'polygon': 'All', 'value': 0.0, 'clip_min': -10.0, 'clip_max': 10.0}
        (data, clip), _ = self._parse([entry])
        self.assertAlmostEqual(clip[0][0], -10.0)
        self.assertAlmostEqual(clip[0][1], 10.0)

    def test_default_clip_is_inf(self):
        (data, clip), _ = self._parse([{'polygon': 'None', 'value': 1.0}])
        self.assertTrue(math.isinf(clip[0][0]) and clip[0][0] < 0)
        self.assertTrue(math.isinf(clip[0][1]) and clip[0][1] > 0)

    def test_multiple_entries_accumulate(self):
        entries = [
            {'polygon': 'All', 'value': 0.0},
            {'polygon': 'sub.csv', 'value': 5.0},
        ]
        (data, clip), _ = self._parse(entries)
        self.assertEqual(len(data), 2)
        self.assertEqual(len(clip), 2)
        self.assertEqual(data[1], ['sub.csv', 5.0])

    def test_special_polygons_not_treated_as_globs(self):
        """'All', 'None', 'Extent' must be passed through unchanged."""
        for special in ('All', 'None', 'Extent'):
            (data, _), _ = self._parse([{'polygon': special, 'value': 0.0}])
            self.assertEqual(data[0][0], special)

    def test_wildcard_more_than_two_files_raises(self):
        """A glob matching >2 files should raise ValueError."""
        with tempfile.TemporaryDirectory() as d:
            for name in ('a.csv', 'b.csv', 'c.csv'):
                _touch(os.path.join(d, name))
            pattern = os.path.join(d, '*.csv')
            with self.assertRaises(ValueError):
                self._parse([{'polygon': pattern, 'value': 1.0}])

    def test_wildcard_two_files_calls_polygon_join(self):
        """A glob matching exactly 2 files should call polygon_from_matching_breaklines."""
        with tempfile.TemporaryDirectory() as d:
            for name in ('left.csv', 'right.csv'):
                _touch(os.path.join(d, name))
            pattern = os.path.join(d, '*.csv')
            fake_polygon = [(0, 0), (1, 0), (1, 1), (0, 1)]
            with patch('anuga.scenario.parse_input_data_toml.su') as mock_su:
                mock_su.read_polygon.return_value = [(0, 0), (1, 1)]
                mock_su.polygon_from_matching_breaklines.return_value = fake_polygon
                (data, _), info = self._parse([{'polygon': pattern, 'value': 2.0}])
            # The polygon in data should be the joined result
            self.assertEqual(data[0][0], fake_polygon)
            # An info message should have been appended
            self.assertTrue(any('Combining' in m for m in info))

    def test_print_info_populated_for_wildcard(self):
        with tempfile.TemporaryDirectory() as d:
            for name in ('l.csv', 'r.csv'):
                _touch(os.path.join(d, name))
            pattern = os.path.join(d, '*.csv')
            with patch('anuga.scenario.parse_input_data_toml.su') as mock_su:
                mock_su.read_polygon.return_value = []
                mock_su.polygon_from_matching_breaklines.return_value = []
                (_, _), info = self._parse([{'polygon': pattern, 'value': 0.0}])
            self.assertEqual(len(info), 1)


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestProjectSection(unittest.TestCase):
    """Tests for ProjectDataTOML._parse_project via the constructor."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.toml_path = os.path.join(self.tmpdir, 'config.toml')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make(self, project_extra='', mesh_extra=''):
        content = textwrap.dedent(f"""\
            [project]
            scenario = "test_scenario"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            {project_extra}

            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 1000000.0
            {mesh_extra}
        """)
        _write_toml(self.toml_path, content)
        return ProjectDataTOML(self.toml_path)

    def test_required_fields_parsed(self):
        p = self._make()
        self.assertEqual(p.scenario, 'test_scenario')
        self.assertEqual(p.output_basedir, 'OUTPUT/')
        self.assertAlmostEqual(p.yieldstep, 60.0)
        self.assertAlmostEqual(p.finaltime, 3600.0)
        self.assertEqual(p.flow_algorithm, 'DE0')

    def _make_with_algorithm(self, alg):
        content = textwrap.dedent(f"""\
            [project]
            scenario = "test_scenario"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "{alg}"
            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 1000000.0
        """)
        _write_toml(self.toml_path, content)
        return ProjectDataTOML(self.toml_path)

    def test_all_supported_flow_algorithms_accepted(self):
        # Must stay in sync with Domain.set_flow_algorithm's accepted list.
        for alg in ('DE0', 'DE1', 'DE2', 'DE0_7', 'DE1_7', 'DE_ader2'):
            p = self._make_with_algorithm(alg)
            self.assertEqual(p.flow_algorithm, alg)

    def test_invalid_flow_algorithm_raises(self):
        with self.assertRaises(ValueError):
            self._make_with_algorithm('DE_nope')

    def test_projection_information_as_int(self):
        p = self._make()
        self.assertEqual(p.projection_information, -55)
        self.assertIsInstance(p.projection_information, int)

    def test_projection_information_as_string(self):
        # Write a TOML where projection_information is a proj4 string (not int)
        content = textwrap.dedent("""\
            [project]
            scenario = "test_scenario"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = "+proj=utm +zone=55 +south"
            flow_algorithm = "DE0"

            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 1000000.0
        """)
        _write_toml(self.toml_path, content)
        p = ProjectDataTOML(self.toml_path)
        self.assertEqual(p.projection_information, '+proj=utm +zone=55 +south')

    def test_optional_defaults(self):
        p = self._make()
        self.assertAlmostEqual(p.output_tif_cellsize, 50.0)
        self.assertIsNone(p.output_tif_bounding_polygon)
        self.assertEqual(p.max_quantity_update_frequency, 1)
        self.assertAlmostEqual(p.max_quantity_collection_start_time, 0.0)
        self.assertFalse(p.store_vertices_uniquely)
        self.assertFalse(p.store_elevation_every_timestep)
        self.assertFalse(p.report_mass_conservation_statistics)
        self.assertFalse(p.report_peak_velocity_statistics)
        self.assertFalse(p.report_smallest_edge_timestep_statistics)
        self.assertFalse(p.report_operator_statistics)
        self.assertIsNone(p.omp_num_threads)
        self.assertEqual(p.compute_mode, 'legacy')
        self.assertEqual(p.multiprocessor_mode, 1)
        self.assertIsNone(p.outputstep)

    def test_optional_overrides(self):
        extra = textwrap.dedent("""\
            output_tif_cellsize = 25.0
            output_tif_bounding_polygon = "clip.shp"
            max_quantity_update_frequency = 5
            store_vertices_uniquely = true
            omp_num_threads = 4
            multiprocessor_mode = 2
            outputstep = 120.0
        """)
        p = self._make(project_extra=extra)
        self.assertAlmostEqual(p.output_tif_cellsize, 25.0)
        self.assertEqual(p.output_tif_bounding_polygon, 'clip.shp')
        self.assertEqual(p.max_quantity_update_frequency, 5)
        self.assertTrue(p.store_vertices_uniquely)
        self.assertEqual(p.omp_num_threads, 4)
        self.assertEqual(p.multiprocessor_mode, 2)
        self.assertAlmostEqual(p.outputstep, 120.0)

    def test_compute_mode_unified(self):
        p = self._make(project_extra='compute_mode = "unified"')
        self.assertEqual(p.compute_mode, 'unified')
        self.assertEqual(p.multiprocessor_mode, 2)   # kept in sync

    def test_compute_mode_legacy(self):
        p = self._make(project_extra='compute_mode = "legacy"')
        self.assertEqual(p.compute_mode, 'legacy')
        self.assertEqual(p.multiprocessor_mode, 1)

    def test_compute_mode_invalid_raises(self):
        with self.assertRaises(ValueError):
            self._make(project_extra='compute_mode = "turbo"')

    def test_multiprocessor_mode_back_compat_maps_to_compute_mode(self):
        """The legacy integer key still works and maps to compute_mode."""
        p = self._make(project_extra='multiprocessor_mode = 2')
        self.assertEqual(p.compute_mode, 'unified')
        self.assertEqual(p.multiprocessor_mode, 2)

    def test_compute_mode_takes_precedence_over_multiprocessor_mode(self):
        p = self._make(project_extra=textwrap.dedent("""\
            compute_mode = "legacy"
            multiprocessor_mode = 2
        """))
        self.assertEqual(p.compute_mode, 'legacy')
        self.assertEqual(p.multiprocessor_mode, 1)

    def test_collection_starttime_exceeds_finaltime_raises(self):
        extra = 'max_quantity_collection_starttime = 7200.0'
        with self.assertRaises(ValueError):
            self._make(project_extra=extra)

    def test_starttime_defaults_to_zero(self):
        self.assertAlmostEqual(self._make().starttime, 0.0)

    def test_negative_starttime_expression(self):
        p = self._make(project_extra='starttime = "-3*600"')
        self.assertAlmostEqual(p.starttime, -1800.0)

    def test_starttime_not_before_finaltime_raises(self):
        with self.assertRaises(ValueError):
            self._make(project_extra='starttime = 3600.0')

    def test_config_filename_attribute(self):
        p = self._make()
        self.assertEqual(p.config_filename, self.toml_path)

    def test_print_info_is_list(self):
        p = self._make()
        self.assertIsInstance(p.print_info, list)
        self.assertTrue(len(p.print_info) > 0)


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestMeshSection(unittest.TestCase):
    """Tests for ProjectDataTOML._parse_mesh via the constructor."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.toml_path = os.path.join(self.tmpdir, 'config.toml')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make(self, mesh_extra=''):
        content = textwrap.dedent(f"""\
            [project]
            scenario = "test"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"

            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 500000.0
            {mesh_extra}
        """)
        _write_toml(self.toml_path, content)
        return ProjectDataTOML(self.toml_path)

    def test_defaults(self):
        p = self._make()
        self.assertFalse(p.use_existing_mesh_pickle)
        self.assertEqual(p.bounding_polygon_and_tags_file, 'extent.shp')
        self.assertAlmostEqual(p.default_res, 500000.0)
        self.assertEqual(p.interior_regions_data, [])
        self.assertEqual(p.interior_holes_data, [])
        self.assertIsNone(p.bounding_polygon_explicit_tags)
        self.assertEqual(p.breakline_files, [])
        self.assertEqual(p.riverwall_csv_files, [])
        self.assertEqual(p.break_line_intersect_point_movement_threshold, 'ignore')
        self.assertIsNone(p.pt_areas)

    def test_interior_regions_parsed(self):
        extra = textwrap.dedent("""\
            [[mesh.interior_regions]]
            polygon = "region_a.shp"
            resolution = 10000.0
            [[mesh.interior_regions]]
            polygon = "region_b.shp"
            resolution = 5000.0
        """)
        p = self._make(mesh_extra=extra)
        self.assertEqual(len(p.interior_regions_data), 2)
        self.assertEqual(p.interior_regions_data[0], ['region_a.shp', 10000.0])
        self.assertEqual(p.interior_regions_data[1], ['region_b.shp', 5000.0])

    def test_interior_holes_parsed(self):
        extra = textwrap.dedent("""\
            [[mesh.interior_holes]]
            polygon = "building_a.csv"
            tag = "building"
            [[mesh.interior_holes]]
            polygon = "building_b.csv"
        """)
        p = self._make(mesh_extra=extra)
        self.assertEqual(len(p.interior_holes_data), 2)
        # tag is optional; absent means None, which pmesh renders as its own
        # 'interior' default rather than a named tag.
        self.assertEqual(p.interior_holes_data[0], ['building_a.csv', 'building'])
        self.assertEqual(p.interior_holes_data[1], ['building_b.csv', None])

    def test_interior_holes_coexist_with_interior_regions(self):
        # Holes and regions are independent: a hole removes triangles, a region
        # only resizes them, so a mesh may legitimately carry both.
        extra = textwrap.dedent("""\
            [[mesh.interior_regions]]
            polygon = "region_a.shp"
            resolution = 10000.0
            [[mesh.interior_holes]]
            polygon = "building_a.csv"
        """)
        p = self._make(mesh_extra=extra)
        self.assertEqual(len(p.interior_regions_data), 1)
        self.assertEqual(len(p.interior_holes_data), 1)

    def test_boundary_tags_parsed(self):
        extra = textwrap.dedent("""\
            [[mesh.boundary_tags]]
            tag = "ocean"
            edges = [0, 1, 2]
            [[mesh.boundary_tags]]
            tag = "land"
            edges = [3]
        """)
        p = self._make(mesh_extra=extra)
        self.assertIsNotNone(p.bounding_polygon_explicit_tags)
        self.assertEqual(len(p.bounding_polygon_explicit_tags), 2)
        self.assertEqual(p.bounding_polygon_explicit_tags[0]['tag'], 'ocean')
        self.assertEqual(p.bounding_polygon_explicit_tags[0]['edges'], [0, 1, 2])

    def test_numeric_breakline_threshold(self):
        p = self._make('breakline_intersection_threshold = 0.5')
        self.assertAlmostEqual(p.break_line_intersect_point_movement_threshold, 0.5)

    def test_interior_regions_and_breaklines_coexist(self):
        """interior_regions may be combined with breaklines/riverwalls.

        setup_mesh passes interior_regions AND breaklines (riverwalls) to
        create_pmesh_from_regions together, so this is a valid combination
        (e.g. a catchment-refined model that also has riverwalls).
        """
        bl = os.path.join(self.tmpdir, 'line.csv')
        _touch(bl)
        extra = textwrap.dedent(f"""\
            breakline_files = ["{_tp(bl)}"]
            [[mesh.interior_regions]]
            polygon = "reg.shp"
            resolution = 1000.0
        """)
        p = self._make(mesh_extra=extra)   # must not raise
        self.assertEqual(len(p.interior_regions_data), 1)
        self.assertEqual(len(p.breakline_files), 1)

    def test_interior_regions_and_region_areas_file_conflict_raises(self):
        """interior_regions and region_areas_file are alternative resolution
        methods and cannot be combined."""
        ra = os.path.join(self.tmpdir, 'areas.csv')
        _touch(ra)
        extra = textwrap.dedent(f"""\
            region_areas_file = "{_tp(ra)}"
            [[mesh.interior_regions]]
            polygon = "reg.shp"
            resolution = 1000.0
        """)
        with self.assertRaises(ValueError):
            self._make(mesh_extra=extra)

    def test_region_areas_type_length(self):
        areas_file = os.path.join(self.tmpdir, 'areas.csv')
        _touch(areas_file)
        extra = textwrap.dedent(f"""\
            region_areas_file = "{_tp(areas_file)}"
            region_areas_type = "length"
        """)
        p = self._make(mesh_extra=extra)
        self.assertTrue(p.region_resolutions_from_length)

    def test_region_areas_type_area(self):
        areas_file = os.path.join(self.tmpdir, 'areas.csv')
        _touch(areas_file)
        extra = textwrap.dedent(f"""\
            region_areas_file = "{_tp(areas_file)}"
            region_areas_type = "area"
        """)
        p = self._make(mesh_extra=extra)
        self.assertFalse(p.region_resolutions_from_length)

    def test_region_areas_type_invalid_raises(self):
        areas_file = os.path.join(self.tmpdir, 'areas.csv')
        _touch(areas_file)
        extra = textwrap.dedent(f"""\
            region_areas_file = "{_tp(areas_file)}"
            region_areas_type = "banana"
        """)
        with self.assertRaises(ValueError):
            self._make(mesh_extra=extra)

    def test_breakline_files_expanded(self):
        for name in ('bl_a.csv', 'bl_b.csv'):
            _touch(os.path.join(self.tmpdir, name))
        pattern = os.path.join(self.tmpdir, 'bl_*.csv')
        extra = f'breakline_files = ["{_tp(pattern)}"]'
        p = self._make(mesh_extra=extra)
        self.assertEqual(len(p.breakline_files), 2)


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestInitialConditions(unittest.TestCase):
    """Tests for ProjectDataTOML._parse_initial_conditions."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.toml_path = os.path.join(self.tmpdir, 'config.toml')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make(self, ic_section='', ic_add_section=''):
        content = textwrap.dedent(f"""\
            [project]
            scenario = "test"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"

            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 1000000.0

            {ic_section}
            {ic_add_section}
        """)
        _write_toml(self.toml_path, content)
        return ProjectDataTOML(self.toml_path)

    def test_all_quantities_initialised_empty(self):
        p = self._make()
        for q in ('elevation', 'friction', 'stage', 'xmomentum', 'ymomentum'):
            self.assertEqual(getattr(p, f'{q}_data'), [])
            self.assertEqual(getattr(p, f'{q}_clip_range'), [])
            self.assertIsNone(getattr(p, f'{q}_mean'))
            self.assertEqual(getattr(p, f'{q}_additions'), [])

    def test_elevation_entry(self):
        ic = textwrap.dedent("""\
            [initial_conditions]
            [[initial_conditions.elevation]]
            polygon = "All"
            value = 0.0
            clip_min = -100.0
            clip_max = 200.0
        """)
        p = self._make(ic_section=ic)
        self.assertEqual(p.elevation_data, [['All', 0.0]])
        self.assertAlmostEqual(p.elevation_clip_range[0][0], -100.0)
        self.assertAlmostEqual(p.elevation_clip_range[0][1], 200.0)

    def test_multiple_elevation_entries_ordered(self):
        ic = textwrap.dedent("""\
            [initial_conditions]
            [[initial_conditions.elevation]]
            polygon = "All"
            value = 0.0
            [[initial_conditions.elevation]]
            polygon = "sub.shp"
            value = "local.asc"
        """)
        p = self._make(ic_section=ic)
        self.assertEqual(len(p.elevation_data), 2)
        self.assertEqual(p.elevation_data[1][0], 'sub.shp')

    def test_spatial_average(self):
        ic = textwrap.dedent("""\
            [initial_conditions]
            elevation_spatial_average = 25.0
        """)
        p = self._make(ic_section=ic)
        self.assertAlmostEqual(p.elevation_mean, 25.0)

    def test_additions_parsed(self):
        add = textwrap.dedent("""\
            [initial_condition_additions]
            [[initial_condition_additions.elevation]]
            polygon = "All"
            value = 1.5
        """)
        p = self._make(ic_add_section=add)
        self.assertEqual(p.elevation_additions, [['All', 1.5]])

    def test_friction_entry(self):
        ic = textwrap.dedent("""\
            [initial_conditions]
            [[initial_conditions.friction]]
            polygon = "All"
            value = 0.03
        """)
        p = self._make(ic_section=ic)
        self.assertEqual(p.friction_data, [['All', 0.03]])


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestBoundaryConditions(unittest.TestCase):
    """Tests for ProjectDataTOML._parse_boundary_conditions."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.toml_path = os.path.join(self.tmpdir, 'config.toml')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make(self, bc_section=''):
        content = textwrap.dedent(f"""\
            [project]
            scenario = "test"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"

            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 1000000.0

            {bc_section}
        """)
        _write_toml(self.toml_path, content)
        return ProjectDataTOML(self.toml_path)

    def test_empty_boundary_conditions(self):
        p = self._make()
        self.assertEqual(p.boundary_data, [])

    def test_default_attribute_name(self):
        p = self._make()
        self.assertEqual(p.boundary_tags_attribute_name, 'bnd_tag')

    def test_custom_attribute_name(self):
        bc = '[boundary_conditions]\nboundary_tags_attribute_name = "TAG"\n'
        p = self._make(bc_section=bc)
        self.assertEqual(p.boundary_tags_attribute_name, 'TAG')

    def test_reflective_boundary(self):
        bc = textwrap.dedent("""\
            [boundary_conditions]
            [[boundary_conditions.boundaries]]
            tag = "land"
            type = "Reflective"
        """)
        p = self._make(bc_section=bc)
        self.assertEqual(len(p.boundary_data), 1)
        self.assertEqual(p.boundary_data[0], ['land', 'Reflective'])

    def test_stage_boundary_with_file(self):
        tfile = os.path.join(self.tmpdir, 'tide.csv')
        _touch(tfile)
        bc = textwrap.dedent(f"""\
            [boundary_conditions]
            [[boundary_conditions.boundaries]]
            tag = "ocean"
            type = "Stage"
            file = "{_tp(tfile)}"
            start_time = 100.0
        """)
        p = self._make(bc_section=bc)
        row = p.boundary_data[0]
        self.assertEqual(row[0], 'ocean')
        self.assertEqual(row[1], 'Stage')
        self.assertEqual(row[2], tfile)
        self.assertAlmostEqual(row[3], 100.0)

    def test_stage_boundary_full_class_name_alias(self):
        """The explicit ANUGA class name is accepted as an alias for "Stage"."""
        tfile = os.path.join(self.tmpdir, 'tide.csv')
        _touch(tfile)
        for type_name in ('Transmissive_n_momentum_zero_t_momentum_set_stage',
                          'Transmissive_n_momentum_zero_t_momentum_set_stage_boundary'):
            bc = textwrap.dedent(f"""\
                [boundary_conditions]
                [[boundary_conditions.boundaries]]
                tag = "east"
                type = "{type_name}"
                file = "{_tp(tfile)}"
            """)
            p = self._make(bc_section=bc)
            self.assertEqual(p.boundary_data[0][1], 'Stage')  # normalised

    def test_flather_stage_boundary(self):
        tfile = os.path.join(self.tmpdir, 'wave.csv')
        _touch(tfile)
        bc = textwrap.dedent(f"""\
            [boundary_conditions]
            [[boundary_conditions.boundaries]]
            tag = "east"
            type = "Flather_Stage"
            file = "{_tp(tfile)}"
        """)
        p = self._make(bc_section=bc)
        self.assertEqual(p.boundary_data[0][1], 'Flather_Stage')

    def test_stage_boundary_default_start_time_zero(self):
        tfile = os.path.join(self.tmpdir, 'tide.csv')
        _touch(tfile)
        bc = textwrap.dedent(f"""\
            [boundary_conditions]
            [[boundary_conditions.boundaries]]
            tag = "ocean"
            type = "Stage"
            file = "{_tp(tfile)}"
        """)
        p = self._make(bc_section=bc)
        self.assertAlmostEqual(p.boundary_data[0][3], 0.0)

    def test_stage_boundary_missing_file_raises(self):
        bc = textwrap.dedent("""\
            [boundary_conditions]
            [[boundary_conditions.boundaries]]
            tag = "ocean"
            type = "Stage"
            file = "nonexistent_tide.csv"
        """)
        with self.assertRaises(FileNotFoundError):
            self._make(bc_section=bc)

    def test_unknown_boundary_type_raises(self):
        bc = textwrap.dedent("""\
            [boundary_conditions]
            [[boundary_conditions.boundaries]]
            tag = "weird"
            type = "Transmissive"
        """)
        with self.assertRaises(ValueError):
            self._make(bc_section=bc)

    def test_multiple_boundaries(self):
        tfile = os.path.join(self.tmpdir, 'tide.csv')
        _touch(tfile)
        bc = textwrap.dedent(f"""\
            [boundary_conditions]
            [[boundary_conditions.boundaries]]
            tag = "ocean"
            type = "Stage"
            file = "{_tp(tfile)}"
            [[boundary_conditions.boundaries]]
            tag = "land"
            type = "Reflective"
        """)
        p = self._make(bc_section=bc)
        self.assertEqual(len(p.boundary_data), 2)


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestInlets(unittest.TestCase):
    """Tests for ProjectDataTOML._parse_inlets."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.toml_path = os.path.join(self.tmpdir, 'config.toml')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make(self, inlet_section=''):
        content = textwrap.dedent(f"""\
            [project]
            scenario = "test"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"

            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 1000000.0

            {inlet_section}
        """)
        _write_toml(self.toml_path, content)
        return ProjectDataTOML(self.toml_path)

    def test_no_inlets(self):
        p = self._make()
        self.assertEqual(p.inlet_data, [])

    def test_single_inlet(self):
        lfile = os.path.join(self.tmpdir, 'inlet_line.csv')
        tfile = os.path.join(self.tmpdir, 'inlet_ts.csv')
        _touch(lfile)
        _touch(tfile)
        section = textwrap.dedent(f"""\
            [[inlets]]
            name = "inlet_1"
            line_file = "{_tp(lfile)}"
            timeseries_file = "{_tp(tfile)}"
            start_time = 30.0
        """)
        p = self._make(inlet_section=section)
        self.assertEqual(len(p.inlet_data), 1)
        row = p.inlet_data[0]
        self.assertEqual(row[0], 'inlet_1')
        self.assertEqual(row[1], lfile)
        self.assertEqual(row[2], tfile)
        self.assertAlmostEqual(row[3], 30.0)

    def test_default_start_time_zero(self):
        lfile = os.path.join(self.tmpdir, 'l.csv')
        tfile = os.path.join(self.tmpdir, 't.csv')
        _touch(lfile)
        _touch(tfile)
        section = textwrap.dedent(f"""\
            [[inlets]]
            name = "q1"
            line_file = "{_tp(lfile)}"
            timeseries_file = "{_tp(tfile)}"
        """)
        p = self._make(inlet_section=section)
        self.assertAlmostEqual(p.inlet_data[0][3], 0.0)

    def test_missing_line_file_raises(self):
        tfile = os.path.join(self.tmpdir, 't.csv')
        _touch(tfile)
        section = textwrap.dedent(f"""\
            [[inlets]]
            name = "q1"
            line_file = "nonexistent.csv"
            timeseries_file = "{_tp(tfile)}"
        """)
        with self.assertRaises(FileNotFoundError):
            self._make(inlet_section=section)

    def test_missing_timeseries_file_raises(self):
        lfile = os.path.join(self.tmpdir, 'l.csv')
        _touch(lfile)
        section = textwrap.dedent(f"""\
            [[inlets]]
            name = "q1"
            line_file = "{_tp(lfile)}"
            timeseries_file = "nonexistent.csv"
        """)
        with self.assertRaises(FileNotFoundError):
            self._make(inlet_section=section)


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestRainfall(unittest.TestCase):
    """Tests for ProjectDataTOML._parse_rainfall."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.toml_path = os.path.join(self.tmpdir, 'config.toml')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make(self, rain_section=''):
        content = textwrap.dedent(f"""\
            [project]
            scenario = "test"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"

            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 1000000.0

            {rain_section}
        """)
        _write_toml(self.toml_path, content)
        return ProjectDataTOML(self.toml_path)

    def test_no_rainfall(self):
        p = self._make()
        self.assertEqual(p.rain_data, [])

    def test_rainfall_polygon_all(self):
        tfile = os.path.join(self.tmpdir, 'rain.csv')
        _touch(tfile)
        section = textwrap.dedent(f"""\
            [[rainfall]]
            timeseries_file = "{_tp(tfile)}"
            polygon = "All"
            multiplier = 2.0
        """)
        p = self._make(rain_section=section)
        row = p.rain_data[0]
        self.assertEqual(row[0], tfile)
        self.assertAlmostEqual(row[1], 0.0)       # start_time default
        self.assertEqual(row[2], 'linear')         # interpolation_type default
        self.assertEqual(row[3], 'All')
        self.assertAlmostEqual(row[4], 2.0)

    def test_rainfall_with_polygon_file(self):
        tfile = os.path.join(self.tmpdir, 'rain.csv')
        pfile = os.path.join(self.tmpdir, 'rain_area.shp')
        _touch(tfile)
        _touch(pfile)
        section = textwrap.dedent(f"""\
            [[rainfall]]
            timeseries_file = "{_tp(tfile)}"
            polygon = "{_tp(pfile)}"
        """)
        p = self._make(rain_section=section)
        self.assertEqual(p.rain_data[0][3], pfile)

    def test_missing_timeseries_raises(self):
        section = textwrap.dedent("""\
            [[rainfall]]
            timeseries_file = "nonexistent.csv"
        """)
        with self.assertRaises(FileNotFoundError):
            self._make(rain_section=section)

    def test_missing_polygon_file_raises(self):
        tfile = os.path.join(self.tmpdir, 'rain.csv')
        _touch(tfile)
        section = textwrap.dedent(f"""\
            [[rainfall]]
            timeseries_file = "{_tp(tfile)}"
            polygon = "missing_area.shp"
        """)
        with self.assertRaises(FileNotFoundError):
            self._make(rain_section=section)

    def test_default_values(self):
        tfile = os.path.join(self.tmpdir, 'rain.csv')
        _touch(tfile)
        section = textwrap.dedent(f"""\
            [[rainfall]]
            timeseries_file = "{_tp(tfile)}"
        """)
        p = self._make(rain_section=section)
        row = p.rain_data[0]
        self.assertAlmostEqual(row[1], 0.0)
        self.assertEqual(row[2], 'linear')
        self.assertEqual(row[3], 'All')
        self.assertAlmostEqual(row[4], 1.0)


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestBridges(unittest.TestCase):
    """Tests for ProjectDataTOML._parse_bridges."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.toml_path = os.path.join(self.tmpdir, 'config.toml')
        # Create stub files for a valid bridge
        for name in ('deck.csv', 'exchange0.csv', 'exchange1.csv', 'curve.csv'):
            _touch(os.path.join(self.tmpdir, name))
        self.deck  = os.path.join(self.tmpdir, 'deck.csv')
        self.ex0   = os.path.join(self.tmpdir, 'exchange0.csv')
        self.ex1   = os.path.join(self.tmpdir, 'exchange1.csv')
        self.curve = os.path.join(self.tmpdir, 'curve.csv')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make(self, bridge_section=''):
        content = textwrap.dedent(f"""\
            [project]
            scenario = "test"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"

            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 1000000.0

            {bridge_section}
        """)
        _write_toml(self.toml_path, content)
        return ProjectDataTOML(self.toml_path)

    def _bridge_toml(self, label='bridge_1', extra=''):
        return textwrap.dedent(f"""\
            [[bridges]]
            label = "{label}"
            deck_file = "{_tp(self.deck)}"
            deck_elevation = 5.0
            exchange_line_0 = "{_tp(self.ex0)}"
            exchange_line_1 = "{_tp(self.ex1)}"
            enquiry_gap = 1.0
            internal_boundary_curve_file = "{_tp(self.curve)}"
            {extra}
        """)

    def test_no_bridges(self):
        p = self._make()
        self.assertEqual(p.bridge_data, [])

    def test_disabled_bridge_skipped(self):
        p = self._make(self._bridge_toml(extra='enabled = false'))
        self.assertEqual(p.bridge_data, [])
        # elevation_data should not have been modified
        self.assertEqual(p.elevation_data, [])

    def test_valid_bridge_row(self):
        p = self._make(self._bridge_toml())
        self.assertEqual(len(p.bridge_data), 1)
        row = p.bridge_data[0]
        self.assertEqual(row[0], 'bridge_1')
        self.assertEqual(row[1], self.deck)
        self.assertAlmostEqual(row[2], 5.0)
        self.assertEqual(row[3], self.ex0)
        self.assertEqual(row[4], self.ex1)
        self.assertAlmostEqual(row[5], 1.0)
        self.assertEqual(row[6], self.curve)
        self.assertAlmostEqual(row[7], 0.0)   # vertical_datum_offset default
        self.assertAlmostEqual(row[8], 0.0)   # smoothing_timescale default

    def test_bridge_prepended_to_breakline_files(self):
        p = self._make(self._bridge_toml())
        # deck file should be first in breakline_files
        self.assertEqual(p.breakline_files[0], self.deck)

    def test_bridge_prepended_to_elevation_data(self):
        p = self._make(self._bridge_toml())
        # deck should be first elevation entry with deck_elevation as value
        self.assertEqual(p.elevation_data[0], [self.deck, 5.0])
        self.assertEqual(p.elevation_clip_range[0], [float('-inf'), float('inf')])

    def test_bridge_optional_fields(self):
        extra = 'vertical_datum_offset = 0.3\nsmoothing_timescale = 60.0'
        p = self._make(self._bridge_toml(extra=extra))
        row = p.bridge_data[0]
        self.assertAlmostEqual(row[7], 0.3)
        self.assertAlmostEqual(row[8], 60.0)

    def test_missing_deck_file_raises(self):
        section = textwrap.dedent(f"""\
            [[bridges]]
            label = "b"
            deck_file = "nonexistent.csv"
            deck_elevation = 5.0
            exchange_line_0 = "{_tp(self.ex0)}"
            exchange_line_1 = "{_tp(self.ex1)}"
            enquiry_gap = 1.0
            internal_boundary_curve_file = "{_tp(self.curve)}"
        """)
        with self.assertRaises(FileNotFoundError):
            self._make(section)

    def test_two_bridges_both_prepended(self):
        deck2 = os.path.join(self.tmpdir, 'deck2.csv')
        _touch(deck2)
        section = self._bridge_toml('bridge_1') + self._bridge_toml('bridge_2')
        # Both bridges use same deck; replace deck2 via extra for bridge_2
        deck2_section = textwrap.dedent(f"""\
            [[bridges]]
            label = "bridge_2"
            deck_file = "{_tp(deck2)}"
            deck_elevation = 3.0
            exchange_line_0 = "{_tp(self.ex0)}"
            exchange_line_1 = "{_tp(self.ex1)}"
            enquiry_gap = 2.0
            internal_boundary_curve_file = "{_tp(self.curve)}"
        """)
        p = self._make(self._bridge_toml('bridge_1') + deck2_section)
        self.assertEqual(len(p.bridge_data), 2)
        # elevation_data: bridge_2 prepended last so is at index 0,
        # bridge_1 at index 1
        self.assertEqual(p.elevation_data[0][0], deck2)
        self.assertEqual(p.elevation_data[1][0], self.deck)


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestPumpingStations(unittest.TestCase):
    """Tests for ProjectDataTOML._parse_pumping_stations."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.toml_path = os.path.join(self.tmpdir, 'config.toml')
        for name in ('basin.shp', 'ex0.csv', 'ex1.csv'):
            _touch(os.path.join(self.tmpdir, name))
        self.basin = os.path.join(self.tmpdir, 'basin.shp')
        self.ex0   = os.path.join(self.tmpdir, 'ex0.csv')
        self.ex1   = os.path.join(self.tmpdir, 'ex1.csv')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make(self, ps_section=''):
        content = textwrap.dedent(f"""\
            [project]
            scenario = "test"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"

            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 1000000.0

            {ps_section}
        """)
        _write_toml(self.toml_path, content)
        return ProjectDataTOML(self.toml_path)

    def _ps_toml(self, label='ps_1', extra=''):
        return textwrap.dedent(f"""\
            [[pumping_stations]]
            label = "{label}"
            pump_capacity = 5.0
            pump_rate_of_increase = 1.0
            pump_rate_of_decrease = 1.0
            hw_to_start_pumping = 2.0
            hw_to_stop_pumping = 1.0
            basin_polygon_file = "{_tp(self.basin)}"
            basin_elevation = -1.5
            exchange_line_0 = "{_tp(self.ex0)}"
            exchange_line_1 = "{_tp(self.ex1)}"
            {extra}
        """)

    def test_no_pumping_stations(self):
        p = self._make()
        self.assertEqual(p.pumping_station_data, [])

    def test_disabled_station_skipped(self):
        p = self._make(self._ps_toml(extra='enabled = false'))
        self.assertEqual(p.pumping_station_data, [])

    def test_valid_station_row(self):
        p = self._make(self._ps_toml())
        self.assertEqual(len(p.pumping_station_data), 1)
        row = p.pumping_station_data[0]
        self.assertEqual(row[0], 'ps_1')
        self.assertAlmostEqual(row[1], 5.0)   # pump_capacity
        self.assertAlmostEqual(row[4], 2.0)   # hw_to_start_pumping
        self.assertAlmostEqual(row[5], 1.0)   # hw_to_stop_pumping
        self.assertEqual(row[6], self.basin)
        self.assertAlmostEqual(row[7], -1.5)  # basin_elevation
        self.assertEqual(row[8], self.ex0)
        self.assertEqual(row[9], self.ex1)
        self.assertAlmostEqual(row[10], 0.0)  # smoothing_timescale default

    def test_station_prepended_to_elevation_data(self):
        p = self._make(self._ps_toml())
        self.assertEqual(p.elevation_data[0], [self.basin, -1.5])

    def test_station_prepended_to_breakline_files(self):
        p = self._make(self._ps_toml())
        self.assertEqual(p.breakline_files[0], self.basin)

    def test_missing_basin_file_raises(self):
        section = textwrap.dedent(f"""\
            [[pumping_stations]]
            label = "ps"
            pump_capacity = 1.0
            pump_rate_of_increase = 1.0
            pump_rate_of_decrease = 1.0
            hw_to_start_pumping = 1.0
            hw_to_stop_pumping = 0.5
            basin_polygon_file = "nonexistent.shp"
            basin_elevation = 0.0
            exchange_line_0 = "{_tp(self.ex0)}"
            exchange_line_1 = "{_tp(self.ex1)}"
        """)
        with self.assertRaises(FileNotFoundError):
            self._make(section)


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestProjectDataWrapper(unittest.TestCase):
    """Tests for the ProjectData wrapper in parse_input_data.py."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.toml_path = os.path.join(self.tmpdir, 'config.toml')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _write_minimal(self):
        _write_toml(self.toml_path, _MINIMAL_PROJECT)

    def test_toml_extension_dispatches_correctly(self):
        self._write_minimal()
        p = ProjectData(self.toml_path)
        self.assertEqual(p.scenario, 'test_scenario')

    def test_unsupported_extension_raises(self):
        bad = os.path.join(self.tmpdir, 'config.json')
        _touch(bad)
        with self.assertRaises(ValueError):
            ProjectData(bad)

    def test_defaults_injected_for_toml_only_attrs(self):
        """Attributes introduced in TOML but absent from Excel get defaults."""
        self._write_minimal()
        p = ProjectData(self.toml_path)
        # These are set by _parse_project but also listed in _defaults
        self.assertIsNotNone(p.multiprocessor_mode)
        self.assertIsNone(p.outputstep)
        self.assertIsNone(p.omp_num_threads)
        self.assertFalse(p.report_operator_statistics)


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestCulverts(unittest.TestCase):
    """Tests for ProjectDataTOML._parse_culverts."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _path(self, name):
        return os.path.join(self.tmp, name)

    def _touch(self, name):
        p = self._path(name)
        _touch(p)
        return p

    def _parse(self, toml_body):
        path = self._path('cfg.toml')
        _write_toml(path, _MINIMAL_PROJECT + toml_body)
        with patch('os.path.exists', return_value=True):
            return ProjectDataTOML(path)

    def test_no_culverts(self):
        p = self._parse('')
        self.assertEqual(p.culvert_data, [])

    def test_disabled_culvert_skipped(self):
        body = """
            [[culverts]]
            enabled = false
            type    = "boyd_box"
            label   = "c1"
            width   = 0.9
            exchange_line_0 = "up.csv"
            exchange_line_1 = "down.csv"
        """
        p = self._parse(body)
        self.assertEqual(p.culvert_data, [])

    def test_losses_dict_preserved(self):
        """losses may be a named-component dict (Boyd operators accept it)."""
        body = """
            [[culverts]]
            type    = "boyd_box"
            label   = "box1"
            width   = 0.9
            exchange_line_0 = "up.csv"
            exchange_line_1 = "down.csv"
            losses  = {inlet = 0.5, outlet = 1.0, pier = 0.0}
        """
        cd = self._parse(body).culvert_data[0]
        self.assertEqual(cd['losses'], {'inlet': 0.5, 'outlet': 1.0, 'pier': 0.0})

    def test_losses_list_preserved(self):
        """losses may be a list (Boyd operators sum the entries)."""
        body = """
            [[culverts]]
            type    = "boyd_box"
            label   = "box1"
            width   = 0.9
            exchange_line_0 = "up.csv"
            exchange_line_1 = "down.csv"
            losses  = [0.5, 1.0]
        """
        cd = self._parse(body).culvert_data[0]
        self.assertEqual(cd['losses'], [0.5, 1.0])

    def test_losses_negative_component_raises(self):
        body = """
            [[culverts]]
            type    = "boyd_box"
            label   = "box1"
            width   = 0.9
            exchange_line_0 = "up.csv"
            exchange_line_1 = "down.csv"
            losses  = {inlet = -0.5}
        """
        with self.assertRaises(ValueError):
            self._parse(body)

    def test_boyd_box_defaults(self):
        body = """
            [[culverts]]
            type    = "boyd_box"
            label   = "box1"
            width   = 0.9
            exchange_line_0 = "up.csv"
            exchange_line_1 = "down.csv"
        """
        p = self._parse(body)
        self.assertEqual(len(p.culvert_data), 1)
        cd = p.culvert_data[0]
        self.assertEqual(cd['type'], 'boyd_box')
        self.assertEqual(cd['label'], 'box1')
        self.assertAlmostEqual(cd['width'], 0.9)
        self.assertIsNone(cd['height'])        # defaults to None → operator uses width
        self.assertAlmostEqual(cd['losses'], 0.0)
        self.assertAlmostEqual(cd['barrels'], 1.0)
        self.assertAlmostEqual(cd['manning'], 0.013)
        self.assertAlmostEqual(cd['enquiry_gap'], 0.2)
        self.assertIsNone(cd['invert_elevations'])
        self.assertIsNone(cd['diameter'])

    def test_boyd_box_explicit_height(self):
        body = """
            [[culverts]]
            type    = "boyd_box"
            label   = "box2"
            width   = 1.2
            height  = 0.6
            exchange_line_0 = "up.csv"
            exchange_line_1 = "down.csv"
        """
        p = self._parse(body)
        self.assertAlmostEqual(p.culvert_data[0]['height'], 0.6)

    def test_boyd_pipe_defaults(self):
        body = """
            [[culverts]]
            type     = "boyd_pipe"
            label    = "pipe1"
            diameter = 0.6
            exchange_line_0 = "up.csv"
            exchange_line_1 = "down.csv"
        """
        p = self._parse(body)
        cd = p.culvert_data[0]
        self.assertEqual(cd['type'], 'boyd_pipe')
        self.assertAlmostEqual(cd['diameter'], 0.6)
        self.assertIsNone(cd['width'])
        self.assertIsNone(cd['height'])

    def test_invalid_type_raises(self):
        body = """
            [[culverts]]
            type    = "unknown_type"
            label   = "bad"
            width   = 1.0
            exchange_line_0 = "up.csv"
            exchange_line_1 = "down.csv"
        """
        with patch('os.path.exists', return_value=True):
            path = self._path('bad.toml')
            _write_toml(path, _MINIMAL_PROJECT + body)
            with self.assertRaises(ValueError):
                ProjectDataTOML(path)

    def test_invert_elevations_parsed(self):
        body = """
            [[culverts]]
            type               = "boyd_box"
            label              = "box3"
            width              = 0.9
            exchange_line_0    = "up.csv"
            exchange_line_1    = "down.csv"
            invert_elevations  = [1.5, 1.2]
        """
        p = self._parse(body)
        self.assertEqual(p.culvert_data[0]['invert_elevations'], [1.5, 1.2])

    def test_end_points_instead_of_exchange_lines(self):
        body = """
            [[culverts]]
            type        = "boyd_box"
            label       = "box4"
            width       = 0.9
            end_point_0 = [100.0, 200.0]
            end_point_1 = [110.0, 200.0]
        """
        p = self._parse(body)
        cd = p.culvert_data[0]
        self.assertIsNone(cd['exchange_line_0'])
        self.assertEqual(cd['end_point_0'], [100.0, 200.0])
        self.assertEqual(cd['end_point_1'], [110.0, 200.0])

    def test_missing_exchange_file_raises(self):
        body = """
            [[culverts]]
            type    = "boyd_box"
            label   = "box5"
            width   = 0.9
            exchange_line_0 = "missing_up.csv"
            exchange_line_1 = "missing_down.csv"
        """
        path = self._path('miss.toml')
        _write_toml(path, _MINIMAL_PROJECT + body)
        with self.assertRaises(FileNotFoundError):
            ProjectDataTOML(path)

    def test_two_culverts_both_parsed(self):
        body = """
            [[culverts]]
            type    = "boyd_box"
            label   = "c1"
            width   = 0.9
            exchange_line_0 = "up1.csv"
            exchange_line_1 = "dn1.csv"

            [[culverts]]
            type     = "boyd_pipe"
            label    = "c2"
            diameter = 0.6
            exchange_line_0 = "up2.csv"
            exchange_line_1 = "dn2.csv"
        """
        p = self._parse(body)
        self.assertEqual(len(p.culvert_data), 2)
        self.assertEqual(p.culvert_data[0]['type'], 'boyd_box')
        self.assertEqual(p.culvert_data[1]['type'], 'boyd_pipe')


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestWeirs(unittest.TestCase):
    """Tests for ProjectDataTOML._parse_weirs."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _path(self, name):
        return os.path.join(self.tmp, name)

    def _parse(self, toml_body):
        path = self._path('cfg.toml')
        _write_toml(path, _MINIMAL_PROJECT + toml_body)
        with patch('os.path.exists', return_value=True):
            return ProjectDataTOML(path)

    def test_no_weirs(self):
        p = self._parse('')
        self.assertEqual(p.weir_data, [])

    def test_disabled_weir_skipped(self):
        body = """
            [[weirs]]
            enabled = false
            label   = "w1"
            width   = 3.0
            exchange_line_0 = "up.csv"
            exchange_line_1 = "down.csv"
        """
        p = self._parse(body)
        self.assertEqual(p.weir_data, [])

    def test_weir_defaults(self):
        body = """
            [[weirs]]
            label = "w1"
            width = 3.0
            exchange_line_0 = "up.csv"
            exchange_line_1 = "down.csv"
        """
        p = self._parse(body)
        self.assertEqual(len(p.weir_data), 1)
        wd = p.weir_data[0]
        self.assertEqual(wd['label'], 'w1')
        self.assertAlmostEqual(wd['width'], 3.0)
        self.assertIsNone(wd['height'])
        self.assertAlmostEqual(wd['losses'], 0.0)
        self.assertAlmostEqual(wd['enquiry_gap'], 0.0)
        self.assertAlmostEqual(wd['manning'], 0.013)
        self.assertIsNone(wd['invert_elevations'])

    def test_weir_explicit_height(self):
        body = """
            [[weirs]]
            label  = "w2"
            width  = 3.0
            height = 1.5
            exchange_line_0 = "up.csv"
            exchange_line_1 = "down.csv"
        """
        p = self._parse(body)
        self.assertAlmostEqual(p.weir_data[0]['height'], 1.5)

    def test_weir_invert_elevations(self):
        body = """
            [[weirs]]
            label              = "w3"
            width              = 2.0
            exchange_line_0    = "up.csv"
            exchange_line_1    = "down.csv"
            invert_elevations  = [0.5, 0.3]
        """
        p = self._parse(body)
        self.assertEqual(p.weir_data[0]['invert_elevations'], [0.5, 0.3])

    def test_weir_end_points(self):
        body = """
            [[weirs]]
            label       = "w4"
            width       = 2.0
            end_point_0 = [50.0, 60.0]
            end_point_1 = [55.0, 60.0]
        """
        p = self._parse(body)
        wd = p.weir_data[0]
        self.assertIsNone(wd['exchange_line_0'])
        self.assertEqual(wd['end_point_0'], [50.0, 60.0])

    def test_missing_exchange_file_raises(self):
        body = """
            [[weirs]]
            label           = "w5"
            width           = 2.0
            exchange_line_0 = "missing.csv"
            exchange_line_1 = "also_missing.csv"
        """
        path = self._path('miss.toml')
        _write_toml(path, _MINIMAL_PROJECT + body)
        with self.assertRaises(FileNotFoundError):
            ProjectDataTOML(path)


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestValidation(unittest.TestCase):
    """Tests for the _Validator-based batch error reporting in ProjectDataTOML."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _path(self, name='config.toml'):
        return os.path.join(self.tmpdir, name)

    def _make(self, content):
        p = self._path()
        _write_toml(p, content)
        return ProjectDataTOML(p)

    def _assert_fails(self, content, *expected_fragments):
        """Assert that parsing raises ValueError and that the message contains
        each fragment in *expected_fragments*."""
        p = self._path()
        _write_toml(p, content)
        with self.assertRaises(ValueError) as ctx:
            ProjectDataTOML(p)
        msg = str(ctx.exception)
        for frag in expected_fragments:
            self.assertIn(frag, msg,
                          f'Expected {frag!r} in error message:\n{msg}')

    # --- project section ---

    def test_missing_scenario_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
        """), 'scenario', 'missing')

    def test_missing_yieldstep_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
        """), 'yieldstep', 'missing')

    def test_negative_yieldstep_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = -10.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
        """), 'yieldstep', '> 0')

    def test_invalid_flow_algorithm_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "BADFLOW"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
        """), 'flow_algorithm', 'BADFLOW')

    def test_invalid_multiprocessor_mode_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            multiprocessor_mode = 5
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
        """), 'multiprocessor_mode', '5')

    def test_outputstep_not_multiple_of_yieldstep_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            outputstep = 91.0
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
        """), 'outputstep', 'yieldstep')

    def test_multiple_project_errors_reported_together(self):
        """All errors in the project section must appear in a single ValueError."""
        path = self._path()
        _write_toml(path, textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = -10.0
            finaltime = -1.0
            projection_information = -55
            flow_algorithm = "BAD"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
        """))
        with self.assertRaises(ValueError) as ctx:
            ProjectDataTOML(path)
        msg = str(ctx.exception)
        self.assertIn('yieldstep', msg)
        self.assertIn('finaltime', msg)
        self.assertIn('flow_algorithm', msg)

    # --- mesh section ---

    def test_missing_default_res_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
        """), 'default_res', 'missing')

    def test_negative_default_res_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = -500.0
        """), 'default_res', '> 0')

    def test_negative_interior_region_resolution_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
            [[mesh.interior_regions]]
            polygon = "reg.shp"
            resolution = -100.0
        """), 'resolution', '> 0')

    # --- culvert physical range checks ---

    def test_culvert_negative_width_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
            [[culverts]]
            label = "c1"
            width = -1.0
            end_point_0 = [0.0, 0.0]
            end_point_1 = [10.0, 0.0]
        """), 'width', '> 0')

    def test_culvert_blockage_above_one_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
            [[culverts]]
            label = "c1"
            width = 1.0
            blockage = 1.5
            end_point_0 = [0.0, 0.0]
            end_point_1 = [10.0, 0.0]
        """), 'blockage', '[0.0, 1.0]')

    def test_culvert_negative_manning_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
            [[culverts]]
            label = "c1"
            width = 1.0
            manning = -0.013
            end_point_0 = [0.0, 0.0]
            end_point_1 = [10.0, 0.0]
        """), 'manning', '> 0')

    def test_culvert_pipe_negative_diameter_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
            [[culverts]]
            label = "c1"
            type = "boyd_pipe"
            diameter = -0.5
            end_point_0 = [0.0, 0.0]
            end_point_1 = [10.0, 0.0]
        """), 'diameter', '> 0')

    # --- weir physical range checks ---

    def test_weir_negative_width_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
            [[weirs]]
            label = "w1"
            width = -2.0
            end_point_0 = [0.0, 0.0]
            end_point_1 = [5.0, 0.0]
        """), 'width', '> 0')

    def test_weir_blockage_negative_raises(self):
        self._assert_fails(textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
            [[weirs]]
            label = "w1"
            width = 2.0
            blockage = -0.1
            end_point_0 = [0.0, 0.0]
            end_point_1 = [5.0, 0.0]
        """), 'blockage', '[0.0, 1.0]')

    def test_error_message_names_file(self):
        """The ValueError message must include the TOML filename."""
        path = self._path('my_scenario.toml')
        _write_toml(path, textwrap.dedent("""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = -1.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 1000.0
        """))
        with self.assertRaises(ValueError) as ctx:
            ProjectDataTOML(path)
        self.assertIn('my_scenario.toml', str(ctx.exception))

    def test_valid_config_no_error(self):
        """A fully valid minimal config must parse without errors."""
        p = self._make(textwrap.dedent("""\
            [project]
            scenario = "ok"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE1"
            outputstep = 120.0
            [mesh]
            bounding_polygon = "e.shp"
            default_res = 10000.0
        """))
        self.assertEqual(p.flow_algorithm, 'DE1')
        self.assertAlmostEqual(p.outputstep, 120.0)


if __name__ == '__main__':
    unittest.main()


class TestErosion(unittest.TestCase):
    """[[erosion]] — five behaviours over seven operator classes."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.toml_path = os.path.join(self.tmp, 'cfg.toml')

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _make(self, erosion_toml='', project_extra=''):
        content = textwrap.dedent("""\
            [project]
            scenario = "test_scenario"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
        """) + textwrap.dedent(project_extra) + textwrap.dedent("""
            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 1000000.0
        """) + textwrap.dedent(erosion_toml)
        _write_toml(self.toml_path, content)
        return ProjectDataTOML(self.toml_path)

    def test_absent_by_default(self):
        p = self._make()
        self.assertEqual(p.erosion_data, [])

    def test_polygon_form_parsed(self):
        p = self._make("""
            [[erosion]]
            type = "simple"
            polygon = "scour.csv"
            threshold = 1.5
            base = -2.0
        """)
        self.assertEqual(len(p.erosion_data), 1)
        e = p.erosion_data[0]
        self.assertEqual(e['type'], 'simple')
        self.assertEqual(e['polygon'], 'scour.csv')
        self.assertAlmostEqual(e['threshold'], 1.5)
        self.assertAlmostEqual(e['base'], -2.0)
        self.assertIsNone(e['center'])

    def test_circle_form_parsed(self):
        p = self._make("""
            [[erosion]]
            type = "simple"
            center = [10.0, 20.0]
            radius = 5.0
        """)
        e = p.erosion_data[0]
        self.assertEqual(e['center'], [10.0, 20.0])
        self.assertAlmostEqual(e['radius'], 5.0)
        self.assertIsNone(e['polygon'])

    def test_type_specific_params_kept(self):
        p = self._make("""
            [[erosion]]
            type = "bed_shear"
            polygon = "a.csv"
            shear_factor = 1234.0
        """)
        self.assertAlmostEqual(p.erosion_data[0]['shear_factor'], 1234.0)

    def test_type_specific_param_on_wrong_type_rejected(self):
        # Silently ignoring it would let a user believe they had configured
        # something they had not.
        with self.assertRaises(ValueError) as cm:
            self._make("""
                [[erosion]]
                type = "flat_slice"
                polygon = "a.csv"
                shear_factor = 1234.0
            """)
        self.assertIn('shear_factor', str(cm.exception))

    def test_unknown_type_rejected(self):
        with self.assertRaises(ValueError) as cm:
            self._make("""
                [[erosion]]
                type = "landslide"
                polygon = "a.csv"
            """)
        self.assertIn('landslide', str(cm.exception))

    def test_region_required(self):
        with self.assertRaises(ValueError) as cm:
            self._make("""
                [[erosion]]
                type = "simple"
            """)
        self.assertIn('region', str(cm.exception))

    def test_both_region_forms_rejected(self):
        with self.assertRaises(ValueError) as cm:
            self._make("""
                [[erosion]]
                type = "simple"
                polygon = "a.csv"
                center = [1.0, 2.0]
                radius = 3.0
            """)
        self.assertIn('not both', str(cm.exception))

    def test_erosion_defaults_to_time_varying_elevation(self):
        # With erosion present and no explicit setting, elevation storage
        # defaults to time-varying so the eroded bed is recorded — silently.
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            p = self._make("""
                [[erosion]]
                type = "simple"
                polygon = "a.csv"
            """)
        self.assertTrue(p.store_elevation_every_timestep)
        msgs = ' '.join(str(x.message) for x in w)
        self.assertNotIn('store_elevation_every_timestep', msgs)

    def test_no_default_flip_without_erosion(self):
        # No erosion -> the static default stands.
        p = self._make()
        self.assertFalse(p.store_elevation_every_timestep)

    def test_warns_when_elevation_explicitly_static_with_erosion(self):
        # If the user deliberately chose static storage, respect it but warn
        # that the eroded bed will not be recorded.
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            p = self._make(
                erosion_toml="""
                    [[erosion]]
                    type = "simple"
                    polygon = "a.csv"
                """,
                project_extra="store_elevation_every_timestep = false\n")
        self.assertFalse(p.store_elevation_every_timestep)
        msgs = ' '.join(str(x.message) for x in w)
        self.assertIn('store_elevation_every_timestep', msgs)



class TestArithmeticValues(unittest.TestCase):
    """Numeric fields accept a quoted simple-arithmetic string (e.g. "5*60")."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.toml_path = os.path.join(self.tmpdir, 'config.toml')

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make(self, yieldstep, finaltime, extra=''):
        _write_toml(self.toml_path, textwrap.dedent(f"""\
            [project]
            scenario = "s"
            output_base_directory = "OUTPUT/"
            yieldstep = {yieldstep}
            finaltime = {finaltime}
            projection_information = -55
            flow_algorithm = "DE0"
            {extra}

            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 1000000.0
        """))
        return ProjectDataTOML(self.toml_path)

    def test_arithmetic_time_fields(self):
        p = self._make(yieldstep='"2*30"', finaltime='"5*60"',
                       extra='outputstep = "10*60"')
        self.assertAlmostEqual(p.yieldstep, 60.0)
        self.assertAlmostEqual(p.finaltime, 300.0)
        self.assertAlmostEqual(p.outputstep, 600.0)

    def test_plain_numbers_still_work(self):
        p = self._make(yieldstep=60.0, finaltime=3600.0)
        self.assertAlmostEqual(p.yieldstep, 60.0)
        self.assertAlmostEqual(p.finaltime, 3600.0)

    def test_parentheses_and_power(self):
        p = self._make(yieldstep='"(1+1)*30"', finaltime='"2**3 * 60"')
        self.assertAlmostEqual(p.yieldstep, 60.0)
        self.assertAlmostEqual(p.finaltime, 480.0)

    def test_bad_expression_reported(self):
        with self.assertRaises(ValueError) as cm:
            self._make(yieldstep=60.0, finaltime='"5 * "')
        self.assertIn('finaltime', str(cm.exception))

    def test_names_and_calls_rejected(self):
        with self.assertRaises(ValueError):
            self._make(yieldstep=60.0, finaltime='"__import__(\'os\')"')


if __name__ == '__main__':
    unittest.main()


class TestSediment(unittest.TestCase):
    """[sediment] and [[sediment.fractions]]."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.toml_path = os.path.join(self.tmp, 'cfg.toml')

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _make(self, sediment_toml='', project_extra=''):
        content = textwrap.dedent("""\
            [project]
            scenario = "test_scenario"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 3600.0
            projection_information = -55
            flow_algorithm = "DE0"
        """) + textwrap.dedent(project_extra) + textwrap.dedent("""
            [mesh]
            bounding_polygon = "extent.shp"
            default_res = 1000000.0
        """) + textwrap.dedent(sediment_toml)
        _write_toml(self.toml_path, content)
        return ProjectDataTOML(self.toml_path)

    ONE_FRACTION = """
        [sediment]
        [[sediment.fractions]]
        name = "sand"
        diameter = 2.0e-4
    """

    def test_absent_by_default(self):
        p = self._make()
        self.assertIsNone(p.sediment_data)
        self.assertFalse(p.store_elevation_every_timestep)

    def test_minimal_parses_and_leaves_defaults_to_the_domain(self):
        p = self._make(self.ONE_FRACTION)
        d = p.sediment_data
        self.assertEqual(d['fractions'], [{'name': 'sand', 'diameter': 2.0e-4}])
        self.assertEqual(d['erodible_regions'], [])
        for key in ('porosity', 'shear_closure', 'bed_material', 'bedload'):
            self.assertNotIn(key, d)

    def test_grass_needs_its_coefficient(self):
        with self.assertRaises(Exception) as cm:
            self._make(self.ONE_FRACTION.replace(
                '[sediment]', '[sediment]\nbedload = "grass"'))
        self.assertIn('bedload_K', str(cm.exception))

    def test_grass_and_open_boundaries_parse(self):
        p = self._make(self.ONE_FRACTION.replace(
            '[sediment]',
            '[sediment]\nbedload = "grass"\nbedload_K = 0.001\n'
            'bedload_open_boundaries = ["left", "right"]'))
        d = p.sediment_data
        self.assertEqual(d['bedload'], 'grass')
        self.assertAlmostEqual(d['bedload_K'], 0.001)
        self.assertEqual(d['bedload_open_boundaries'], ['left', 'right'])

    def test_open_boundaries_must_be_tags(self):
        with self.assertRaises(Exception) as cm:
            self._make(self.ONE_FRACTION.replace(
                '[sediment]', '[sediment]\nbedload_open_boundaries = 3'))
        self.assertIn('bedload_open_boundaries', str(cm.exception))

    def test_de_leeuw_entrainment_parses(self):
        p = self._make(self.ONE_FRACTION.replace(
            '[sediment]',
            '[sediment]\nentrainment = "de_leeuw"\nde_leeuw_fit = "nghiem_2022"\n'
            'skin_roughness = 0.002'))
        d = p.sediment_data
        self.assertEqual(d['entrainment'], 'de_leeuw')
        self.assertEqual(d['de_leeuw_fit'], 'nghiem_2022')
        self.assertAlmostEqual(d['skin_roughness'], 0.002)

    def test_de_leeuw_needs_a_noncohesive_bed(self):
        with self.assertRaises(Exception) as cm:
            self._make(self.ONE_FRACTION.replace(
                '[sediment]', '[sediment]\nentrainment = "de_leeuw"\nbed_material = "cohesive"'))
        self.assertIn('de_leeuw', str(cm.exception))

    def test_full_table_parses(self):
        p = self._make("""
            [sediment]
            porosity = 0.28
            c_max = 0.25
            bed_evolution = true
            rho_w = 1025.0
            shear_closure = "energy_slope"
            bed_material = "cohesive"
            tau_crit = 0.088
            deposition_law = "threshold"
            tau_d = 0.1
            near_bed = "rouse"
            friction_mode = "wilson"
            bed = "gravel"
            grain_size = 0.02
            bedload = "wong_parker_eq24"
            bedload_K = 3.97
            angle_of_repose = 34.0
            erodible_base_depth = "3/2"

            [[sediment.fractions]]
            name = "sand"
            diameter = 2.0e-4
            tau_c_star = 0.045
            initial_concentration = 0.001
            [sediment.fractions.boundary]
            left = 0.01
            right = "1/100"

            [[sediment.fractions]]
            name = "silt"
            diameter = 2.0e-5

            [[sediment.erodible_regions]]
            polygon = "channel.csv"
            [[sediment.erodible_regions]]
            center = [10.0, 20.0]
            radius = 5.0
            erodible = false
        """)
        d = p.sediment_data
        self.assertAlmostEqual(d['porosity'], 0.28)
        self.assertEqual(d['shear_closure'], 'energy_slope')
        self.assertEqual(d['bed_material'], 'cohesive')
        self.assertAlmostEqual(d['tau_crit'], 0.088)
        self.assertEqual(d['deposition_law'], 'threshold')
        self.assertEqual(d['friction_mode'], 'wilson')
        self.assertEqual(d['bed'], 'gravel')
        self.assertEqual(d['bedload'], 'wong_parker_eq24')
        self.assertAlmostEqual(d['bedload_K'], 3.97)
        self.assertAlmostEqual(d['angle_of_repose'], 34.0)
        self.assertAlmostEqual(d['erodible_base_depth'], 1.5)   # arithmetic string
        self.assertEqual(len(d['fractions']), 2)
        sand = d['fractions'][0]
        self.assertAlmostEqual(sand['tau_c_star'], 0.045)
        self.assertEqual(sand['boundary'], {'left': 0.01, 'right': 0.01})
        self.assertEqual(d['fractions'][1], {'name': 'silt', 'diameter': 2.0e-5})
        regions = d['erodible_regions']
        self.assertEqual(regions[0]['polygon'], 'channel.csv')
        self.assertTrue(regions[0]['erodible'])
        self.assertEqual(regions[1]['center'], [10.0, 20.0])
        self.assertFalse(regions[1]['erodible'])
        self.assertTrue(p.store_elevation_every_timestep)

    def test_needs_a_fraction(self):
        with self.assertRaises(ValueError) as cm:
            self._make("""
                [sediment]
                porosity = 0.3
            """)
        self.assertIn('at least one [[sediment.fractions]]', str(cm.exception))

    def test_bad_values_are_all_reported(self):
        with self.assertRaises(ValueError) as cm:
            self._make("""
                [sediment]
                porosity = 1.5
                shear_closure = "bed_slope"
                friction_mode = "larsen_lamb"
                tau_crit = 0.1
                erodible_base_elevation = 0.0
                erodible_base_depth = 1.0
                nonsense = 3
                [[sediment.fractions]]
                name = "sand"
                diameter = -1.0
                [[sediment.fractions]]
                name = "sand"
                diameter = 1.0e-4
            """)
        msg = str(cm.exception)
        for fragment in ("'porosity'", "'shear_closure'", "larsen_lamb",
                         "'tau_crit' / 'K_e'", "not both", "unknown key 'nonsense'",
                         "'diameter'", "duplicate fraction name"):
            self.assertIn(fragment, msg)

    def test_slope_bounds_parse_and_are_checked(self):
        p = self._make("""
            [sediment]
            shear_closure = "depth_slope"
            max_slope = 0.05
            freeze_slope = true
            [[sediment.fractions]]
            name = "sand"
            diameter = 2.0e-4
        """)
        self.assertAlmostEqual(p.sediment_data['max_slope'], 0.05)
        self.assertTrue(p.sediment_data['freeze_slope'])
        with self.assertRaises(ValueError) as cm:
            self._make("""
                [sediment]
                shear_closure = "energy_slope"
                freeze_slope = true
                max_slope = -1.0
                [[sediment.fractions]]
                name = "sand"
                diameter = 2.0e-4
            """)
        self.assertIn("'freeze_slope' applies", str(cm.exception))
        self.assertIn("'max_slope'", str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            self._make("""
                [sediment]
                max_slope = 0.1
                [[sediment.fractions]]
                name = "sand"
                diameter = 2.0e-4
            """)
        self.assertIn("quadratic_drag uses no slope", str(cm.exception))

    def test_fixed_bed_does_not_flip_elevation_storage(self):
        p = self._make("""
            [sediment]
            bed_evolution = false
            [[sediment.fractions]]
            name = "sand"
            diameter = 2.0e-4
        """)
        self.assertFalse(p.store_elevation_every_timestep)

    def test_warns_when_elevation_explicitly_static_with_bed_evolution(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            p = self._make(self.ONE_FRACTION,
                           project_extra='store_elevation_every_timestep = false\n')
        self.assertFalse(p.store_elevation_every_timestep)
        self.assertTrue(any('sediment' in str(x.message) for x in w))

