"""Unit tests for anuga.scenario.scenario_summary (the --dry-run HTML summary)."""
import os
import tempfile
import shutil
import textwrap
import unittest

try:
    from anuga.scenario.scenario_summary import (
        build_summary_html, write_scenario_summary,
        _read_timeseries, _aggregate_rainfall, _friction_tier)
    HAS_MODULE = True
    SKIP = ''
except ImportError as _e:  # pragma: no cover
    HAS_MODULE = False
    SKIP = str(_e)


_BASE = """\
[project]
scenario = "demo_scenario"
output_base_directory = "OUTPUT/"
yieldstep = 60.0
finaltime = 3600.0
projection_information = "EPSG:32756"
flow_algorithm = "DE0"
compute_mode = "unified"

[mesh]
bounding_polygon = "extent.csv"
default_res = 500.0
[[mesh.interior_regions]]
polygon = "refine.csv"
resolution = 50.0

[boundary_conditions]
[[boundary_conditions.boundaries]]
tag = "west"
type = "Reflective"
[[boundary_conditions.boundaries]]
tag = "east"
type = "Reflective"

[initial_conditions]
[[initial_conditions.elevation]]
polygon = "All"
value = 0.0
[[initial_conditions.friction]]
polygon = "creek.csv"
value = 0.03
[[initial_conditions.friction]]
polygon = "houses.csv"
value = 10.0
[[initial_conditions.friction]]
polygon = "All"
value = 0.04
"""


@unittest.skipUnless(HAS_MODULE, SKIP)
class TestScenarioSummary(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, body, name='scenario.toml'):
        path = os.path.join(self.tmp, name)
        with open(path, 'w') as fh:
            fh.write(textwrap.dedent(body))
        return path

    # ---- core rendering ----------------------------------------------------

    def test_build_html_core(self):
        html = build_summary_html(self._write(_BASE))
        self.assertIn('<title>', html)
        self.assertIn('demo_scenario'.replace('_', ' '), html)
        self.assertIn('class="stats"', html)
        # badges reflect project settings
        self.assertIn('EPSG:32756', html)
        self.assertIn('unified', html)
        # boundary tags appear
        self.assertIn('west', html)
        self.assertIn('east', html)
        # theme-aware + self-contained (inline style, no external asset URLs)
        self.assertIn('prefers-color-scheme', html)
        self.assertNotIn('http://', html)
        self.assertNotIn('https://', html)

    def test_sediment_fractions_counted(self):
        html = build_summary_html(self._write(_BASE + """
            [sediment]
            porosity = 0.3
            [[sediment.fractions]]
            name = "sand"
            diameter = 2.0e-4
            [[sediment.fractions]]
            name = "silt"
            diameter = 2.0e-5
        """))
        self.assertIn('sediment fractions', html)
        self.assertIn('sand, silt', html)

    def test_friction_breakdown_present(self):
        html = build_summary_html(self._write(_BASE))
        self.assertIn("Friction", html)
        # three distinct n-values -> three bars
        self.assertEqual(html.count('class="fbar"'), 3)

    def test_quoted_times_and_starttime_badge(self):
        body = _BASE.replace('finaltime = 3600.0',
                             'finaltime = "3600"\nstarttime = "-3*600"')
        html = build_summary_html(self._write(body))
        self.assertIn('-1800 s', html)
        self.assertIn('3600 s', html)

    def test_multiprocessor_mode_maps_to_compute_label(self):
        body = _BASE.replace('compute_mode = "unified"', 'multiprocessor_mode = 2')
        html = build_summary_html(self._write(body))
        self.assertIn('unified', html)  # 2 -> unified

    def test_missing_optional_sections_degrade(self):
        minimal = textwrap.dedent("""\
            [project]
            scenario = "bare"
            output_base_directory = "OUTPUT/"
            yieldstep = 60.0
            finaltime = 600.0
            projection_information = -56
            flow_algorithm = "DE0"
            [mesh]
            bounding_polygon = "extent.csv"
            default_res = 1000.0
        """)
        html = build_summary_html(self._write(minimal))
        self.assertIn('<title>', html)
        # no rainfall / structures sections, but must not error. (The CSS always
        # defines .hyeto; the chart is only present as an SVG class="hyeto".)
        self.assertNotIn('<h2>Structures', html)
        self.assertNotIn('class="hyeto"', html)

    # ---- rainfall / hyetograph --------------------------------------------

    def test_csv_rainfall_hyetograph(self):
        # a simple triangular hyetograph in mm/hr
        rain = os.path.join(self.tmp, 'rain.csv')
        with open(rain, 'w') as fh:
            fh.write('time,rate_mm_hr\n0,0\n1800,20\n3600,0\n')
        body = _BASE + textwrap.dedent("""\
            [[rainfall]]
            timeseries_file = "rain.csv"
            polygon = "All"
            start_time = 0.0
        """)
        html = build_summary_html(self._write(body))
        self.assertIn('class="hyeto"', html)
        self.assertIn('peak mm/hr', html)

    def test_read_timeseries_csv_units(self):
        rain = os.path.join(self.tmp, 'r.csv')
        with open(rain, 'w') as fh:
            fh.write('time,mmhr\n0,36\n60,36\n')   # 36 mm/hr -> 0.01 mm/s
        t, r = _read_timeseries(rain)
        self.assertAlmostEqual(r[0], 0.01, places=6)

    def test_aggregate_rainfall_stats(self):
        rain = os.path.join(self.tmp, 'r.csv')
        with open(rain, 'w') as fh:
            fh.write('time,mmhr\n0,0\n1800,60\n3600,0\n')  # peak 60 mm/hr at 0.5 h
        agg = _aggregate_rainfall([{'timeseries_file': 'r.csv'}], self.tmp)
        self.assertIsNotNone(agg)
        self.assertAlmostEqual(agg['peak'], 60.0, delta=1.0)
        self.assertAlmostEqual(agg['peak_t_h'], 0.5, delta=0.1)
        # triangle area: 0.5 * base(1h=3600s) * height(60mm/hr=0.01667mm/s) ~ 30 mm
        self.assertAlmostEqual(agg['depth_mm'], 30.0, delta=1.5)

    def test_no_rainfall_returns_none(self):
        self.assertIsNone(_aggregate_rainfall([], self.tmp))

    # ---- inlets: constant rate vs time-varying hydrograph ------------------

    def _inlet_body(self, name='hydro.csv'):
        return _BASE + textwrap.dedent(f"""\
            [[inlets]]
            name = "creek"
            line_file = "line.csv"
            timeseries_file = "{name}"
            start_time = 0.0
        """)

    def test_inlet_constant_rate_reported(self):
        with open(os.path.join(self.tmp, 'hydro.csv'), 'w') as fh:
            fh.write('time,discharge\n0,20\n1000000,20\n')
        html = build_summary_html(self._write(self._inlet_body()))
        self.assertIn('constant 20 m³/s', html)
        self.assertNotIn('class="ln-reed"', html)   # no chart element

    def test_inlet_varying_gets_chart(self):
        with open(os.path.join(self.tmp, 'hydro.csv'), 'w') as fh:
            fh.write('time,discharge\n0,0\n1800,50\n3600,10\n')
        html = build_summary_html(self._write(self._inlet_body()))
        self.assertIn('class="ln-reed"', html)       # chart element present
        self.assertIn('peak 50', html)
        self.assertIn('mean', html)

    def test_inlet_missing_timeseries_falls_back(self):
        html = build_summary_html(self._write(self._inlet_body('nope.csv')))
        self.assertIn('creek', html)
        self.assertIn('line source', html)           # generic fallback, no crash

    def test_is_constant_helper(self):
        from anuga.scenario.scenario_summary import _is_constant
        self.assertTrue(_is_constant([5.0, 5.0, 5.0]))
        self.assertFalse(_is_constant([0.0, 5.0, 1.0]))

    # ---- friction tiers ----------------------------------------------------

    def test_friction_tiers(self):
        self.assertEqual(_friction_tier(0.01)[0], 'water')
        self.assertEqual(_friction_tier(0.03)[0], 'stone')
        self.assertEqual(_friction_tier(0.15)[0], 'reed')
        self.assertEqual(_friction_tier(10.0)[0], 'silt')

    # ---- file output -------------------------------------------------------

    def test_write_creates_file(self):
        cfg = self._write(_BASE)
        out = write_scenario_summary(cfg, open_browser=False)
        self.assertTrue(os.path.exists(out))
        self.assertTrue(out.endswith('scenario_summary.html'))
        self.assertGreater(os.path.getsize(out), 1000)

    def test_write_custom_output_path(self):
        cfg = self._write(_BASE)
        dest = os.path.join(self.tmp, 'custom.html')
        out = write_scenario_summary(cfg, output_html=dest, open_browser=False)
        self.assertEqual(out, dest)
        self.assertTrue(os.path.exists(dest))


if __name__ == '__main__':
    unittest.main()
