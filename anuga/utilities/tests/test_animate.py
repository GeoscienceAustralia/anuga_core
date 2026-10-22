"""
Tests for anuga.utilities.animate

Covers:
  - _face_to_vertex     (pure numpy helper)
  - _nice_contour_levels (pure math helper)
  - _draw_elev_contours  (matplotlib helper)
  - SWW_plotter          (init, save_*_frame, plot_mesh, counters, absolute mode)
  - Domain_plotter       (init with mock domain, save_depth_frame)

All matplotlib work uses the Agg backend to avoid any display requirement.
"""

import matplotlib
matplotlib.use('Agg')

import os
import tempfile
import shutil
import unittest
import numpy as np

try:
    from anuga.utilities.animate import (
        SWW_plotter, Domain_plotter,
        _face_to_vertex, _nice_contour_levels, _draw_elev_contours,
        BASEMAP_PROVIDERS, BASEMAP_DEFAULT,
    )
    HAS_MODULE = True
    SKIP_REASON = ''
except ImportError as _e:
    HAS_MODULE = False
    SKIP_REASON = str(_e)


# ---------------------------------------------------------------------------
# SWW fixture factory
# ---------------------------------------------------------------------------

def _write_minimal_sww(path, n_timesteps=3, with_epsg=False):
    """Write a tiny synthetic SWW file (2 triangles, 4 nodes)."""
    from anuga.file.netcdf import NetCDFFile
    from anuga.config import netcdf_mode_w

    n_nodes = 4
    n_tris = 2

    # Mesh: axis-aligned square split on the diagonal
    x    = np.array([0.0, 2.0, 2.0, 0.0], dtype=np.float32)
    y    = np.array([0.0, 0.0, 1.0, 1.0], dtype=np.float32)
    vols = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int32)
    t    = np.array([0.0, 10.0, 20.0], dtype=np.float64)

    elev  = np.array([-1.0, -0.5], dtype=np.float32)
    stage = np.array([[0.5, 0.3],
                      [0.6, 0.4],
                      [0.7, 0.5]], dtype=np.float32)
    xmom  = np.zeros((n_timesteps, n_tris), dtype=np.float32)
    ymom  = np.zeros((n_timesteps, n_tris), dtype=np.float32)

    fid = NetCDFFile(path, netcdf_mode_w)
    fid.starttime = 0.0
    fid.timezone  = 'UTC'
    fid.xllcorner = 300000.0
    fid.yllcorner = 6000000.0
    fid.zone      = 55

    if with_epsg:
        fid.epsg = 32755

    fid.createDimension('number_of_volumes',   n_tris)
    fid.createDimension('number_of_vertices',  3)
    fid.createDimension('number_of_points',    n_nodes)
    fid.createDimension('number_of_timesteps', n_timesteps)

    fid.createVariable('x',            'f', ('number_of_points',))
    fid.createVariable('y',            'f', ('number_of_points',))
    fid.createVariable('volumes',      'i', ('number_of_volumes', 'number_of_vertices'))
    fid.createVariable('time',         'd', ('number_of_timesteps',))
    fid.createVariable('elevation_c',  'f', ('number_of_volumes',))
    fid.createVariable('stage_c',      'f', ('number_of_timesteps', 'number_of_volumes'))
    fid.createVariable('xmomentum_c',  'f', ('number_of_timesteps', 'number_of_volumes'))
    fid.createVariable('ymomentum_c',  'f', ('number_of_timesteps', 'number_of_volumes'))

    fid.variables['x'][:]           = x
    fid.variables['y'][:]           = y
    fid.variables['volumes'][:]     = vols
    fid.variables['time'][:]        = t
    fid.variables['elevation_c'][:] = elev
    fid.variables['stage_c'][:]     = stage
    fid.variables['xmomentum_c'][:] = xmom
    fid.variables['ymomentum_c'][:] = ymom
    fid.close()


def _make_mock_domain(plot_dir):
    """Return a minimal MagicMock that satisfies Domain_plotter.__init__."""
    from unittest.mock import MagicMock

    nodes    = np.array([[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]])
    tris     = np.array([[0, 1, 2], [0, 2, 3]])
    centroids= np.array([[4/3, 1/3], [2/3, 2/3]])

    geo = MagicMock()
    geo.zone      = 55
    geo.xllcorner = 300000.0
    geo.yllcorner = 6000000.0
    geo.epsg      = None

    def _qty(vals):
        q = MagicMock()
        q.centroid_values = np.array(vals, dtype=float)
        return q

    domain = MagicMock()
    domain.geo_reference        = geo
    domain.nodes                = nodes
    domain.triangles            = tris
    domain.centroid_coordinates = centroids
    domain.quantities = {
        'elevation':  _qty([-1.0, -0.5]),
        'stage':      _qty([0.5,   0.3]),
        'xmomentum':  _qty([0.0,   0.0]),
        'ymomentum':  _qty([0.0,   0.0]),
        'friction':   _qty([0.025, 0.025]),
    }
    domain.get_name.return_value     = os.path.join(plot_dir, 'test_dom')
    domain.get_time.return_value     = 0.0
    domain.get_timestep.return_value = 1.0
    return domain


# ---------------------------------------------------------------------------
# _face_to_vertex
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestFaceToVertex(unittest.TestCase):

    def _tri(self, x, y, triangles):
        import matplotlib.tri as mtri
        return mtri.Triangulation(x, y, triangles)

    def test_all_same_value(self):
        x = np.array([0., 1., 0., 1.])
        y = np.array([0., 0., 1., 1.])
        t = np.array([[0, 1, 2], [1, 3, 2]])
        triang = self._tri(x, y, t)
        result = _face_to_vertex(triang, [5.0, 5.0])
        np.testing.assert_allclose(result, 5.0)

    def test_different_values_shared_vertex(self):
        # Triangle 0 uses vertices 0,1,2; triangle 1 uses 0,2,3.
        # Vertex 2 is shared → its average should be (v0 + v1) / 2.
        x = np.array([0., 2., 2., 0.])
        y = np.array([0., 0., 1., 1.])
        t = np.array([[0, 1, 2], [0, 2, 3]])
        triang = self._tri(x, y, t)
        face_values = np.array([4.0, 2.0])
        result = _face_to_vertex(triang, face_values)
        # Vertex 0: contributes tri0(4) and tri1(2) → mean 3.0
        self.assertAlmostEqual(result[0], 3.0)
        # Vertex 1: only tri0 → 4.0
        self.assertAlmostEqual(result[1], 4.0)
        # Vertex 2: tri0 and tri1 → mean 3.0
        self.assertAlmostEqual(result[2], 3.0)
        # Vertex 3: only tri1 → 2.0
        self.assertAlmostEqual(result[3], 2.0)

    def test_face_mask_excludes_triangles(self):
        x = np.array([0., 2., 2., 0.])
        y = np.array([0., 0., 1., 1.])
        t = np.array([[0, 1, 2], [0, 2, 3]])
        triang = self._tri(x, y, t)
        # Mask first triangle → only tri1 (value 2.0) contributes
        mask = np.array([True, False])
        result = _face_to_vertex(triang, [4.0, 2.0], face_mask=mask)
        # Vertices touched only by tri0 (0→masked) fall back to 0.0
        self.assertAlmostEqual(result[1], 0.0)   # only in tri0
        self.assertAlmostEqual(result[3], 2.0)   # only in tri1

    def test_all_masked_returns_zeros(self):
        x = np.array([0., 1., 0.])
        y = np.array([0., 0., 1.])
        t = np.array([[0, 1, 2]])
        triang = self._tri(x, y, t)
        result = _face_to_vertex(triang, [7.0], face_mask=[True])
        np.testing.assert_array_equal(result, 0.0)

    def test_output_length_equals_n_vertices(self):
        x = np.array([0., 1., 0., 1.])
        y = np.array([0., 0., 1., 1.])
        t = np.array([[0, 1, 2], [1, 3, 2]])
        triang = self._tri(x, y, t)
        result = _face_to_vertex(triang, [1.0, 2.0])
        self.assertEqual(len(result), 4)


# ---------------------------------------------------------------------------
# _nice_contour_levels
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestNiceContourLevels(unittest.TestCase):

    def test_zero_span_returns_fallback(self):
        result = _nice_contour_levels(5.0, 5.0, 10)
        self.assertEqual(result, 10)

    def test_round_levels_in_range(self):
        levels = _nice_contour_levels(0.0, 100.0, 5)
        self.assertIsInstance(levels, np.ndarray)
        self.assertGreaterEqual(len(levels), 2)
        # All levels should be within [vmin, vmax]
        self.assertTrue(np.all(levels >= 0.0))
        self.assertTrue(np.all(levels <= 100.0 + 1e-6))

    def test_levels_are_multiples_of_step(self):
        levels = _nice_contour_levels(0.0, 50.0, 5)
        if isinstance(levels, np.ndarray) and len(levels) >= 2:
            step = levels[1] - levels[0]
            # All levels should be (approx) multiples of step
            remainder = np.mod(levels, step)
            np.testing.assert_allclose(remainder, 0.0, atol=1e-6)

    def test_small_range(self):
        levels = _nice_contour_levels(0.0, 2.0, 5)
        if isinstance(levels, np.ndarray):
            self.assertGreaterEqual(len(levels), 2)

    def test_large_range(self):
        levels = _nice_contour_levels(0.0, 10000.0, 8)
        if isinstance(levels, np.ndarray):
            self.assertGreaterEqual(len(levels), 2)
            # Step should be at least 100 for this range
            self.assertGreater(levels[1] - levels[0], 50.0)


# ---------------------------------------------------------------------------
# SWW_plotter
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestSWWPlotter(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.swwfile = os.path.join(self.tmpdir, 'test.sww')
        _write_minimal_sww(self.swwfile)
        self.plot_dir = os.path.join(self.tmpdir, 'frames')

    def tearDown(self):
        import matplotlib.pyplot as plt
        plt.close('all')
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make(self, **kwargs):
        return SWW_plotter(self.swwfile,
                           plot_dir=self.plot_dir,
                           **kwargs)

    # --- initialisation ---

    def test_init_loads_geometry(self):
        sp = self._make()
        self.assertEqual(len(sp.x), 4)
        self.assertEqual(len(sp.y), 4)
        self.assertEqual(sp.triangles.shape, (2, 3))

    def test_init_loads_time(self):
        sp = self._make()
        np.testing.assert_allclose(sp.time, [0.0, 10.0, 20.0])

    def test_init_loads_stage_and_elev(self):
        sp = self._make()
        self.assertEqual(sp.stage.shape, (3, 2))
        # elevation is static (1-D)
        self.assertEqual(sp.elev.shape, (2,))

    def test_init_depth_computed(self):
        sp = self._make()
        # depth[0,:] = stage[0,:] - elev = [0.5-(-1), 0.3-(-0.5)] = [1.5, 0.8]
        np.testing.assert_allclose(sp.depth[0], [1.5, 0.8], atol=1e-5)

    def test_init_speed_non_negative(self):
        sp = self._make()
        self.assertTrue(np.all(sp.speed >= 0.0))

    def test_init_counters_zero(self):
        sp = self._make()
        self.assertEqual(sp._depth_frame_count,     0)
        self.assertEqual(sp._stage_frame_count,     0)
        self.assertEqual(sp._speed_frame_count,     0)
        self.assertEqual(sp._elev_frame_count,      0)
        self.assertEqual(sp._max_depth_frame_count, 0)

    def test_init_geo_attributes(self):
        sp = self._make()
        self.assertAlmostEqual(sp.xllcorner, 300000.0)
        self.assertAlmostEqual(sp.yllcorner, 6000000.0)
        self.assertEqual(sp.zone, 55)
        self.assertEqual(sp.timezone, 'UTC')

    def test_init_without_timezone_defaults_to_utc(self):
        """Older SWW files have no timezone attribute."""
        from anuga.file.netcdf import NetCDFFile
        fid = NetCDFFile(self.swwfile, 'a')
        fid.delncattr('timezone')
        fid.close()
        sp = self._make()
        self.assertEqual(sp.timezone, 'UTC')
        self.assertEqual(sp.depth.shape, (3, 2))

    def test_init_absolute_mode_shifts_coords(self):
        sp_rel = self._make(absolute=False)
        sp_abs = self._make(absolute=True)
        np.testing.assert_allclose(sp_abs.x, sp_rel.x + sp_rel.xllcorner, atol=1.0)
        np.testing.assert_allclose(sp_abs.y, sp_rel.y + sp_rel.yllcorner, atol=1.0)

    def test_init_absolute_mode_triangulation_is_absolute(self):
        # Regression: triang was built before the absolute shift, and
        # Triangulation copies its coordinate arrays, so triang.x stayed
        # relative while x/xc were absolute -- anything plotted against
        # absolute axes landed off the map.
        sp = self._make(absolute=True)
        np.testing.assert_allclose(sp.triang.x, sp.x)
        np.testing.assert_allclose(sp.triang.y, sp.y)

    def test_name_attribute(self):
        sp = self._make()
        self.assertEqual(sp.name, 'test')

    # --- save_depth_frame ---

    def test_save_depth_frame_creates_png(self):
        sp = self._make()
        sp.save_depth_frame(frame=0, dpi=20)
        png = os.path.join(self.plot_dir, 'test_depth_0000000000.png')
        self.assertTrue(os.path.exists(png), f'Expected {png}')

    def test_save_depth_frame_increments_counter(self):
        sp = self._make()
        sp.save_depth_frame(frame=0, dpi=20)
        sp.save_depth_frame(frame=1, dpi=20)
        self.assertEqual(sp._depth_frame_count, 2)
        png1 = os.path.join(self.plot_dir, 'test_depth_0000000001.png')
        self.assertTrue(os.path.exists(png1))

    def test_save_depth_frame_last(self):
        sp = self._make()
        sp.save_depth_frame(frame=-1, dpi=20)
        png = os.path.join(self.plot_dir, 'test_depth_0000000000.png')
        self.assertTrue(os.path.exists(png))

    # --- save_stage_frame ---

    def test_save_stage_frame_creates_png(self):
        sp = self._make()
        sp.save_stage_frame(frame=0, dpi=20)
        png = os.path.join(self.plot_dir, 'test_stage_0000000000.png')
        self.assertTrue(os.path.exists(png))

    def test_save_stage_frame_increments_counter(self):
        sp = self._make()
        sp.save_stage_frame(frame=0, dpi=20)
        self.assertEqual(sp._stage_frame_count, 1)

    # --- save_speed_frame ---

    def test_save_speed_frame_creates_png(self):
        sp = self._make()
        sp.save_speed_frame(frame=0, dpi=20)
        png = os.path.join(self.plot_dir, 'test_speed_0000000000.png')
        self.assertTrue(os.path.exists(png))

    def test_save_speed_frame_increments_counter(self):
        sp = self._make()
        sp.save_speed_frame(frame=0, dpi=20)
        self.assertEqual(sp._speed_frame_count, 1)

    # --- save_elev_frame ---

    def test_save_elev_frame_static_creates_png(self):
        sp = self._make()
        sp.save_elev_frame(frame=0, dpi=20)
        png = os.path.join(self.plot_dir, 'test_elev_0000000000.png')
        self.assertTrue(os.path.exists(png))

    def test_save_elev_frame_increments_counter(self):
        sp = self._make()
        sp.save_elev_frame(frame=0, dpi=20)
        self.assertEqual(sp._elev_frame_count, 1)

    # --- save_max_depth_frame ---

    def test_save_max_depth_frame_creates_png(self):
        sp = self._make()
        sp.save_max_depth_frame(dpi=20)
        png = os.path.join(self.plot_dir, 'test_max_depth_0000000000.png')
        self.assertTrue(os.path.exists(png))

    def test_save_max_depth_frame_increments_counter(self):
        sp = self._make()
        sp.save_max_depth_frame(dpi=20)
        self.assertEqual(sp._max_depth_frame_count, 1)

    # --- plot helpers ---

    def test_plot_mesh_returns_fig_ax_im(self):
        sp = self._make()
        result = sp.plot_mesh()
        self.assertEqual(len(result), 3)
        fig, ax, im = result
        self.assertIsNotNone(fig)
        self.assertIsNotNone(ax)

    def test_plot_depth_frame_returns_fig_ax(self):
        sp = self._make()
        result = sp.plot_depth_frame(frame=0, dpi=20)
        self.assertEqual(len(result), 2)
        fig, ax = result
        self.assertIsNotNone(fig)

    def test_plot_stage_frame_returns_fig_ax(self):
        sp = self._make()
        result = sp.plot_stage_frame(frame=0, dpi=20)
        self.assertEqual(len(result), 2)

    def test_plot_elev_frame_returns_fig_ax(self):
        sp = self._make()
        result = sp.plot_elev_frame(frame=0, dpi=20)
        self.assertEqual(len(result), 2)

    # --- xlim / ylim zoom ---

    def test_save_depth_frame_with_zoom(self):
        sp = self._make()
        sp.save_depth_frame(frame=0, dpi=20, xlim=(0, 1), ylim=(0, 0.5))
        png = os.path.join(self.plot_dir, 'test_depth_0000000000.png')
        self.assertTrue(os.path.exists(png))

    # --- show_mesh overlay ---

    def test_save_depth_frame_show_mesh(self):
        sp = self._make()
        sp.save_depth_frame(frame=0, dpi=20, show_mesh=True)
        png = os.path.join(self.plot_dir, 'test_depth_0000000000.png')
        self.assertTrue(os.path.exists(png))

    # --- show_elev contours ---

    def test_save_depth_frame_show_elev_contours(self):
        sp = self._make()
        sp.save_depth_frame(frame=0, dpi=20, show_elev=True, elev_levels=3)
        png = os.path.join(self.plot_dir, 'test_depth_0000000000.png')
        self.assertTrue(os.path.exists(png))

    # --- figure cache ---

    def test_figure_cache_reuses_figure(self):
        sp = self._make()
        sp.save_depth_frame(frame=0, dpi=20)
        sp.save_depth_frame(frame=1, dpi=20)
        # Both frames with identical kwargs → same cached figure object
        self.assertEqual(len(sp._figure_cache), 1)

    def test_clear_figure_cache(self):
        sp = self._make()
        sp.save_depth_frame(frame=0, dpi=20)
        self.assertGreater(len(sp._figure_cache), 0)
        sp._clear_figure_cache()
        self.assertEqual(len(sp._figure_cache), 0)

    def test_closed_figure_is_not_reused(self):
        # save_* closes its figure; drawing on it again fails under TkAgg
        # (anuga_sww_gui), and a notebook would not show it.
        import matplotlib.pyplot as plt
        sp = self._make()
        sp.save_depth_frame(frame=0, dpi=20)
        fig, _ = sp.plot_depth_frame(frame=1, dpi=20)
        self.assertIn(fig.number, plt.get_fignums())
        self.assertIs(plt.figure(fig.number), fig)
        plt.close(fig)

    def test_save_max_frames_close_their_figure(self):
        import matplotlib.pyplot as plt
        sp = self._make()
        before = plt.get_fignums()
        for qty in ('depth', 'stage', 'speed', 'speed_depth'):
            getattr(sp, 'save_max_%s_frame' % qty)(dpi=20)
            self.assertEqual(plt.get_fignums(), before)

    def test_frames_save_repeatedly_with_tkagg(self):
        # What anuga_sww_gui does: max_* and static elev frames are drawn in
        # the GUI process, where pyplot figures get a Tk toolbar.
        if os.name != 'nt' and not os.environ.get('DISPLAY'):
            self.skipTest('needs a display for TkAgg')
        try:
            import tkinter
            tkinter.Tk().destroy()
        except Exception as e:
            self.skipTest('Tk is not usable: %s' % e)
        sww = self.swwfile
        script = (
            "import matplotlib; matplotlib.use('TkAgg')\n"
            "from anuga.utilities.animate import SWW_plotter\n"
            "sp = SWW_plotter(%r, plot_dir=%r)\n"
            "for _ in range(2):\n"
            "    sp.save_max_depth_frame(dpi=20)\n"
            "    sp.save_elev_frame(dpi=20)\n"
            "    sp.save_depth_frame(frame=0, dpi=20)\n"
            % (sww, self.plot_dir))
        import subprocess
        import sys
        result = subprocess.run([sys.executable, '-c', script],
                                capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stderr)

    # --- plot_dir=None writes to cwd ---

    def test_plot_dir_none_writes_to_cwd(self):
        sp = SWW_plotter(self.swwfile, plot_dir=None)
        old_cwd = os.getcwd()
        try:
            os.chdir(self.tmpdir)
            sp.save_depth_frame(frame=0, dpi=20)
            png = os.path.join(self.tmpdir, 'test_depth_0000000000.png')
            self.assertTrue(os.path.exists(png))
        finally:
            os.chdir(old_cwd)


# ---------------------------------------------------------------------------
# SWW_plotter — elev_delta (time-varying elevation)
# ---------------------------------------------------------------------------

def _write_sww_with_time_varying_elev(path, n_timesteps=3):
    """Write a SWW file where elevation_c is 2-D (time-varying, e.g. erosion)."""
    from anuga.file.netcdf import NetCDFFile
    from anuga.config import netcdf_mode_w

    n_nodes = 4
    n_tris = 2

    x    = np.array([0.0, 2.0, 2.0, 0.0], dtype=np.float32)
    y    = np.array([0.0, 0.0, 1.0, 1.0], dtype=np.float32)
    vols = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int32)
    t    = np.array([0.0, 10.0, 20.0], dtype=np.float64)[:n_timesteps]

    # elevation changes over time: erodes by 0.1 m/step on first triangle
    elev_t0 = np.array([-1.0, -0.5], dtype=np.float32)
    elev = np.zeros((n_timesteps, n_tris), dtype=np.float32)
    for i in range(n_timesteps):
        elev[i, 0] = elev_t0[0] - 0.1 * i   # first tri erodes
        elev[i, 1] = elev_t0[1]              # second tri static

    stage = np.zeros((n_timesteps, n_tris), dtype=np.float32)
    xmom  = np.zeros((n_timesteps, n_tris), dtype=np.float32)
    ymom  = np.zeros((n_timesteps, n_tris), dtype=np.float32)

    fid = NetCDFFile(path, netcdf_mode_w)
    fid.starttime = 0.0
    fid.timezone  = 'UTC'
    fid.xllcorner = 300000.0
    fid.yllcorner = 6000000.0
    fid.zone      = 55

    fid.createDimension('number_of_volumes',   n_tris)
    fid.createDimension('number_of_vertices',  3)
    fid.createDimension('number_of_points',    n_nodes)
    fid.createDimension('number_of_timesteps', n_timesteps)

    fid.createVariable('x',            'f', ('number_of_points',))
    fid.createVariable('y',            'f', ('number_of_points',))
    fid.createVariable('volumes',      'i', ('number_of_volumes', 'number_of_vertices'))
    fid.createVariable('time',         'd', ('number_of_timesteps',))
    fid.createVariable('elevation_c',  'f', ('number_of_timesteps', 'number_of_volumes'))
    fid.createVariable('stage_c',      'f', ('number_of_timesteps', 'number_of_volumes'))
    fid.createVariable('xmomentum_c',  'f', ('number_of_timesteps', 'number_of_volumes'))
    fid.createVariable('ymomentum_c',  'f', ('number_of_timesteps', 'number_of_volumes'))

    fid.variables['x'][:]           = x
    fid.variables['y'][:]           = y
    fid.variables['volumes'][:]     = vols
    fid.variables['time'][:]        = t
    fid.variables['elevation_c'][:] = elev
    fid.variables['stage_c'][:]     = stage
    fid.variables['xmomentum_c'][:] = xmom
    fid.variables['ymomentum_c'][:] = ymom
    fid.close()

    return elev   # return for assertion use


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestElevDelta(unittest.TestCase):
    """Tests for the elev_delta property and save_elev_delta_frame method."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.plot_dir = os.path.join(self.tmpdir, 'frames')
        os.makedirs(self.plot_dir, exist_ok=True)
        # static-elev SWW
        self.static_sww = os.path.join(self.tmpdir, 'static.sww')
        _write_minimal_sww(self.static_sww)
        # time-varying-elev SWW
        self.dynamic_sww = os.path.join(self.tmpdir, 'dynamic.sww')
        self._elev_array = _write_sww_with_time_varying_elev(self.dynamic_sww)

    def tearDown(self):
        import matplotlib.pyplot as plt
        plt.close('all')
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_elev_delta_is_none_for_static_elev(self):
        sp = SWW_plotter(self.static_sww, plot_dir=self.plot_dir)
        self.assertIsNone(sp.elev_delta)

    def test_elev_delta_shape_for_time_varying_elev(self):
        sp = SWW_plotter(self.dynamic_sww, plot_dir=self.plot_dir)
        self.assertEqual(sp.elev.ndim, 2)
        delta = sp.elev_delta
        self.assertIsNotNone(delta)
        self.assertEqual(delta.shape, sp.elev.shape)

    def test_elev_delta_values_correct(self):
        sp = SWW_plotter(self.dynamic_sww, plot_dir=self.plot_dir)
        delta = sp.elev_delta
        # frame 0: delta should be zero everywhere
        np.testing.assert_allclose(delta[0, :], 0.0, atol=1e-6)
        # frame 1: first tri eroded by 0.1
        self.assertAlmostEqual(float(delta[1, 0]), -0.1, places=5)
        self.assertAlmostEqual(float(delta[1, 1]),  0.0, places=5)
        # frame 2: first tri eroded by 0.2
        self.assertAlmostEqual(float(delta[2, 0]), -0.2, places=5)

    def test_save_elev_delta_frame_creates_png(self):
        sp = SWW_plotter(self.dynamic_sww, plot_dir=self.plot_dir)
        sp.save_elev_delta_frame(frame=1, dpi=20)
        png = os.path.join(self.plot_dir, 'dynamic_elev_delta_0000000000.png')
        self.assertTrue(os.path.exists(png))

    def test_save_elev_delta_frame_increments_counter(self):
        sp = SWW_plotter(self.dynamic_sww, plot_dir=self.plot_dir)
        sp.save_elev_delta_frame(frame=0, dpi=20)
        sp.save_elev_delta_frame(frame=1, dpi=20)
        self.assertEqual(sp._elev_delta_frame_count, 2)

    def test_save_elev_delta_frame_static_elev_creates_png(self):
        # For static elevation, delta is always zero — should still produce a file
        sp = SWW_plotter(self.static_sww, plot_dir=self.plot_dir)
        sp.save_elev_delta_frame(frame=0, dpi=20)
        png = os.path.join(self.plot_dir, 'static_elev_delta_0000000000.png')
        self.assertTrue(os.path.exists(png))


# ---------------------------------------------------------------------------
# Domain_plotter
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestDomainPlotter(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.plot_dir = os.path.join(self.tmpdir, 'dp_frames')

    def tearDown(self):
        import matplotlib.pyplot as plt
        plt.close('all')
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make(self, **kwargs):
        domain = _make_mock_domain(self.tmpdir)
        return Domain_plotter(domain, plot_dir=self.plot_dir, **kwargs)

    def test_init_loads_geometry(self):
        dp = self._make()
        self.assertEqual(len(dp.x), 4)
        self.assertEqual(len(dp.y), 4)
        self.assertEqual(dp.triangles.shape, (2, 3))

    def test_init_elev_and_stage(self):
        dp = self._make()
        np.testing.assert_allclose(dp.elev,  [-1.0, -0.5])
        np.testing.assert_allclose(dp.stage, [ 0.5,  0.3])

    def test_init_depth_computed(self):
        dp = self._make()
        np.testing.assert_allclose(dp.depth, [1.5, 0.8])

    def test_init_speed_non_negative(self):
        dp = self._make()
        self.assertTrue(np.all(dp.speed >= 0.0))

    def test_init_counters_zero(self):
        dp = self._make()
        self.assertEqual(dp._depth_frame_count, 0)
        self.assertEqual(dp._stage_frame_count, 0)
        self.assertEqual(dp._speed_frame_count, 0)

    def test_derived_quantities_track_domain(self):
        # Regression: depth/xvel/yvel/speed/speed_depth were cached in
        # __init__, so sampling them inside an evolve loop returned the t = 0
        # values for the whole run.
        domain = _make_mock_domain(self.tmpdir)
        dp = Domain_plotter(domain, plot_dir=self.plot_dir)

        np.testing.assert_allclose(dp.depth, [1.5, 0.8])

        # Advance the domain the way evolve() does: in place, on the centroids.
        domain.quantities['stage'].centroid_values[:] = [2.5, 1.3]
        domain.quantities['xmomentum'].centroid_values[:] = [7.0, 0.0]
        domain.quantities['ymomentum'].centroid_values[:] = [0.0, 0.0]

        np.testing.assert_allclose(dp.depth, [3.5, 1.8])
        np.testing.assert_allclose(dp.xvel[0], 2.0)
        np.testing.assert_allclose(dp.speed[0], 2.0)
        np.testing.assert_allclose(dp.speed_depth[0], 7.0)

    def test_init_absolute_mode(self):
        domain = _make_mock_domain(self.tmpdir)
        dp = Domain_plotter(domain, plot_dir=self.plot_dir, absolute=True)
        # In absolute mode, x/y are offset by xllcorner/yllcorner
        self.assertAlmostEqual(dp.x[0], 300000.0, places=0)
        self.assertAlmostEqual(dp.y[0], 6000000.0, places=0)

    def test_save_depth_frame_creates_png(self):
        dp = self._make()
        dp.save_depth_frame(dpi=20)
        domain_name = os.path.basename(dp.domain.get_name())
        png = os.path.join(self.plot_dir,
                           domain_name + '_depth_0000000000.png')
        self.assertTrue(os.path.exists(png),
                        f'Expected {png}; dir={os.listdir(self.plot_dir)}')

    def test_save_depth_frame_increments_counter(self):
        dp = self._make()
        dp.save_depth_frame(dpi=20)
        dp.save_depth_frame(dpi=20)
        self.assertEqual(dp._depth_frame_count, 2)

    def test_plot_mesh_returns_fig_ax(self):
        dp = self._make()
        result = dp.plot_mesh(dpi=20)
        self.assertEqual(len(result), 2)
        fig, ax = result
        self.assertIsNotNone(fig)

    def test_save_stage_frame_creates_png(self):
        dp = self._make()
        dp.save_stage_frame(dpi=20)
        domain_name = os.path.basename(dp.domain.get_name())
        png = os.path.join(self.plot_dir,
                           domain_name + '_stage_0000000000.png')
        self.assertTrue(os.path.exists(png))

    def test_save_speed_frame_creates_png(self):
        dp = self._make()
        dp.save_speed_frame(dpi=20)
        domain_name = os.path.basename(dp.domain.get_name())
        png = os.path.join(self.plot_dir,
                           domain_name + '_speed_0000000000.png')
        self.assertTrue(os.path.exists(png))


class Test_tile_client_config(unittest.TestCase):
    """ANUGA must identify itself to tile servers (no network needed).

    contextily defaults to USER_AGENT = "contextily-" + uuid4().hex, i.e. a
    fresh random agent per process.  OpenStreetMap's tile usage policy requires
    a stable agent identifying the application and blocks rotating ones -- and
    the block arrives as an HTTP 200 tile whose image reads "Access blocked",
    so it silently renders onto the map instead of raising.
    """

    def setUp(self):
        try:
            import contextily  # noqa: F401
        except ImportError:
            self.skipTest('contextily is not installed')

    def test_contextily_default_agent_is_replaced(self):
        import contextily as cx
        import contextily.tile as tile
        from anuga.utilities import animate

        original = tile.USER_AGENT
        try:
            tile.USER_AGENT = 'contextily-0123456789abcdef'
            animate._configure_tile_client(cx)
            self.assertIn('ANUGA', tile.USER_AGENT)
            self.assertFalse(tile.USER_AGENT.startswith('contextily-'))
        finally:
            tile.USER_AGENT = original

    def test_user_supplied_agent_is_left_alone(self):
        """Only contextily's own default is replaced, never a deliberate one."""
        import contextily as cx
        import contextily.tile as tile
        from anuga.utilities import animate

        original = tile.USER_AGENT
        try:
            tile.USER_AGENT = 'MyApp/1.0 (+https://example.org)'
            animate._configure_tile_client(cx)
            self.assertEqual(tile.USER_AGENT, 'MyApp/1.0 (+https://example.org)')
        finally:
            tile.USER_AGENT = original


if __name__ == '__main__':
    unittest.main()


class Test_set_plotter_returns_the_plotter(unittest.TestCase):
    """`domain.set_plotter()` hands back the plotter it attaches, so the
    common two-line form collapses to one."""

    def _domain(self):
        import anuga
        d = anuga.rectangular_cross_domain(6, 6, len1=10.0, len2=10.0)
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 1.0)
        return d

    def test_it_returns_the_plotter_it_attached(self):
        import anuga
        d = self._domain()
        plotter = d.set_plotter(absolute=True)
        assert isinstance(plotter, anuga.Domain_plotter)
        assert plotter is d.dplotter

    def test_a_second_call_returns_the_same_plotter(self):
        """The plotter is built once; later calls hand back the same one."""
        d = self._domain()
        first = d.set_plotter()
        assert d.set_plotter() is first

    def test_the_domain_still_mirrors_the_plotter(self):
        """Returning it does not stop the arrays and methods being copied
        onto the domain, which existing scripts rely on."""
        d = self._domain()
        plotter = d.set_plotter()
        assert d.triang is plotter.triang
        assert d.plot_depth_frame == plotter.plot_depth_frame


class Test_domain_derived_plot_quantities(unittest.TestCase):
    """depth, xvel, yvel, speed and speed_depth are properties of the domain
    that read through to the plotter, so they follow the run. Copying their
    values in set_plotter froze them at the moment it was called."""

    def _domain(self):
        import anuga
        d = anuga.rectangular_cross_domain(6, 6, len1=10.0, len2=10.0)
        d.store = False
        d.set_quantity('elevation', 0.0)
        d.set_quantity('stage', 1.0)
        d.set_quantity('xmomentum', 0.5)
        d.set_boundary({t: anuga.Reflective_boundary(d)
                        for t in d.get_boundary_tags()})
        return d

    def test_they_exist_and_match_the_plotter(self):
        d = self._domain()
        p = d.set_plotter()
        for name in ('depth', 'xvel', 'yvel', 'speed', 'speed_depth'):
            assert np.allclose(getattr(d, name), getattr(p, name)), name

    def test_they_follow_the_run(self):
        d = self._domain()
        p = d.set_plotter()
        before = {n: getattr(d, n).copy()
                  for n in ('depth', 'xvel', 'speed')}
        for _ in d.evolve(yieldstep=1.0, finaltime=2.0):
            pass
        for n, was in before.items():
            assert not np.allclose(getattr(d, n), was), \
                '%s did not change over the run' % n
            assert np.allclose(getattr(d, n), getattr(p, n)), n

    def test_without_a_plotter_they_say_so(self):
        d = self._domain()
        with self.assertRaises(AttributeError) as cm:
            d.depth
        assert 'set_plotter' in str(cm.exception)
