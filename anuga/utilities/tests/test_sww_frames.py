"""SWW_plotter's quantities, read from the file as they are used.

SWW_plotter used to hold every timestep of eight (time, triangle) arrays in
memory.  They are now FrameArrays that read a timestep, or a few triangles'
timeseries, at a time.  These tests check every way the plotter and
anuga_sww_gui use them gives the same values as the arrays did.
"""

import multiprocessing
import platform

import numpy as np
import pytest

pytest.importorskip('matplotlib')
pytest.importorskip('netCDF4')

from anuga.utilities import _sww_frames
from anuga.utilities.animate import SWW_plotter

NTIME, NTRI = 7, 40
MIN_DEPTH = 0.001


def _write_sww(path, animated_elevation=False):
    """A random SWW with wet and dry triangles, and a tracer."""
    from anuga.file.netcdf import NetCDFFile
    from anuga.config import netcdf_mode_w

    rng = np.random.default_rng(1)
    npts = NTRI + 2
    x = rng.uniform(0, 100, npts).astype(np.float32)
    y = rng.uniform(0, 100, npts).astype(np.float32)
    vols = np.array([[i, i + 1, i + 2] for i in range(NTRI)], dtype=np.int32)
    elev = rng.uniform(-1, 1, (NTIME, NTRI) if animated_elevation else NTRI).astype(np.float32)
    depth = np.where(rng.uniform(size=(NTIME, NTRI)) < 0.3, 0.0,
                     rng.uniform(0, 2, (NTIME, NTRI))).astype(np.float32)
    data = {
        'elevation_c': elev,
        'stage_c': (elev + depth).astype(np.float32),
        'xmomentum_c': (depth * rng.uniform(-1, 1, (NTIME, NTRI))).astype(np.float32),
        'ymomentum_c': (depth * rng.uniform(-1, 1, (NTIME, NTRI))).astype(np.float32),
        'dye_c': rng.uniform(0, 1, (NTIME, NTRI)).astype(np.float32),
    }

    fid = NetCDFFile(path, netcdf_mode_w)
    fid.starttime = 0.0
    fid.timezone = 'UTC'
    fid.xllcorner = 0.0
    fid.yllcorner = 0.0
    fid.zone = 55
    fid.tracer_names = 'dye'
    fid.createDimension('number_of_volumes', NTRI)
    fid.createDimension('number_of_vertices', 3)
    fid.createDimension('number_of_points', npts)
    fid.createDimension('number_of_timesteps', NTIME)
    fid.createVariable('x', 'f', ('number_of_points',))[:] = x
    fid.createVariable('y', 'f', ('number_of_points',))[:] = y
    fid.createVariable('volumes', 'i', ('number_of_volumes', 'number_of_vertices'))[:] = vols
    fid.createVariable('time', 'd', ('number_of_timesteps',))[:] = 10.0 * np.arange(NTIME)
    for name, values in data.items():
        dims = ('number_of_timesteps', 'number_of_volumes') if values.ndim == 2 \
            else ('number_of_volumes',)
        fid.createVariable(name, 'f', dims)[:] = values
    fid.close()
    return data


def _reference(data):
    """The arrays as SWW_plotter used to compute them."""
    elev, stage = data['elevation_c'], data['stage_c']
    xmom, ymom = data['xmomentum_c'], data['ymomentum_c']
    depth = stage - elev
    with np.errstate(invalid='ignore', divide='ignore'):
        xvel = np.where(depth > MIN_DEPTH, xmom / depth, 0.0)
        yvel = np.where(depth > MIN_DEPTH, ymom / depth, 0.0)
    speed = np.sqrt(xvel**2 + yvel**2)
    return dict(stage=stage, xmom=xmom, ymom=ymom, depth=depth, xvel=xvel,
                yvel=yvel, speed=speed, speed_depth=speed * depth)


@pytest.fixture(params=[False, True], ids=['static-elev', 'animated-elev'])
def plotter(request, tmp_path):
    path = str(tmp_path / 'run.sww')
    data = _write_sww(path, animated_elevation=request.param)
    sp = SWW_plotter(path, plot_dir=None, min_depth=MIN_DEPTH)
    yield sp, data, _reference(data)
    sp.close()


QUANTITIES = ['stage', 'xmom', 'ymom', 'depth', 'xvel', 'yvel', 'speed', 'speed_depth']


@pytest.mark.parametrize('name', QUANTITIES)
def test_a_timestep_matches(plotter, name):
    sp, _, ref = plotter
    for frame in (0, 3, -1):
        np.testing.assert_array_equal(getattr(sp, name)[frame, :], ref[name][frame, :])
        np.testing.assert_array_equal(getattr(sp, name)[frame], ref[name][frame])


@pytest.mark.parametrize('name', QUANTITIES)
def test_a_few_triangles_over_time_match(plotter, name):
    """As the timeseries and cross-section discharge read them."""
    sp, _, ref = plotter
    a = getattr(sp, name)
    for triangles in ([5, 2, 30], [7, 7], np.array([1, 39]), 4, slice(3, 9)):
        np.testing.assert_array_equal(a[:, triangles], ref[name][:, triangles])
    np.testing.assert_array_equal(a[2, [4, 1]], ref[name][2, [4, 1]])


@pytest.mark.parametrize('name', QUANTITIES)
def test_ranges_and_maxima_match(plotter, name):
    sp, _, ref = plotter
    a = getattr(sp, name)
    assert a.min() == ref[name].min()
    assert a.max() == ref[name].max()
    np.testing.assert_array_equal(np.max(a, axis=0), ref[name].max(axis=0))
    np.testing.assert_array_equal(a.min(axis=0), ref[name].min(axis=0))
    np.testing.assert_array_equal(a.max(axis=1), ref[name].max(axis=1))


def test_the_shape_is_time_by_triangles(plotter):
    sp, _, _ = plotter
    assert sp.depth.shape == (NTIME, NTRI)
    assert sp.depth.ndim == 2
    assert len(sp.depth) == NTIME
    assert sp.depth.dtype == np.float32


def test_anything_else_gets_the_whole_array(plotter):
    """Numpy functions, arithmetic and ndarray methods still work."""
    sp, _, ref = plotter
    np.testing.assert_array_equal(np.asarray(sp.depth), ref['depth'])
    np.testing.assert_array_equal(sp.depth[1:4], ref['depth'][1:4])
    np.testing.assert_array_equal(sp.stage - sp.depth, ref['stage'] - ref['depth'])
    np.testing.assert_array_equal(sp.speed * 2.0, ref['speed'] * 2.0)
    np.testing.assert_allclose(sp.depth.mean(), ref['depth'].mean())
    np.testing.assert_array_equal(np.percentile(sp.speed, 90), np.percentile(ref['speed'], 90))


def test_the_maximum_cannot_be_changed_by_its_caller(plotter):
    sp, _, ref = plotter
    m = sp.depth.max(axis=0)
    m[:] = -99.0
    np.testing.assert_array_equal(sp.max_depth, ref['depth'].max(axis=0))


def test_the_plotter_max_properties_match(plotter):
    sp, _, ref = plotter
    np.testing.assert_array_equal(sp.max_depth, ref['depth'].max(axis=0))
    np.testing.assert_array_equal(sp.max_speed, ref['speed'].max(axis=0))
    np.testing.assert_array_equal(sp.max_speed_depth, ref['speed_depth'].max(axis=0))
    np.testing.assert_array_equal(sp.max_stage, ref['stage'].max(axis=0))


def test_tracers_and_elevation_match(plotter):
    sp, data, _ = plotter
    np.testing.assert_array_equal(sp.tracers['dye'][2], data['dye_c'][2])
    np.testing.assert_array_equal(sp.tracer_dye[:, [3, 1]], data['dye_c'][:, [3, 1]])
    assert sp.tracer_dye.max() == data['dye_c'].max()
    elev = data['elevation_c']
    if elev.ndim == 1:
        assert isinstance(sp.elev, np.ndarray)
        np.testing.assert_array_equal(sp.elev, elev)
        assert sp.elev_delta is None
    else:
        np.testing.assert_array_equal(sp.elev[2], elev[2])
        np.testing.assert_array_equal(sp.elev[:, [6]], elev[:, [6]])
        np.testing.assert_array_equal(sp.elev_delta[3], elev[3] - elev[0])
        assert sp.elev_delta.min() == (elev - elev[0]).min()


def test_water_volume_and_discharge_still_work(plotter):
    sp, data, ref = plotter

    class Segment:
        def __init__(self, triangle_id):
            self.triangle_id, self.normal, self.length = triangle_id, (0.6, 0.8), 2.0

    class Mesh:
        """Unit areas, and a cross-section through some triangles, one
        twice: the file has no vertex quantities to build a mesh from."""
        areas = np.ones(NTRI)

        def get_intersecting_segments(self, polyline, verbose=False):
            return [Segment(t) for t in (9, 3, 9, 20)]

    sp.mesh = Mesh()
    stage, elev = data['stage_c'], data['elevation_c']
    expected = ((stage - elev) * (stage > elev)).sum(axis=1)
    np.testing.assert_allclose(sp.water_volume(), expected, rtol=1e-6)

    _, Q = sp.get_flow_through_cross_section([[0.0, 50.0], [100.0, 50.0]])
    ids = [9, 3, 9, 20]
    expected = ((0.6 * ref['xmom'][:, ids] + 0.8 * ref['ymom'][:, ids]) * 2.0).sum(axis=1)
    np.testing.assert_allclose(Q, expected, rtol=1e-6)


def test_opening_the_plotter_reads_no_timesteps(tmp_path, monkeypatch):
    """The point of the change: nothing time-varying is read up front."""
    path = str(tmp_path / 'run.sww')
    _write_sww(path)
    reads = []
    for method in ('frame', 'columns'):
        original = getattr(_sww_frames.SWWSource, method)
        monkeypatch.setattr(_sww_frames.SWWSource, method,
                            lambda self, *a, _o=original: reads.append(a) or _o(self, *a))
    sp = SWW_plotter(path, plot_dir=None)
    assert reads == []
    sp.depth[2]
    assert reads, 'reading a timestep should go through SWWSource'
    sp.close()


def test_a_closed_plotter_opens_the_file_again(plotter):
    sp, _, ref = plotter
    sp.depth[1]
    sp.close()
    np.testing.assert_array_equal(sp.speed[5], ref['speed'][5])


def _depth_in_child(sp, queue):
    queue.put(np.asarray(sp.depth[4]).tolist())


@pytest.mark.skipif(platform.system() != 'Linux', reason='fork start method')
def test_a_forked_process_opens_its_own_file(plotter):
    """anuga_sww_gui's animation workers are forked on Linux."""
    sp, _, ref = plotter
    sp.depth[4]                         # the parent has the file open
    ctx = multiprocessing.get_context('fork')
    queue = ctx.Queue()
    child = ctx.Process(target=_depth_in_child, args=(sp, queue))
    child.start()
    result = queue.get(timeout=30)
    child.join(timeout=30)
    np.testing.assert_array_equal(np.array(result, dtype=np.float32), ref['depth'][4])
    np.testing.assert_array_equal(sp.depth[4], ref['depth'][4])
