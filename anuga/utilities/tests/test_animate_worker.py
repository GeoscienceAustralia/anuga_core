"""The parallel frame workers must render every quantity the GUI offers.

anuga_sww_gui renders a multi-frame animation in worker processes, which
number each frame by setting the plotter's frame counter before saving.
Tracer counters live in a dict rather than an attribute per quantity, so
the workers need to set those differently.
"""

import glob
import os

import pytest

import anuga
from anuga.utilities import _animate_worker

pytest.importorskip('matplotlib')
pytest.importorskip('netCDF4')


def _tracer_run(tmp_path, name='tw'):
    d = anuga.rectangular_cross_domain(4, 4, len1=40.0, len2=40.0)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    d.add_tracer('dye', initial_value=0.5)
    d.set_name(name)
    d.set_datadir(str(tmp_path))
    for _ in d.evolve(yieldstep=0.5, finaltime=1.0):
        pass
    return str(tmp_path / (name + '.sww'))


def _render(sww, plot_dir, qty, out_pos, frame=1):
    _animate_worker.worker_init(sww, str(plot_dir), 0.001, None)
    try:
        return _animate_worker.worker_frame(
            frame, out_pos, 'save_%s_frame' % qty, qty,
            dpi=40, vmin=0.0, vmax=1.0, cmap='viridis', basemap=False,
            alpha=1.0, basemap_provider=None)
    finally:
        _animate_worker._wp = None


@pytest.mark.parametrize('qty', ['depth', 'tracer_dye'])
def test_a_worker_numbers_the_frame_it_is_given(tmp_path, qty):
    sww = _tracer_run(tmp_path)
    plot_dir = tmp_path / 'frames'
    plot_dir.mkdir()
    assert _render(sww, plot_dir, qty, out_pos=3) == 3
    frames = glob.glob(os.path.join(str(plot_dir), '*.png'))
    assert [os.path.basename(f)[-14:] for f in frames] == ['0000000003.png'], \
        'expected one frame numbered 3, got %s' % frames


def test_a_quantity_without_a_counter_is_reported(tmp_path):
    sww = _tracer_run(tmp_path)
    with pytest.raises(KeyError, match='no_such_quantity'):
        _render(sww, tmp_path, 'no_such_quantity', out_pos=0)
